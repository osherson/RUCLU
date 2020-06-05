# RU-Reclusterizer: should run on miniAOD
# Run with:
#> python RUCLU.py pathtofiles(will use all .roots it finds) TAG(Process label from edmDump) N(number of events to process, or -1 if you want to run them all) name (name this run)
# Standard Includes: ----------------------------------------#
import os
import array
from array import *
import glob
import math
import ROOT
from ROOT import *
import sys
import itertools
from itertools import *
from optparse import OptionParser
# Obviously can only be run in a CMSSW Framework
from DataFormats.FWLite import * 
from HLTrigger import *
# -----------------------------------------------------------#

ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Tools:-----------------------------------------------------#
def HardGet(e, L, H): # shorthand def for getting collection from event
	e.getByLabel(L, H)
	if H.isValid() and len(H.product()) > 0: 
		return H.product()
	return False
def RtoCXY(eta, phi):
	if phi >= -0.171: cryphi = (phi)*57.2598
	else: cryphi = (phi+2*math.pi)*57.2598
	cryeta = 85 - (1.479 - math.fabs(eta))*85/1.479
	cryeta = math.copysign(cryeta, eta)
	return (cryeta, cryphi+11)
def GetColType(i):
	i = math.fabs(i)
	if i == 11: return 38
	if i == 13: return 8
	if i == 22: return 33
	if i == 211: return 2
	if i == 130: return 6
	return 18
# -----------------------------------------------------------#

# Objects:---------------------------------------------------#
class CryEE: #individual crystal
	def __init__(self, x, y, E):
		self.x = x
		self.y = y
		self.E = E
class ProtoCluster: #CMS cluster collection, with a location and a set of crystals
	def __init__(self, V4, Crystals):
		self.V4 = V4
		self.Crystals = Crystals	
	def __lt__(self, other): # need this for ordering stuff with .sort()
		return self.V4.E() < other.V4.E()
class Cluster:
	def __init__(self, Cs):
		self.V4 = TLorentzVector()
		self.Cs = Cs
	def addC(self, C):
		self.Cs.append(C)
	def Set(self):
		for c in self.Cs:
			self.V4 += c.V4
	def Center(self):
		maxE = 0.
		cEta = 0
		cPhi = 0
		for c in self.Cs:
			for cc in c.Crystals:
				if cc.E > maxE:
					maxE = cc.E
					cEta = cc.x
					cPhi = cc.y
		return cEta, cPhi
	def GetCrystals(self):
		C = []
		for c in self.Cs:
			for cc in c.Crystals:
				C.append([cc.x,cc.y,cc.E])
		return C
		
# -----------------------------------------------------------#

##############################################################
files = glob.glob( sys.argv[1]+"*root" )
events = Events (files)
print files
PH = Handle("vector<pat::Photon>")
PL = ("slimmedPhotons", "", sys.argv[2])
JH = Handle("vector<pat::Jet>")
JL = ("slimmedJets", "", sys.argv[2])  
GH = Handle("vector<reco::GenParticle>")
GL = ("prunedGenParticles", "", sys.argv[2])
CH = Handle("std::vector<reco::CaloCluster>")
CL = ('reducedEgamma', 'reducedEBEEClusters', sys.argv[2])
HHb = Handle("edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >")
HLb = ("reducedEgamma", "reducedEBRecHits", sys.argv[2])
n = 0
for event in events:
	n+=1
	if not sys.argv[3] == "-1":
		if n > int(sys.argv[3]): break
	Hits = HardGet(event, HLb, HHb)
	Clus = HardGet(event, CL, CH)
	if Clus == False: continue
	if Hits == False: continue
	CLUS = []
	for c in Clus:
		cV = TLorentzVector()
		R = math.sqrt(c.x()**2 + c.y()**2 + c.z()**2)
		px = c.x()/R*c.energy()
		py = c.y()/R*c.energy()
		pz = c.z()/R*c.energy()
		cV.SetPxPyPzE(px,py,pz,c.energy())
		if math.fabs(cV.Eta()) < 1.44:
			pair = c.hitsAndFractions()
			Cry = []
			for i in pair:
				if math.fabs(cV.Eta()) < 1.44:
					eeDI = EBDetId(i[0])
					Cry.append(CryEE(eeDI.ieta(), eeDI.iphi(), Hits.find(i[0]).energy()*i[1]))
			CLUS.append(ProtoCluster(cV, Cry))
	if len(CLUS) < 1: continue
	Clusters = []
	while len(CLUS) > 0:
		CLUS.sort()
		NewCluster = Cluster([CLUS[-1]])
		del CLUS[-1]
		drop = []
		for cl in CLUS:
			if cl.V4.DeltaR(NewCluster.Cs[0].V4) < 0.2625:
				NewCluster.addC(cl)
				drop.append(cl)
		for d in drop: CLUS.remove(d)
		NewCluster.Set()
		Clusters.append(NewCluster)
	Jets = HardGet(event, JL, JH)
	if Jets == False: continue
	m = 0
	VL = TLine(0,-15.,0,15)
	HL = TLine(-15,0,15,0)
	for l in [VL, HL]: l.SetLineStyle(2)
	for C in Clusters:
		print "-----------"
		m+=1
		ceta, cphi = C.Center()
		H = TH2F("H"+str(n)+"_"+str(m), "(#eta = "+str(C.V4.Eta())+", #phi = "+str(C.V4.Phi())+");#Delta i#eta;#Delta i#phi",31, -15.5, 15.5, 31, -15.5, 15.5)
		H.SetStats(0)
		H2 = TH2F("H2"+str(n)+"_"+str(m), "(#eta = "+str(C.V4.Eta())+", #phi = "+str(C.V4.Phi())+");#Delta i#eta;#Delta i#phi",31, -15.5, 15.5, 31, -15.5, 15.5)
		H2.SetStats(0)
		for c in C.GetCrystals():
			H.Fill(c[0]-ceta, c[1]-cphi, c[2])
		Circles = []
		for j in Jets:
			allconstituents = []
			for constituent in j.daughterPtrVector() :
				if constituent.numberOfDaughters() > 0 :
					allconstituents.extend([x for x in constituent.daughterPtrVector()])
				else:
					allconstituents.append(constituent)
			for x in allconstituents:
				jV = TLorentzVector()
				jV.SetPtEtaPhiE(x.pt(), x.eta(), x.phi(), x.energy())
				if jV.DeltaR(C.V4) < 0.2625:
						ieta, iphi = RtoCXY(x.eta(), x.phi())
						H2.Fill(ieta-ceta, iphi-cphi, x.energy())
						circ = TEllipse(ieta-ceta, iphi-cphi, 1.33, 1.33)
						circ.SetFillStyle(0)
						circ.SetLineWidth(2)
						circ.SetLineColor(GetColType(x.pdgId()))
						Circles.append(circ)

		C = TCanvas("C"+str(n)+"_"+str(m), "C"+str(n)+"_"+str(m), 1200, 700)
		C.Divide(2,1)
		C.cd(1)
		H.Draw("colz")
		HL.Draw("same")
		VL.Draw("same")
		C.cd(2)
		H2.Draw("colz")
		VL.Draw("same")
		HL.Draw("same")
		for c in Circles: c.Draw("same")
		C.Print("qcd_"+str(n)+"_"+str(m)+".png")
		
		
		
		
		
		
		
		
