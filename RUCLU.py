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
	def Match(self, GP):
		self.Ms = []
		for g in GP:
			if g.numberOfDaughters() == 0:
				gL = TLorentzVector()
				gL.SetPtEtaPhiE(g.pt(), g.eta(), g.phi(), g.energy())
				if gL.DeltaR(self.V4) < 0.2625:
					self.Ms.append(TRANSLATE(g.pdgId()))
# -----------------------------------------------------------#

# Plotting:--------------------------------------------------#
class ECALRollout:
	def __init__(self, name):
		self.name = name
		self.Frame = TH2F(name, ";crystal #eta; crystal #phi", 195, -97, 97,  374, -16, 378)
		self.Frame.SetStats(0)
		self.Elps = []
		self.Pairs = []
	def AddCrystals(self, C):
		for c in C:
			ceta, cphi = RtoCXY(c.V4.Eta(), c.V4.Phi())
			E = TEllipse(ceta, cphi, 5, 5)
			E.SetFillColor(33)
			E.SetLineColor(33)
			self.Elps.append(E)
			for cc in c.Crystals:
				self.Frame.Fill(cc.x,cc.y,cc.E)
	def AddRUCLUs(self, C):
		for c in C:
			ceta, cphi = RtoCXY(c.V4.Eta(), c.V4.Phi())
			E = TEllipse(ceta, cphi, 15, 15)
			E.SetFillStyle(0)
			E.SetLineColor(kBlue)
			self.Elps.append(E)
	def AddJets(self, v4s):
		for j in v4s:
			ceta, cphi = RtoCXY(j.Eta(), j.Phi())
			E = TEllipse(ceta, cphi, 0.4/0.0175, 0.4/0.0175)
			E.SetFillStyle(0)
			E.SetLineColor(kBlack)
			E.SetLineStyle(2)
			self.Elps.append(E)
	def AddPairs(self, P):
		for p in P:
			etas = []
			phis = []
			for i in p:
				ceta, cphi = RtoCXY(i.Eta(), i.Phi())
				E = TEllipse(ceta, cphi, 2, 2)
				E.SetFillStyle(0)
				E.SetLineColor(kRed)
				self.Pairs.append(E)
				etas.append(ceta)
				phis.append(cphi)
			L = TLine(etas[0], phis[0], etas[1], phis[1])
			L.SetLineColor(kRed)
			L.SetLineStyle(2)
			self.Pairs.append(L)
				
				
				
			
			
	def Print(self):
		Sides = TBox(-85,1,85,361)
		Sides.SetLineColor(33)
		Sides.SetLineStyle(2)
		Sides.SetFillStyle(0)
		C = TCanvas()
		C.cd()
		self.Frame.Draw("col")
		Sides.Draw("same")
		for E in self.Elps: E.Draw("same")
		self.Frame.Draw("colsame")
		for p in self.Pairs: p.Draw("same")
		C.Print(self.name+".png")
		C.Print(self.name+".root")
		
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
	tmpEROL = ECALRollout(sys.argv[4]+"Evt"+str(n))
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
	tmpEROL.AddCrystals(CLUS)
	Clusters = []
	while len(CLUS) > 0:
		CLUS.sort()
		print "working on cluster centered at: "+str(CLUS[-1].V4.Eta())+", "+str(CLUS[-1].V4.Phi())+" (e = "+str(CLUS[-1].V4.E())+")" 
		NewCluster = Cluster([CLUS[-1]])
		del CLUS[-1]
		drop = []
		for cl in CLUS:
			print "checking: " +str(cl.V4.Eta())+", "+str(cl.V4.Phi())+" (e = "+str(cl.V4.E())+"), DeltaR = " + str(cl.V4.DeltaR(NewCluster.Cs[0].V4))
			if cl.V4.DeltaR(NewCluster.Cs[0].V4) < 0.2625:
				NewCluster.addC(cl)
				drop.append(cl)
		for d in drop: CLUS.remove(d)
		NewCluster.Set()
		Clusters.append(NewCluster)
	tmpEROL.AddRUCLUs(Clusters)
	Jets = HardGet(event, JL, JH)
	Photons = HardGet(event, PL, PH)
	if Jets == False: continue
	if Photons == False: continue
	myJets = []
	for j in Jets:
		jV = TLorentzVector()
		if math.fabs(j.eta()) < 1.44:
			jV.SetPtEtaPhiE(j.pt(), j.eta(), j.phi(), j.energy())
			keep = True
			for p in Photons:
				pV = TLorentzVector()
				pV.SetPtEtaPhiE(p.pt(), p.eta(), p.phi(), p.energy())
				if pV.DeltaR(jV) < 0.4 and p.energy() > j.energy()*0.75:
					keep = False
					break
			if keep: myJets.append(jV)
	tmpEROL.AddJets(myJets)
	DauPairs = []
	GP = HardGet(event, GL, GH)
	for ig in GP:
		if math.fabs(ig.eta()) < 1.44 and ig.numberOfDaughters() == 2:
			if ig.daughter(0).pdgId() == 22 and ig.daughter(1).pdgId() == 22:
				D1 = TLorentzVector()
				D1.SetPtEtaPhiE(ig.daughter(0).pt(), ig.daughter(0).eta(), ig.daughter(0).phi(), ig.daughter(0).energy())
				D2 = TLorentzVector()
				D2.SetPtEtaPhiE(ig.daughter(1).pt(), ig.daughter(1).eta(), ig.daughter(1).phi(), ig.daughter(1).energy())
				DauPairs.append([D1,D2])
	tmpEROL.AddPairs(DauPairs)
	tmpEROL.Print()
		
		
		
		
		
		
		
