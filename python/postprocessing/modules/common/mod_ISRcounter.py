import ROOT
import os
import numpy as np
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

# inspired by https://github.com/CERN-PH-CMG/cmgtools-lite/blob/759bdc4213c50db48cb695ae498f7a97794a1410/TTHAnalysis/python/analyzers/nIsrAnalyzer.py

class ISRcounter(Module):

    def __init__(self, isData=False):
        self.isData             = isData

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nISR", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def deltaPhi(self, phi1, phi2):
        dphi = phi2-phi1
        if  dphi > math.pi:
            dphi -= 2.0*math.pi
        if dphi <= -math.pi:
            dphi += 2.0*math.pi
        return abs(dphi)

    def deltaR2(self, l1, l2):
        return self.deltaPhi(l1.phi, l2.phi)**2 + (l1.eta - l2.eta)**2

    def deltaR(self, l1, l2):
        return math.sqrt(self.deltaR2(l1,l2))

    def ptRatio(self, l1, l2):
	return float( l1.pt) / float(l2.pt)


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets        = Collection(event, "Jet")
        genParts    = Collection(event, "GenPart")
        electrons   = Collection(event, "Electron")
        muons       = Collection(event, "Muon")

        # get the daughter indices for each generated particle. needed later on
        motherIndices   = []
        for genPart in genParts:
            if genPart.genPartIdxMother > 0:
                genPart.motherId = genParts[genPart.genPartIdxMother].pdgId
            else: genPart.motherId = None
            motherIndices += [genPart.genPartIdxMother]
            genPart.daughterIds = []
            genPart.daughterPdgIds = []

        for i, motherIndex in enumerate(motherIndices):
            if motherIndex > 0:
                genParts[motherIndex].daughterIds += [i]
                genParts[motherIndex].daughterPdgIds += [genParts[i].pdgId]

        # clean jets from leptons

        # use loosely isolated veto electrons
        selectedElectrons = []
        for e in electrons:
		if e.pt>5 and e.cutBased >= 1 and e.dxy < 0.02 and e.dz < 0.1 and ( e.pfRelIso03_all*e.pt < 5 if e.pt<=25 else e.pfRelIso03_all<0.2): selectedElectrons.append(e)

        # use loosely isolated loose muons
        selectedMuons = []
        for m in muons:
            if m.pt>3.5 and m.isPFcand and (m.isTracker or m.isGlobal) and m.dxy < 0.02 and m.dz < 0.1 and ( m.pfRelIso03_all*m.pt < 5 if m.pt<=25 else m.pfRelIso03_all<0.2): selectedMuons.append(m)

        cleanJets = []
        for j in jets:
            clean = True
            if j.pt < 25. or abs(j.eta)>4.7: clean = False
            for coll in [selectedElectrons,selectedMuons]:
                for l in coll:
                    if self.ptRatio(j, l) < 2 and self.deltaR(j, l) < 0.4: clean = False

            if clean:
                cleanJets += [j]

        nISR = 0
        for jet in cleanJets:
            if jet.pt < 30.0: continue
            if abs(jet.eta)>2.4: continue
            matched = False
            for gen in genParts:
                if gen.status != 23 or abs(gen.pdgId)>5: continue
                absMomId = abs(gen.motherId) if gen.motherId is not None else 0
                if not (absMomId==6 or absMomId==23 or absMomId==24 or absMomId==25 or absMomId>1e6): continue
                for daughter in gen.daughterIds:
                    dR = self.deltaR(jet, genParts[daughter])
                    if dR < 0.3:
                        matched = True
                        break
            if not matched:
                nISR += 1

        self.out.fillBranch("nISR", nISR)
        
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

ISRCounter = lambda : ISRcounter( )

