import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ElectronScaleSmear(Module):
    def __init__(self,fileName,doApplyCor=False,doSysVar=True,isData=False,run=None):
        self.doApplyCor = doApplyCor
        self.doSysVar = doSysVar
        self.fileName = fileName
        self.run = run
        self.isData = isData
        ROOT.gROOT.cd()

        # Load the class EnergyScaleCorrection
        #pathToLib = '%s/src/PhysicsTools/NanoAODTools/src/' % os.environ['CMSSW_BASE']
        pathToLib = '%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/elecES/' % os.environ['CMSSW_BASE']
        try:
          ROOT.gSystem.Load("libPhysicsToolsNanoAODTools")
          dummy = ROOT.EnergyScaleCorrection
        #Load it via ROOT ACLIC. NB: this creates the object file in the CMSSW directory,
        #causing problems if many jobs are working from the same CMSSW directory
        except Exception as e:
          print "Could not load module via python, trying via ROOT", e
          if "/EnergyScaleCorrection_cc.so" not in ROOT.gSystem.GetLibraries():
            plib = '%s/EnergyScaleCorrection.cc+' % pathToLib
            print 'Loading C++ library from ' + plib
            ROOT.gROOT.ProcessLine('.L ' + plib)
            dummy = ROOT.EnergyScaleCorrection
        self.eleCorr = ROOT.EnergyScaleCorrection(self.fileName, ROOT.EnergyScaleCorrection.ECALELF) 
  

    def loadHisto(self,filename,hname):
      pass

    def beginJob(self):
      pass

    def endJob(self):
      pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
      if outputFile: 
        outputFile.cd()
      self.out = wrappedOutputTree
      if self.doApplyCor:
        self.out.branch("Electron_pt_nom","F", lenVar='nElectron')
        self.out.branch("Electron_mass_nom","F", lenVar='nElectron')
      if self.doSysVar:
        self.out.branch("Electron_pt_up","F", lenVar='nElectron')
        self.out.branch("Electron_pt_down","F", lenVar='nElectron')
        self.out.branch("Electron_mass_up","F", lenVar='nElectron')
        self.out.branch("Electron_mass_down","F", lenVar='nElectron')
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
      """process event, return True (go to next module) or False (fail, go to next event)"""
      run = long(event.run) if self.run == None else self.run
      if self.doApplyCor:
        elecPtNom = []
        elecMassNom = []
      if self.doSysVar:
        elecPtUp = []
        elecPtDo = []
        elecMassUp = []
        elecMassDo = []
      elec = Collection(event, 'Electron')
      nrandom = ROOT.gRandom.Gaus(0,1)
      for el in elec:
        p = ROOT.TLorentzVector()
        p.SetPtEtaPhiM(el.pt, el.eta, el.phi, el.mass)
        r9 = el.r9
        seed = el.seedGain if hasattr(el, 'seedGain') else 12
        ecor = el.eCorr if hasattr(el, 'eCorr') else 1
        et = p.Et()
        abseta = abs(p.Eta()+el.deltaEtaSC)
        if self.isData:
          if not self.doApplyCor: return True
          # Apply scale to data
          praw = p*(1./ecor) # Correct back
          et = praw.Et()
          abseta = abs(praw.Eta()+el.deltaEtaSC)
          escale = self.eleCorr.scaleCorr(run, et, abseta, r9)
          vEle = praw*escale
          elecPtNom.append(vEle.Pt())
          elecMassNom.append(vEle.M())
        else: # MC
          if self.doApplyCor: # Apply smear to MC
            praw = p*(1./ecor) # Correct back
            et = praw.Et()
            abseta = abs(praw.Eta()+el.deltaEtaSC)
            eleSmear = self.eleCorr.smearingSigma(run, et, abseta, r9, 12, 0, 0.)
            vEle = praw*(1+eleSmear*nrandom)
            elecPtNom.append(vEle.Pt())
            elecMassNom.append(vEle.M())
            et = vEle.Et()
          if self.doSysVar:
            eleSmear    = self.eleCorr.smearingSigma(run, et, abseta, r9, 12, 0, 0.)
            escaleErr   = self.eleCorr.scaleCorrUncert(run, et, abseta, r9)
            eleSmearUp  = self.eleCorr.smearingSigma(run, et, abseta, r9, 12,  1, 0.)
            eleSmearDo  = self.eleCorr.smearingSigma(run, et, abseta, r9, 12, -1, 0.)
            eleSmearUnc = nrandom*np.sqrt( (eleSmearUp-eleSmear)*(eleSmearUp-eleSmear) + (eleSmearDo-eleSmear)*(eleSmearDo-eleSmear) )
            vEleUp = p*(1+escaleErr+eleSmearUnc)
            vEleDo = p*(1-escaleErr-eleSmearUnc)
            elecPtUp.append(vEleUp.Pt())
            elecPtDo.append(vEleDo.Pt())
            elecMassUp.append(vEleUp.M())
            elecMassDo.append(vEleDo.M())
      if self.doApplyCor:
        self.out.fillBranch("Electron_pt_nom", elecPtNom)
        self.out.fillBranch("Electron_mass_nom", elecMassNom)
      if self.doSysVar:
        self.out.fillBranch("Electron_pt_up",     elecPtUp)
        self.out.fillBranch("Electron_pt_down",   elecPtDo)
        self.out.fillBranch("Electron_mass_up",   elecMassUp)
        self.out.fillBranch("Electron_mass_down", elecMassDo)
      return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

fname2016="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/elecES/Legacy2016_07Aug2017_FineEtaR9_v3_ele_unc" % os.environ['CMSSW_BASE']
fname2017="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/elecES/Run2017_17Nov2017_v1_ele_unc" % os.environ['CMSSW_BASE']
fname2018="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/elecES/Run2018_Step2Closure_CoarseEtaR9Gain_v2" % os.environ['CMSSW_BASE']
fname2017LowPU="%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/elecES/Run2017_LowPU_v2" % os.environ['CMSSW_BASE']

elecUnc2016MC    = lambda : ElectronScaleSmear(fname2016, doApplyCor=False, doSysVar=True, isData=False)
elecUnc2017MC    = lambda : ElectronScaleSmear(fname2017, doApplyCor=False, doSysVar=True, isData=False)
elecUnc2018MC    = lambda : ElectronScaleSmear(fname2018, doApplyCor=False, doSysVar=True, isData=False)
elecScale2017Data= lambda : ElectronScaleSmear(fname2017, doApplyCor=True,  doSysVar=False,isData=True)

elecScale5TeVData     = lambda : ElectronScaleSmear(fname2017LowPU, doApplyCor=True, doSysVar=False, isData=True, run=306936)
elecScale2017LowPuData= lambda : ElectronScaleSmear(fname2017LowPU, doApplyCor=True, doSysVar=False, isData=True)
elecScale5TeVMC       = lambda : ElectronScaleSmear(fname2017LowPU, doApplyCor=True, doSysVar=True, isData=False, run=306936)
elecScale2017LowPuMC  = lambda : ElectronScaleSmear(fname2017LowPU, doApplyCor=True, doSysVar=True, isData=False)
