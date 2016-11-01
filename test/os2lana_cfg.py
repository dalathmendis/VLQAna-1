import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('isData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is data?"
    )
options.register('zdecaymode', 'zmumu',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Z->mumu or Z->elel? Choose: 'zmumu' or 'zelel'"
    )
options.register('lepID', 'TIGHT',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "lepton ID? Choose: 'TIGHT' or 'LOOSE'"
    )
options.register('outFileName', 'os2lana.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('doPUReweightingOfficial', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting using official recipe"
    )
options.register('filterSignal', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Select only tZtt or bZbZ modes"
    )
options.register('signalType', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Select one of EvtType_MC_tZtZ, EvtType_MC_tZtH, EvtType_MC_tZbW, EvtType_MC_tHtH, EvtType_MC_tHbW, EvtType_MC_bWbW, EvtType_MC_bZbZ, EvtType_MC_bZbH, EvtType_MC_bZtW, EvtType_MC_bHbH, EvtType_MC_bHtW, EvtType_MC_tWtW" 
    )
options.register('applyLeptonSFs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply lepton SFs to the MC"
    )
options.register('applyBTagSFs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tagging SFs to the MC"
    )
options.register('applyDYNLOCorr', True, ### Set to true only for DY process ### Only EWK NLO k-factor is applied
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply DY EWK k-factor to DY MC"
    )
options.register('FileNames', 'FileNames_QCD_HT1000to1500',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of list of input files"
    )
options.register('optimizeReco', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Optimize mass reconstruction"
    )

options.register('applyHtCorr', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "applyig ht correction"
    )


options.register('doSkim', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Produce skim 1 or 0"
     )


options.setDefault('maxEvents', -1)
options.parseArguments()
print options

hltpaths = []
if options.isData:
  options.filterSignal = False 
  options.signalType = "" 
  options.optimizeReco = False
  options.applyLeptonSFs = False 
  options.applyBTagSFs   = False 
  options.applyDYNLOCorr = False 
  options.applyHtCorr = False
  if options.doSkim == False:
    if options.zdecaymode == "zmumu":
      hltpaths = [
        "HLT_DoubleIsoMu17_eta2p1_v", 
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        #"HLT_DoubleMu8_Mass8_PFHT300_v",
        ]
    elif options.zdecaymode == "zelel":
      hltpaths = [
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v"
        ]
    else:
      sys.exit("!!!Error: Wrong Z decay mode option chosen. Choose either 'zmumu' or 'zelel'!!!") 

if options.filterSignal == True and options.doSkim == False and len(options.signalType) == 0:
  sys.exit("!!!Error: Cannot keep signalType empty when filterSignal switched on!!!")  

process = cms.Process("OS2LAna")

from inputFiles_cfi import * 

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
      #FileNames[options.FileNames]
   # 'root://cms-xrd-global.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160811_185354/0000/B2GEDMNtuple_1.root',     
   # 'root://cms-xrd-global.cern.ch//store/group/phys_b2g/skhi/B2GAnaFW_80X_V2p0_RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/DoubleEG/RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/160914_162605/0001/B2GEDMNtuple_1345.root',
#'root://cms-xrd-global.cern.ch//store/group/phys_b2g/skhi/B2GAnaFW_80X_V2p0_RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/DoubleMuon/RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/160914_162636/0000/B2GEDMNtuple_230.root',

    'root://cms-xrd-global.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0_1/160830_162627/0000/B2GEDMNtuple_1.root',

#'root://cms-xrd-global.cern.ch//store/group/phys_b2g/skhi/B2GAnaFW_80X_V2p0_RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/DoubleMuon/RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/160914_163007/0005/B2GEDMNtuple_5946.root',


#    'root://cms-xrd-global.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160825_160530/0000/B2GEDMNtuple_1.root',
#'root://cms-xrd-global.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160825_162354/0000/B2GEDMNtuple_13.root',
#'root://eoscms.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160825_162354/0000/B2GEDMNtuple_13.root',
#'root://eoscms.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160811_185354/0000/B2GEDMNtuple_1.root',

#'root://eoscms.cern.ch//store/group/phys_b2g/skhi/B2GAnaFW_80X_V2p0_RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/DoubleMuon/RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/160914_162636/0000/B2GEDMNtuple_230.root',   
#'root://cms-xrd-global.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160811_185354/0000/B2GEDMNtuple_1.root',

   # 'root://eoscms.cern.ch//store/group/phys_b2g/B2GAnaFW_80X_V2p0/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0/160811_185354/0000/B2GEDMNtuple_1.root',

#'root://cms-xrd-global.cern.ch//store/group/phys_b2g/skhi/B2GAnaFW_80X_V2p0_RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/DoubleMuon/RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0/160914_162636/0003/B2GEDMNtuple_3452.root',


    #'file:/afs/cern.ch/user/d/devdatta/eos/cms/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160410_205229/0000/B2GEDMNtuple_1.root', 
      #'file:/afs/cern.ch/work/d/devdatta/CMSREL/B2GAnaFW_76X/CMSSW_7_6_3_patch2/src/Analysis/B2GAnaFW/test/B2GEDMNtuple_DYJets_M50Madgraph.root',
      #'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_76X_V1p1_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160401_102909/0000/B2GEDMNtuple_1.root'
      #'root://cms-xrd-global.cern.ch//store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/TT_TuneCUETP8M1noCR_13TeV-powheg-pythia8/B2GAnaFW_76X_V1p1_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160401_092930/0000/B2GEDMNtuple_13.root'
    ) 
    )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Output Report                                                                          
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Analysis.VLQAna.EventCleaner_cff") 
process.evtcleaner.isData = options.isData 
process.evtcleaner.hltPaths = cms.vstring (hltpaths)  
process.evtcleaner.DoPUReweightingOfficial = cms.bool(options.doPUReweightingOfficial)  
#process.evtcleaner.storeLHEWts = options.storeLHEWts

from Analysis.VLQAna.OS2LAna_cfi import * 

### Z candidate and jet selections 
process.ana = ana.clone(
    filterSignal = cms.bool(options.filterSignal),
    signalType = cms.string(options.signalType),
    zdecayMode = cms.string(options.zdecaymode),
    applyLeptonSFs = cms.bool(options.applyLeptonSFs),
    applyBTagSFs = cms.bool(options.applyBTagSFs),
    applyHtCorr = cms.bool(options.applyHtCorr),
    applyDYNLOCorr = cms.bool(options.applyDYNLOCorr),
    optimizeReco = cms.bool(options.optimizeReco),
    doSkim       = cms.bool(options.doSkim),
    )
process.ana.elselParams.elidtype = cms.string(options.lepID)
process.ana.muselParams.muidtype = cms.string(options.lepID)
process.ana.muselParams.muIsoMax = cms.double(0.15)
process.ana.lepsfsParams.lepidtype = cms.string(options.lepID)
process.ana.lepsfsParams.zdecayMode = cms.string(options.zdecaymode)
#process.ana.BoostedZCandParams.ptMin = cms.double(80.)
process.ana.jetAK8selParams.jetPtMin = cms.double(200) 
process.ana.jetAK4BTaggedselParams.jetPtMin = cms.double(50) 
process.ana.STMin = cms.double(1000.)
process.ana.vlqMass = cms.double(1000.) #M=1000
process.ana.bosonMass = cms.double(91.2) #Z

process.TFileService = cms.Service("TFileService",
       fileName = cms.string(
         options.outFileName
         )
       )


#outCommand = ['keep *', 'drop *_evtcleaner_*_*', 'drop *_photons_*_*', 'drop *_photonjets_*_*', 'drop *_TriggerResults_*_*']          

outCommand = ['keep *', 'drop *_evtcleaner_*_*', 'drop *_photons_*_*', 'drop *_photonjets_*_*', 'drop *_*Puppi_*_*', 'drop *_TriggerResults_*_*']

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFileName.split('.',1)[0]+'_skim.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    dropMetaData = cms.untracked.string('DROPPED'),#'type_label_instance_process'                                                      
    outputCommands = cms.untracked.vstring(outCommand )
    )


## Event counters
from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)

process.load("Analysis.VLQAna.VLQCandProducer_cff")

process.p = cms.Path(
    process.allEvents
    *process.evtcleaner
    *process.cleanedEvents
    *cms.ignore(process.ana)
    * process.finalEvents
    )

#process.schedule = cms.Schedule(process.p)

#open('dump.py','w').write(process.dumpPython())
if options.doSkim:
  process.outpath = cms.EndPath(process.out)
