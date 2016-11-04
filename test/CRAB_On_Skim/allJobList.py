#! /bin/python

def list_Zelel(isData):
    if isData:  
        jobList = [
           # ['/DoubleEG/devdatta-Run2015D-05Oct2015-v1_B2GAnaFW_v74x_v8p4-ff3168b63d5cee365f49bf7ea3ba6ef3/USER', 'DoubleEG-Run2015D-05Oct2015-v1', '40', ''],
           # ['/DoubleEG/devdatta-Run2015D-PromptReco-v4_B2GAnaFW_v74x_v8p4-5daaa7fbf157b0642c1d3dfb260fbeff/USER', 'DoubleEG-Run2015D-PromptReco-v4', '40', ''],
            ['/DoubleEG/skhi-RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0-642b48c1c707cb5cfa93fc168fde448b/USER', 'DoubleEG-Run2016', '100', ''],
           
            ]
        return jobList
    else:
        jobList = [  
            ['/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_Run2Spring15MiniAODv2_25ns_v74x_v84-50153fb607659b6f9fb41d9f35391d0e/USER', 'DY_amcatnlo', '10',''],
            ['/TT_TuneCUETP8M1_13TeV-powheg-pythia8/jkarancs-B2GAnaFW_v74x_V8p4_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1-7b09a84a5c42d0a63b01d8e8e63c7a89/USER', 'TT_powheg', '20', ''],
            ]
        return jobList

def list_Zmumu(isData):
    if isData:    
        jobList = [
            ['/DoubleMuon/skhi-RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0-642b48c1c707cb5cfa93fc168fde448b/USER', 'DoubleMuon', '1', ''], 
            ]
        return jobList
    else:
        jobList = [ 
            #['/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_100-200', '20',''],          
           # ['/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_200-400', '20',''], 
           # ['/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0_1-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_400-600', '20',''],
           # ['/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_600-800', '15',''],
           # ['/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_800-1200', '10',''],
           # ['/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_1200-2500', '2',''],
           # ['/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'DY_2500-Inf', '2',''],
           # ['/TT_TuneCUETP8M1_13TeV-powheg-pythia8/asparker-RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0-9c09e10dd1f806cf9fdf5818b1c7d288/USER', 'TTbar', '40',''],
           
            #TPrimesamples
           # ['/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZtZ800', '2','EvtType_MC_tZtZ'],
           # ['/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZtH800', '2','EvtType_MC_tZtH'],
           # ['/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZbW800', '2','EvtType_MC_tZbW'],

           # ['/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZtZ1000', '2','EvtType_MC_tZtZ'],
           # ['/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZtH1000', '2','EvtType_MC_tZtH'],
           # ['/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZbW1000', '2','EvtType_MC_tZbW'],

           # ['/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZtZ1200', '2','EvtType_MC_tZtZ'],
           # ['/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZtH1200', '2','EvtType_MC_tZtH'],
           # ['/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p0-a5c23f9700079f8bc7b9d2b4fb46cf81/USER', 'tZbW1200', '2','EvtType_MC_tZbW'],
            ]
        return jobList
