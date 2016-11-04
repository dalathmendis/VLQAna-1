from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.General.requestName = 'DY_amcatnlo_os2lana_v4'
config.General.workArea = 'os2lana_v4'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'os2lana_cfg.py'
config.JobType.pyCfgParams = ['isData=False', 'zdecaymode=zelel', 'doPUReweightingOfficial=True', 'applyLeptonSFs=True']
config.JobType.inputFiles = ['hnpv_data_Run2015D_mc_RunIISpring15DR74-Asympt25ns_pvwt.root','RunII2015_25ns_PUXsec69000nb.root','RunII2015_25ns_PUXsec72450nb.root', 'RunII2015_25ns_PUXsec65550nb.root', 'PUDistMC_2015_25ns_Startup_PoissonOOTPU.root','Summer15_25nsV7_MC_L1FastJet_AK4PFchs.txt', 'Summer15_25nsV7_MC_L2Relative_AK4PFchs.txt', 'Summer15_25nsV7_MC_L3Absolute_AK4PFchs.txt', 'Summer15_25nsV7_MC_Uncertainty_AK4PFchs.txt', 'Summer15_25nsV7_MC_L1FastJet_AK8PFchs.txt', 'Summer15_25nsV7_MC_L2Relative_AK8PFchs.txt', 'Summer15_25nsV7_MC_L3Absolute_AK8PFchs.txt', 'Summer15_25nsV7_MC_Uncertainty_AK8PFchs.txt']
#config.JobType.inputFiles = ['Summer15_25nsV6_MC_L2L3Residual_AK8PFchs.txt', 'Summer15_25nsV6_MC_L3Absolute_AK8PFchs.txt', 'hnpv_data_Run2015D_mc_RunIISpring15DR74-Asympt25ns_pvwt.root', 'PUDistData_Run2015ABCD.root', 'PUDistMC_2015_25ns_Startup_PoissonOOTPU.root']

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_Run2Spring15MiniAODv2_25ns_v74x_v84-50153fb607659b6f9fb41d9f35391d0e/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.ignoreLocality = False
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/lpcbprime/noreplica/dmendis/HistonewtodayInc1/Zelel'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.section_('User')

