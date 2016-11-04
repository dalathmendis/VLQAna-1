from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.General.requestName = 'DoubleMuon_os2lanav1'
config.General.workArea = 'os2lanav1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'os2lana_cfg.py'
config.JobType.pyCfgParams = ['isData=True', 'zdecaymode=zmumu', 'doPUReweightingOfficial=False', 'applyLeptonSFs=False','applyBTagSFs=False','applyDYNLOCorr=False','optimizeReco=False','doSkim=True']
#config.JobType.inputFiles = ['hnpv_data_Run2015D_mc_RunIISpring15DR74-Asympt25ns_pvwt.root','RunII2015_25ns_PUXsec69000nb.root','RunII2015_25ns_PUXsec72450nb.root', 'RunII2015_25ns_PUXsec65550nb.root', 'PUDistMC_2015_25ns_Startup_PoissonOOTPU.root','Summer15_25nsV7_MC_L1FastJet_AK4PFchs.txt', 'Summer15_25nsV7_MC_L2Relative_AK4PFchs.txt', 'Summer15_25nsV7_MC_L3Absolute_AK4PFchs.txt', 'Summer15_25nsV7_MC_Uncertainty_AK4PFchs.txt', 'Summer15_25nsV7_MC_L1FastJet_AK8PFchs.txt', 'Summer15_25nsV7_MC_L2Relative_AK8PFchs.txt', 'Summer15_25nsV7_MC_L3Absolute_AK8PFchs.txt', 'Summer15_25nsV7_MC_Uncertainty_AK8PFchs.txt']
#config.JobType.inputFiles = ['Summer15_25nsV6_MC_L2L3Residual_AK8PFchs.txt', 'Summer15_25nsV6_MC_L3Absolute_AK8PFchs.txt', 'hnpv_data_Run2015D_mc_RunIISpring15DR74-Asympt25ns_pvwt.root', 'PUDistData_Run2015ABCD.root', 'PUDistMC_2015_25ns_Startup_PoissonOOTPU.root']

config.JobType.inputFiles = ['2016_25ns_Spring_PUXsec65740nb50.root','PUDistMC_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU.root','scalefactors_v4.root','2016_25ns_Spring_PUXsec69200nb50.root', 'Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt','2016_25ns_Spring_PUXsec72660nb50.root','Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt','CSVv2_ichep.csv' , 'Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt','Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt']


config.section_("Data")
config.Data.inputDataset = '/DoubleMuon/skhi-RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p0-642b48c1c707cb5cfa93fc168fde448b/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = False
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/lpcbprime/noreplica/dmendis/CrabonSkim/Zmumu'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.section_('User')

