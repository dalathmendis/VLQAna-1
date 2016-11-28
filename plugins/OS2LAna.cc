// -*- C++ -*-
//
// Package:    Analysis/VLQAna
// Class:      OS2LAna
// 
/**\class VLQAna OS2LAna.cc Analysis/VLQAna/plugins/OS2LAna.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Fri, 27 Feb 2015 16:09:10 GMT
// Modified: Sadia Khalil
//           25 Mar 2016 17:11 CDT
//

#include <iostream>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"

#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/DileptonCandsProducer.h"
#include "Analysis/VLQAna/interface/CandidateFilter.h"
#include "Analysis/VLQAna/interface/MuonMaker.h"
#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/ApplyLeptonSFs.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"
#include "Analysis/VLQAna/interface/METMaker.h"
#include "Analysis/VLQAna/interface/PickGenPart.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "Analysis/VLQAna/interface/MassReco.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"

#include "Analysis/VLQAna/interface/TopCandsProducer.h"
#include "Analysis/VLQAna/interface/ZCandsProducer.h"
#include "Analysis/VLQAna/interface/HCandsProducer.h"
#include "Analysis/VLQAna/interface/ZHCandsProducer.h"
#include "Analysis/VLQAna/interface/ApplyTriggerSFs.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

//
// class declaration
//

class OS2LAna : public edm::EDFilter {
  public:
    explicit OS2LAna(const edm::ParameterSet&);
    ~OS2LAna();

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt);
    double GetDYNLOCorr(const double dileppt);
  double htCorr(double ht, double p0, double p1); 
  double htCorr(double ht, double p0, double p1, double p2, double p3);
  
  // ----------member data ---------------------------
    edm::EDGetTokenT<string>   t_evttype         ;
    edm::EDGetTokenT<double>   t_evtwtGen        ;
    edm::EDGetTokenT<double>   t_evtwtPV         ;
    edm::EDGetTokenT<unsigned> t_npv             ;
    edm::EDGetTokenT<bool>     t_hltdecision     ;
    edm::ParameterSet DilepCandParams_           ; 
    edm::ParameterSet ZCandParams_               ; 
  //    edm::ParameterSet BoostedZCandParams_        ; 
    edm::ParameterSet GenHSelParams_             ;
    edm::ParameterSet genParams_                 ;
    const double HTMin_                          ;
    const double STMin_                          ; 
    const bool filterSignal_                     ;
    const bool additionalPlots_                  ;
    const std::string signalType_                ;
    const std::string zdecayMode_                ;
    const bool optimizeReco_                     ;
    const bool doSkim_                           ;
    const bool categorize_                           ;
    const double vlqMass_                        ;
    const double bosonMass_                      ;
    const bool applyLeptonSFs_                   ;
    const bool applyTriggerSFs_                   ;
    const bool applyBTagSFs_                     ;
    const bool applyHtCorr_                     ;
    const bool applyDYNLOCorr_                   ;
    const std::string fname_DYNLOCorr_           ; 
    const std::string funname_DYNLOCorr_         ; 
    ApplyLeptonSFs lepsfs                        ;
   ApplyTriggerSFs triggersfs                        ;

    METMaker metmaker                            ;
    MuonMaker muonmaker                          ; 
    ElectronMaker electronmaker                  ; 
    JetMaker jetAK4maker                         ; 
    JetMaker jetAK4BTaggedmaker                  ; 
    JetMaker jetAK8maker                         ; 
    JetMaker jetHTaggedmaker                     ; 
    JetMaker jetWTaggedmaker                     ; 
    JetMaker jetTopTaggedmaker                   ; 
    edm::Service<TFileService> fs                ; 
    std::map<std::string, TH1D*> h1_             ; 
    std::map<std::string, TH2D*> h2_             ; 
    std::string lep; 
    PickGenPart genpart                          ;
    const std::string fnamebtagSF_               ;
    std::unique_ptr<BTagSFUtils> btagsfutils_    ; 

};

using namespace std;

// static data member definitions
void OS2LAna::fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt){

  for  (unsigned int iele=0; iele<goodElectrons.size(); ++iele){
    float scEta = goodElectrons.at(iele).getscEta();
    if(fabs(scEta) <= 1.479){
      h1_["Eta_EB_el_pre"]-> Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EB_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaIn_EB_el_pre"]->Fill(goodElectrons.at(iele).getdEtaIn(), evtwt);
      h1_["dPhiIn_EB_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EB_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["D0_EB_el_pre"]->Fill(goodElectrons.at(iele).getD0(), evtwt);
      h1_["SCETA_EB_el_pre"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
      h1_["Full5x5siee_EB_el_pre"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
      h1_["HoE_EB_el_pre"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
      h1_["ooEmooP_EB_el_pre"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
      h1_["missHits_EB_el_pre"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
      h1_["conveto_EB_el_pre"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);
    }
    else if  (fabs(scEta) > 1.479 && fabs(scEta) < 2.5){
      h1_["Eta_EE_el_pre"]->Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EE_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdEtaIn(), evtwt);
      h1_["dPhiIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EE_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["D0_EE_el_pre"]->Fill(goodElectrons.at(iele).getD0(), evtwt);
      h1_["SCETA_EE_el_pre"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
      h1_["Full5x5siee_EE_el_pre"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
      h1_["HoE_EE_el_pre"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
      h1_["ooEmooP_EE_el_pre"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
      h1_["missHits_EE_el_pre"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
      h1_["conveto_EE_el_pre"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);
    }
  }
}

// constructors and destructor
OS2LAna::OS2LAna(const edm::ParameterSet& iConfig) :
  t_evttype               (consumes<string>  (iConfig.getParameter<edm::InputTag>("evttype"))),
  t_evtwtGen              (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtGen"))),
  t_evtwtPV               (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPV"))),
  t_npv                   (consumes<unsigned>(iConfig.getParameter<edm::InputTag>("npv"))),
  t_hltdecision           (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltdecision"))),
  DilepCandParams_        (iConfig.getParameter<edm::ParameterSet> ("DilepCandParams")),
  ZCandParams_            (iConfig.getParameter<edm::ParameterSet> ("ZCandParams")),
  //BoostedZCandParams_     (iConfig.getParameter<edm::ParameterSet> ("BoostedZCandParams")),
  GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
  genParams_              (iConfig.getParameter<edm::ParameterSet> ("genParams")),
  HTMin_                  (iConfig.getParameter<double>            ("HTMin")),
  STMin_                  (iConfig.getParameter<double>            ("STMin")), 
  filterSignal_           (iConfig.getParameter<bool>              ("filterSignal")), 
  additionalPlots_        (iConfig.getParameter<bool>              ("additionalPlots")), 
  signalType_             (iConfig.getParameter<std::string>       ("signalType")), 
  zdecayMode_             (iConfig.getParameter<std::string>       ("zdecayMode")),
  optimizeReco_           (iConfig.getParameter<bool>              ("optimizeReco")),
  doSkim_                 (iConfig.getParameter<bool>              ("doSkim")),
  categorize_             (iConfig.getParameter<bool>              ("categorize")),
  vlqMass_                (iConfig.getParameter<double>            ("vlqMass")),
  bosonMass_              (iConfig.getParameter<double>            ("bosonMass")),
  applyLeptonSFs_         (iConfig.getParameter<bool>              ("applyLeptonSFs")), 
  applyTriggerSFs_         (iConfig.getParameter<bool>              ("applyTriggerSFs")),
  applyBTagSFs_           (iConfig.getParameter<bool>              ("applyBTagSFs")), 
  applyHtCorr_            (iConfig.getParameter<bool>              ("applyHtCorr")),
  applyDYNLOCorr_         (iConfig.getParameter<bool>              ("applyDYNLOCorr")), 
  fname_DYNLOCorr_        (iConfig.getParameter<std::string>       ("File_DYNLOCorr")),
  funname_DYNLOCorr_      (iConfig.getParameter<std::string>       ("Fun_DYNLOCorr")),
  lepsfs                  (iConfig.getParameter<edm::ParameterSet> ("lepsfsParams")),
  metmaker                (iConfig.getParameter<edm::ParameterSet> ("metselParams"),consumesCollector()),
  muonmaker               (iConfig.getParameter<edm::ParameterSet> ("muselParams"),consumesCollector()),
  electronmaker           (iConfig.getParameter<edm::ParameterSet> ("elselParams"),consumesCollector()),
  jetAK4maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK4selParams"),consumesCollector()),
  jetAK4BTaggedmaker      (iConfig.getParameter<edm::ParameterSet> ("jetAK4BTaggedselParams"),consumesCollector()),
  jetAK8maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK8selParams"),consumesCollector()),
  jetHTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetHTaggedselParams"),consumesCollector()),
  jetWTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetWTaggedselParams"),consumesCollector()),
  jetTopTaggedmaker       (iConfig.getParameter<edm::ParameterSet> ("jetTopTaggedselParams"),consumesCollector()),   
  lep                     (iConfig.getParameter<std::string>       ("lep")), 
  genpart                 (genParams_, consumesCollector()),
  fnamebtagSF_            (iConfig.getParameter<std::string>       ("fnamebtagSF")),
  btagsfutils_            (new BTagSFUtils(fnamebtagSF_,BTagEntry::OP_MEDIUM,30., 670., 30., 670., 20., 1000.))
  {

   produces<vlq::JetCollection>("tjets") ; 
   produces<vlq::JetCollection>("wjets") ; 
   produces<vlq::JetCollection>("bjets") ; 
   produces<vlq::JetCollection>("jets") ; 
   produces<vlq::CandidateCollection>("zllcands") ; 
  }



OS2LAna::~OS2LAna() {}

bool OS2LAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<string>   h_evttype     ; evt.getByToken(t_evttype    , h_evttype    ) ; 
  Handle<double>   h_evtwtGen    ; evt.getByToken(t_evtwtGen   , h_evtwtGen   ) ; 
  Handle<double>   h_evtwtPV     ; evt.getByToken(t_evtwtPV    , h_evtwtPV    ) ; 
  Handle<unsigned> h_npv         ; evt.getByToken(t_npv        , h_npv        ) ; 
  Handle<bool>     h_hltdecision ; evt.getByToken(t_hltdecision, h_hltdecision) ; 

  //  unsigned npv(*h_npv.product()) ; 
  h1_["checkPU"]->Fill(*h_npv.product(), *h_evtwtGen.product());
  if(zdecayMode_ == "zmumu") {lep = "mu";}
  else if ( zdecayMode_ == "zelel") {lep = "el";}
  else edm::LogError("OS2LAna::filter") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;

  if (filterSignal_) {
    if(doSkim_){
      if( signalType_.empty()){
        //cout << *h_evttype.product() << endl;                                                                                              
        if      (*h_evttype.product() == "EvtType_MC_bZbZ"){h1_["signalEvts_all"] -> Fill(1);}
        else if (*h_evttype.product() == "EvtType_MC_bZbH"){h1_["signalEvts_all"] -> Fill(2);}
        else if (*h_evttype.product() == "EvtType_MC_bZtW"){h1_["signalEvts_all"] -> Fill(3);}
        else if (*h_evttype.product() == "EvtType_MC_bHbH"){h1_["signalEvts_all"] -> Fill(4);}
        else if (*h_evttype.product() == "EvtType_MC_bHtW"){h1_["signalEvts_all"] -> Fill(5);}
        else if (*h_evttype.product() == "EvtType_MC_tWtW"){h1_["signalEvts_all"] -> Fill(6);}
        else if (*h_evttype.product() == "EvtType_MC_tZtZ"){h1_["signalEvts_all"] -> Fill(7);}
        else if (*h_evttype.product() == "EvtType_MC_tZtH"){h1_["signalEvts_all"] -> Fill(8);}
        else if (*h_evttype.product() == "EvtType_MC_tZbW"){h1_["signalEvts_all"] -> Fill(9);}
        else if (*h_evttype.product() == "EvtType_MC_tHtH"){h1_["signalEvts_all"] -> Fill(10);}
        else if (*h_evttype.product() == "EvtType_MC_tHbW"){h1_["signalEvts_all"] -> Fill(11);}
        else if (*h_evttype.product() == "EvtType_MC_bWbW"){h1_["signalEvts_all"] -> Fill(12);}
      }
      else{edm::LogError(">>>>ERROR>>>>OS2LAna >>>>  Please do not specify any signal type when skimming the signal" );}
    }
    else{
      if (*h_evttype.product()!=signalType_) return false ;
      else  h1_["signalEvts"] -> Fill(1) ;
    }
  }

  const bool hltdecision(*h_hltdecision.product()) ; 
  if ( !hltdecision ) return false;

  //double evtwtgen(*h_evtwtGen.product());
  double evtwt((*h_evtwtGen.product()) * (*h_evtwtPV.product())) ; 
  double prewt=1.0;
  // cout << " evtwt , prewt = "<< evtwt << " , " <<prewt <<endl;
  vlq::MuonCollection goodMuons; 
  muonmaker(evt, goodMuons) ; 

  vlq::ElectronCollection goodElectrons; 
  electronmaker(evt, goodElectrons) ;

  vlq::MetCollection goodMet;
  metmaker(evt, goodMet) ;

  vlq::CandidateCollection dimuons, dielectrons, dileptons,L1,L2, el1,el2;   
  vlq::CandidateCollection zll,zllBoosted; //generic collection
  vlq::CandidateCollection tops, W, B, BC, D, Z,ZB, H,Hb,ZH, ZHb; 

  // dilepton properties: M > 50, lead pt > 45, second pt > 25
  DileptonCandsProducer dileptonsprod(DilepCandParams_) ; 
  dileptonsprod.operator()<vlq::MuonCollection>(dimuons, goodMuons); 
  dileptonsprod.operator()<vlq::ElectronCollection>(dielectrons, goodElectrons) ; 
  dileptonsprod.operator()<vlq::ElectronCollection>(goodElectrons,L1,L2) ;


  //================================================================
  //First pre-selection: 1) 2 OS dileptons from boosted Z, >=3 jets
  //================================================================

  //dilepton candidate
  if (zdecayMode_ == "zmumu") {dileptons = dimuons; }
  else if (zdecayMode_ == "zelel") {dileptons = dielectrons;}
  if (dileptons.size()  < 1) return false;

  //// Get Dy EWK correction
  if ( applyDYNLOCorr_ &&  *h_evttype.product() != "EvtType_Data" ) {
    double EWKNLOkfact(GetDYNLOCorr(dileptons.at(0).getPt())) ; 
    evtwt *= EWKNLOkfact ;
  }

  //get lepton ID and Iso SF
  if (applyLeptonSFs_ && *h_evttype.product() != "EvtType_Data") {
    if ( zdecayMode_ == "zmumu" ){
      evtwt *= lepsfs(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) * lepsfs(goodMuons.at(1).getPt(), goodMuons.at(1).getEta() ) ;}
    else if ( zdecayMode_ == "zelel" ){ 
      evtwt *= lepsfs(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) * lepsfs(goodElectrons.at(1).getPt(), goodElectrons.at(1).getEta() ) ;}
  }

  //apply flat trigger scale factors
  if (applyTriggerSFs_ && *h_evttype.product() != "EvtType_Data") {                                                     
    if ( zdecayMode_ == "zmumu" ){                                                                                        
      evtwt *= triggersfs.TrigSFMu1(goodMuons.at(0).getPt(),goodMuons.at(0).getEta())*triggersfs.TrigSFMu2(goodMuons.at(1).getPt(),goodMuons.at(1).getEta());}                                                                                                 
    else if ( zdecayMode_ == "zelel" ){                                                                                   
      evtwt *=0.968256;}                                                                                                  
  }


  h1_["cutflow"] -> Fill(1, evtwt) ;

  //Z mass candidate filter: 75 < M < 105, lead pt > 45, 2nd pt > 25, Z pt > 0
  CandidateFilter zllfilter(ZCandParams_) ; 
  zllfilter(dileptons, zll);
  zllfilter(L1,L2,el1,el2);
 
  //electron energies before preslection
  if (el1.size()>0) h1_["le1E"] -> Fill(el1.at(0).getEnergy(), evtwt) ;
  if (el2.size()>0) h1_["le2E"] -> Fill(el2.at(0).getEnergy(), evtwt) ;
  if (zll.size()>0) h1_["zllE"] -> Fill(zll.at(0).getEnergy(), evtwt) ;
  if (goodElectrons.size()>0) h1_["el1E"] -> Fill(goodElectrons.at(0).getEnergy(), evtwt) ;  
  if (goodElectrons.size()>0) h1_["el2E"] -> Fill(goodElectrons.at(1).getEnergy(), evtwt) ;
  

  // jets
  vlq::JetCollection goodAK4Jets;
  jetAK4maker(evt, goodAK4Jets) ;

  CandidateCleaner cleanjets(0.4);
  if (zdecayMode_ == "zmumu") {cleanjets(goodAK4Jets, goodMuons);}
  else if (zdecayMode_ == "zelel") {cleanjets(goodAK4Jets, goodElectrons);} 

  HT htak4(goodAK4Jets) ; 
  double ht = htak4.getHT();

  if (applyHtCorr_ && *h_evttype.product() != "EvtType_Data"){
    double corr(1.);
    if (zdecayMode_ == "zmumu"){
      if (ht <1300){
	//corr = htCorr(ht, 1.17876, -0.000576353);//for norebin
	corr = htCorr(ht, 1.09227, -0.000473877);//for rebin =4
	//corr = htCorr(ht, 1.06501 , 0.000273571, -0.00000149895,0.000000000688488);
      }
      else if (ht >= 1300){
	//corr = 0.602407;//for no rebining
	corr = 0.47623;//for rebin = 4
      }
    }
    
    else if (zdecayMode_ == "zelel"){
      if (ht <1200){
        corr = htCorr(ht, 1.88713, -0.000739714);//for rebin =4                                                                                                                        
      }
      else if (ht >= 1200){
        corr = 0.9994732;//for rebin = 4                                                                                                                                               
      }




    }
    if (corr > 0.) evtwt *= corr;
  }
  
  ////////////////////////////////////////////////////////// 
  //Fill N-1 selected plots for dilepton mass, Ht, ad Njets
  //////////////////////////////////////////////////////////
  //at least one Z cand in event
  if(zll.size() > 0) {h1_["cutflow"] -> Fill(2, evtwt) ;}
  else return false ;

  // at least HT > 200 in event                                                                                                     
  if ( htak4.getHT() > HTMin_ ) h1_["cutflow"] -> Fill(3, evtwt) ;                                                                
  else return false ;   


  //at least 3 AK4 jets in event                                                                                                        
  if (goodAK4Jets.size() > 2 ) {h1_["cutflow"] -> Fill(4, evtwt) ;}
  else return false;

  if (doSkim_){
    return true;
  }
  //deltaPhi(MET,lep)
  double ST = htak4.getHT() ;
  ST += zll.at(0).getPt() + goodMet.at(0).getFullPt();

  ///////////////////////////////////////////// 
  // fill rest of the plots after pre-selection
  /////////////////////////////////////////////
 
  h1_["ht_pre"] -> Fill(htak4.getHT(), evtwt);

  for (auto izll : zll){
    h1_["mass_z"+lep+lep+"_pre"] -> Fill(izll.getMass(), evtwt) ;
    h1_["mass_Z"+lep+lep+"_pre"] -> Fill(izll.getMass(), evtwt) ;
  }

  //lepton energies after preselection
  if (el1.size()>0) h1_["le1E_pre"] -> Fill(el1.at(0).getEnergy(), evtwt) ;
  if (el2.size()>0) h1_["le2E_pre"] -> Fill(el2.at(0).getEnergy(), evtwt) ;
  if (zll.size()>0) h1_["zllE_pre"] -> Fill(zll.at(0).getEnergy(), evtwt) ;
  if (goodElectrons.size()>0) h1_["el1E_pre"] -> Fill(goodElectrons.at(0).getEnergy(), evtwt) ;
  if (goodElectrons.size()>0) h1_["el2E_pre"] -> Fill(goodElectrons.at(1).getEnergy(), evtwt) ;


  h1_["nak4_pre"] -> Fill(goodAK4Jets.size(), evtwt) ;
  h1_["st_pre"] -> Fill(ST, evtwt) ;


  //ak4 jet plots
  for(int j=0; j<2; ++j){
    h1_[Form("ptak4jet%d_pre", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
    h1_[Form("etaak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
    h1_[Form("cvsak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
  }
  h1_["phi_jet1MET_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

  // npv
  h1_["npv_noweight_pre"] -> Fill(*h_npv.product(), *h_evtwtGen.product()); 
  h1_["npv_pre"] -> Fill(*h_npv.product(), evtwt);

  // met
  h1_["met_pre"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
  h1_["metPhi_pre"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

  //lepton specfic properties
  if ( zdecayMode_ == "zmumu" ){       
    for(int l=0; l<2; ++l){
      h1_["pt_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
      h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
    } 

    for  (unsigned int imu=0; imu<goodMuons.size(); ++imu){
      h1_["Iso04_mu_pre"]->Fill(goodMuons.at(imu).getIso04(), evtwt);
      h1_["Dz_mu_pre"]->Fill(goodMuons.at(imu).getDz(), evtwt);
      h1_["D0_mu_pre"]->Fill(goodMuons.at(imu).getD0(), evtwt);
    }

    h1_["dr_mumu_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
  }
  else if (zdecayMode_ == "zelel" ) {
    for(int l=0; l<2; ++l){
      h1_["pt_"+lep+Form("%d_pre", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
      h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
    } 

    h1_["dr_elel_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
    if(additionalPlots_) fillAdditionalPlots(goodElectrons,evtwt);
  }
   
  // Z pt
  for (auto izll : zll) h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), evtwt) ;
  //========================================================
  // Preselection done, proceeding with control selections
  //========================================================
  
  //Flavor Check                                                                                                            
  for (vlq::Jet jet : goodAK4Jets) {
    if ( abs(jet.getHadronFlavour()) == 5){
      for (auto izll : zll) h1_["pt_zb_pre"] -> Fill(jet.getPt()+izll.getPt(), evtwt) ;
    }
    if ( abs(jet.getHadronFlavour()) == 4){
      for (auto izll : zll) h1_["pt_zc_pre"] -> Fill(jet.getPt()+izll.getPt(), evtwt) ;
    }
    if ( abs(jet.getHadronFlavour()) < 4){
      for (auto izll : zll) h1_["pt_zlight_pre"] -> Fill(jet.getPt()+izll.getPt(), evtwt) ;
    }
  }

  
  //b-tagging:
  vlq::JetCollection goodBTaggedAK4Jets;
  jetAK4BTaggedmaker(evt, goodBTaggedAK4Jets) ; 


  if (zdecayMode_ == "zmumu") {cleanjets(goodBTaggedAK4Jets, goodMuons); }
  else if (zdecayMode_ == "zelel") {cleanjets(goodBTaggedAK4Jets, goodElectrons); }  


  double btagsf(1) ;
  double btagsf_bcUp(1) ; 
  double btagsf_bcDown(1) ; 
  double btagsf_lUp(1) ; 
  double btagsf_lDown(1) ; 
  if ( applyBTagSFs_  && *h_evttype.product() != "EvtType_Data") {
     std::vector<double>csvs;
     std::vector<double>pts;
     std::vector<double>etas;
     std::vector<int>   flhads;

     for (vlq::Jet jet : goodAK4Jets) {
       csvs.push_back(jet.getCSV()) ; 
       pts.push_back(jet.getPt()) ; 
       etas.push_back(jet.getEta()) ; 
       flhads.push_back(jet.getHadronFlavour()) ; 
     }

     btagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetAK4maker.idxjetCSVDiscMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ; 


  }
  MassReco reco;
  TLorentzVector Leptons = zll.at(0).getP4();
  // Mass reconstruction                                                                                                                    
  pair<double, double> chi2_result_cnt;

  //fill the control plots related to mass reconstruction                                                                                   
  if ( goodBTaggedAK4Jets.size() == 0) {
    for (auto izll : zll) {
      h1_["nob_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
    }
    h1_["nob_ht"] ->Fill(htak4.getHT(), evtwt);
    h1_["nob_st"] ->Fill(ST, evtwt);
    if (goodMet.at(0).getFullPt()<60){ 
      h1_["nbjets_met_0btagcnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      h1_["ht_met_0btagcnt"] ->Fill(htak4.getHT(), evtwt);
      h1_["st_met_0btagcnt"] ->Fill(ST, evtwt);
    }

    h1_["btagSF_0btag"] ->Fill(btagsf, evtwt);
  }
  
  else if (goodBTaggedAK4Jets.size() > 0){ 
    for (auto izll : zll) {
      h1_["b_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
      h1_["b_st"] ->Fill(ST, evtwt);
    }
    h1_["1b_ht"] ->Fill(htak4.getHT(), evtwt);
    h1_["1b_st"] ->Fill(ST, evtwt);
    if (goodMet.at(0).getFullPt()<60){
      h1_["nbjets_met_1btagcnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      h1_["ht_met_1btagcnt"] ->Fill(htak4.getHT(), evtwt);
      h1_["st_met_1btagcnt"] ->Fill(ST, evtwt);
    
    }
    h1_["btagSF_1btag"] ->Fill(btagsf, evtwt);
  }
  
  if (goodMet.at(0).getFullPt()<60){ 
    h1_["nbjets_met_cnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
    h1_["lowmet_ht"] ->Fill(htak4.getHT(), evtwt);
    h1_["lowmet_st"] ->Fill(ST, evtwt);
  }
  
  h1_["nbjets_cnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
  
  //fill control plots
  if ( goodBTaggedAK4Jets.size() > 0 && ST < 700) {
    for (auto izll : zll) {
      h1_["mass_z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;  
      h1_["pt_z"+lep+lep+"_cnt"] -> Fill(izll.getPt(), evtwt) ; 
    }
    h1_["nak4_cnt"] -> Fill(goodAK4Jets.size(), evtwt) ;
    h1_["ht_cnt"] -> Fill(htak4.getHT(), evtwt) ;
    h1_["st_cnt"] -> Fill(ST, evtwt) ;   
    h1_["npv_noweight_cnt"] -> Fill(*h_npv.product(), *h_evtwtGen.product()); 
    h1_["npv_cnt"] -> Fill(*h_npv.product(), evtwt);
    h1_["met_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
    h1_["metPhi_cnt"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

    //lepton specfic properties
    if ( zdecayMode_ == "zmumu" ){       
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
        h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
      } 
      h1_["dr_mumu_cnt"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
    }
    else if (zdecayMode_ == "zelel" ) {
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_cnt", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
        h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
      } 
      h1_["dr_elel_cnt"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
    }

    //ak4 jet plots
    for(int j=0; j<2; ++j){
      h1_[Form("ptak4jet%d_cnt", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
      h1_[Form("etaak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
      h1_[Form("cvsak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
    }
    h1_["phi_jet1MET_cnt"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
  }
  if (categorize_){
    if ( goodBTaggedAK4Jets.size() >=1  && ST < 1000) {
      // HCandidate cut(boosted)                                                                             
      for (unsigned i=0; i<goodAK4Jets.size(); i++) { 
	if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450 && goodAK4Jets.at(i).getCSV()>0.800){
	  TLorentzVector h;
	  h= goodAK4Jets.at(i).getP4();
	  vlq::Candidate h1(h);
	  Hb.push_back(h1);
	}
      }
      //HCandidates (nonboosted)
      HCandsProducer h;
      h.operator()(goodAK4Jets.size(), 2, goodAK4Jets,H)  ;
      
      //Zcandidate cut (boosted)                                     
      for (unsigned i=0; i<goodAK4Jets.size(); i++) {                                                  
	if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
	  TLorentzVector zb;    
	  zb= goodAK4Jets.at(i).getP4();
	  vlq::Candidate zb1(zb);       
	  ZB.push_back(zb1);
	  
	}                                                                                                                                       
      }              

      //Z candidate cut (non boosted)                                                                                                            
      ZCandsProducer z;
      z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;
      
      //Top Candidates (Category D)    
      for (unsigned i=0; i<goodAK4Jets.size(); i++) {
	if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
	  TLorentzVector d;
	  d= goodAK4Jets.at(i).getP4();
	  vlq::Candidate d1(d);
	  D.push_back(d1);
	}
      }
      
      TopCandsProducer top,w;
      // Category BC
      w.operator()(goodAK4Jets, W,B) ;
      for (unsigned i=0; i<W.size(); i++) {
	for (unsigned j=0; j<B.size(); j++) {
	  TLorentzVector bc1;
	  bc1= W.at(i).getP4()+B.at(j).getP4();
	  if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
	    vlq::Candidate bc2(bc1);
	    BC.push_back(bc2); 
	  }  
	}
      }
      // category A
      top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;
  
      // HCandidate cut(boosted)                                                                                                                 
      
      for (unsigned i=0; i<Hb.size(); i++) {
	h1_["H_mass_b_cnt"] -> Fill(Hb.at(i).getMass(), evtwt) ;
	h1_["H_Pt_b_cnt"] -> Fill(Hb.at(i).getPt(), evtwt) ;
      }
      h1_["nHcandidatejets_b_cnt"] -> Fill(Hb.size(), evtwt) ;
      
      //Hcandidate cut(nonboosted)
      for (unsigned i=0; i<H.size(); i++) {
	h1_["H_mass_nb_cnt"] -> Fill(H.at(i).getMass(), evtwt) ;
	h1_["H_Pt_nb_cnt"] -> Fill(H.at(i).getPt(), evtwt) ;
      }
      
      h1_["nHcandidatejets_nb_cnt"] -> Fill(H.size(), evtwt) ;
      
      
      double  nHcandidates=0.0;
      double  nHcandidates1=0.0;
      if(Hb.size()>0 ||H.size()>0){
	nHcandidates = Hb.size()+H.size();
	
	h1_["nHcandidatejets_cnt"] -> Fill(nHcandidates, evtwt) ;
	
      }
      
      nHcandidates1 = Hb.size()+H.size();
      h1_["nHcandidatejets1_cnt"] -> Fill(nHcandidates1, evtwt) ;
      
      
      //Zcandidate cut (boosted)
      for (unsigned i=0; i<ZB.size(); i++) {     
	h1_["Z_mass_a_cnt"] -> Fill(ZB.at(i).getMass(), evtwt) ;
	h1_["Z_Pt_a_cnt"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
      }                                                                                                                                                                                       
      h1_["nzcandidatejets_a_cnt"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
      
      //Z candidate cut (non boosted)                                                                                                                                                         
      for (unsigned i=0; i<Z.size(); i++) {
	h1_["Z_mass_b_cnt"] -> Fill(Z.at(i).getMass(), evtwt) ;
	h1_["Z_Pt_b_cnt"] -> Fill(Z.at(i).getPt(), evtwt) ;
      }
      h1_["nzcandidatejets_b_cnt"] -> Fill(Z.size(), evtwt) ;
      
      double nzcandidates=0.0;
      double nzcandidates1=0.0;
      if ( ZB.size() || Z.size()){
	h1_["cutflow"] -> Fill(12, evtwt);
	nzcandidates = ZB.size()+ Z.size();
	h1_["nzcandidatejets_tot_cnt"] -> Fill(nzcandidates, evtwt) ;
	
      }
      nzcandidates1 = ZB.size()+ Z.size();
      h1_["nzcandidatejets1_tot_cnt"] -> Fill(nzcandidates1, evtwt) ;
      
      
      // Category D    
      for (unsigned i=0; i<D.size(); i++) {
	h1_["top_mass_d_cnt"] -> Fill(D.at(i).getMass(), evtwt) ;
	h1_["top_Pt_d_cnt"] -> Fill(D.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_d_cnt"] -> Fill(D.size(), evtwt) ;
      
      // Category BC
      for (unsigned i=0; i<W.size(); i++) {
	h1_["W_mass_bc_cnt"] -> Fill(W.at(i).getMass(), evtwt) ;
      }
      h1_["nWcandidatejets_bc_cnt"] -> Fill(W.size(), evtwt) ;
      
      for (unsigned i=0; i<B.size(); i++) {
	h1_["lightjet_mass_bc_cnt"] -> Fill(B.at(i).getMass(), evtwt) ; 
      }
      h1_["nlightjetcandidatejets_bc_cnt"] -> Fill(B.size(), evtwt) ;
      
      
      
      for (unsigned i=0; i<BC.size(); i++) {
	h1_["top_mass_bc_cnt"] -> Fill(BC.at(i).getMass(), evtwt) ;
	h1_["top_Pt_bc_cnt"] -> Fill(BC.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_bc_cnt"] -> Fill(BC.size(), evtwt) ;
     
      for (unsigned i=0; i<tops.size(); i++) {
	h1_["top_mass_a_cnt"] -> Fill(tops.at(i).getMass(), evtwt) ;
	h1_["top_Pt_a_cnt"] -> Fill(tops.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_a_cnt"] -> Fill(tops.size(), evtwt) ;
      
      double  ntopcandidates=0.0;
      double  ntopcandidates1=0.0;
      if(D.size() || BC.size() ||tops.size()){
	ntopcandidates = D.size()+BC.size()+tops.size();
	
	h1_["ntopcandidatejets_cnt"] -> Fill(ntopcandidates, evtwt) ;  
	
      } 
      ntopcandidates1 = D.size()+BC.size()+tops.size();
      h1_["ntopcandidatejets1_cnt"] -> Fill(ntopcandidates1, evtwt) ;
      
      //Z and top corelations and ST tempelates
      
      //h1_["st_cnt"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(1, evtwt) ;
      if (goodBTaggedAK4Jets.size() == 1){
	h1_["cutflow1"] -> Fill(2, evtwt) ;   
	//h1_["st_cnt1b"] -> Fill(ST,evtwt);
      }
      
      if (goodBTaggedAK4Jets.size() >=2){
	//h1_["st_cnt2b"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(3, evtwt) ;
      }
      
      //n,Z,H,B
      //(1)
      if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
	//h1_["st_cntT1Z1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(4, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(8, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(1, evtwt) ;
	    // h1_["st_cntT1Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(2, evtwt) ;
	    //h1_["st_cntT1Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	  
	}
	else if (nHcandidates == 0.0){
	  //h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(9, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(3, evtwt) ;
	    // h1_["st_cntT1Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(4, evtwt) ;
	    //h1_["st_cntT1Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}  
      }
      //(2)                                                                                                                                                                                        
      if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
	//h1_["st_cntT0Z1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(5, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(10, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(5, evtwt) ;
	    // h1_["st_cntT1Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(6, evtwt) ;
	    // h1_["st_cntT1Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  // h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(11, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(7, evtwt) ;
	    // h1_["st_cntT1Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(8, evtwt) ;
	    // h1_["st_cntT1Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
      //(3)                                                                                                                                                                                        
      if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
	//h1_["st_cntT1Z0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(6, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(12, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(9, evtwt) ;
	    // h1_["st_cntT1Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(10, evtwt) ;
	    //h1_["st_cntT1Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  // h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(13, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(11, evtwt) ;
	    // h1_["st_cntT1Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(12, evtwt) ;
	    // h1_["st_cntT1Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
      //(4)                                                                                                                                                                                        
      if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
	//h1_["st_cntT0Z0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(7, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(14, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(13, evtwt) ;
	    // h1_["st_cntT1Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(14, evtwt) ;
	    // h1_["st_cntT1Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  //h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(15, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(15, evtwt) ;
	    // h1_["st_cntT1Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(16, evtwt) ;
	    // h1_["st_cntT1Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
      
    }
  }
  
 Hb.clear();
 H.clear();
 ZB.clear();
 Z.clear();
 D.clear();
 BC.clear();
 tops.clear();
 h1_["ptak4led"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
 h1_["ptak4sub"] -> Fill(goodAK4Jets.at(1).getEta(), evtwt) ;
 h1_["httest"] -> Fill(htak4.getHT(), evtwt) ;


  //===========================================================
  // Control selection done, proceeding with signal selections
  //===========================================================

  //Boosted Z candidates: Z pt > 150 GeV
  // CandidateFilter boostedzllfilter(BoostedZCandParams_) ; 
  //boostedzllfilter(dileptons, zllBoosted) ;    
  //if(zllBoosted.size() > 0) h1_["cutflow"] -> Fill(5, evtwt) ;
  //else return false ; 

  // leading jet pt > 100 GeV
  if (goodAK4Jets.at(0).getPt() > 100){h1_["cutflow"] -> Fill(5, prewt) ; }
  else return false;

  // 2nd laeding jet pt > 50 GeV
  if (goodAK4Jets.at(1).getPt() > 50){h1_["cutflow"] -> Fill(6, prewt) ; }
  else return false;

  // at least one b-jet 
  if ( goodBTaggedAK4Jets.size() > 0 ) { h1_["cutflow"] -> Fill(7, evtwt) ;}
  else return false;

  // ST > 1000 GeV
  if ( ST > STMin_ ) h1_["cutflow"] -> Fill(8, evtwt) ;  
  else return false ; 

  // get AK8, top, b, Z, and W jets
  vlq::JetCollection goodAK8Jets, goodHTaggedJets, goodWTaggedJets, goodTopTaggedJets;
  jetAK8maker(evt, goodAK8Jets); 
  cleanjets(goodAK8Jets, goodMuons); 
  cleanjets(goodAK8Jets, goodElectrons); 

  jetWTaggedmaker(evt, goodWTaggedJets);
  cleanjets(goodWTaggedJets, goodMuons); 
  cleanjets(goodWTaggedJets, goodElectrons); 

  jetHTaggedmaker(evt, goodHTaggedJets);
  cleanjets(goodHTaggedJets, goodMuons); 
  cleanjets(goodHTaggedJets, goodElectrons); 

  jetTopTaggedmaker(evt, goodTopTaggedJets);
  cleanjets(goodTopTaggedJets, goodMuons); 
  cleanjets(goodTopTaggedJets, goodElectrons); 

  //if (matchedAK8Jets.size() > 0) h1_["cutflow"] -> Fill(9, evtwt) ;
  //else return false 

  ///////////////////////////////////////
  // fill all the plots in signal region
  ///////////////////////////////////////
  for (auto izll : zll) {
    h1_["mass_z"+lep+lep] -> Fill(izll.getMass(), evtwt) ;  
    h1_["pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ; 
  }

  h1_["nak4"] -> Fill(goodAK4Jets.size(), evtwt) ;
  h1_["ht"] -> Fill(htak4.getHT(), evtwt) ;
  h1_["st"] -> Fill(ST, evtwt) ;   
  h1_["npv_noweight"] -> Fill(*h_npv.product(), *h_evtwtGen.product()); 
  h1_["npv"] -> Fill(*h_npv.product(), evtwt);
  h1_["met"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
  h1_["metPhi"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

  //lepton specfic properties
  if ( zdecayMode_ == "zmumu" ){       
    for(int l=0; l<2; ++l){
      h1_["pt_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
      h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
    } 
    h1_["dr_mumu"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
  }
  else if (zdecayMode_ == "zelel" ) {
    for(int l=0; l<2; ++l){
      h1_["pt_"+lep+Form("%d", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
      h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
    } 
    h1_["dr_elel"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
  }
  //ak4 jet plots
  for(int j=0; j<3; ++j){
    h1_[Form("ptak4jet%d", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
    h1_[Form("etaak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
    h1_[Form("cvsak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
  }
  h1_["phi_jet1MET"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

  // fill the b-tagging plots and efficiency maps
  h1_["nbjets"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
  if (goodMet.at(0).getFullPt()<60){ h1_["nbjets_met_sig"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;}
  h1_["ptbjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
  h1_["etabjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;

  for (vlq::Jet jet : goodAK4Jets) {
    if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
    else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
    else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
  }

  for (vlq::Jet jet : goodBTaggedAK4Jets) {
    if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
    else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
    else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
  }

  // fill the additional plots
  h1_["nak8"] -> Fill(goodAK8Jets.size(), evtwt) ;
  h1_["nwjet"] -> Fill(goodWTaggedJets.size(), evtwt) ; 
  h1_["nhjet"] -> Fill(goodHTaggedJets.size(), evtwt) ; 
  h1_["ntjet"] -> Fill(goodTopTaggedJets.size(), evtwt) ; 

  if (goodAK8Jets.size() > 0) {
    h1_["ptak8leading"] -> Fill((goodAK8Jets.at(0)).getPt(), evtwt) ; 
    h1_["etaak8leading"] -> Fill((goodAK8Jets.at(0)).getEta(), evtwt) ;
    h1_["mak8leading"] -> Fill((goodAK8Jets.at(0)).getMass(), evtwt) ; 
    h1_["trimmedmak8leading"] -> Fill((goodAK8Jets.at(0)).getTrimmedMass(), evtwt) ;
    h1_["prunedmak8leading"] -> Fill((goodAK8Jets.at(0)).getPrunedMass(), evtwt) ;
    h1_["softdropmak8leading"] -> Fill((goodAK8Jets.at(0)).getSoftDropMass(), evtwt) ;
  }
  if (goodAK8Jets.size() > 1) {
    h1_["ptak82nd"] -> Fill((goodAK8Jets.at(1)).getPt(), evtwt) ; 
    h1_["etaak82nd"] -> Fill((goodAK8Jets.at(1)).getEta(), evtwt) ;
    h1_["mak82nd"] -> Fill((goodAK8Jets.at(1)).getMass(), evtwt) ; 
    h1_["trimmedmak82nd"] -> Fill((goodAK8Jets.at(1)).getTrimmedMass(), evtwt) ;
    h1_["prunedmak82nd"] -> Fill((goodAK8Jets.at(1)).getPrunedMass(), evtwt) ;
    h1_["softdropmak82nd"] -> Fill((goodAK8Jets.at(1)).getSoftDropMass(), evtwt) ;
  }

  //// Make B->bZ and T->tZ->bWZ candidates
  TLorentzVector tp_p4, bp_p4;
  tp_p4.SetPtEtaPhiM(0,0,0,0);
  bp_p4.SetPtEtaPhiM(0,0,0,0);

  if (goodTopTaggedJets.size() > 0 && zll.size()>0) {
    tp_p4 = zll.at(0).getP4() + goodTopTaggedJets.at(0).getP4() ;
  }
  else if ( goodWTaggedJets.size() > 0 && goodBTaggedAK4Jets.size() > 0  && zll.size()>0 ) { 
    tp_p4 = zll.at(0).getP4() + goodBTaggedAK4Jets.at(0).getP4() + goodWTaggedJets.at(0).getP4() ;
  }
  else if ( goodBTaggedAK4Jets.size() > 0 && zll.size()>0 && goodAK8Jets.size()>0) {    
    tp_p4 = zll.at(0).getP4() + goodBTaggedAK4Jets.at(0).getP4() + goodAK8Jets.at(0).getP4() ; 
  }
  else if ( goodBTaggedAK4Jets.size() > 0 && zll.size()>0 ){
    bp_p4 = zll.at(0).getP4() + goodBTaggedAK4Jets.at(0).getP4() ; 
  }

  h1_["ptTprime"]->Fill(tp_p4.Pt(), evtwt) ; 
  h1_["yTprime"] ->Fill(tp_p4.Rapidity(), evtwt) ; 
  h1_["mTprime"] ->Fill(tp_p4.Mag(), evtwt) ; 

  h1_["ptBprime"]->Fill(bp_p4.Pt(), evtwt) ; 
  h1_["yBprime"] ->Fill(bp_p4.Rapidity(), evtwt) ; 
  h1_["mBprime"] ->Fill(bp_p4.Mag(), evtwt) ;

  if (categorize_){

    // HCandidate cut(boosted)                                                                             
    for (unsigned i=0; i<goodAK4Jets.size(); i++) { 
      if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450 && goodAK4Jets.at(i).getCSV()>0.800){
	TLorentzVector h;
	h= goodAK4Jets.at(i).getP4();
	vlq::Candidate h1(h);
	Hb.push_back(h1);
      }
    }
    //HCandidates (nonboosted)
    HCandsProducer h;
    h.operator()(goodAK4Jets.size(), 2, goodAK4Jets,H)  ;

    //Zcandidate cut (boosted)                                     
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {                                                  
      if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
	TLorentzVector zb;    
	zb= goodAK4Jets.at(i).getP4();
	vlq::Candidate zb1(zb);       
	ZB.push_back(zb1);
	
      }                                                                                                                                       
    }              
    
    //Z candidate cut (non boosted)                                                                                                            
    ZCandsProducer z;
    z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;
    
    //Top Candidates (Category D)    
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {
      if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
	TLorentzVector d;
	d= goodAK4Jets.at(i).getP4();
	vlq::Candidate d1(d);
	D.push_back(d1);
      }
    }
    
    TopCandsProducer top,w;
    // Category BC
    w.operator()(goodAK4Jets, W,B) ;
    for (unsigned i=0; i<W.size(); i++) {
      for (unsigned j=0; j<B.size(); j++) {
	TLorentzVector bc1;
	bc1= W.at(i).getP4()+B.at(j).getP4();
	if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
	  vlq::Candidate bc2(bc1);
	  BC.push_back(bc2); 
	}  
      }
    }
    // category A
    top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;
    
    for (unsigned i=0; i<Hb.size(); i++) {
      h1_["H_mass_b_sig"] -> Fill(Hb.at(i).getMass(), evtwt) ;
      h1_["H_Pt_b_sig"] -> Fill(Hb.at(i).getPt(), evtwt) ;
    }
    h1_["nHcandidatejets_b_sig"] -> Fill(Hb.size(), evtwt) ;
    
    //Hcandidate cut(nonboosted)
    for (unsigned i=0; i<H.size(); i++) {
      h1_["H_mass_nb_sig"] -> Fill(H.at(i).getMass(), evtwt) ;
      h1_["H_Pt_nb_sig"] -> Fill(H.at(i).getPt(), evtwt) ;
    }
    
    h1_["nHcandidatejets_nb_sig"] -> Fill(H.size(), evtwt) ;

    double  nHcandidates=0.0;
    double  nHcandidates1=0.0;
    if(Hb.size()>0 ||H.size()>0){
      nHcandidates = Hb.size()+H.size();
      
      h1_["nHcandidatejets_sig"] -> Fill(nHcandidates, evtwt) ;
      
    }
    
    nHcandidates1 = Hb.size()+H.size();
    h1_["nHcandidatejets1_sig"] -> Fill(nHcandidates1, evtwt) ;
    
    
    //Zcandidate cut (boosted)
    for (unsigned i=0; i<ZB.size(); i++) {     
      h1_["Z_mass_a_sig"] -> Fill(ZB.at(i).getMass(), evtwt) ;
      h1_["Z_Pt_a_sig"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
    }                                                                                                                                                                                       
    h1_["nzcandidatejets_a_sig"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
      
    //Z candidate cut (non boosted)                                                                                                                                                         
    for (unsigned i=0; i<Z.size(); i++) {
      h1_["Z_mass_b_sig"] -> Fill(Z.at(i).getMass(), evtwt) ;
      h1_["Z_Pt_b_sig"] -> Fill(Z.at(i).getPt(), evtwt) ;
    }
    h1_["nzcandidatejets_b_sig"] -> Fill(Z.size(), evtwt) ;
    double nzcandidates=0.0;
    double nzcandidates1=0.0;
    if ( ZB.size() || Z.size()){
      nzcandidates = ZB.size()+ Z.size();
      h1_["nzcandidatejets_tot_sig"] -> Fill(nzcandidates, evtwt) ;
      
    }
    nzcandidates1 = ZB.size()+ Z.size();
    h1_["nzcandidatejets1_tot_sig"] -> Fill(nzcandidates1, evtwt) ;
    
    
    // Category D    
    for (unsigned i=0; i<D.size(); i++) {
      h1_["top_mass_d_sig"] -> Fill(D.at(i).getMass(), evtwt) ;
      h1_["top_Pt_d_sig"] -> Fill(D.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_d_sig"] -> Fill(D.size(), evtwt) ;
    
    // Category BC
    for (unsigned i=0; i<W.size(); i++) {
      h1_["W_mass_bc_sig"] -> Fill(W.at(i).getMass(), evtwt) ;
    }
      h1_["nWcandidatejets_bc_sig"] -> Fill(W.size(), evtwt) ;
      
      for (unsigned i=0; i<B.size(); i++) {
	h1_["lightjet_mass_bc_sig"] -> Fill(B.at(i).getMass(), evtwt) ; 
      }
      h1_["nlightjetcandidatejets_bc_sig"] -> Fill(B.size(), evtwt) ;
      


      for (unsigned i=0; i<BC.size(); i++) {
	h1_["top_mass_bc_sig"] -> Fill(BC.at(i).getMass(), evtwt) ;
	h1_["top_Pt_bc_sig"] -> Fill(BC.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_bc_sig"] -> Fill(BC.size(), evtwt) ;
      
      for (unsigned i=0; i<tops.size(); i++) {
	h1_["top_mass_a_sig"] -> Fill(tops.at(i).getMass(), evtwt) ;
	h1_["top_Pt_a_sig"] -> Fill(tops.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_a_sig"] -> Fill(tops.size(), evtwt) ;
      
      double  ntopcandidates=0.0;
      double  ntopcandidates1=0.0;
      if(D.size() || BC.size() ||tops.size()){
	ntopcandidates = D.size()+BC.size()+tops.size();
	
	h1_["ntopcandidatejets_sig"] -> Fill(ntopcandidates, evtwt) ;  
	
      } 
      ntopcandidates1 = D.size()+BC.size()+tops.size();
      h1_["ntopcandidatejets1_sig"] -> Fill(ntopcandidates1, evtwt) ;
      //Z and top corelations and ST tempelates
      h1_["st_sig"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(1, evtwt) ;
      if (goodBTaggedAK4Jets.size() == 1){
	h1_["cutflow3"] -> Fill(2, evtwt) ;   
	h1_["st_sig1b"] -> Fill(ST,evtwt);
      }
      
      if (goodBTaggedAK4Jets.size() >=2){
	h1_["st_sig2b"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(3, evtwt) ;
      }
      //n,Z,H,B
      //(1)
      if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
	h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(4, evtwt) ;
	if (nHcandidates >= 1.0){
	  h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(8, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(1, evtwt) ;
	    h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(2, evtwt) ;
	    h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	  
	}
	else if (nHcandidates == 0.0){
	  h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(9, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(3, evtwt) ;
	    h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(4, evtwt) ;
	    h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}  
      }
      //(2)                                                                                                                                                                                        
      if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
	h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(5, evtwt) ;
	if (nHcandidates >= 1.0){
	  h1_["st_sigT0Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(10, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(5, evtwt) ;
	    h1_["st_sigT0Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(6, evtwt) ;
	    h1_["st_sigT0Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  h1_["st_sigT0Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(11, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(7, evtwt) ;
	    h1_["st_sigT0Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(8, evtwt) ;
	    h1_["st_sigT0Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
      //(3)                                                                                                                                                                                        
      if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
	h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(6, evtwt) ;
	if (nHcandidates >= 1.0){
	  h1_["st_sigT1Z0H1"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(12, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(9, evtwt) ;
	    h1_["st_sigT1Z0H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(10, evtwt) ;
	    h1_["st_sigT1Z0H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  h1_["st_sigT1Z0H0"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(13, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(11, evtwt) ;
	    h1_["st_sigT1Z0H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(12, evtwt) ;
	    h1_["st_sigT1Z0H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
      //(4)                                                                                                                                                                                        
      if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
	h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(7, evtwt) ;
	if (nHcandidates >= 1.0){
	  h1_["st_sigT0Z0H1"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(14, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(13, evtwt) ;
	    h1_["st_sigT0Z0H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(14, evtwt) ;
	    h1_["st_sigT0Z0H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  h1_["st_sigT0Z0H0"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(15, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow4"] -> Fill(15, evtwt) ;
	    h1_["st_sigT0Z0H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow4"] -> Fill(16, evtwt) ;
	    h1_["st_sigT0Z0H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
  }


  //Do mass reconstruction
  //  MassReco reco;

  TLorentzVector lep1, lep2;
  if (zdecayMode_ == "zelel"){
    lep1 = goodElectrons.at(0).getP4();
    lep2 = goodElectrons.at(1).getP4();
  }
  else if (zdecayMode_ == "zmumu"){
    lep1 = goodMuons.at(0).getP4();
    lep2 = goodMuons.at(1).getP4();
  }
  Leptons = lep1 + lep2;
  //TLorentzVector Leptons = zll.at(0).getP4();

  if (optimizeReco_ && *h_evttype.product() != "EvtType_Data"){

    GenParticleCollection genPartsInfo;
    genPartsInfo = genpart(evt) ;
    // Declare TLorentzVectors to fill with genParticles
    TLorentzVector bGen, bbarGen, q1, q2;// Z1, Z2;
    TLorentzVector qJet, qbarJet, bJet, bbarJet;
    TLorentzVector had_bjet, lep_bjet, had_bGen, lep_bGen;
    bGen = reco.getGen(genPartsInfo, 5, 8000002);
    bbarGen = reco.getGen(genPartsInfo, -5, 8000002);
    q1 = reco.getGen(genPartsInfo, 1, 5, 8000002);
    q2 = reco.getGen(genPartsInfo, -5, -1, 8000002);

    qJet = reco.getMatchedJet(q1, goodAK4Jets, 0.3);
    qbarJet = reco.getMatchedJet(q2, goodAK4Jets, 0.3);
    bJet = reco.getMatchedJet(bGen, goodAK4Jets, 0.3);
    bbarJet = reco.getMatchedJet(bbarGen, goodAK4Jets, 0.3);

    //Choose charge of b (Change based on B' mass)
    double bcheck = abs((bGen + q1 + q2).M() - vlqMass_);
    double bbarcheck = abs((bbarGen + q1 + q2).M() - vlqMass_);
    if (bcheck < bbarcheck){
      had_bGen = bGen;
      lep_bGen = bbarGen;
      had_bjet = bJet;
      lep_bjet = bbarJet;
    }
    else {
      had_bGen = bbarGen;
      lep_bGen = bGen;
      had_bjet = bbarJet;
      lep_bjet = bJet;
    }
    double genZ = reco.findInvMass(q1, q2);
    double genB = reco.findInvMass(q1, q2, had_bGen);
    double ZJet = reco.findInvMass(qJet, qbarJet);
    double hadBJet = reco.findInvMass(qJet, qbarJet, had_bjet);
    double lepBJet = reco.findInvMass(lep1, lep2, lep_bjet);

    h1_["genZ"]->Fill(genZ, evtwt);
    h1_["genBMass"] -> Fill(genB, evtwt);
    h1_["ZJetMass"]->Fill(ZJet,evtwt);
    h1_["hadBJetMass"] ->Fill(hadBJet, evtwt);
    h1_["lepBJetMass"] ->Fill(lepBJet, evtwt);
  }
  pair<double, double> chi2_result;
  if (goodAK4Jets.size() > 4)
    chi2_result = reco.doReco(goodAK4Jets, bosonMass_, Leptons);
  //else if (goodAK8Jets.size() > 0)
  //chi2_result = reco.doReco(goodAK4Jets, goodAK8Jets.at(0).getP4(), bosonMass_, Leptons);
  else{
    chi2_result.first = -999;
    chi2_result.second = -999;
  }

  //Fill Histograms
  h1_["ZJetMasslep"] ->Fill(Leptons.M(), evtwt);
  h1_["chi2_chi"] ->Fill(chi2_result.first, evtwt);
  h1_["sqrtChi2"] ->Fill(sqrt(chi2_result.first), evtwt);
  h1_["chi2_mass"] ->Fill(chi2_result.second, evtwt);

  std::auto_ptr<vlq::JetCollection> ptr_tjets( new vlq::JetCollection(goodTopTaggedJets) ) ; 
  std::auto_ptr<vlq::JetCollection> ptr_wjets( new vlq::JetCollection(goodWTaggedJets) ) ; 
  std::auto_ptr<vlq::JetCollection> ptr_bjets( new vlq::JetCollection(goodBTaggedAK4Jets ) ) ; 
  std::auto_ptr<vlq::JetCollection> ptr_jets ( new vlq::JetCollection(goodAK4Jets ) ) ; 
  std::auto_ptr<vlq::CandidateCollection> ptr_zllcands ( new vlq::CandidateCollection(zll) ) ; 

  evt.put(ptr_tjets, "tjets") ; 
  evt.put(ptr_wjets, "wjets") ; 
  evt.put(ptr_bjets, "bjets") ; 
  evt.put(ptr_jets , "jets")  ; 
  evt.put(ptr_zllcands , "zllcands")  ; 

  return true ; 
}

// ------------ method called once each job just before starting event loop  ------------
void OS2LAna::beginJob() {
  if(zdecayMode_ == "zmumu") {lep = "mu";}
  else if ( zdecayMode_ == "zelel") {lep = "el";}
  else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;

  h1_["checkPU"] = fs->make<TH1D>("checkPU", "Initial NPV", 51, -0.5, 50.5);

  if (filterSignal_){
    if(doSkim_){
      const int nCh = 12;
      const char *channel[nCh] = {"bZbZ", "bZbH", "bZtW", "bHbH", "bHtW", "tWtW",
                                  "tZtZ", "tZtH", "tZbW", "tHtH", "tHbW", "bWbW"};
      h1_["signalEvts_all"] = fs->make<TH1D>("signalEvts_all", "All signal events", 12, 0.5, 12.5) ;
      for (int i=1;i<=nCh;i++) h1_["signalEvts_all"]->GetXaxis()->SetBinLabel(i,channel[i-1]);
    }
    else{
      h1_["signalEvts"] = fs->make<TH1D>("signalEvts", "signal events", 2, 0.5, 2.5) ;
    }
  }

  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 8, 0.5, 8.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1, "Trig.+l^{+}l^{-}") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2, "75 #lt M(l^{+}l^{-}) #lt 105") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4, "N(AK4) #geq 3") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3, "H_{T} #geq 200") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5, "leading jet pt > 100") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6, "2nd jet pt > 50") ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7, "N(b jet) #geq 1") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8, "S_{T} #geq 1000") ; 
  //h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9, "N(AK8) #geq 1") ; //this is not the event cut


  if (!doSkim_){
    TFileDirectory pre = fs->mkdir ("pre");
    TFileDirectory sig = fs->mkdir ("sig");
    TFileDirectory cnt = fs->mkdir ("cnt");
    TFileDirectory cat = fs->mkdir ("cat");
    TFileDirectory *bookDir[4]; bookDir[0] = &pre; bookDir[1] = &cnt; bookDir[2] = &sig; bookDir[3] = &cat;
    std::vector<string> suffix = {"_pre", "_cnt", ""};
    
    for (int i=0; i<3; i++){
      h1_[("npv_noweight"+suffix[i]).c_str()] = bookDir[i]->make<TH1D>( ("npv_noweight"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("npv"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("npv"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("nak4"+suffix[i]).c_str()] =  bookDir[i]->make<TH1D>( ("nak4"+suffix[i]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
      h1_[("ht"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("ht"+suffix[i]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
      h1_[("st"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("st"+suffix[i]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
      h1_[("met"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met"+suffix[i]).c_str(), "MET [GeV]", 100, 0., 1000.);
      h1_[("metPhi"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("metPhi"+suffix[i]).c_str(), "#Phi(MET)", 20, -5., 5.);
      
      //jets
      for(int j=1; j<4; ++j){
	string jetPtName = Form("ptak4jet%d", j)+suffix[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
	h1_[jetPtName.c_str()] = bookDir[i]->make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
	string jetEtaName = Form("etaak4jet%d", j)+suffix[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
	h1_[jetEtaName.c_str()] = bookDir[i]->make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-4. ,4.) ;
	string jetCVSName = Form("cvsak4jet%d", j)+suffix[i]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j); 
	h1_[jetCVSName.c_str()] = bookDir[i]->make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
      }
      string jet1METPhiName = "phi_jet1MET"+suffix[i];
      h1_[jet1METPhiName.c_str()] = bookDir[i]->make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;
      
      //leptons
      string mass_Z = "mass_z"+lep+lep+suffix[i];
      h1_[mass_Z.c_str()] = bookDir[i]->make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
      
      //shrinked mass window
      string mass_Z1 = "mass_Z"+lep+lep+suffix[i];
      h1_[mass_Z1.c_str()] = bookDir[i]->make<TH1D>(mass_Z1.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 70., 120.) ;
      
      string dr_ll = "dr_"+lep+lep+suffix[i];  
      h1_[dr_ll.c_str()] = bookDir[i]->make<TH1D>(dr_ll.c_str(), ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
      string pt_Z = "pt_z"+lep+lep+suffix[i];
      h1_[pt_Z.c_str()] = bookDir[i]->make<TH1D>(pt_Z.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ; 
      for(int l=1; l<3; ++l){
	string lepPtName = "pt_"+lep+Form("%d",l)+suffix[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
	h1_[lepPtName.c_str()] = bookDir[i]->make<TH1D>(lepPtName.c_str(), lepPtTitle.c_str(), 50, 0., 500.) ;
	string lepEtaName = "eta_"+lep+Form("%d",l)+suffix[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
	h1_[lepEtaName.c_str()] = bookDir[i]->make<TH1D>(lepEtaName.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
      }
    }
    h1_["test1"] = pre.make<TH1D>("test1", ";Mzll;;" , 100,20.,220.) ;
    h1_["test2"] = pre.make<TH1D>("test2", ";Mzll;;" , 100,20.,220.) ;
    h1_["test3"] = pre.make<TH1D>("test3", ";Mzll;;" , 100,20.,220.) ;
    

    //additional plots
    h1_["nbjets"] = sig.make<TH1D>("nbjets", ";N(b jets);;" , 11, -0.5,10.5) ; 
    //addition nb plots
    h1_["nbjets_met_sig"] = sig.make<TH1D>("nbjets_met_sig", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_cnt"] = sig.make<TH1D>("nbjets_met_cnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_0btagcnt"] = sig.make<TH1D>("nbjets_met_0btagcnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_1btagcnt"] = sig.make<TH1D>("nbjets_met_1btagcnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    
    h1_["ht_met_0btagcnt"]   =  sig.make<TH1D>( "ht_met_0btagcnt", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["1b_ht"]   =  sig.make<TH1D>( "1b_ht", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["ht_met_1btagcnt"]   =  sig.make<TH1D>( "ht_met_1btagcnt", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["lowmet_ht"]   =  sig.make<TH1D>("lowmet_ht", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    
    h1_["st_met_0btagcnt"]   =  sig.make<TH1D>( "st_met_0btagcnt", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["1b_st"]   =  sig.make<TH1D>( "1b_st", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["st_met_1btagcnt"]   =  sig.make<TH1D>( "st_met_1btagcnt", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["lowmet_st"]   =  sig.make<TH1D>("lowmet_st", ";S_{T} [GeV]", 100, 0., 4000.) ;

    h1_["btagSF_0btag"] = sig.make<TH1D>("btagSF_0btag", ";btag SF (0 Btag region);;" ,100 ,0.,1.0) ;
    h1_["btagSF_1btag"] = sig.make<TH1D>("btagSF_1btag", ";btag SF (>=1 Btag region);;" ,100 ,0.,1.0) ;

    h1_["nbjets_cnt"] = sig.make<TH1D>("nbjets_cnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    
    h1_["ptbjetleading"]  = sig.make<TH1D>("ptbjetleading", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etabjetleading"] = sig.make<TH1D>("etabjetleading", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;
    
    h1_["nak8"] = sig.make<TH1D>("nak8", ";N(AK8 jets);;" , 11, -0.5,10.5) ; 
    h1_["nwjet"] = sig.make<TH1D>("nwjet", ";N(W jets );;" , 6, -0.5,5.5) ; 
    h1_["nhjet"] = sig.make<TH1D>("nhjet", ";N(H jets );;" , 6, -0.5,5.5) ; 
    h1_["ntjet"] = sig.make<TH1D>("ntjet", ";N(top jets);;" , 6, -0.5,5.5) ; 
    
    h1_["ptak8leading"]  = sig.make<TH1D>("ptak8leading", ";p_{T}(leading AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak8leading"] = sig.make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["mak8leading"] = sig.make<TH1D>("mak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["prunedmak8leading"] = sig.make<TH1D>("prunedmak8leading", ";M(pruned leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["trimmedmak8leading"] = sig.make<TH1D>("trimmedmak8leading", ";M(trimmed leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["softdropmak8leading"] = sig.make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["ptak82nd"]  = sig.make<TH1D>("ptak82nd", ";p_{T}(2nd AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak82nd"] = sig.make<TH1D>("etaak82nd", ";#eta(2nd AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["mak82nd"] = sig.make<TH1D>("mak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["prunedmak82nd"] = sig.make<TH1D>("prunedmak82nd", ";M(pruned 2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["trimmedmak82nd"] = sig.make<TH1D>("trimmedmak82nd", ";M(trimmed 2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["softdropmak82nd"] = sig.make<TH1D>("softdropmak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ;

    h1_["ptTprime"]  = sig.make<TH1D>("ptTprime", ";p_{T}(T quark) [GeV];;" , 100, 0., 2000.) ; 
    h1_["yTprime"] = sig.make<TH1D>("yTprime", ";y(T quark);;" , 40 ,-4. ,4.) ; 
    h1_["mTprime"] = sig.make<TH1D>("mTprime", ";M(T quark) [GeV];;" ,100 ,0., 2000.) ; 
    
    h1_["ptBprime"]  = sig.make<TH1D>("ptBprime", ";p_{T}(B quark) [GeV];;" , 100, 0., 2000.) ; 
    h1_["yBprime"] = sig.make<TH1D>("yBprime", ";y(B quark);;" , 40 ,-4. ,4.) ; 
    h1_["mBprime"] = sig.make<TH1D>("mBprime", ";M(B quark) [GeV];;" ,100 ,0., 2000.) ; 
    
    // b-tagging efficiency maps:
    h2_["pt_eta_b_all"] = fs->make<TH2D>("pt_eta_b_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_all"] = fs->make<TH2D>("pt_eta_c_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_all"] = fs->make<TH2D>("pt_eta_l_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    
    h2_["pt_eta_b_btagged"] = fs->make<TH2D>("pt_eta_b_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_btagged"] = fs->make<TH2D>("pt_eta_c_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_btagged"] = fs->make<TH2D>("pt_eta_l_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    
    // plots to study B mass recontruction
    if (optimizeReco_){
      h1_["genZ"] = sig.make<TH1D>("genZ", ";M (Gen Z Boson) [GeV];;", 20, 0., 200.);
      h1_["genBMass"] = sig.make<TH1D>("genBMass", ";M(Gen B quark) [GeV];;", 100, 800., 1200);
      h1_["ZJetMass"] = sig.make<TH1D>("ZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
      h1_["hadBJetMass"] = sig.make<TH1D>("BJetMass", ";JetM (Hadronic B quark) [GeV];;", 100, 500., 1500.);
      h1_["lepBJetMass"] = sig.make<TH1D>("BJetMasslep", ";M (B Jet Leptonic);;", 100, 500., 1500.);
    }
    
    h1_["ZJetMasslep"] = sig.make<TH1D>("ZJetMasslep", ";M (Z Jet Leptonic);;", 20, 0., 200.);
    h1_["chi2_chi"] = sig.make<TH1D>("chi2_chi", ";#chi^{2};;", 100, 0., 500.);
    h1_["sqrtChi2"] = sig.make<TH1D>("sqrtChi2", ";#chi;;", 100, 0., 500.);
    h1_["chi2_mass"] = sig.make<TH1D>("chi_mass", ";M_{#chi^{2}}(B);;", 180, 200., 2000.);
  
    //electrons specific varaibles in EE and EB at preselection level
  if (zdecayMode_ == "zelel" && additionalPlots_){
    h1_["Eta_EB_el_pre"] = pre.make<TH1D>("Eta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
    h1_["Eta_EE_el_pre"] = pre.make<TH1D>("Eta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ;
    h1_["Iso03_EB_el_pre"] = pre.make<TH1D>("Iso03_EB_el_pre", ";Iso03 (EB);;", 100,0,0.3) ;
    h1_["Iso03_EE_el_pre"] = pre.make<TH1D>("Iso03_EE_el_pre", ";Iso03 (EE);;", 100,0,0.3) ;
    h1_["dEtaIn_EB_el_pre"] = pre.make<TH1D>("dEtaIn_EB_el_pre", ";dEtaIn (EB);;", 200,-0.05,0.05) ;
    h1_["dEtaIn_EE_el_pre"] = pre.make<TH1D>("dEtaIn_EE_el_pre", ";dEtaIn (EE);;", 200,-0.05,0.05) ;
    h1_["dPhiIn_EB_el_pre"] = pre.make<TH1D>("dPhiIn_EB_el_pre", ";dPhiIn (EB);;", 100,-0.2,0.2) ;
    h1_["dPhiIn_EE_el_pre"] = pre.make<TH1D>("dPhiIn_EE_el_pre", ";dPhiIn (EE);;", 100,-0.2,0.2);
    h1_["Dz_EB_el_pre"] = pre.make<TH1D>("Dz_EB_el_pre",";dZ (EB);;", 200,-0.1,0.1) ;
    h1_["Dz_EE_el_pre"] = pre.make<TH1D>("Dz_EE_el_pre", ";dZ (EE);;", 200,-0.4,0.4) ;
    h1_["D0_EB_el_pre"] = pre.make<TH1D>("D0_EB_el_pre", ";d0 (EB);;", 100,-0.1,0.1) ;
    h1_["D0_EE_el_pre"] = pre.make<TH1D>("D0_EE_el_pre", ";d0 (EE);;", 100,-0.1,0.1) ;
    h1_["SCETA_EB_el_pre"] = pre.make<TH1D>("scEta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
    h1_["SCETA_EE_el_pre"] = pre.make<TH1D>("scEta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ; 
    h1_["Full5x5siee_EB_el_pre"] = pre.make<TH1D>("Full5x5siee_EB_el_pre", ";Full5X5SigmaIEtaIEta (EB);;", 200,0,0.01) ;
    h1_["Full5x5siee_EE_el_pre"] = pre.make<TH1D>("Full5x5siee_EE_el_pre", ";Full5X5SigmaIEtaIEta (EE);;", 100,0,0.03) ;
    h1_["HoE_EB_el_pre"] = pre.make<TH1D>("HoE_EB_el_pre", ";H/E (EB);;", 200,0,0.05) ;
    h1_["HoE_EE_el_pre"] = pre.make<TH1D>("HoE_EE_el_pre", ";H/E (EE);;", 200,0,0.1) ;
    h1_["ooEmooP_EB_el_pre"] = pre.make<TH1D>("ooEmooP_EB_el_pre", ";(1/E - 1/P) (EB);;", 200,0,0.02) ;
    h1_["ooEmooP_EE_el_pre"] = pre.make<TH1D>("ooEmooP_EE_el_pre", ";(1/E - 1/P) (EE);;", 200,0,0.02) ;
    h1_["missHits_EB_el_pre"] = pre.make<TH1D>("missHits_EB_el_pre", ";Expected missing Hits (EB);;", 4,-0.5,3.5) ;
    h1_["missHits_EE_el_pre"] = pre.make<TH1D>("missHits_EE_el_pre", ";Expected missing Hits (EE);;", 4,-0.5,3.5) ;
    h1_["conveto_EB_el_pre"] = pre.make<TH1D>("conveto_EB_El_pre", ";has matched Conveto (EB);;", 4,-0.5,3.5) ;
    h1_["conveto_EE_el_pre"] = pre.make<TH1D>("conveto_EE_el_pre", ";has matched Conveto (EE);;", 4,-0.5,3.5) ;


 }

  h1_["Iso04_mu_pre"] = pre.make<TH1D>("Iso04_mu_pre", ";Iso04 ;;", 100,0,0.3) ;
  h1_["Dz_mu_pre"] = pre.make<TH1D>("Dz_mu_pre",";dZ ;;", 200,-0.1,0.1) ;
  h1_["D0_mu_pre"] = pre.make<TH1D>("D0_mu_pre", ";d0;;", 100,-0.1,0.1) ;
  //o batg region plots
  h1_["nob_ht"]= cnt.make<TH1D>("nob_ht", "HT no b tags", 100, 0., 3000.);
  h1_["nob_st"] = fs->make<TH1D>("nob_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["nob_pt_z"+lep+lep] = fs->make<TH1D>(("nob_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["b_pt_z"+lep+lep] = cnt.make<TH1D>(("b_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["b_st"] = cnt.make<TH1D>("b_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["pt_zlight_pre"] = fs->make<TH1D>("pt_zlight_pre", "p_{T} (Z + q_{light}) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zb_pre"] = fs->make<TH1D>("pt_zb_pre", "p_{T} (Z + b) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zc_pre"] = fs->make<TH1D>("pt_zc_pre", "p_{T} (Z + c) [GeV]", 100, 0., 2000.) ;
  
  //test plota
  
  
  h1_["httest"]   = sig.make<TH1D>( "httest",";H_{T}(AK4 jets) [GeV];;", 100, 0., 4000.) ;
  h1_["ptak4led"]  =sig.make<TH1D>("ptak4led", ";p_{T}(leading AK4 jet) [GeV];;" , 50, 0., 1000.) ;
  h1_["ptak4sub"]  =sig.make<TH1D>("ptak4sub", ";p_{T}(leading AK4 jet) [GeV];;" , 50, 0., 1000.) ;
  
  h1_["le1E"]   = pre.make<TH1D>( "le1E",";Energy (lepton 1) [GeV];;", 50, 0., 2000.) ;
  h1_["le2E"]   = pre.make<TH1D>( "le2E",";Energy (lepton 2) [GeV];;", 50, 0., 1500.) ;
  h1_["zllE"]   = pre.make<TH1D>( "zllE",";Energy (Z Candidate) [GeV];;", 50, 0., 1500.) ;
  h1_["el1E"]   = pre.make<TH1D>( "el1E",";Energy (lepton 1) [GeV];;", 50, 0., 2000.) ;
  h1_["el2E"]   = pre.make<TH1D>( "el2E",";Energy (lepton 2) [GeV];;", 50, 0., 1500.) ;

  h1_["le1E_pre"]   = pre.make<TH1D>( "le1E_pre",";Energy (lepton 1) [GeV];;", 50, 0., 2000.) ;
  h1_["le2E_pre"]   = pre.make<TH1D>( "le2E_pre",";Energy (lepton 2) [GeV];;", 50, 0., 1500.) ;
  h1_["zllE_pre"]   = pre.make<TH1D>( "zllE_pre",";Energy (Z Candidate) [GeV];;", 50, 0., 1500.) ;
  h1_["el1E_pre"]   = pre.make<TH1D>( "el1E_pre",";Energy (lepton 1) [GeV];;", 50, 0., 2000.) ;
  h1_["el2E_pre"]   = pre.make<TH1D>( "el2E_pre",";Energy (lepton 2) [GeV];;", 50, 0., 1500.) ;


  h1_["nle1"] = pre.make<TH1D>("nle1", ";N;;" , 20, -0.5,20.5) ;
  h1_["nle2"] = pre.make<TH1D>("nle2", ";N;;" , 20, -0.5,20.5) ;
  h1_["nzll1"] = pre.make<TH1D>("nzll1", ";N;;" , 20, -0.5,20.5) ;
  h1_["nel1"] = pre.make<TH1D>("nel1", ";N;;" , 20, -0.5,20.5) ;
  h1_["sumle"] = pre.make<TH1D>("sumle", ";N;;" , 20, -0.5,20.5) ;
  
  if (categorize_){
    
    h1_["cutflow1"] = cat.make<TH1D>("cutflow1", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(1, "no T,Z,b ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(4, "T1Z1 ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(5, "T0Z1") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(6, "T1Z0 ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(7, "T0Z0") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(8, "T1Z1H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(9, "T1Z1H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(10, "T0Z1H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(11, "T0Z1H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(12, "T1Z0H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(13, "T1Z0H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(14, "T0Z0H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(15, "T0Z0H0 " ) ;


    h1_["cutflow2"] = cat.make<TH1D>("cutflow2", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(1, "T1Z1H1b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(2, "T1Z1H1b2") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(3, "T1Z1H0b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(4, "T1Z1H0b2  ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(5, "T0Z1H1b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(6, "T0Z1H1b2 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(7, "T0Z1H0b1") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(8, "T0Z1H0b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(9, "T1Z0H1b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(10, "T1Z0H1b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(11, "T1Z0H0b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(12, "T1Z0H0b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(13, "T0Z0H1b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(14, "T0Z0H1b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(15, "T0Z0H0b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(16, "T0Z0H0b2 " ) ;


    h1_["cutflow3"] = cat.make<TH1D>("cutflow3", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(1, "no T,Z,b ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(4, "T1Z1 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(5, "T0Z1") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(6, "T1Z0 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(7, "T0Z0") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(8, "T1Z1H1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(9, "T1Z1H0 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(10, "T0Z1H1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(11, "T0Z1H0 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(12, "T1Z0H1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(13, "T1Z0H0 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(14, "T0Z0H1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(15, "T0Z0H0 " ) ;


    h1_["cutflow4"] = cat.make<TH1D>("cutflow4", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(1, "T1Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(2, "T1Z1H1b2") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(3, "T1Z1H0b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(4, "T1Z1H0b2  ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(5, "T0Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(6, "T0Z1H1b2 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(7, "T0Z1H0b1") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(8, "T0Z1H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(9, "T1Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(10, "T1Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(11, "T1Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(12, "T1Z0H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(13, "T0Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(14, "T0Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(15, "T0Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(16, "T0Z0H0b2 " ) ;


    //H candidates                                                                                                           
    h1_["H_mass_b_cnt"]  = sig.make<TH1D>("Hmass-boosted-cnt", ";M(H-boosted) [GeV];;", 100, 0., 400);
    h1_["H_Pt_b_cnt"]  = sig.make<TH1D>("HPt-boosted-cnt", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_b_cnt"] = sig.make<TH1D>("nHcandidate-boosted-cnt", ";N(H jets-boosted);;" , 21, -0.5, 20.5);

    h1_["H_mass_nb_cnt"]  = sig.make<TH1D>("Hmassnb-cnt", ";M(H) [GeV];;", 100, 0., 400);
    h1_["H_Pt_nb_cnt"]  = sig.make<TH1D>("HPtnb-cnt", ";Pt(H) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_nb_cnt"] = sig.make<TH1D>("nHcandidatesnb-cnt", ";N(H jets);;" , 21, -0.5, 20.5);

    h1_["nHcandidatejets_cnt"] = sig.make<TH1D>("nHcandidates-tot-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
    h1_["nHcandidatejets1_cnt"] = sig.make<TH1D>("nHcandidates1-tot-cnt", ";N(H jets);;" , 21, -0.5, 20.5);

    // Z candidates
    h1_["Z_mass_a_cnt"]  = sig.make<TH1D>("Zmass-boosted-cnt", ";M(Z-boosted) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_a_cnt"]  = sig.make<TH1D>("ZPt-boosted-cnt", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_a_cnt"] = sig.make<TH1D>("nzcandidate-boosted-cnt", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);

    h1_["Z_mass_b_cnt"]  = sig.make<TH1D>("Zmass-cnt", ";M(Z) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_b_cnt"]  = sig.make<TH1D>("ZPt-cnt", ";Pt(Z) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_b_cnt"] = sig.make<TH1D>("nzcandidates-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);
  
    h1_["nzcandidatejets_tot_cnt"] = sig.make<TH1D>("nzcandidates-tot-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);
    h1_["nzcandidatejets1_tot_cnt"] = sig.make<TH1D>("nzcandidates1-tot-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);  
    // cat A
    h1_["top_mass_a_cnt"]  = sig.make<TH1D>("topmas-A-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_a_cnt"]  = sig.make<TH1D>("topPt-A-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_a_cnt"] = sig.make<TH1D>("ntopcandidate-A-cnt", ";N(top jets);;" , 21, -0.5, 20.5);

    h1_["top_mass_bc_cnt"]  = sig.make<TH1D>("topmass-Bc-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_bc_cnt"]  = sig.make<TH1D>("topPt-BC-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_bc_cnt"] = sig.make<TH1D>("ntopcandidate-BC-cnt", ";N(top jets);;" , 21, -0.5, 20.5);

    // cat D
    h1_["top_mass_d_cnt"]  = sig.make<TH1D>("topmass-D-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_d_cnt"]  = sig.make<TH1D>("topPt-D-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_d_cnt"] = sig.make<TH1D>("ntopcandidate-D-cnt", ";N(top jets);;" , 21, -0.5, 20.5);

    //W and light jet(BC)
    h1_["W_mass_bc_cnt"]  = sig.make<TH1D>("Wmass-BC-cnt", ";M( W boson) [GeV];;", 100, 0., 400);
    h1_["nWcandidatejets_bc_cnt"] = sig.make<TH1D>("nWcandidate-BC-cnt", ";N(W candidate jets);;" , 21, -0.5, 20.5);

    h1_["lightjet_mass_bc_cnt"]  = sig.make<TH1D>("lightjetmass-BC-cnt", ";M( light jet) [GeV];;", 100, 0., 400);
    h1_["nlightjetcandidatejets_bc_cnt"] = sig.make<TH1D>("nlightjetcandidate-cnt", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
    //total top ( A+ BC+D)
    h1_["ntopcandidatejets_cnt"] = sig.make<TH1D>("ntopcandidate-tot-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
    h1_["ntopcandidatejets1_cnt"] = sig.make<TH1D>("ntopcandidate1-tot-cnt", ";N(top jets);;" , 21, -0.5, 20.5);

    //signal region

    //H candidates                                                                                                                                                              
    h1_["H_mass_b_sig"]  = sig.make<TH1D>("Hmass-boosted-sig", ";M(H-boosted) [GeV];;", 100, 0., 400);
    h1_["H_Pt_b_sig"]  = sig.make<TH1D>("HPt-boosted-sig", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_b_sig"] = sig.make<TH1D>("nHcandidate-boosted-sig", ";N(H jets-boosted);;" , 21, -0.5, 20.5);

    h1_["H_mass_nb_sig"]  = sig.make<TH1D>("Hmassnb-sig", ";M(H) [GeV];;", 100, 0., 400);
    h1_["H_Pt_nb_sig"]  = sig.make<TH1D>("HPtnb-sig", ";Pt(H) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_nb_sig"] = sig.make<TH1D>("nHcandidatesnb-sig", ";N(H jets);;" , 21, -0.5, 20.5);

    h1_["nHcandidatejets_sig"] = sig.make<TH1D>("nHcandidates-tot-sig", ";N(H jets);;" , 21, -0.5, 20.5);
    h1_["nHcandidatejets1_sig"] = sig.make<TH1D>("nHcandidates1-tot-sig", ";N(H jets);;" , 21, -0.5, 20.5);

    // Z candidates                                                                                                                                                             
    h1_["Z_mass_a_sig"]  = sig.make<TH1D>("Zmass-boosted-sig", ";M(Z-boosted) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_a_sig"]  = sig.make<TH1D>("ZPt-boosted-sig", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_a_sig"] = sig.make<TH1D>("nzcandidate-boosted-sig", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);

    h1_["Z_mass_b_sig"]  = sig.make<TH1D>("Zmass-sig", ";M(Z) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_b_sig"]  = sig.make<TH1D>("ZPt-sig", ";Pt(Z) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_b_sig"] = sig.make<TH1D>("nzcandidates-sig", ";N(Z jets);;" , 21, -0.5, 20.5);

    h1_["nzcandidatejets_tot_sig"] = sig.make<TH1D>("nzcandidates-tot-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
    h1_["nzcandidatejets1_tot_sig"] = sig.make<TH1D>("nzcandidates1-tot-sig", ";N(Z jets);;" , 21, -0.5, 20.5);

    // cat A                                                                                                                                                                    
    h1_["top_mass_a_sig"]  = sig.make<TH1D>("topmas-A-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_a_sig"]  = sig.make<TH1D>("topPt-A-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_a_sig"] = sig.make<TH1D>("ntopcandidate-A-sig", ";N(top jets);;" , 21, -0.5, 20.5);

    h1_["top_mass_bc_sig"]  = sig.make<TH1D>("topmass-Bc-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_bc_sig"]  = sig.make<TH1D>("topPt-BC-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_bc_sig"] = sig.make<TH1D>("ntopcandidate-BC-sig", ";N(top jets);;" , 21, -0.5, 20.5);

    // cat D                                                                                                                                                                    
    h1_["top_mass_d_sig"]  = sig.make<TH1D>("topmass-D-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_d_sig"]  = sig.make<TH1D>("topPt-D-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_d_sig"] = sig.make<TH1D>("ntopcandidate-D-sig", ";N(top jets);;" , 21, -0.5, 20.5);

    //W and light jet(BC)                                                                                                                                                       
    h1_["W_mass_bc_sig"]  = sig.make<TH1D>("Wmass-BC-sig", ";M( W boson) [GeV];;", 100, 0., 400);
    h1_["nWcandidatejets_bc_sig"] = sig.make<TH1D>("nWcandidate-BC-sig", ";N(W candidate jets);;" , 21, -0.5, 20.5);

    h1_["lightjet_mass_bc_sig"]  = sig.make<TH1D>("lightjetmass-BC-sig", ";M( light jet) [GeV];;", 100, 0., 400);
    h1_["nlightjetcandidatejets_bc_sig"] = sig.make<TH1D>("nlightjetcandidate-sig", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
    //total top ( A+ BC+D)                                                                                                                                                       
    h1_["ntopcandidatejets_sig"] = sig.make<TH1D>("ntopcandidate-tot-sig", ";N(top jets);;" , 21, -0.5, 20.5);
    h1_["ntopcandidatejets1_sig"] = sig.make<TH1D>("ntopcandidate1-tot-sig", ";N(top jets);;" , 21, -0.5, 20.5);

    //ST plots
    //top,Z, Higgs

    h1_["st_sig"] =cat.make<TH1D>("ST_sig", ";S_{T} [Gev];;" , 100,0.,3000.);

    h1_["st_sig1b"] =cat.make<TH1D>("ST_sig1b", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sig2b"] =cat.make<TH1D>("ST_sig2b", ";S_{T} [Gev];;" , 100,0.,3000.);

    h1_["st_sigT1Z1"] =cat.make<TH1D>("ST_sigT1Z1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z1"] =cat.make<TH1D>("ST_sigT0Z1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z0"] =cat.make<TH1D>("ST_sigT1Z0", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z0"] =cat.make<TH1D>("ST_sigT0Z0", ";S_{T} [Gev];;" , 100,0.,3000.);

    h1_["st_sigT1Z1H1"] =cat.make<TH1D>("ST_sigT1Z1H1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z1H0"] =cat.make<TH1D>("ST_sigT1Z1H0", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z1H1"] =cat.make<TH1D>("ST_sigT0Z1H1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z1H0"] =cat.make<TH1D>("ST_sigT0Z1H0", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z0H1"] =cat.make<TH1D>("ST_sigT1Z0H1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z0H0"] =cat.make<TH1D>("ST_sigT1Z0H0", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z0H1"] =cat.make<TH1D>("ST_sigT0Z0H1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z0H0"] =cat.make<TH1D>("ST_sigT0Z0H0", ";S_{T} [Gev];;" , 100,0.,3000.);


    h1_["st_sigT1Z1H1b1"] =cat.make<TH1D>("ST_sigT1Z1H1b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z1H1b2"] =cat.make<TH1D>("ST_sigT1Z1H1b2", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z1H0b1"] =cat.make<TH1D>("ST_sigT1Z1H0b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z1H0b2"] =cat.make<TH1D>("ST_sigT1Z1H0b2", ";S_{T} [Gev];;" , 100,0.,3000.);
  
    h1_["st_sigT0Z1H1b1"] =cat.make<TH1D>("ST_sigT0Z1H1b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z1H1b2"] =cat.make<TH1D>("ST_sigT0Z1H1b2", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z1H0b1"] =cat.make<TH1D>("ST_sigT0Z1H0b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z1H0b2"] =cat.make<TH1D>("ST_sigT0Z1H0b2", ";S_{T} [Gev];;" , 100,0.,3000.);

    h1_["st_sigT1Z0H1b1"] =cat.make<TH1D>("ST_sigT1Z0H1b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z0H1b2"] =cat.make<TH1D>("ST_sigT1Z0H1b2", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z0H0b1"] =cat.make<TH1D>("ST_sigT1Z0H0b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT1Z0H0b2"] =cat.make<TH1D>("ST_sigT1Z0H0b2", ";S_{T} [Gev];;" , 100,0.,3000.);

    h1_["st_sigT0Z0H1b1"] =cat.make<TH1D>("ST_sigT0Z0H1b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z0H1b2"] =cat.make<TH1D>("ST_sigT0Z0H1b2", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z0H0b1"] =cat.make<TH1D>("ST_sigT0Z0H0b1", ";S_{T} [Gev];;" , 100,0.,3000.);
    h1_["st_sigT0Z0H0b2"] =cat.make<TH1D>("ST_sigT0Z0H0b2", ";S_{T} [Gev];;" , 100,0.,3000.);
    
  }
  }
}
double OS2LAna::GetDYNLOCorr(const double dileppt) {

  std::unique_ptr<TFile >file_DYNLOCorr = std::unique_ptr<TFile>(new TFile(fname_DYNLOCorr_.c_str())) ;
  std::unique_ptr<TF1>fun_DYNLOCorr = std::unique_ptr<TF1>(dynamic_cast<TF1*>(file_DYNLOCorr->Get(funname_DYNLOCorr_.c_str()))) ; 
  double EWKNLOkfact = fun_DYNLOCorr->Eval(dileppt);
  return EWKNLOkfact; 

}

double OS2LAna::htCorr(double ht, double p0, double p1){
  return(p1*ht + p0);
}

double OS2LAna::htCorr(double ht, double p0, double p1, double p2, double p3){
  return((p3*ht*ht*ht) +(p2*ht*ht) +(p1*ht) + p0);
}

void OS2LAna::endJob() {

  return ; 
}

DEFINE_FWK_MODULE(OS2LAna);
