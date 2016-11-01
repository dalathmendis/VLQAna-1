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



#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TFitResult.h>

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
    double ZptCorr(vlq::Candidate, double, double);
  double htCorr(double ht, double p0, double p1);

    // ----------member data ---------------------------
    edm::EDGetTokenT<string>   t_evttype         ;
    edm::EDGetTokenT<double>   t_evtwtGen        ;
    edm::EDGetTokenT<double>   t_evtwtPV         ;
    edm::EDGetTokenT<unsigned> t_npv             ;
    edm::EDGetTokenT<bool>     t_hltdecision     ;
    edm::ParameterSet DilepCandParams_           ; 
    edm::ParameterSet ZCandParams_               ; 
    edm::ParameterSet BoostedZCandParams_        ; 
    edm::ParameterSet GenHSelParams_             ;
    edm::ParameterSet genParams_                 ;
    const double HTMin_                          ;
    const double STMin_                          ; 
    const double PtZMin_                         ;



    const bool filterSignal_                     ;
    const bool additionalPlots_                  ;
    const bool doEWKcorr_                        ;
    const std::string signalType_                ;
    const std::string zdecayMode_                ;
    const bool optimizeReco_                     ;
    const double vlqMass_                        ;
    const double bosonMass_                      ;
    const bool applyLeptonSFs_                   ;
    const bool applyBTagSFs_                     ;
    const bool applyZptCorr_                     ;
    const bool applyDYNLOCorr_                   ;
    const std::string fname_DYNLOCorr_           ; 
    const std::string funname_DYNLOCorr_         ; 
    ApplyLeptonSFs lepsfs                        ;
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
    const std::string file_EWK_                            ;
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
      }
      else if  (fabs(scEta > 1.479) && fabs(scEta < 2.5)){
         h1_["Eta_EE_el_pre"]->Fill(goodElectrons.at(iele).getEta(), evtwt);
         h1_["Iso03_EE_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
         h1_["dEtaIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdEtaIn(), evtwt);
         h1_["dPhiIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
         h1_["Dz_EE_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
         h1_["D0_EE_el_pre"]->Fill(goodElectrons.at(iele).getD0(), evtwt);
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
  BoostedZCandParams_     (iConfig.getParameter<edm::ParameterSet> ("BoostedZCandParams")),
  GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
  genParams_              (iConfig.getParameter<edm::ParameterSet> ("genParams")),
  HTMin_                  (iConfig.getParameter<double>            ("HTMin")),
  STMin_                  (iConfig.getParameter<double>            ("STMin")), 
  PtZMin_                 (iConfig.getParameter<double>            ("PtZMin")),

  filterSignal_           (iConfig.getParameter<bool>              ("filterSignal")), 
  additionalPlots_        (iConfig.getParameter<bool>              ("additionalPlots")), 
  doEWKcorr_              (iConfig.getParameter<bool>              ("doEWKcorr")),
  signalType_             (iConfig.getParameter<std::string>       ("signalType")), 
  zdecayMode_             (iConfig.getParameter<std::string>       ("zdecayMode")),
  optimizeReco_           (iConfig.getParameter<bool>              ("optimizeReco")),
  vlqMass_                (iConfig.getParameter<double>            ("vlqMass")),
  bosonMass_              (iConfig.getParameter<double>            ("bosonMass")),
  applyLeptonSFs_         (iConfig.getParameter<bool>              ("applyLeptonSFs")), 
  applyBTagSFs_           (iConfig.getParameter<bool>              ("applyBTagSFs")), 
  applyZptCorr_           (iConfig.getParameter<bool>              ("applyZptCorr")),
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
  //  file_EWK_               (iConfig.getParameter<std::string>              ("File_EWK"))
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

  h1_["checkPU"]->Fill(*h_npv.product(), *h_evtwtGen.product());

  //  unsigned npv(*h_npv.product()) ; 

  if(zdecayMode_ == "zmumu") {lep = "mu";}
  else if ( zdecayMode_ == "zelel") {lep = "el";}
  else edm::LogError("OS2LAna::filter") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;
 
  if (filterSignal_) {
     if (*h_evttype.product()!=signalType_) return false ;
     else  h1_["signalEvts"] -> Fill(1) ;
  }

  const bool hltdecision(*h_hltdecision.product()) ; 
  if ( !hltdecision ) return false;

  //double evtwtgen(*h_evtwtGen.product());
  double evtwt((*h_evtwtGen.product()) * (*h_evtwtPV.product())) ; 

  vlq::MuonCollection goodMuons; 
  muonmaker(evt, goodMuons) ; 

  vlq::ElectronCollection goodElectrons; 
  electronmaker(evt, goodElectrons) ;

  vlq::MetCollection goodMet;
  metmaker(evt, goodMet) ;

  vlq::CandidateCollection dimuons, dielectrons, dileptons, tops, W, B, BC, D, Z,ZB, H,Hb,ZH, ZHb;   
  vlq::CandidateCollection zll, zllBoosted , lepton1, lepton2, zlep1,zlep2; //generic collection

  // dilepton properties: M > 50, lead pt > 45, second pt > 25
  DileptonCandsProducer dileptonsprod(DilepCandParams_) ; 
  // dileptonsprod.operator()<vlq::MuonCollection>(dimuons, goodMuons, lepton1,lepton2); 
  // dileptonsprod.operator()<vlq::ElectronCollection>(dielectrons, goodElectrons, lepton1,lepton2) ;

  dileptonsprod.operator()<vlq::MuonCollection>(dimuons, goodMuons);
  dileptonsprod.operator()<vlq::ElectronCollection>(dielectrons, goodElectrons) ;

  //================================================================
  //First pre-selection: 1) 2 OS dileptons from boosted Z, >=3 jets
  //================================================================
  //dilepton candidate
  if (zdecayMode_ == "zmumu") {dileptons = dimuons; }
  else if (zdecayMode_ == "zelel") {dileptons = dielectrons;}
  if (dileptons.size()  < 1) return false;
  
  //// Get Dy EWK correction
  if ( applyDYNLOCorr_ ) {
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
  
  h1_["cutflow"] -> Fill(1, evtwt) ;

  //Z mass candidate filter: 75 < M < 105, lead pt > 45, 2nd pt > 25, Z pt > 0
  CandidateFilter zllfilter(ZCandParams_) ; 
  zllfilter(dileptons, zll);  
  /* zllfilter(dileptons, zll,lepton1,lepton2, zlep1,zlep2);

  if (zll.size()>0){

    cout <<"zll mass is = "<< zll.at(0).getMass()<<endl;
  }
   TLorentzVector p4l1,p4l2;
   if ( zlep1.size()>0 &&zlep2.size()>0){
     //  TLorentzVector p4l1,p4l2;
    p4l1= zlep1.at(0).getP4();
    p4l2= zlep2.at(0).getP4();
    cout << " deliepton mass is " <<(p4l1+p4l2).Mag()<<endl;
    }*/
   //  zllfilter(dileptons, zll);  // Do Z pt correction
  // if ( applyZptCorr_  && zll.size() > 0) {
  //   evtwt *= ZptCorr(zll.at(0), 1.26117, -0.000903805) ;
  //    }
  
  // jets
  vlq::JetCollection goodAK4Jets;
  jetAK4maker(evt, goodAK4Jets) ;

  CandidateCleaner cleanjets(0.4);
  if (zdecayMode_ == "zmumu") {cleanjets(goodAK4Jets, goodMuons);}
  else if (zdecayMode_ == "zelel") {cleanjets(goodAK4Jets, goodElectrons);} 

  HT htak4(goodAK4Jets) ; 
  double ht = htak4.getHT();
  // cout << "ht: " << ht << endl;
  // cout << evtwt << endl;
  if (applyZptCorr_){
     double corr = htCorr(ht,1.456,-0.000695945); //for electrons
     // double corr = htCorr(ht,1.16045,-0.000492635); // for muons
    if (corr > 0){
      evtwt *= corr;
    }
  }

    // evtwt *= htCorr(ht, 1.19587, -0.000547296); //for muons
  // cout << "ht corr: " << htCorr(ht,  1.19587, -0.000547296) << " evtwt: " << evtwt;

 ////////////////////////////////////////////////////////// 
  //Fill N-1 selected plots for dilepton mass, Ht, ad Njets
  //////////////////////////////////////////////////////////

  //Fill Z mass: >=3 jets, HT > 300
  if (goodAK4Jets.size() > 2 && htak4.getHT() > HTMin_){
     for (auto idilepton : dileptons) h1_["mass_z"+lep+lep+"_pre"] -> Fill(idilepton.getMass(), evtwt) ; 
  }
  //at least one Z cand in event
  if(zll.size() > 0) {h1_["cutflow"] -> Fill(2, evtwt) ;}
  else return false ;

  //Fill Njets: HT > 300 && Z mass
  if (htak4.getHT() > HTMin_){h1_["nak4_pre"] -> Fill(goodAK4Jets.size(), evtwt) ;} 

  //at least 3 AK4 jets in event
  if (goodAK4Jets.size() > 2 ) {h1_["cutflow"] -> Fill(3, evtwt) ;} 
  else return false;

  //Fill HT: >=3 jets, && Z mass
  h1_["ht_pre"] -> Fill(htak4.getHT(), evtwt) ; 
   
  // at least HT > 150 in event
  if ( htak4.getHT() > HTMin_ ) h1_["cutflow"] -> Fill(4, evtwt) ;  
  else return false ; 

  //deltaPhi(MET,lep)
  
  ///////////////////////////////////////////// 
  // fill rest of the plots after pre-selection
  /////////////////////////////////////////////
  /* cout << " zll.at(0).getPt()                       = " <<zll.at(0).getPt()<<endl;
  cout << " zll.at(0).getMass()                     = " <<zll.at(0).getMass()<<endl;
  cout << " (p4l1+ p4l2).getPt()                    = " <<(p4l1+p4l2).Pt()<<endl;
  cout << " p4l1.getPt() + p4l2.getPt()             =  "<< p4l1.Pt() + p4l2.Pt()<<endl;
  cout << " zlep1.at(0).getPt()                     = "   << zlep1.at(0).getPt()<<endl;
  cout << " zlep2.at(0).getPt()                     = " << zlep2.at(0).getPt()<<endl;
  cout << "zlep2.at(0).getPt()+ zlep1.at(0).getPt() = " << zlep1.at(0).getPt()+ zlep2.at(0).getPt()<<endl;
  */
  double ST = htak4.getHT() ;
  ST += zll.at(0).getPt() + goodMet.at(0).getFullPt(); 
  //cout <<"zll PT is :   "<< zll.at(0).getPt();
  // cout <<"lep 1 lep2 pt sum is " << zlep1.at(0).getPt()+zlep2.at(0).getPt();

  // if ( zdecayMode_ == "zmumu" ){
  //  ST += goodMuons.at(0).getPt() + goodMuons.at(1).getPt() + goodMet.at(0).getFullPt();}
  // else if ( zdecayMode_ == "zelel" ){
  //  ST += goodElectrons.at(0).getPt() + goodElectrons.at(1).getPt() + goodMet.at(0).getFullPt();}


  h1_["st_pre"] -> Fill(ST, evtwt) ;   

  //ak4 jet plots
  for(int j=0; j<3; ++j){
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
  vlq::JetCollection goodBTaggedAK4Jets;// btag,nonbtag;
  vlq::CandidateCollection btag,nonbtag;
  jetAK4BTaggedmaker(evt, goodBTaggedAK4Jets) ; 
  /*
  cout << " size of good AK4 jet collection = " << goodAK4Jets.size()<<endl;
  cout << " size of BTagged AK4 jet collection = "<< goodBTaggedAK4Jets.size()<<endl;


  vector<TLorentzVector> A;
  for (unsigned i=0; i<goodAK4Jets.size(); i++){

    //cout << i <<"th good AK4Jet CSV is = " << goodAK4Jets.at(i).getCSV()<<endl;
    // cout << i <<"th good AK4Jet Pt is = " << goodAK4Jets.at(i).getPt()<<endl;
  
    //vector<TLorentzVector> A;
    TLorentzVector jetinfo = goodAK4Jets.at(i).getP4();
    if (goodAK4Jets.at(i).getCSV()>0.800 && goodAK4Jets.at(i).getPt()>50){
      A.push_back(jetinfo);
      cout<< "CSV afetr selection is = " << goodAK4Jets.at(i).getCSV()<<endl;
      cout<< "Pt afetr selection is = " << goodAK4Jets.at(i).getPt()<<endl;
    }

  }
  cout <<"size of selected btags from good AK 4 jets is = " << A.size()<<endl;
  

  cout << " *************************************"<< endl;
  for (unsigned i=0; i<goodBTaggedAK4Jets.size(); i++){

    cout << i <<"th BTagged AK4Jet CSV is = " << goodBTaggedAK4Jets.at(i).getCSV()<<endl;
    cout << i <<"th BTagged AK4Jet Pt is = " << goodBTaggedAK4Jets.at(i).getPt()<<endl;  
}

  */

  if (zdecayMode_ == "zmumu") {cleanjets(goodBTaggedAK4Jets, goodMuons); }
  else if (zdecayMode_ == "zelel") {cleanjets(goodBTaggedAK4Jets, goodElectrons); }  

  MassReco reco;
  TLorentzVector Leptons = zll.at(0).getP4();
  pair<double, double> chi2_result_cnt;

  //fill control plots

  if ( goodBTaggedAK4Jets.size() == 0) {
    /*

 //Zcandidate cut (boosted)
 //cout << " good AK4 jet size 1 = " << goodAK4Jets.size()<<endl;                                                                                                                        
 for (unsigned i=0; i<goodAK4Jets.size(); i++) {                  
   // cout <<i<<"th  AK4 jet mass  mass is " << goodAK4Jets.at(i).getMass()<< endl;                                                                                                   
   if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
     TLorentzVector zb;    
     zb= goodAK4Jets.at(i).getP4();
     vlq::Candidate zb1(zb);       
     ZB.push_back(zb1);
     
   }                                                                                                                                                                                   
 }                                                                                                                                                                                       
 //cout << " good AK4 jet size 2 = " << goodAK4Jets.size()<<endl;                                                                                                                                            
                                   
 //cout<< " ZB size is " << ZB.size()<<endl;
 for (unsigned i=0; i<ZB.size(); i++) {     
   h1_["Z_mass_a"] -> Fill(ZB.at(i).getMass(), evtwt) ;
   h1_["Z_Pt_a"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
 }                                                                                                                                                                                       
 h1_["nzcandidatejets_a"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
 //if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     

//cout << " good AK4 jet size 3 = " << goodAK4Jets.size()<<endl;
//Z candidate cut (non boosted)                                                                                                                                                         
 ZCandsProducer z;
 z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;

 for (unsigned i=0; i<Z.size(); i++) {
   h1_["Z_mass_b"] -> Fill(Z.at(i).getMass(), evtwt) ;
   h1_["Z_Pt_b"] -> Fill(Z.at(i).getPt(), evtwt) ;
 }
 h1_["nzcandidatejets_b"] -> Fill(Z.size(), evtwt) ;

 // if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
 
 double nzcandidates=0.0;
 if ( ZB.size()>0 || Z.size()>0){

  h1_["cutflow"] -> Fill(12, evtwt);
   nzcandidates = ZB.size()+ Z.size();
  h1_["nzcandidatejets_tot"] -> Fill(nzcandidates, evtwt) ;

 }
 
 //cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;
 // else return false ;                                                                                                                                                          
 // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                                                                

 
 
 // Category D    
 
 //cout << " good AK4 jet size 1 = " << goodAK4Jets.size()<<endl;
 for (unsigned i=0; i<goodAK4Jets.size(); i++) {
   // cout <<i<<"th  AK4 jet mass  mass is " << goodAK4Jets.at(i).getMass()<< endl;
   if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
     TLorentzVector d;
     d= goodAK4Jets.at(i).getP4();
     vlq::Candidate d1(d);
     D.push_back(d1);
     }
 }
 
 //cout<< " D size is " << D.size()<<endl;
     
 for (unsigned i=0; i<D.size(); i++) {
   h1_["top_mass_d"] -> Fill(D.at(i).getMass(), evtwt) ;
   h1_["top_Pt_d"] -> Fill(D.at(i).getPt(), evtwt) ;
 }
 h1_["ntopcandidatejets_d"] -> Fill(D.size(), evtwt) ;
 
 //if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
 // else return false ;
 
 //cout << " good AK4 jet size 2 = " << goodAK4Jets.size()<<endl;

 
 TopCandsProducer top,w;
 
 // Category BC
 w.operator()(goodAK4Jets, W,B) ;
 for (unsigned i=0; i<W.size(); i++) {
   //  cout <<i<<"th  W mass is " << W.at(i).getMass()<< endl;
   h1_["W_mass_bc"] -> Fill(W.at(i).getMass(), evtwt) ;
 }
 h1_["nWcandidatejets_bc"] -> Fill(W.size(), evtwt) ;

 for (unsigned i=0; i<B.size(); i++) {
   //cout <<i<<"th  B mass is " << B.at(i).getMass()<< endl;
   h1_["lightjet_mass_bc"] -> Fill(B.at(i).getMass(), evtwt) ; 

}
 h1_["nlightjetcandidatejets_bc"] -> Fill(B.size(), evtwt) ;

 for (unsigned i=0; i<W.size(); i++) {
   for (unsigned j=0; j<B.size(); j++) {
     //  cout <<"mass W/b =[ " <<W.at(i).getMass()<<","<<B.at(j).getMass() <<"]"<< endl;
     TLorentzVector bc1;
     bc1= W.at(i).getP4()+B.at(j).getP4();
     if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
       // cout << "top candidate is " << bc1.M() << endl;
       //cout << "top candidate pt " << bc1.Pt() << endl;  
       vlq::Candidate bc2(bc1);
       BC.push_back(bc2); 
       
     }  
   }
   
   
 }
 for (unsigned i=0; i<BC.size(); i++) {
   // cout <<i<<"th  BC mass is " << BC.at(i).getMass()<< endl;
   // cout <<i<<"th  BC pt is " << BC.at(i).getPt()<< endl;
   h1_["top_mass_bc"] -> Fill(BC.at(i).getMass(), evtwt) ;
   h1_["top_Pt_bc"] -> Fill(BC.at(i).getPt(), evtwt) ;
 }
 h1_["ntopcandidatejets_bc"] -> Fill(BC.size(), evtwt) ;
 
 //if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}
 //else return false ;
 
 
 //cout << " BC category size is " << BC.size()<<endl;
 
 
 // cout << " W size is " << W.size()<< endl;
 // cout << " B size is " << B.size()<< endl;
 
 //cout << " good AK4 jet size 3 = " << goodAK4Jets.size()<<endl;
 
 // category A
 top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;
 
 for (unsigned i=0; i<tops.size(); i++) {
   h1_["top_mass_a"] -> Fill(tops.at(i).getMass(), evtwt) ;
   h1_["top_Pt_a"] -> Fill(tops.at(i).getPt(), evtwt) ;
 }
 h1_["ntopcandidatejets_a"] -> Fill(tops.size(), evtwt) ;
 
 //if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
 // else return false ;
 // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;
  
 double  ntopcandidates=0.0;
 if(D.size()>0 || BC.size()>0 ||tops.size()>0){
   h1_["cutflow"] -> Fill(16, evtwt) ;
   ntopcandidates = D.size()+BC.size()+tops.size();
   
   h1_["ntopcandidatejets"] -> Fill(ntopcandidates, evtwt) ;  
   
 } 
 //Z and top corelations
 //(1)
 if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
   
   h1_["cutflow1"] -> Fill(1, evtwt) ;
   h1_["st_Cat1"] -> Fill(ST, evtwt) ;


 
}
 
 //2)                                                                                                                                                                                     
 else if (ntopcandidates==0.0    &&   nzcandidates >= 1.0){
   
   h1_["cutflow1"] -> Fill(2, evtwt) ;
 
   h1_["st_Cat2"] -> Fill(ST, evtwt) ;
}
 
 //(3)                                                                                                                                                                                     
 else if (ntopcandidates>=1.0    &&   nzcandidates == 0.0){
   
   h1_["cutflow1"] -> Fill(3, evtwt) ;
   h1_["st_Cat3"] -> Fill(ST, evtwt) ;
 }
 
 //(4)                                                                                                                                                                                     
 else if (ntopcandidates == 0.0    &&   nzcandidates== 0.0){
     
     h1_["cutflow1"] -> Fill(4, evtwt) ;
     h1_["st_Cat4"] -> Fill(ST, evtwt) ;   

}
   


 */






















    for (auto izll : zll) {
      h1_["nob_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
    }
    h1_["nob_ht"] ->Fill(htak4.getHT(), evtwt);
    h1_["nob_st"] ->Fill(ST, evtwt);    
    if (goodAK4Jets.size() > 3){
      chi2_result_cnt = reco.doReco(goodAK4Jets, bosonMass_, Leptons);
    }
    else if (goodAK4Jets.size() == 3){
      chi2_result_cnt.first = -998;
      chi2_result_cnt.second = -998;
  }
    else{
      chi2_result_cnt.first = -999;
      chi2_result_cnt.second = -999;
  }
  }
  double btagsf(1) ;
  double btagsf_bcUp(1) ; 
  double btagsf_bcDown(1) ; 
  double btagsf_lUp(1) ; 
  double btagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
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

     btagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetAK4BTaggedmaker.idxjetCSVDiscMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ; 

  }
  //cout << "btag SF in OS2LAna: " << btagsf << endl;
  // apply b tag scale factors
  evtwt *= btagsf;
  if (chi2_result_cnt.second == -998)
     h1_["3jets_cnt"] ->Fill(1, evtwt);
  else if (chi2_result_cnt.second > 0){
    h1_["chi2_chi_cnt"] ->Fill(chi2_result_cnt.first, evtwt);
    h1_["chi_mass_cnt"] ->Fill(chi2_result_cnt.second, evtwt);
  }

  if (goodBTaggedAK4Jets.size() > 0){
   
    
    
    for (auto izll : zll) {
      h1_["b_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
      h1_["b_st"] ->Fill(ST, evtwt);
    }
  }

  //fill control plots
  if ( goodBTaggedAK4Jets.size() > 0 && ST < 700) {
    for (auto izll : zll) {
        h1_["mass_z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;  
        h1_["pt_z"+lep+lep+"_cnt"] -> Fill(izll.getPt(), evtwt) ; 
     }
     h1_["nak4_cnt"] -> Fill(goodAK4Jets.size(), evtwt) ;
     h1_["ht_cnt"] -> Fill(htak4.getHT(), evtwt) ;
     h1_["ht_cnt1"] -> Fill(htak4.getHT(), evtwt) ;
     h1_["st_cnt"] -> Fill(ST, evtwt) ;   
     h1_["st_cnt1"] -> Fill(ST, evtwt) ;
     h1_["npv_noweight_cnt"] -> Fill(*h_npv.product(), *h_evtwtGen.product()); 
     h1_["npv_cnt"] -> Fill(*h_npv.product(), evtwt);
     h1_["met_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
     h1_["met_cnt1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
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
     for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d_cnt", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
     }
     h1_["phi_jet1MET_cnt"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
  }
  /*
  if ( goodBTaggedAK4Jets.size() >=1  && ST < 1000) {
    // cout<< " ak 4 jeet size is " << goodAK4Jets.size() <<endl;
    
    //ZH combined category, only difference from Z is that an extended mass window

    
    // ZHCandidate cut(boosted)                                              
    vector<TLorentzVector> A;
    for (unsigned i=0; i<goodAK4Jets.size(); i++){

      //cout << i <<"th good AK4Jet CSV is = " << goodAK4Jets.at(i).getCSV()<<endl;                                                                            
      // cout << i <<"th good AK4Jet Pt is = " << goodAK4Jets.at(i).getPt()<<endl;                                                                             

      //vector<TLorentzVector> A;                                                                                                                              
      TLorentzVector jetinfo = goodAK4Jets.at(i).getP4();
      if (goodAK4Jets.at(i).getCSV()>0.800 && goodAK4Jets.at(i).getPt()>50){
	A.push_back(jetinfo);
	cout<< "CSV afetr selection is = " << goodAK4Jets.at(i).getCSV()<<endl;
	cout<< "Pt afetr selection is = " << goodAK4Jets.at(i).getPt()<<endl;
	// cout<< "Eta afetr selection is = " << goodAK4Jets.at(i).getEta()<<endl;                                                                             
	// cout<< "Index afetr selection is = " << goodAK4Jets.at(i).getIndex()<<endl;                                                                         
	// cout<< "Hadron flavour afetr selection is = " << goodAK4Jets.at(i).getHadronFlavour()<<endl;                                                        
	//cout<< "Parton flavour afetr selection is = " << goodAK4Jets.at(i).getPartonFlavour()<<endl;                                                         
      }
      // cout <<"size of selected btags from good AK 4 jets is = " << A.size()<<endl;                                                                          

    }
    cout <<"size of selected btags from good AK 4 jets is = " << A.size()<<endl;


    cout << " *************************************"<< endl;
    for (unsigned i=0; i<goodBTaggedAK4Jets.size(); i++){

      cout << i <<"th BTagged AK4Jet CSV is = " << goodBTaggedAK4Jets.at(i).getCSV()<<endl;
      cout << i <<"th BTagged AK4Jet Pt is = " << goodBTaggedAK4Jets.at(i).getPt()<<endl;
      // cout<< "Eta afetr selection is = " << goodBTaggedAK4Jets.at(i).getEta()<<endl;                                                                        
      //cout<< "Index afetr selection is = " << goodBTaggedAK4Jets.at(i).getIndex()<<endl;                                                                     
      //cout<< "Hadron flavour afetr selection is = " << goodBTaggedAK4Jets.at(i).getHadronFlavour()<<endl;                                                    
      //cout<< "Parton flavour afetr selection is = " << goodBTaggedAK4Jets.at(i).getPartonFlavour()<<endl;                                                    
    }
    



                                                               
    cout <<"size of selected btags from good AK 4 jets is ****2 ***** = " << A.size()<<endl;
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {
      if (A.size()>0){
	if (goodAK4Jets.at(i).getMass()>=70 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 300){
	  TLorentzVector zhb;
	  zhb= goodAK4Jets.at(i).getP4();
	  vlq::Candidate zhb1(zhb);
	  ZHb.push_back(zhb1);
	  cout <<"size of selected btags from good AK 4 jets is ****3 ***** = " << A.size()<<endl;
	}
      }
    }
    for (unsigned i=0; i<ZHb.size(); i++) {
      h1_["ZH_mass_b"] -> Fill(ZHb.at(i).getMass(), evtwt) ;
      h1_["ZH_Pt_b"] -> Fill(ZHb.at(i).getPt(), evtwt) ;
    }
    h1_["nZHcandidatejets_b"] -> Fill(ZHb.size(), evtwt) ;

    // if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
    // else return false ;                                   
 
    //ZHcandidate cut(nonboosted)                                                                                                                  
    ZHCandsProducer zh;
    zh.operator()(goodAK4Jets.size(), 2, goodAK4Jets,ZH)  ;

    for (unsigned i=0; i<ZH.size(); i++) {
      h1_["ZH_mass_nb"] -> Fill(ZH.at(i).getMass(), evtwt) ;
      h1_["ZH_Pt_nb"] -> Fill(ZH.at(i).getPt(), evtwt) ;
    }

    h1_["nZHcandidatejets_nb"] -> Fill(ZH.size(), evtwt) ;

    //if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    // else return false ;                                                                                                                        
    // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                             

    double  nZHcandidates=0.0;
    double  nZHcandidates1=0.0;
    if(ZHb.size()>0 ||ZH.size()>0){
      //h1_["cutflow"] -> Fill(16, evtwt) ;
      nZHcandidates = ZHb.size()+ZH.size();

      h1_["nZHcandidatejets"] -> Fill(nZHcandidates, evtwt) ;

    }

    nZHcandidates1 = ZHb.size()+ZH.size();
    h1_["nZHcandidatejets1"] -> Fill(nZHcandidates1, evtwt) ;


    */



    /*
     
    // HCandidate cut(boosted)                                                                                                                 
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {
      
      if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450){
	TLorentzVector h;
	h= goodAK4Jets.at(i).getP4();
	vlq::Candidate h1(h);
	Hb.push_back(h1);
      }
    }
    
    for (unsigned i=0; i<Hb.size(); i++) {
      h1_["H_mass_b"] -> Fill(Hb.at(i).getMass(), evtwt) ;
      h1_["H_Pt_b"] -> Fill(Hb.at(i).getPt(), evtwt) ;
    }
    h1_["nHcandidatejets_b"] -> Fill(Hb.size(), evtwt) ;
    
    // if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
    // else return false ;                                   
    
    //Hcandidate cut(nonboosted)                                                                                                                  
    HCandsProducer h;
    h.operator()(goodAK4Jets.size(), 2, goodAK4Jets,H)  ;
    
    for (unsigned i=0; i<H.size(); i++) {
      h1_["H_mass_nb"] -> Fill(H.at(i).getMass(), evtwt) ;
      h1_["H_Pt_nb"] -> Fill(H.at(i).getPt(), evtwt) ;
    }
    
    h1_["nHcandidatejets_nb"] -> Fill(H.size(), evtwt) ;
    
    //if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    // else return false ;                                                                                                                        
    // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                             
    
    double  nHcandidates=0.0;
    double  nHcandidates1=0.0;
    if(Hb.size()>0 ||H.size()>0){
      //h1_["cutflow"] -> Fill(16, evtwt) ;
      nHcandidates = Hb.size()+H.size();
      
      h1_["nHcandidatejets"] -> Fill(nHcandidates, evtwt) ;
      
    }
    
    nHcandidates1 = Hb.size()+H.size();
    h1_["nHcandidatejets1"] -> Fill(nHcandidates1, evtwt) ;
    

    
    
    
    
    //Zcandidate cut (boosted)
    //cout << " good AK4 jet size 1 = " << goodAK4Jets.size()<<endl;                                                                                                                        
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {                  
      // cout <<i<<"th  AK4 jet mass  mass is " << goodAK4Jets.at(i).getMass()<< endl;                                                                                                   
      if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
	TLorentzVector zb;    
	zb= goodAK4Jets.at(i).getP4();
	vlq::Candidate zb1(zb);       
	ZB.push_back(zb1);
	
      }                                                                                                                                                                                   
    }                                                                                                                                                                                       
    //cout << " good AK4 jet size 2 = " << goodAK4Jets.size()<<endl;                                                                                                                                            
    
    //cout<< " ZB size is " << ZB.size()<<endl;
    for (unsigned i=0; i<ZB.size(); i++) {     
      h1_["Z_mass_a"] -> Fill(ZB.at(i).getMass(), evtwt) ;
      h1_["Z_Pt_a"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
    }                                                                                                                                                                                       
    h1_["nzcandidatejets_a"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
    if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     

    //cout << " good AK4 jet size 3 = " << goodAK4Jets.size()<<endl;
    //Z candidate cut (non boosted)                                                                                                                                                         
    ZCandsProducer z;
    z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;
    
    for (unsigned i=0; i<Z.size(); i++) {
      h1_["Z_mass_b"] -> Fill(Z.at(i).getMass(), evtwt) ;
      h1_["Z_Pt_b"] -> Fill(Z.at(i).getPt(), evtwt) ;
    }
    h1_["nzcandidatejets_b"] -> Fill(Z.size(), evtwt) ;

    if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
 
    double nzcandidates=0.0;
    double nzcandidates1=0.0;
    if ( ZB.size() || Z.size()){

      h1_["cutflow"] -> Fill(12, evtwt);
      nzcandidates = ZB.size()+ Z.size();
      h1_["nzcandidatejets_tot"] -> Fill(nzcandidates, evtwt) ;
      
    }
    nzcandidates1 = ZB.size()+ Z.size();
    h1_["nzcandidatejets1_tot"] -> Fill(nzcandidates1, evtwt) ;
 
    //cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;
    // else return false ;                                                                                                                                                          
    // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                                                                
   
 
    
    // Category D    
 
    //cout << " good AK4 jet size 1 = " << goodAK4Jets.size()<<endl;
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {
      // cout <<i<<"th  AK4 jet mass  mass is " << goodAK4Jets.at(i).getMass()<< endl;
      if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
	TLorentzVector d;
	d= goodAK4Jets.at(i).getP4();
	vlq::Candidate d1(d);
	D.push_back(d1);
      }
    }
 
    //cout<< " D size is " << D.size()<<endl;
     
    for (unsigned i=0; i<D.size(); i++) {
      h1_["top_mass_d"] -> Fill(D.at(i).getMass(), evtwt) ;
      h1_["top_Pt_d"] -> Fill(D.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_d"] -> Fill(D.size(), evtwt) ;
    
    if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
    // else return false ;
    
    //cout << " good AK4 jet size 2 = " << goodAK4Jets.size()<<endl;

 
    TopCandsProducer top,w;
 
    // Category BC
    w.operator()(goodAK4Jets, W,B) ;
    for (unsigned i=0; i<W.size(); i++) {
      //  cout <<i<<"th  W mass is " << W.at(i).getMass()<< endl;
      h1_["W_mass_bc"] -> Fill(W.at(i).getMass(), evtwt) ;
    }
    h1_["nWcandidatejets_bc"] -> Fill(W.size(), evtwt) ;

    for (unsigned i=0; i<B.size(); i++) {
      //cout <<i<<"th  B mass is " << B.at(i).getMass()<< endl;
      h1_["lightjet_mass_bc"] -> Fill(B.at(i).getMass(), evtwt) ; 

    }
    h1_["nlightjetcandidatejets_bc"] -> Fill(B.size(), evtwt) ;

    for (unsigned i=0; i<W.size(); i++) {
      for (unsigned j=0; j<B.size(); j++) {
	//  cout <<"mass W/b =[ " <<W.at(i).getMass()<<","<<B.at(j).getMass() <<"]"<< endl;
	TLorentzVector bc1;
	bc1= W.at(i).getP4()+B.at(j).getP4();
	if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
	  // cout << "top candidate is " << bc1.M() << endl;
	  //cout << "top candidate pt " << bc1.Pt() << endl;  
	  vlq::Candidate bc2(bc1);
	  BC.push_back(bc2); 
       
	}  
      }
      
   
    }
    for (unsigned i=0; i<BC.size(); i++) {
      // cout <<i<<"th  BC mass is " << BC.at(i).getMass()<< endl;
      // cout <<i<<"th  BC pt is " << BC.at(i).getPt()<< endl;
      h1_["top_mass_bc"] -> Fill(BC.at(i).getMass(), evtwt) ;
      h1_["top_Pt_bc"] -> Fill(BC.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_bc"] -> Fill(BC.size(), evtwt) ;
    
    if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}
    //else return false ;
 
 
    //cout << " BC category size is " << BC.size()<<endl;
 
 
    // cout << " W size is " << W.size()<< endl;
    // cout << " B size is " << B.size()<< endl;
 
    //cout << " good AK4 jet size 3 = " << goodAK4Jets.size()<<endl;
 
    // category A
    top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;
    
    for (unsigned i=0; i<tops.size(); i++) {
      h1_["top_mass_a"] -> Fill(tops.at(i).getMass(), evtwt) ;
      h1_["top_Pt_a"] -> Fill(tops.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_a"] -> Fill(tops.size(), evtwt) ;
    
    if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    // else return false ;
    // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;
 
    double  ntopcandidates=0.0;
    double  ntopcandidates1=0.0;
    if(D.size() || BC.size() ||tops.size()){
      h1_["cutflow"] -> Fill(16, evtwt) ;
      ntopcandidates = D.size()+BC.size()+tops.size();
   
      h1_["ntopcandidatejets"] -> Fill(ntopcandidates, evtwt) ;  
   
    } 
    ntopcandidates1 = D.size()+BC.size()+tops.size();
    h1_["ntopcandidatejets1"] -> Fill(ntopcandidates1, evtwt) ;
    //Z and top corelations and ST tempelates
    
    h1_["st_sig"] -> Fill(ST,evtwt);
    h1_["cutflow1"] -> Fill(1, evtwt) ;
    if (goodBTaggedAK4Jets.size() == 1){
      h1_["cutflow1"] -> Fill(2, evtwt) ;   
      h1_["st_sig1b"] -> Fill(ST,evtwt);
    }
 
    if (goodBTaggedAK4Jets.size() >=2){
      h1_["st_sig2b"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(3, evtwt) ;
    }
    
    
    //n,Z,H,B

    //(1)
    if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
      h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(4, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(8, evtwt) ;

	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow2"] -> Fill(1, evtwt) ;
	  h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow2"] -> Fill(2, evtwt) ;
	  h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
	}
	
      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(9, evtwt) ;

        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(3, evtwt) ;
	  h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
        }
	
        else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow2"] -> Fill(4, evtwt) ;
	  h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
        }	
      }  
    }
    //(2)                                                                                                                                                                                        
    if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
      h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(5, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(10, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(5, evtwt) ;
          h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
        }
        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(6, evtwt) ;
          h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
        }

      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(11, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(7, evtwt) ;
          h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
        }

        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(8, evtwt) ;
          h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }

    //(3)                                                                                                                                                                                        
    if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
      h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(6, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(12, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(9, evtwt) ;
          h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
        }
        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(10, evtwt) ;
          h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
        }

      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(13, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(11, evtwt) ;
          h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
        }

        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(12, evtwt) ;
          h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }

    //(4)                                                                                                                                                                                        
    if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
      h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(7, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(14, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(13, evtwt) ;
          h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
        }
        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(14, evtwt) ;
          h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
        }

      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(15, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(15, evtwt) ;
          h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
        }

        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(16, evtwt) ;
          h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }
    
    
    //Top ,ZH
    //(1)
    if (ntopcandidates >=1.0    &&   nZHcandidates>=1.0){
      h1_["st_sigT1ZH1"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(4, evtwt) ;

      if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow3"] -> Fill(8, evtwt) ;
	h1_["st_sigT1ZH1b1"] -> Fill(ST, evtwt) ;
      }
   
  

      else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow3"] -> Fill(9, evtwt) ;
	h1_["st_sigT1ZH1b2"] -> Fill(ST, evtwt) ;
      }
 


    }
   
 
    //(2)                                                                                                                                                       
    if (ntopcandidates == 0.0    &&   nZHcandidates >=1.0){
      h1_["st_sigT0ZH1"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(5, evtwt) ;

      if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow3"] -> Fill(10, evtwt) ;
	h1_["st_sigT0ZH1b1"] -> Fill(ST, evtwt) ;
      }



      else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow3"] -> Fill(11, evtwt) ;
	h1_["st_sigT0ZH1b2"] -> Fill(ST, evtwt) ;
      }



    }



    //(3)                                                                                                                                                     
    if (ntopcandidates >=1.0    &&   nZHcandidates==0.0){
      h1_["st_sigT1ZH0"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(6, evtwt) ;

      if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow3"] -> Fill(12, evtwt) ;
	h1_["st_sigT1ZH0b1"] -> Fill(ST, evtwt) ;
      }



      else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow3"] -> Fill(13, evtwt) ;
	h1_["st_sigT1ZH0b2"] -> Fill(ST, evtwt) ;
      }

    }


    //(4)                                                                                                                                                       
    if (ntopcandidates ==0.0    &&   nZHcandidates==0.0){
      h1_["st_sigT0ZH0"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(7, evtwt) ;

      if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow3"] -> Fill(14, evtwt) ;
	h1_["st_sigT0ZH0b1"] -> Fill(ST, evtwt) ;
      }



      else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow1"] -> Fill(15, evtwt) ;
	h1_["st_sigT0ZH0b2"] -> Fill(ST, evtwt) ;
      }



    }



    

  }

    */





  //===========================================================
  // Control selection done, proceeding with signal selections
  //===========================================================

  //Boosted Z candidates: Z pt > 150 GeV
  CandidateFilter boostedzllfilter(BoostedZCandParams_) ; 
  boostedzllfilter(dileptons, zllBoosted) ;    
  if(zllBoosted.size() > 0) h1_["cutflow"] -> Fill(5, evtwt) ;
  else return false ; 

  // leading jet pt > 100 GeV
  if (goodAK4Jets.at(0).getPt() > 100){h1_["cutflow"] -> Fill(6, evtwt) ; }
  else return false;

  // 2nd laeding jet pt > 50 GeV
  if (goodAK4Jets.at(1).getPt() > 50){h1_["cutflow"] -> Fill(7, evtwt) ; }
  else return false;

  // at least one b-jet 
  if ( goodBTaggedAK4Jets.size() > 0 ) { h1_["cutflow"] -> Fill(8, evtwt) ;}
  else return false;

  // ST > 1000 GeV
  if ( ST > STMin_ ) h1_["cutflow"] -> Fill(9, evtwt) ;  
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
  h1_["htsig1"] -> Fill(htak4.getHT(), evtwt) ;
  h1_["st"] -> Fill(ST, evtwt) ;   
  h1_["stsig1"] -> Fill(ST, evtwt) ;
  h1_["npv_noweight"] -> Fill(*h_npv.product(), *h_evtwtGen.product()); 
  h1_["npv"] -> Fill(*h_npv.product(), evtwt);
  h1_["met"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
  h1_["metsig1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
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

  if ( goodAK4Jets.size()==1){

    h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
  }
  
  
  
  else if ( goodAK4Jets.size()==2){
    
    h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
    h1_["ptak4jet2_sig"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
  }
  
  
  else if ( goodAK4Jets.size()==3){
    
    h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
    h1_["ptak4jet2_sig"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
    h1_["ptak4jet3_sig"]  -> Fill(goodAK4Jets.at(2).getPt(), evtwt) ;
  }
  
  
  
  
  else if ( goodAK4Jets.size()==4){
    h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
    h1_["ptak4jet2_sig"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
    h1_["ptak4jet3_sig"]  -> Fill(goodAK4Jets.at(2).getPt(), evtwt) ;
    h1_["ptak4jet4_sig"]  -> Fill(goodAK4Jets.at(3).getPt(), evtwt) ;
  }
  
  
  else if ( goodAK4Jets.size()==5){
   
    h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
    h1_["ptak4jet2_sig"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
    h1_["ptak4jet3_sig"]  -> Fill(goodAK4Jets.at(2).getPt(), evtwt) ;
    h1_["ptak4jet4_sig"]  -> Fill(goodAK4Jets.at(3).getPt(), evtwt) ;
    h1_["ptak4jet5_sig"]  -> Fill(goodAK4Jets.at(4).getPt(), evtwt) ;
  }
  
  else if ( goodAK4Jets.size()==6){
  
  h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
  h1_["ptak4jet2_sig"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
  h1_["ptak4jet3_sig"]  -> Fill(goodAK4Jets.at(2).getPt(), evtwt) ;
  h1_["ptak4jet4_sig"]  -> Fill(goodAK4Jets.at(3).getPt(), evtwt) ;
  h1_["ptak4jet5_sig"]  -> Fill(goodAK4Jets.at(4).getPt(), evtwt) ;
  h1_["ptak4jet6_sig"]  -> Fill(goodAK4Jets.at(5).getPt(), evtwt) ;
  }
  
  
/*
  if (goodAK4Jets.size()>0)  h1_["ptak4jet1_sig"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
  else if (goodAK4Jets.size()>=1) h1_["ptak4jet2_sig"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
  else if (goodAK4Jets.size()>=2) h1_["ptak4jet3_sig"]  -> Fill(goodAK4Jets.at(2).getPt(), evtwt) ;
  else if (goodAK4Jets.size()>=3) h1_["ptak4jet4_sig"]  -> Fill(goodAK4Jets.at(3).getPt(), evtwt) ;
  else if (goodAK4Jets.size()>=4) h1_["ptak4jet5_sig"]  -> Fill(goodAK4Jets.at(4).getPt(), evtwt) ;
  else if (goodAK4Jets.size()>=5) h1_["ptak4jet6_sig"]  -> Fill(goodAK4Jets.at(5).getPt(), evtwt) ;
  */


  h1_["nak4jet_sig"]  -> Fill(goodAK4Jets.size(), evtwt) ;

  h1_["phi_jet1MET"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

  // fill the b-tagging plots and efficiency maps
  h1_["nbjets"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
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


  /*
  //TprimeTprime categoriztion                                                                                                                                                                     
  //apply a cut of ptZ>150 GeV                                                                                                                                                                     
  for  (unsigned int izll=0; izll<zll.size(); ++izll){
    float zllpt = zll.at(izll).getPt();


    if( zllpt > PtZMin_){

      h1_["nZll"] -> Fill(zll.size(), evtwt) ;
      h1_["st_Cat"] -> Fill(ST, evtwt) ;


      //will only use one btag and two btag requirment and will look at hte ST plots 
      //Cat A                                                                                                                                                       
      if (zll.size()>0  &&  goodBTaggedAK4Jets.size() == 0 ){
        h1_["cutflow1"] -> Fill(1, evtwt) ;
        h1_["st_CatA"] -> Fill(ST, evtwt) ;
        h1_["nak4_A"] -> Fill(goodAK4Jets.size(), evtwt) ;
      }
      //Cat A
      
      if (zll.size()>0  &&  goodBTaggedAK4Jets.size() == 1 ){
        h1_["cutflow1"] -> Fill(2, evtwt) ;
        h1_["st_CatB"] -> Fill(ST, evtwt) ;
        h1_["nak4_B"] -> Fill(goodAK4Jets.size(), evtwt) ;
	}
      //Cat B                                                                                                                                                      
      else if (zll.size()>0  && goodBTaggedAK4Jets.size() >= 2 ){
        h1_["cutflow1"] -> Fill(3, evtwt);
	h1_["st_CatC"] -> Fill(ST, evtwt) ;
        h1_["nak4_C"] -> Fill(goodAK4Jets.size(), evtwt) ;
      }
    }
  }


  */

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
  //Do mass reconstruction TPrime


  //  CandidateFilter boostedzllfilter(BoostedZCandParams_) ;

  TLorentzVector lep1, lep2;
  if (zdecayMode_ == "zelel"){
    lep1 = goodElectrons.at(0).getP4();
    lep2 = goodElectrons.at(1).getP4();
  }
  else if (zdecayMode_ == "zmumu"){
    lep1 = goodMuons.at(0).getP4();
    lep2 = goodMuons.at(1).getP4();
  }
  //  cout << " mass of lepton 1 = " << lep1.M()<<endl;
  // cout << " mass of lepton 2 = " << lep2.M()<<endl;
  // cout << " mass of Z boson " << (lep1+lep2).M()<<endl; 
 if (optimizeReco_ && *h_evttype.product() != "EvtType_Data"){

    GenParticleCollection genPartsInfo;
    genPartsInfo = genpart(evt) ;

    // Declare TLorentzVectors to fill with genParticles                                                                                                                                      
     TLorentzVector tGen, tbarGen, q1, q2;// Z1, Z2;                                                                                    
     TLorentzVector qJet, qbarJet, tJet, tbarJet , WJet , WbarJet, bJet;
     TLorentzVector W1Jet,Wbar1Jet,W2Jet,Wbar2Jet,b1Jet,b2Jet;  
     TLorentzVector had_tjet, lep_tjet, had_tGen, lep_tGen;

    if (bosonMass_ == 91.2){
      q1 = reco.getGen(genPartsInfo, 1, 5, 23);
      q2 = reco.getGen(genPartsInfo, -5, -1, 23);
    }
    else{
      q1 = reco.getGen(genPartsInfo, 1, 5, 25);
      q2 = reco.getGen(genPartsInfo, -5, -1, 25);
    }
    qJet = reco.getMatchedJet(q1, goodAK4Jets, 0.3);
    qbarJet = reco.getMatchedJet(q2, goodAK4Jets, 0.3);
   
    //top reconstruction from W and b
    TLorentzVector q3, q4,q5,q6, b1,b2;
    for (auto& gen1 : genPartsInfo){
      
      
      if (gen1.getPdgID() == 24 && gen1.getMom0PdgID()==6){
	for (auto& gen : genPartsInfo){	
	  if (gen.getPdgID() >=1 && gen.getPdgID() <= 4 && gen.getMom0PdgID() == 24){q3 = gen.getP4();}
	  else  if (gen.getPdgID() >=-4 && gen.getPdgID() <= -1 && gen.getMom0PdgID() == 24){q4 = gen.getP4();}
	  else if (gen.getPdgID() == 5 && gen.getMom0PdgID() == 6){ b1 = gen.getP4();}
	  W1Jet = reco.getMatchedJet(q3, goodAK4Jets, 0.3);
	  Wbar1Jet = reco.getMatchedJet(q4, goodAK4Jets, 0.3);
	  b1Jet = reco.getMatchedJet(b1, goodAK4Jets, 0.3);
	  
	  if (b1Jet.M() == 0 && b1Jet.Pt()==0 && W1Jet.M() == 0 && W1Jet.Pt()==0 && Wbar1Jet.M() == 0 && Wbar1Jet.Pt()==0) {continue;}	
	  // else if (W1Jet.M() == Wbar1Jet.M()) {continue;}
	  tGen = q3+q4+b1;
	  if (W1Jet.M() != Wbar1Jet.M()){tJet =W1Jet+Wbar1Jet+ b1Jet;}
	  if (tGen.M() == 0 && tGen.Pt()==0 && tJet.M() == 0 && tJet.Pt()==0) {continue;}
	  
	}
      }
      
      else if (gen1.getPdgID() == -24 && gen1.getMom0PdgID() == -6){
	for (auto& gen : genPartsInfo){
	  if (gen.getPdgID() >=1 && gen.getPdgID() <= 4 && gen.getMom0PdgID() == -24){ q5 = gen.getP4();}
	  else if (gen.getPdgID() >=-4 && gen.getPdgID() <= -1 && gen.getMom0PdgID() ==- 24){ q6 = gen.getP4();}
	  else if (gen.getPdgID() == -5 && gen.getMom0PdgID() == -6){ b2 = gen.getP4();}
	  W2Jet = reco.getMatchedJet(q5, goodAK4Jets, 0.3);
	  Wbar2Jet = reco.getMatchedJet(q6, goodAK4Jets, 0.3);
	  b2Jet = reco.getMatchedJet(b2, goodAK4Jets, 0.3);
	  
	   if (b2Jet.M() == 0 && b2Jet.Pt()==0 && W2Jet.M() == 0 && W2Jet.Pt()==0 && Wbar2Jet.M() == 0 && Wbar2Jet.Pt()==0) {continue;}
	   // else if (W2Jet.M() == Wbar2Jet.M()) {continue;}
	   tbarGen = q5+q6+b2;
	   if (W2Jet.M() != Wbar2Jet.M()){tbarJet =W2Jet+Wbar2Jet+ b2Jet;}
	   
	   if (tbarGen.M() == 0 && tbarGen.Pt()==0 && tbarJet.M() == 0 && tbarJet.Pt()==0 && tGen.M() == 0 && tGen.Pt()==0 && tJet.M() == 0 && tJet.Pt()==0) {continue;}
	   
	   
	}
      }
    } 
    //Higgs and Z candidates combines , only difference from Z candidate is that we have expand the mass window



    //Higgs candidates
    //case (1)non boosted

    if (q1.M() >0 && q2.M()> 0  && q1 != q2 &&  q1.M() != q2.M() && q1.Pt() > 0 && q2.Pt()>0  && (q1+q2).M()> 0 && qJet.M()> 0 && qJet.Pt()>0 && qbarJet.M()> 0 && qbarJet.Pt()>0 && qJet.M() != qbarJet.M()){ 
	 

    if (bosonMass_ == 91.2){

      double hadgenZ = reco.findInvMass(q1, q2);
      double hadZJet = reco.findInvMass(qJet, qbarJet);
      double hadZJetPt = reco.findPt(qJet, qbarJet);
      
      
      h1_["hadgenZ_nonb"]->Fill(hadgenZ, evtwt);
      h1_["hadZJetMass_nonb"]->Fill(hadZJet,evtwt);
      h1_["hadZJetPt_nonb"]->Fill(hadZJetPt,evtwt);
    }
    
    else {
     
      
      double hadgenH = reco.findInvMass(q1, q2);
      double hadHJet = reco.findInvMass(qJet, qbarJet);
      double hadHJetPt =reco.findPt(qJet, qbarJet);
      h1_["hadgenH_nonb"]->Fill(hadgenH, evtwt);
      h1_["hadHJetMass_nonb"]->Fill(hadHJet,evtwt);
      h1_["hadHJetPt_nonb"]->Fill(hadHJetPt,evtwt);  
    }


 
    
    }

 //case (2)boosted

    else if (q1.M() >0 && q2.M()> 0  && q1.Pt() > 0 && q2.Pt()>0  && (q1+q2).M()> 0 && qJet.M()> 0 && qJet.Pt()>0 && qbarJet.M()> 0 && qbarJet.Pt()>0 && qJet.M() == qbarJet.M()){ 
	 

    if (bosonMass_ == 91.2){

      double hadgenZ = reco.findInvMass(q1, q2);
      double hadZJet = reco.findInvMass(qJet);
      double hadZJetPt = reco.findPt(qJet);
      
      
      h1_["hadgenZ_b"]->Fill(hadgenZ, evtwt);
      h1_["hadZJetMass_b"]->Fill(hadZJet,evtwt);
      h1_["hadZJetPt_b"]->Fill(hadZJetPt,evtwt);
    }
    
    else {
     
      
      double hadgenH = reco.findInvMass(q1, q2);
      double hadHJet = reco.findInvMass(qJet);
      double hadHJetPt =reco.findPt(qJet);
      h1_["hadgenH_b"]->Fill(hadgenH, evtwt);
      h1_["hadHJetMass_b"]->Fill(hadHJet,evtwt);
      h1_["hadHJetPt_b"]->Fill(hadHJetPt,evtwt);  
    }


 
    
    }
   

   

    //top candidates
   
    //  if (q3.M() >0 && q4.M()> 0 && b1.M()>0 && q5.M() >0 && q6.M()> 0 && b2.M()>0 && q3 != q4 && q5 != q6 && q3.M()!= q4.M() && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()>0 &&  q5.M() != q6.M() && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0 && Wbar2Jet.Pt()>0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0 && Wbar1Jet.Pt()>0){
    if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5 != q6 &&  q5.M() != q6.M() && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0 && W2Jet.M() != Wbar2Jet.M()  && W2Jet.M() != b2Jet.M() &&  Wbar2Jet.M() != b2Jet.M()){ 
	 
      /*
	  cout << " mass of gen particle q5     =" << q5.M()<<endl;
          cout << " mass of gen particle q6     =" << q6.M()<<endl;
          cout << " mass of gen particle b2     = " << b2.M()<<endl;

          cout << " mass of  q5 jet (W2jet)     ="  << W2Jet.M()<<endl;
          cout << " mass of q6 jet (Wbar2jet)   ="  << Wbar2Jet.M()<<endl;
          cout << " mass of b2 jet              ="  << b2Jet.M()<<endl;
	 
	  cout << " mass of gen particle tbargen   = " << tbarGen.M()<<endl;
          cout << " mass of tbarJet                = " << tbarJet.M()<<endl;

	  cout<< "   " << endl;
      */	  
	  had_tGen = tbarGen;
	  had_tjet = tbarJet;    
	  /*
	  if (W2Jet.M() == Wbar2Jet.M()){

	    had_tjet =W2Jet+ b2Jet;
	  }
	  else {had_tjet=W2Jet+Wbar2Jet+ b2Jet;}
	  */
	  // cout << " mass of  had_tGen  ************  =" << had_tGen.M()<<endl;
	  // cout << " mass of  had_tjet  ************  =" << had_tjet.M()<<endl;
     
      
      
    // double genZ = reco.findInvMass(q1, q2);
    double hadgenT = reco.findInvMass(q1, q2, had_tGen);
    // double gent = reco.findInvMass(lep_tGen);   
    double hadgent = reco.findInvMass(had_tGen);
    // double ZJet = reco.findInvMass(qJet, qbarJet);
    double hadTJet = reco.findInvMass(qJet, qbarJet, had_tjet);
    // double lepTJet = reco.findInvMass(lep1, lep2, lep_tjet);
    double hadtJet = reco.findInvMass(had_tjet);
    // double leptJet = reco.findInvMass(lep_tjet);

    double hadtJetPt = reco.findPt(had_tjet);
    // double leptJetPt = reco.findPt(lep_tjet);

    if (bosonMass_ == 91.2){

      double hadgenZ = reco.findInvMass(q1, q2);
      double hadZJet = reco.findInvMass(qJet, qbarJet);
      double hadZJetPt = reco.findPt(qJet, qbarJet);
      double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);
      //double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);    
   

      h1_["hadgenZ"]->Fill(hadgenZ, evtwt);
      h1_["hadZJetMass"]->Fill(hadZJet,evtwt);
      h1_["hadZJetPt"]->Fill(hadZJetPt,evtwt);
      h1_["hadTJetPt"]->Fill(hadTJetPt,evtwt);
      // h1_["lepTJetPt"]->Fill(lepTJetPt,evtwt);
    }

    else {

      double hadgenH = reco.findInvMass(q1, q2);
      double hadHJet = reco.findInvMass(qJet, qbarJet);
      double hadHJetPt =reco.findPt(qJet, qbarJet);
      double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);
      //  double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);
      h1_["hadgenH"]->Fill(hadgenH, evtwt);
      h1_["hadHJetMass"]->Fill(hadHJet,evtwt);
      h1_["hadHJetPt"]->Fill(hadHJetPt,evtwt);  
      h1_["hadTJetPt"]->Fill(hadTJetPt,evtwt);
      //   h1_["lepTJetPt"]->Fill(lepTJetPt,evtwt);
  }



    // h1_["genZ"]->Fill(genZ, evtwt);
    h1_["hadgenTMass"] -> Fill(hadgenT, evtwt);
    // h1_["ZJetMass"]->Fill(ZJet,evtwt);
    h1_["hadgentMass"] -> Fill(hadgent, evtwt);   
    h1_["hadTJetMass"] ->Fill(hadTJet, evtwt);
    // h1_["lepTJetMass"] ->Fill(lepTJet, evtwt);
   
    h1_["hadtJetMass"] ->Fill(hadtJet, evtwt);
    // h1_["leptJetMass"] ->Fill(leptJet, evtwt);
    h1_["hadtJetPt"] ->Fill(hadtJetPt, evtwt);
    // h1_["leptJetPt"] ->Fill(leptJetPt, evtwt);

    h2_["genhadtJetMasshad"] ->Fill(hadgent,hadtJet, evtwt);
 
    }
    
    else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3 != q4 &&  q3.M() != q4.M() && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0 && W1Jet.M() != Wbar1Jet.M() && W1Jet.M() != b1Jet.M() &&  Wbar1Jet.M() != b1Jet.M()){
      /*
      cout << " mass of gen particle q3     =" << q3.M()<<endl;
      cout << " mass of gen particle q4     =" << q4.M()<<endl;
      cout << " mass of gen particle b1     = " << b1.M()<<endl;

      cout << " mass of  q3 jet (W1jet)     ="  << W1Jet.M()<<endl;
      cout << " mass of q4 jet (Wbar1jet)   ="  << Wbar1Jet.M()<<endl;
      cout << " mass of b jet               ="  << b1Jet.M()<<endl;

      cout << " mass of gen particle tgen   = " << tGen.M()<<endl;
      cout << " mass of tJet                = " << tJet.M()<<endl;

      cout << "   " <<endl;


      cout<< "   " << endl;
      */
      lep_tGen = tGen;
      lep_tjet = tJet;
      
      //  cout << " mass of  lep_tGen  ************   =" << lep_tGen.M()<<endl;
      // cout << " mass of  lep_tjet  ************  =" << lep_tjet.M()<<endl;
      
      double lepgenT = reco.findInvMass(lep1, lep2, lep_tGen);
      double lepgent = reco.findInvMass(lep_tGen);
      double lepTJet = reco.findInvMass(lep1, lep2, lep_tjet);
      double leptJet = reco.findInvMass(lep_tjet);

      double leptJetPt = reco.findPt(lep_tjet);

      if (bosonMass_ == 91.2){

	double lepZ = reco.findInvMass(lep1, lep2);
	double lepZPt = reco.findPt(lep1, lep2);
	double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);


	h1_["lepZ"]->Fill(lepZ, evtwt);
	h1_["lepZPt"]->Fill(lepZPt,evtwt);
      	h1_["lepTJetPt"]->Fill(lepTJetPt,evtwt);
      }

      else {

	double lepH = reco.findInvMass(lep1, lep2);
	double lepHPt =reco.findPt(lep1, lep2);
	double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);
	h1_["lepH"]->Fill(lepH, evtwt);
	h1_["lepHPt"]->Fill(lepHPt, evtwt);
	h1_["lepTJetPt"]->Fill(lepTJetPt,evtwt);
      }

      
      
      
      h1_["lepgenTMass"] -> Fill(lepgenT, evtwt);
      
      h1_["lepgentMass"] -> Fill(lepgent, evtwt);
      h1_["lepTJetMass"] ->Fill(lepTJet, evtwt);
      
      h1_["leptJetMass"] ->Fill(leptJet, evtwt);
      h1_["leptJetPt"] ->Fill(leptJetPt, evtwt);
      
      h2_["genhadtJetMasslep"] ->Fill(lepgent,leptJet, evtwt);

    }
    

    // W jets ( Category B)

    else if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0  && W2Jet.M() == Wbar2Jet.M()  &&  W2Jet.M() != b2Jet.M() &&  Wbar2Jet.M() != b2Jet.M()){ 
      
      /*
	  cout << " mass of gen particle q5  wwww   =" << q5.M()<<endl;
          cout << " mass of gen particle q6  wwww   =" << q6.M()<<endl;
          cout << " mass of gen particle b2  wwww   = " << b2.M()<<endl;

          cout << " mass of  q5 jet (W2jet)  wwww   ="  << W2Jet.M()<<endl;
          cout << " mass of q6 jet (Wbar2jet) wwww  ="  << Wbar2Jet.M()<<endl;
          cout << " mass of b2 jet            www w ="  << b2Jet.M()<<endl;
	 
	  cout << " mass of gen particle tbargen www  = " << tbarGen.M()<<endl;
          cout << " mass of tbarJet             www   = " << tbarJet.M()<<endl;

	  cout<< "   " << endl;
      */
          tbarJet =W2Jet+ b2Jet;
     	  had_tGen = tbarGen;
	  had_tjet = tbarJet;    
   
	  //  cout << " mass of  had_tGen  ************  WWWW =" << had_tGen.M()<<endl;
	  // cout << " mass of  had_tjet  ************  WWWW =" << had_tjet.M()<<endl;
	  // cout << " mass of  had_Wjet  ************ WWWWw  =" << W2Jet.M()<<endl;
	 
      
      
    double hadgenT = reco.findInvMass(q1, q2, had_tGen);
    double hadgent = reco.findInvMass(had_tGen);
    double hadTJet = reco.findInvMass(qJet, qbarJet, had_tjet);
    double hadtJet = reco.findInvMass(had_tjet);
    double hadWJetmass = W2Jet.M();
    double hadWJetpt = W2Jet.Pt();
    double hadtJetPt = reco.findPt(had_tjet);

    if (bosonMass_ == 91.2){

      double hadgenZ = reco.findInvMass(q1, q2);
      double hadZJet = reco.findInvMass(qJet, qbarJet);
      double hadZJetPt = reco.findPt(qJet, qbarJet);
      double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);
    
      h1_["hadgenZ_B"]->Fill(hadgenZ, evtwt);
      h1_["hadZJetMass_B"]->Fill(hadZJet,evtwt);
      h1_["hadZJetPt_B"]->Fill(hadZJetPt,evtwt);
      h1_["hadTJetPt_B"]->Fill(hadTJetPt,evtwt);
    }

    else {

      double hadgenH = reco.findInvMass(q1, q2);
      double hadHJet = reco.findInvMass(qJet, qbarJet);
      double hadHJetPt =reco.findPt(qJet, qbarJet);
      double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);
      h1_["hadgenH_B"]->Fill(hadgenH, evtwt);
      h1_["hadHJetMass_B"]->Fill(hadHJet,evtwt);
      h1_["hadHJetPt_B"]->Fill(hadHJetPt,evtwt);  
      h1_["hadTJetPt_B"]->Fill(hadTJetPt,evtwt);
  }



    h1_["hadgenTMass_B"] -> Fill(hadgenT, evtwt);
    h1_["hadgentMass_B"] -> Fill(hadgent, evtwt);   
    h1_["hadTJetMass_B"] ->Fill(hadTJet, evtwt);
   
    h1_["hadtJetMass_B"] ->Fill(hadtJet, evtwt);
    h1_["hadWJetMass_B"] ->Fill(hadWJetmass, evtwt);
    h1_["hadWJetPt_B"] ->Fill(hadWJetpt, evtwt);
    h1_["hadtJetPt_B"] ->Fill(hadtJetPt, evtwt);
    h2_["genhadtJetMasshad_B"] ->Fill(hadgent,hadtJet, evtwt);
 
    }
    
    else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0 && W1Jet.M() == Wbar1Jet.M()  &&  W1Jet.M() != b1Jet.M() &&  Wbar1Jet.M() != b1Jet.M()){
      /*
      cout << " mass of gen particle q3 WWW    =" << q3.M()<<endl;
      cout << " mass of gen particle q4 WWW    =" << q4.M()<<endl;
      cout << " mass of gen particle b1 WWW    = " << b1.M()<<endl;

      cout << " mass of  q3 jet (W1jet)   WWW  ="  << W1Jet.M()<<endl;
      cout << " mass of q4 jet (Wbar1jet) WWW  ="  << Wbar1Jet.M()<<endl;
      cout << " mass of b jet             www  ="  << b1Jet.M()<<endl;

      cout << " mass of gen particle tgen WWW  = " << tGen.M()<<endl;
      cout << " mass of tJet           www     = " << tJet.M()<<endl;

      cout << "   " <<endl;

      
      cout<< "   " << endl;
      */
      // tGen = q3+b1;
      
      tbarJet =W1Jet+ b1Jet;     
      lep_tGen = tGen;
      lep_tjet = tJet;
      
      // cout << " mass of  lep_tGen  ************ wwww  =" << lep_tGen.M()<<endl;
      // cout << " mass of  lep_tjet  ************ WWWWw  =" << lep_tjet.M()<<endl;
      // cout << " mass of  lep_Wjet  ************ WWWWw  =" << W1Jet.M()<<endl;
      
      double lepgenT = reco.findInvMass(lep1, lep2, lep_tGen);
      double lepgent = reco.findInvMass(lep_tGen);
      double lepTJet = reco.findInvMass(lep1, lep2, lep_tjet);
      double leptJet = reco.findInvMass(lep_tjet);
      double lepWJetmass = W1Jet.M();
      double lepWJetpt = W1Jet.Pt();
      double leptJetPt = reco.findPt(lep_tjet);
       
      if (bosonMass_ == 91.2){
	
	double lepZ = reco.findInvMass(lep1, lep2);
	double lepZPt = reco.findPt(lep1, lep2);
	double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);
	
	
	h1_["lepZ_B"]->Fill(lepZ, evtwt);
	h1_["lepZPt_B"]->Fill(lepZPt,evtwt);
      	h1_["lepTJetPt_B"]->Fill(lepTJetPt,evtwt);
      }

      else {
	
	double lepH = reco.findInvMass(lep1, lep2);
	double lepHPt =reco.findPt(lep1, lep2);
	double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);
	h1_["lepH_B"]->Fill(lepH, evtwt);
	h1_["lepHPt_B"]->Fill(lepHPt, evtwt);
	h1_["lepTJetPt_B"]->Fill(lepTJetPt,evtwt);
      }
      
      
      
      
      h1_["lepgenTMass_B"] -> Fill(lepgenT, evtwt);
      h1_["lepgentMass_B"] -> Fill(lepgent, evtwt);
      h1_["lepTJetMass_B"] ->Fill(lepTJet, evtwt);
      h1_["lepWJetMass_B"] ->Fill(lepWJetmass, evtwt);
      h1_["lepWJetPt_B"] ->Fill(lepWJetpt, evtwt);
      h1_["leptJetMass_B"] ->Fill(leptJet, evtwt);
      h1_["leptJetPt_B"] ->Fill(leptJetPt, evtwt);
      h2_["genhadtJetMasslep_B"] ->Fill(lepgent,leptJet, evtwt);
      
    }
    // Category C ( one of q jet from W overlaps with b jet, no W jet category here)

    else if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0  &&  W2Jet.M() != Wbar2Jet.M() && ( W2Jet.M() == b2Jet.M() ||  Wbar2Jet.M() == b2Jet.M())){

      tbarJet =W2Jet+ Wbar2Jet; 
      /*
	cout << " mass of gen particle q5  (cat C)   =" << q5.M()<<endl;
	cout << " mass of gen particle q6  (cat C)   =" << q6.M()<<endl;
	cout << " mass of gen particle b2  (cat C)   = " << b2.M()<<endl;

	cout << " mass of  q5 jet (W2jet)  (cat C)   ="  << W2Jet.M()<<endl;
	cout << " mass of q6 jet (Wbar2jet) (cat C)  ="  << Wbar2Jet.M()<<endl;
	cout << " mass of b2 jet            (cat C) ="  << b2Jet.M()<<endl;

	cout << " mass of gen particle tbargen (cat C)  = " << tbarGen.M()<<endl;
	cout << " mass of tbarJet             (cat C)   = " << tbarJet.M()<<endl;

	cout<< "   " << endl;
      */
	//	tbarJet =W2Jet+ Wbar2Jet;
	//	if ( W2Jet.M() == b2Jet.M() ){ tbarJet =W2Jet+ Wbar2Jet;}
	//	else if ( Wbar2Jet.M() == b2Jet.M()){ tbarJet =W2Jet+ Wbar2Jet;}
	had_tGen = tbarGen;
	had_tjet = tbarJet;


	//	cout << " mass of  had_tGen  ************  (cat C) =" << had_tGen.M()<<endl;
	//	cout << " mass of  had_tjet  ************  (cat C) =" << had_tjet.M()<<endl;



	double hadgenT = reco.findInvMass(q1, q2, had_tGen);
	double hadgent = reco.findInvMass(had_tGen);
	double hadTJet = reco.findInvMass(qJet, qbarJet, had_tjet);
	double hadtJet = reco.findInvMass(had_tjet);
	double hadtJetPt = reco.findPt(had_tjet);
	double hadmergedjetmass = b2Jet.M();
	double hadmergedjetpt = b2Jet.Pt();
	if (bosonMass_ == 91.2){

	  double hadgenZ = reco.findInvMass(q1, q2);
	  double hadZJet = reco.findInvMass(qJet, qbarJet);
	  double hadZJetPt = reco.findPt(qJet, qbarJet);
	  double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);

	  h1_["hadgenZ_C"]->Fill(hadgenZ, evtwt);
	  h1_["hadZJetMass_C"]->Fill(hadZJet,evtwt);
	  h1_["hadZJetPt_C"]->Fill(hadZJetPt,evtwt);
	  h1_["hadTJetPt_C"]->Fill(hadTJetPt,evtwt);
	}

	else {

	  double hadgenH = reco.findInvMass(q1, q2);
	  double hadHJet = reco.findInvMass(qJet, qbarJet);
	  double hadHJetPt =reco.findPt(qJet, qbarJet);
	  double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);
	  h1_["hadgenH_C"]->Fill(hadgenH, evtwt);
	  h1_["hadHJetMass_C"]->Fill(hadHJet,evtwt);
	  h1_["hadHJetPt_C"]->Fill(hadHJetPt,evtwt);
	  h1_["hadTJetPt_C"]->Fill(hadTJetPt,evtwt);
	}

	h1_["hadgenTMass_C"] -> Fill(hadgenT, evtwt);
	h1_["hadgentMass_C"] -> Fill(hadgent, evtwt);
	h1_["hadTJetMass_C"] ->Fill(hadTJet, evtwt);

	h1_["hadtJetMass_C"] ->Fill(hadtJet, evtwt);
	h1_["hadtJetPt_C"] ->Fill(hadtJetPt, evtwt);
	h1_["hadmergedJetMass_C"] ->Fill(hadmergedjetmass, evtwt);
	h1_["hadmergedJetPt_C"] ->Fill(hadmergedjetpt , evtwt);

	h2_["genhadtJetMasshad_C"] ->Fill(hadgent,hadtJet, evtwt);

    }
 

    else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0  &&  W1Jet.M()!= Wbar1Jet.M() && ( W1Jet.M() == b1Jet.M() ||  Wbar1Jet.M() == b1Jet.M())){


      tJet =W1Jet+ Wbar1Jet; 
      /*
      cout << " mass of gen particle q3 (cat C)    =" << q3.M()<<endl;
      cout << " mass of gen particle q4 (cat C)    =" << q4.M()<<endl;
      cout << " mass of gen particle b1 (cat C)    = " << b1.M()<<endl;

      cout << " mass of  q3 jet (W1jet)   (cat C)  ="  << W1Jet.M()<<endl;
      cout << " mass of q4 jet (Wbar1jet) (cat C)  ="  << Wbar1Jet.M()<<endl;
      cout << " mass of b jet             (cat C)  ="  << b1Jet.M()<<endl;

      cout << " mass of gen particle tgen (cat C)  = " << tGen.M()<<endl;
      cout << " mass of tJet           (cat C)     = " << tJet.M()<<endl;

      cout << "   " <<endl;


      cout<< "   " << endl;
      */
      // tJet =W1Jet+ Wbar1Jet;
      lep_tGen = tGen;
      lep_tjet = tJet;

      // cout << " mass of  lep_tGen  ************ (cat C)  =" << lep_tGen.M()<<endl;
      //cout << " mass of  lep_tjet  ************ (cat C)  =" << lep_tjet.M()<<endl;
      // cout << " mass of  lep_Wjet  ************ (cat C)  =" << W1Jet.M()<<endl;

      double lepgenT = reco.findInvMass(lep1, lep2, lep_tGen);
      double lepgent = reco.findInvMass(lep_tGen);
      double lepTJet = reco.findInvMass(lep1, lep2, lep_tjet);
      double leptJet = reco.findInvMass(lep_tjet);
      double leptJetPt = reco.findPt(lep_tjet);
      double lepmergedjetmass = b1Jet.M();
      double lepmergedjetpt = b1Jet.Pt();

      if (bosonMass_ == 91.2){

        double lepZ = reco.findInvMass(lep1, lep2);
        double lepZPt = reco.findPt(lep1, lep2);
        double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);


        h1_["lepZ_C"]->Fill(lepZ, evtwt);
        h1_["lepZPt_C"]->Fill(lepZPt,evtwt);
        h1_["lepTJetPt_C"]->Fill(lepTJetPt,evtwt);
      }

      else {

        double lepH = reco.findInvMass(lep1, lep2);
        double lepHPt =reco.findPt(lep1, lep2);
        double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);
        h1_["lepH_C"]->Fill(lepH, evtwt);
        h1_["lepHPt_C"]->Fill(lepHPt, evtwt);
        h1_["lepTJetPt_C"]->Fill(lepTJetPt,evtwt);
      }

      h1_["lepgenTMass_C"] -> Fill(lepgenT, evtwt);
      h1_["lepgentMass_C"] -> Fill(lepgent, evtwt);
      h1_["lepTJetMass_C"] ->Fill(lepTJet, evtwt);
      h1_["leptJetMass_C"] ->Fill(leptJet, evtwt);
      h1_["leptJetPt_C"] ->Fill(leptJetPt, evtwt);
     
      h1_["lepmergedJetMass_C"] ->Fill(lepmergedjetmass, evtwt);
      h1_["lepmergedJetPt_C"] ->Fill(lepmergedjetpt , evtwt);


      h2_["genhadtJetMasslep_C"] ->Fill(lepgent,leptJet, evtwt);

    }
 

    // Category D ( All three jet masses are equal , i.e. fully merged case)

    else if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0  &&  W2Jet.M() == Wbar2Jet.M() &&  W2Jet.M() == b2Jet.M() &&  Wbar2Jet.M() == b2Jet.M()){
    
      tbarJet =W2Jet;
      /*
      cout << " mass of gen particle q5  (cat D)   =" << q5.M()<<endl;
      cout << " mass of gen particle q6  (cat D)   =" << q6.M()<<endl;
      cout << " mass of gen particle b2  (cat D)   = " << b2.M()<<endl;

      cout << " mass of  q5 jet (W2jet)  (cat D)   ="  << W2Jet.M()<<endl;
      cout << " mass of q6 jet (Wbar2jet) (cat D)  ="  << Wbar2Jet.M()<<endl;
      cout << " mass of b2 jet            (cat D) ="  << b2Jet.M()<<endl;

      cout << " mass of gen particle tbargen (cat D)  = " << tbarGen.M()<<endl;
      cout << " mass of tbarJet             (cat D)   = " << tbarJet.M()<<endl;

      cout<< "   " << endl;
      */
      had_tGen = tbarGen;
      had_tjet = tbarJet;


      //cout << " mass of  had_tGen  ************  (cat D) =" << had_tGen.M()<<endl;
      //cout << " mass of  had_tjet  ************  (cat D) =" << had_tjet.M()<<endl;



      double hadgenT = reco.findInvMass(q1, q2, had_tGen);
      double hadgent = reco.findInvMass(had_tGen);
      double hadTJet = reco.findInvMass(qJet, qbarJet, had_tjet);
      double hadtJet = reco.findInvMass(had_tjet);
      double hadtJetPt = reco.findPt(had_tjet);

      if (bosonMass_ == 91.2){

	double hadgenZ = reco.findInvMass(q1, q2);
	double hadZJet = reco.findInvMass(qJet, qbarJet);
	double hadZJetPt = reco.findPt(qJet, qbarJet);
	double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);

	h1_["hadgenZ_D"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass_D"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt_D"]->Fill(hadZJetPt,evtwt);
	h1_["hadTJetPt_D"]->Fill(hadTJetPt,evtwt);
      }

      else {

	double hadgenH = reco.findInvMass(q1, q2);
	double hadHJet = reco.findInvMass(qJet, qbarJet);
	double hadHJetPt =reco.findPt(qJet, qbarJet);
	double hadTJetPt = reco.findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenH_D"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass_D"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt_D"]->Fill(hadHJetPt,evtwt);
	h1_["hadTJetPt_D"]->Fill(hadTJetPt,evtwt);
      }

      h1_["hadgenTMass_D"] -> Fill(hadgenT, evtwt);
      h1_["hadgentMass_D"] -> Fill(hadgent, evtwt);
      h1_["hadTJetMass_D"] ->Fill(hadTJet, evtwt);

      h1_["hadtJetMass_D"] ->Fill(hadtJet, evtwt);
      h1_["hadtJetPt_D"] ->Fill(hadtJetPt, evtwt);
      h2_["genhadtJetMasshad_D"] ->Fill(hadgent,hadtJet, evtwt);

    }
 
 
    else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0  &&  W1Jet.M()== Wbar1Jet.M() &&  W1Jet.M() == b1Jet.M() &&  Wbar1Jet.M() == b1Jet.M()){
         
      tJet =W1Jet;
      /*
      cout << " mass of gen particle q3 (cat D)    =" << q3.M()<<endl;
      cout << " mass of gen particle q4 (cat D)    =" << q4.M()<<endl;
      cout << " mass of gen particle b1 (cat D)    = " << b1.M()<<endl;

      cout << " mass of  q3 jet (W1jet)   (cat D)  ="  << W1Jet.M()<<endl;
      cout << " mass of q4 jet (Wbar1jet) (cat D)  ="  << Wbar1Jet.M()<<endl;
      cout << " mass of b jet             (cat D)  ="  << b1Jet.M()<<endl;

      cout << " mass of gen particle tgen (cat D)  = " << tGen.M()<<endl;
      cout << " mass of tJet           (cat D)     = " << tJet.M()<<endl;

      cout << "   " <<endl;


      cout<< "   " << endl;
      */
      // tJet =W1Jet+ Wbar1Jet;                                                                                                              
      lep_tGen = tGen;
      lep_tjet = tJet;
      // cout << " mass of  lep_tGen  ************ (cat D)  =" << lep_tGen.M()<<endl;
      // cout << " mass of  lep_tjet  ************ (cat D)  =" << lep_tjet.M()<<endl;
      // cout << " mass of  lep_Wjet  ************ (cat D)  =" << W1Jet.M()<<endl;                                                           

      double lepgenT = reco.findInvMass(lep1, lep2, lep_tGen);
      double lepgent = reco.findInvMass(lep_tGen);
      double lepTJet = reco.findInvMass(lep1, lep2, lep_tjet);
      double leptJet = reco.findInvMass(lep_tjet);
      double leptJetPt = reco.findPt(lep_tjet);


      if (bosonMass_ == 91.2){

        double lepZ = reco.findInvMass(lep1, lep2);
        double lepZPt = reco.findPt(lep1, lep2);
        double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);

        h1_["lepZ_D"]->Fill(lepZ, evtwt);
        h1_["lepZPt_D"]->Fill(lepZPt,evtwt);
        h1_["lepTJetPt_D"]->Fill(lepTJetPt,evtwt);
      }
   
      else {

        double lepH = reco.findInvMass(lep1, lep2);
        double lepHPt =reco.findPt(lep1, lep2);
        double lepTJetPt = reco.findPt(lep1, lep2, lep_tjet);
        h1_["lepH_D"]->Fill(lepH, evtwt);
        h1_["lepHPt_D"]->Fill(lepHPt, evtwt);
        h1_["lepTJetPt_D"]->Fill(lepTJetPt,evtwt);
      }
    
                                                                                                                                                
      h1_["lepgenTMass_D"] -> Fill(lepgenT, evtwt);
      h1_["lepgentMass_D"] -> Fill(lepgent, evtwt);
      h1_["lepTJetMass_D"] ->Fill(lepTJet, evtwt);
      h1_["leptJetMass_D"] ->Fill(leptJet, evtwt);
      h1_["leptJetPt_D"] ->Fill(leptJetPt, evtwt);
      h2_["genhadtJetMasslep_D"] ->Fill(lepgent,leptJet, evtwt);
    }


    else if (q1.M() >0 && q2.M()> 0 && q1.Pt() > 0 && q2.Pt()>0 && (q1+q2).M()> 0  && qJet.M()==qbarJet.M()){
      //  cout<< "q1 jet m is " << qJet.M()<<endl;
      // cout <<"q2 jet m is " << qbarJet.M()<<endl;
      if (bosonMass_ == 91.2){

        double hadgenZ = reco.findInvMass(q1, q2);
        double hadZJet = reco.findInvMass(qJet);
        double hadZJetPt = reco.findPt(qJet);
        double hadTJetPt = reco.findPt(qJet, had_tjet);
     
        h1_["hadgenZ_E"]->Fill(hadgenZ, evtwt);
        h1_["hadZJetMass_E"]->Fill(hadZJet,evtwt);
        h1_["hadZJetPt_E"]->Fill(hadZJetPt,evtwt);
        h1_["hadTJetPt_E"]->Fill(hadTJetPt,evtwt);
      }



    }






 }

 
 //ZH combined category, only difference from Z is that an extended mass window
 /*
 //start ( choosing Btaged jets from goodAK4Jets using CSV>0.8 and Pt>50
 cout << " size of good AK4 jet collection = " << goodAK4Jets.size()<<endl;                                                                                 
 cout << " size of BTagged AK4 jet collection = "<< goodBTaggedAK4Jets.size()<<endl;                                                                       
                                                                                                                                                            
                                                                                                                                                            
  vector<TLorentzVector> A;                                                                                                                                 
  for (unsigned i=0; i<goodAK4Jets.size(); i++){                                                                                                            
                                                                                                                                                            
  //cout << i <<"th good AK4Jet CSV is = " << goodAK4Jets.at(i).getCSV()<<endl;                                                                            
  // cout << i <<"th good AK4Jet Pt is = " << goodAK4Jets.at(i).getPt()<<endl;                                                                            
                                                                                                                                                            
    //vector<TLorentzVector> A;                                                                                                                             
    TLorentzVector jetinfo = goodAK4Jets.at(i).getP4();                                                                                                     
    if (goodAK4Jets.at(i).getCSV()>0.800 && goodAK4Jets.at(i).getPt()>50){                                                                                  
      A.push_back(jetinfo);                                                                                                                                 
      cout<< "CSV afetr selection is = " << goodAK4Jets.at(i).getCSV()<<endl;                                                                               
      cout<< "Pt afetr selection is = " << goodAK4Jets.at(i).getPt()<<endl;                                                                               
    }                                                                                                                                                       
                                                                                                                                                            
  }                                                                                                                                                         
  cout <<"size of selected btags from good AK 4 jets is = " << A.size()<<endl;                                                                         
                                                                                                                                                       
                                                                                                                                                          
  cout << " *************************************"<< endl;                                                                                       
  for (unsigned i=0; i<goodBTaggedAK4Jets.size(); i++){                                                                     
  cout << i <<"th BTagged AK4Jet CSV is = " << goodBTaggedAK4Jets.at(i).getCSV()<<endl;                                                                    
  cout << i <<"th BTagged AK4Jet Pt is = " << goodBTaggedAK4Jets.at(i).getPt()<<endl;    

  }
 
 // ZHCandidate cut(boosted)                                                                                                                 
  
 for (unsigned i=0; i<goodAK4Jets.size(); i++) {
   //cout <<"size of selected btags from good AK 4 jets is ****2**** = " << A.size()<<endl;
   if (goodAK4Jets.at(i).getMass()>=70 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 300){
     TLorentzVector zhb;
     zhb= goodAK4Jets.at(i).getP4();
     vlq::Candidate zhb1(zhb);
     ZHb.push_back(zhb1);
     //  cout <<"size of selected btags from good AK 4 jets is ****3**** = " << A.size()<<endl;
   }
 }

 //cout <<"size of selected btags from good AK 4 jets is ****4**** = " << A.size()<<endl;
 for (unsigned i=0; i<ZHb.size(); i++) {
   h1_["ZH_mass_b"] -> Fill(ZHb.at(i).getMass(), evtwt) ;
   h1_["ZH_Pt_b"] -> Fill(ZHb.at(i).getPt(), evtwt) ;
 }
 h1_["nZHcandidatejets_b"] -> Fill(ZHb.size(), evtwt) ;

 // if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
 // else return false ;                                   
 
 //ZHcandidate cut(nonboosted)                                                                                                                  
 ZHCandsProducer zh;
 zh.operator()(goodAK4Jets.size(), 2, goodAK4Jets,ZH)  ;

 for (unsigned i=0; i<ZH.size(); i++) {
   h1_["ZH_mass_nb"] -> Fill(ZH.at(i).getMass(), evtwt) ;
   h1_["ZH_Pt_nb"] -> Fill(ZH.at(i).getPt(), evtwt) ;
 }

 h1_["nZHcandidatejets_nb"] -> Fill(ZH.size(), evtwt) ;

 //if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
 // else return false ;                                                                                                                        
 // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                             

 double  nZHcandidates=0.0;
 double  nZHcandidates1=0.0;
 if(ZHb.size()>0 ||ZH.size()>0){
   //h1_["cutflow"] -> Fill(16, evtwt) ;
   nZHcandidates = ZHb.size()+ZH.size();

   h1_["nZHcandidatejets"] -> Fill(nZHcandidates, evtwt) ;

 }

 nZHcandidates1 = ZHb.size()+ZH.size();
 h1_["nZHcandidatejets1"] -> Fill(nZHcandidates1, evtwt) ;


 
 */

 
 
  
    
 // HCandidate cut(boosted)                                                                                                                 
 for (unsigned i=0; i<goodAK4Jets.size(); i++) {
   
   if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450 && goodAK4Jets.at(i).getCSV()>0.800){
     TLorentzVector h;
     h= goodAK4Jets.at(i).getP4();
     vlq::Candidate h1(h);
     Hb.push_back(h1);
   }
 }
 
 for (unsigned i=0; i<Hb.size(); i++) {
   h1_["H_mass_b"] -> Fill(Hb.at(i).getMass(), evtwt) ;
   h1_["H_Pt_b"] -> Fill(Hb.at(i).getPt(), evtwt) ;
 }
 h1_["nHcandidatejets_b"] -> Fill(Hb.size(), evtwt) ;

 // if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
 // else return false ;                                   
 
 //Hcandidate cut(nonboosted)                                                                                                                  
 HCandsProducer h;
 h.operator()(goodAK4Jets.size(), 2, goodAK4Jets,H)  ;

 for (unsigned i=0; i<H.size(); i++) {
   h1_["H_mass_nb"] -> Fill(H.at(i).getMass(), evtwt) ;
   h1_["H_Pt_nb"] -> Fill(H.at(i).getPt(), evtwt) ;
 }

 h1_["nHcandidatejets_nb"] -> Fill(H.size(), evtwt) ;

 //if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
 // else return false ;                                                                                                                        
 // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                             

 double  nHcandidates=0.0;
 double  nHcandidates1=0.0;
 if(Hb.size()>0 ||H.size()>0){
   //h1_["cutflow"] -> Fill(16, evtwt) ;
   nHcandidates = Hb.size()+H.size();

   h1_["nHcandidatejets"] -> Fill(nHcandidates, evtwt) ;

 }

 nHcandidates1 = Hb.size()+H.size();
 h1_["nHcandidatejets1"] -> Fill(nHcandidates1, evtwt) ;


 

 //Zcandidate cut (boosted)
 //cout << " good AK4 jet size 1 = " << goodAK4Jets.size()<<endl;                                                                                                                        
 for (unsigned i=0; i<goodAK4Jets.size(); i++) {                  
   // cout <<i<<"th  AK4 jet mass  mass is " << goodAK4Jets.at(i).getMass()<< endl;                                                                                                   
   if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
     TLorentzVector zb;    
     zb= goodAK4Jets.at(i).getP4();
     vlq::Candidate zb1(zb);       
     ZB.push_back(zb1);
     
   }                                                                                                                                                                                   
 }                                                                                                                                                                                       
 //cout << " good AK4 jet size 2 = " << goodAK4Jets.size()<<endl;                                                                                                                                            
                                   
 //cout<< " ZB size is " << ZB.size()<<endl;
 for (unsigned i=0; i<ZB.size(); i++) {     
   h1_["Z_mass_a"] -> Fill(ZB.at(i).getMass(), evtwt) ;
   h1_["Z_Pt_a"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
 }                                                                                                                                                                                       
 h1_["nzcandidatejets_a"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     

//cout << " good AK4 jet size 3 = " << goodAK4Jets.size()<<endl;
//Z candidate cut (non boosted)                                                                                                                                                         
 ZCandsProducer z;
 z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;

 for (unsigned i=0; i<Z.size(); i++) {
   h1_["Z_mass_b"] -> Fill(Z.at(i).getMass(), evtwt) ;
   h1_["Z_Pt_b"] -> Fill(Z.at(i).getPt(), evtwt) ;
 }
 h1_["nzcandidatejets_b"] -> Fill(Z.size(), evtwt) ;

 if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
 
 double nzcandidates=0.0;
 double nzcandidates1=0.0; 
 if ( ZB.size()>0 || Z.size()>0){

  h1_["cutflow"] -> Fill(12, evtwt);
   nzcandidates = ZB.size()+ Z.size();
  h1_["nzcandidatejets_tot"] -> Fill(nzcandidates, evtwt) ;

 }
 nzcandidates1 = ZB.size()+ Z.size();                                                                                                       
 h1_["nzcandidatejets1_tot"] -> Fill(nzcandidates1, evtwt) ;   
 
 //cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;
 // else return false ;                                                                                                                                                          
 // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;                                                                                                                

 
 
 // Category D    
 
 //cout << " good AK4 jet size 1 = " << goodAK4Jets.size()<<endl;
 for (unsigned i=0; i<goodAK4Jets.size(); i++) {
   // cout <<i<<"th  AK4 jet mass  mass is " << goodAK4Jets.at(i).getMass()<< endl;
   if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
     TLorentzVector d;
     d= goodAK4Jets.at(i).getP4();
     vlq::Candidate d1(d);
     D.push_back(d1);
     }
 }
 
 //cout<< " D size is " << D.size()<<endl;
     
 for (unsigned i=0; i<D.size(); i++) {
   h1_["top_mass_d"] -> Fill(D.at(i).getMass(), evtwt) ;
   h1_["top_Pt_d"] -> Fill(D.at(i).getPt(), evtwt) ;
 }
 h1_["ntopcandidatejets_d"] -> Fill(D.size(), evtwt) ;
 
 if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
 // else return false ;
 
 //cout << " good AK4 jet size 2 = " << goodAK4Jets.size()<<endl;

 
 TopCandsProducer top,w;
 
 // Category BC
 w.operator()(goodAK4Jets, W,B) ;
 for (unsigned i=0; i<W.size(); i++) {
   //  cout <<i<<"th  W mass is " << W.at(i).getMass()<< endl;
   h1_["W_mass_bc"] -> Fill(W.at(i).getMass(), evtwt) ;
 }
 h1_["nWcandidatejets_bc"] -> Fill(W.size(), evtwt) ;

 for (unsigned i=0; i<B.size(); i++) {
   //cout <<i<<"th  B mass is " << B.at(i).getMass()<< endl;
   h1_["lightjet_mass_bc"] -> Fill(B.at(i).getMass(), evtwt) ; 

}
 h1_["nlightjetcandidatejets_bc"] -> Fill(B.size(), evtwt) ;

 for (unsigned i=0; i<W.size(); i++) {
   for (unsigned j=0; j<B.size(); j++) {
     //  cout <<"mass W/b =[ " <<W.at(i).getMass()<<","<<B.at(j).getMass() <<"]"<< endl;
     TLorentzVector bc1;
     bc1= W.at(i).getP4()+B.at(j).getP4();
     if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
       // cout << "top candidate is " << bc1.M() << endl;
       //cout << "top candidate pt " << bc1.Pt() << endl;  
       vlq::Candidate bc2(bc1);
       BC.push_back(bc2); 
       
     }  
   }
   
   
 }
 for (unsigned i=0; i<BC.size(); i++) {
   // cout <<i<<"th  BC mass is " << BC.at(i).getMass()<< endl;
   // cout <<i<<"th  BC pt is " << BC.at(i).getPt()<< endl;
   h1_["top_mass_bc"] -> Fill(BC.at(i).getMass(), evtwt) ;
   h1_["top_Pt_bc"] -> Fill(BC.at(i).getPt(), evtwt) ;
 }
 h1_["ntopcandidatejets_bc"] -> Fill(BC.size(), evtwt) ;
 
if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}
 //else return false ;
 
 
 //cout << " BC category size is " << BC.size()<<endl;
 
 
 // cout << " W size is " << W.size()<< endl;
 // cout << " B size is " << B.size()<< endl;
 
 //cout << " good AK4 jet size 3 = " << goodAK4Jets.size()<<endl;
 
 // category A
 top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;
 
 for (unsigned i=0; i<tops.size(); i++) {
   h1_["top_mass_a"] -> Fill(tops.at(i).getMass(), evtwt) ;
   h1_["top_Pt_a"] -> Fill(tops.at(i).getPt(), evtwt) ;
 }
 h1_["ntopcandidatejets_a"] -> Fill(tops.size(), evtwt) ;
 
 if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
 // else return false ;
 // cout << " good AK4 jet size 4 = " << goodAK4Jets.size()<<endl;
 
 double  ntopcandidates=0.0;
 double  ntopcandidates1=0.0;
 if(D.size()>0 || BC.size()>0 ||tops.size()>0){
   h1_["cutflow"] -> Fill(16, evtwt) ;
   ntopcandidates = D.size()+BC.size()+tops.size();
   
   h1_["ntopcandidatejets"] -> Fill(ntopcandidates, evtwt) ;  
   
 } 

 ntopcandidates1 = D.size()+BC.size()+tops.size();                                                                                          
 h1_["ntopcandidatejets1"] -> Fill(ntopcandidates1, evtwt) ; 
 //Z and top corelations and ST tempelates

 h1_["st_sig"] -> Fill(ST,evtwt);
 h1_["cutflow1"] -> Fill(1, evtwt) ;
 if (goodBTaggedAK4Jets.size() == 1){
   h1_["cutflow1"] -> Fill(2, evtwt) ;   
   h1_["st_sig1b"] -> Fill(ST,evtwt);
 }
 
 if (goodBTaggedAK4Jets.size() >=2){
   h1_["st_sig2b"] -> Fill(ST,evtwt);
   h1_["cutflow1"] -> Fill(3, evtwt) ;
 }

 /*
 //only top,Z,b
 
 //(1)
 if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
   h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
   h1_["cutflow1"] -> Fill(4, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow1"] -> Fill(8, evtwt) ;
     h1_["st_sigT1Z1b1"] -> Fill(ST, evtwt) ;
   }
   
  

   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow1"] -> Fill(9, evtwt) ;
     h1_["st_sigT1Z1b2"] -> Fill(ST, evtwt) ;
   }
 


 }
   
 
 //(2)                                                                                                                                                       
 if (ntopcandidates == 0.0    &&   nzcandidates >=1.0){
   h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
   h1_["cutflow1"] -> Fill(5, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow1"] -> Fill(10, evtwt) ;
     h1_["st_sigT0Z1b1"] -> Fill(ST, evtwt) ;
   }



   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow1"] -> Fill(11, evtwt) ;
     h1_["st_sigT0Z1b2"] -> Fill(ST, evtwt) ;
   }



 }



 //(3)                                                                                                                                                     
 if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
   h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
   h1_["cutflow1"] -> Fill(6, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow1"] -> Fill(12, evtwt) ;
     h1_["st_sigT1Z0b1"] -> Fill(ST, evtwt) ;
   }



   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow1"] -> Fill(13, evtwt) ;
     h1_["st_sigT1Z0b2"] -> Fill(ST, evtwt) ;
   }

 }


 //(4)                                                                                                                                                       
 if (ntopcandidates ==0.0    &&   nzcandidates==0.0){
   h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
   h1_["cutflow1"] -> Fill(7, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow1"] -> Fill(14, evtwt) ;
     h1_["st_sigT0Z0b1"] -> Fill(ST, evtwt) ;
   }



   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow1"] -> Fill(15, evtwt) ;
     h1_["st_sigT0Z0b2"] -> Fill(ST, evtwt) ;
   }



 }

 */

 //top,Z,H,b
    //(1)
    if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
      h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(4, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(8, evtwt) ;

	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow2"] -> Fill(1, evtwt) ;
	  h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow2"] -> Fill(2, evtwt) ;
	  h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
	}
	
      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(9, evtwt) ;

        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(3, evtwt) ;
	  h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
        }
	
        else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow2"] -> Fill(4, evtwt) ;
	  h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
        }	
      }  
    }
    //(2)                                                                                                                                                                                        
    if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
      h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(5, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT0Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(10, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(5, evtwt) ;
          h1_["st_sigT0Z1H1b1"] -> Fill(ST, evtwt) ;
        }
        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(6, evtwt) ;
          h1_["st_sigT0Z1H1b2"] -> Fill(ST, evtwt) ;
        }

      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT0Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(11, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(7, evtwt) ;
          h1_["st_sigT0Z1H0b1"] -> Fill(ST, evtwt) ;
        }

        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(8, evtwt) ;
          h1_["st_sigT0Z1H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }

    //(3)                                                                                                                                                                                        
    if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
      h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(6, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z0H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(12, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(9, evtwt) ;
          h1_["st_sigT1Z0H1b1"] -> Fill(ST, evtwt) ;
        }
        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(10, evtwt) ;
          h1_["st_sigT1Z0H1b2"] -> Fill(ST, evtwt) ;
        }

      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z0H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(13, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(11, evtwt) ;
          h1_["st_sigT1Z0H0b1"] -> Fill(ST, evtwt) ;
        }

        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(12, evtwt) ;
          h1_["st_sigT1Z0H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }

    //(4)                                                                                                                                                                                        
    if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
      h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(7, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT0Z0H1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(14, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(13, evtwt) ;
          h1_["st_sigT0Z0H1b1"] -> Fill(ST, evtwt) ;
        }
        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(14, evtwt) ;
          h1_["st_sigT0Z0H1b2"] -> Fill(ST, evtwt) ;
        }

      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT0Z0H0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(15, evtwt) ;
        if( goodBTaggedAK4Jets.size() == 1 ){
          h1_["cutflow2"] -> Fill(15, evtwt) ;
          h1_["st_sigT0Z0H0b1"] -> Fill(ST, evtwt) ;
        }

        else if( goodBTaggedAK4Jets.size() >= 2 ){
          h1_["cutflow2"] -> Fill(16, evtwt) ;
          h1_["st_sigT0Z0H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }
    /*
    
   //Top ,ZH
 //(1)
 if (ntopcandidates >=1.0    &&   nZHcandidates>=1.0){
   h1_["st_sigT1ZH1"] -> Fill(ST,evtwt);
   h1_["cutflow3"] -> Fill(4, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow3"] -> Fill(8, evtwt) ;
     h1_["st_sigT1ZH1b1"] -> Fill(ST, evtwt) ;
   }
   
  

   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow3"] -> Fill(9, evtwt) ;
     h1_["st_sigT1ZH1b2"] -> Fill(ST, evtwt) ;
   }
 


 }
   
 
 //(2)                                                                                                                                                       
 if (ntopcandidates == 0.0    &&   nZHcandidates >=1.0){
   h1_["st_sigT0ZH1"] -> Fill(ST,evtwt);
   h1_["cutflow3"] -> Fill(5, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow3"] -> Fill(10, evtwt) ;
     h1_["st_sigT0ZH1b1"] -> Fill(ST, evtwt) ;
   }



   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow3"] -> Fill(11, evtwt) ;
     h1_["st_sigT0ZH1b2"] -> Fill(ST, evtwt) ;
   }



 }



 //(3)                                                                                                                                                     
 if (ntopcandidates >=1.0    &&   nZHcandidates==0.0){
   h1_["st_sigT1ZH0"] -> Fill(ST,evtwt);
   h1_["cutflow3"] -> Fill(6, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow3"] -> Fill(12, evtwt) ;
     h1_["st_sigT1ZH0b1"] -> Fill(ST, evtwt) ;
   }



   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow3"] -> Fill(13, evtwt) ;
     h1_["st_sigT1ZH0b2"] -> Fill(ST, evtwt) ;
   }

 }


 //(4)                                                                                                                                                       
 if (ntopcandidates ==0.0    &&   nZHcandidates==0.0){
   h1_["st_sigT0ZH0"] -> Fill(ST,evtwt);
   h1_["cutflow3"] -> Fill(7, evtwt) ;

   if( goodBTaggedAK4Jets.size() == 1 ){
     h1_["cutflow3"] -> Fill(14, evtwt) ;
     h1_["st_sigT0ZH0b1"] -> Fill(ST, evtwt) ;
   }



   else if( goodBTaggedAK4Jets.size() >= 2 ){
     h1_["cutflow1"] -> Fill(15, evtwt) ;
     h1_["st_sigT0ZH0b2"] -> Fill(ST, evtwt) ;
   }



 }

    */




 
 // h1_["b_st"] ->Fill(ST, evtwt)
 
 /*
 if ( goodAK4Jets.at(1).getMass() != goodAK4Jets.at(2).getMass() && goodAK4Jets.at(1).getMass() != goodAK4Jets.at(3).getMass() && goodAK4Jets.at(2).getMass() != goodAK4Jets.at(3).getMass()){ 
   TopCandsProducer top;
   top.operator()(tops, goodAK4Jets) ;
   
   double topmass = tops.at(1).getMass()+ tops.at(2).getMass()+tops.at(3).getMass();
   if (topmass>120 && topmass <240){
     cout<< " top mass = " << topmass << endl;
     h1_["topmass"] ->Fill(topmass, evtwt);
   }
 }
 */
 /*
  //Do mass reconstruction
  TLorentzVector lep1, lep2;
  if (zdecayMode_ == "zelel"){
     lep1 = goodElectrons.at(0).getP4();
     lep2 = goodElectrons.at(1).getP4();
  }
  else if (zdecayMode_ == "zmumu"){
     lep1 = goodMuons.at(0).getP4();
     lep2 = goodMuons.at(1).getP4();
}
  if (optimizeReco_ && *h_evttype.product() != "EvtType_Data"){

     GenParticleCollection genPartsInfo;
     genPartsInfo = genpart(evt) ;
     // Declare TLorentzVectors to fill with genParticles
     TLorentzVector bGen, bbarGen, q1, q2;// Z1, Z2;
     TLorentzVector qJet, qbarJet, bJet, bbarJet;
     TLorentzVector had_bjet, lep_bjet, had_bGen, lep_bGen;
     bGen = reco.getGen(genPartsInfo, 5, 8000002);
     bbarGen = reco.getGen(genPartsInfo, -5, 8000002);
     if (bosonMass_ == 91.2){
       q1 = reco.getGen(genPartsInfo, 1, 5, 23);
       q2 = reco.getGen(genPartsInfo, -5, -1, 23);
     }
     else{
       q1 = reco.getGen(genPartsInfo, 1, 5, 25);
       q2 = reco.getGen(genPartsInfo, -5, -1, 25);
     }

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
  if (goodAK4Jets.size() >= 4){
     chi2_result = reco.doReco(goodAK4Jets, bosonMass_, Leptons);
  }
  else if (goodAK4Jets.size() == 3){
     chi2_result.first = -998;
     chi2_result.second = -998;
  }
  else{
     chi2_result.first = -999;
     chi2_result.second = -999;
  }
  //Fill Histograms
  h1_["ZJetMasslep"] ->Fill(Leptons.M(), evtwt);
  h1_["chi2_chi"] ->Fill(chi2_result.first, evtwt);
  h1_["sqrtChi2"] ->Fill(sqrt(chi2_result.first), evtwt);
  if (chi2_result.second == -998)
    h1_["3jets"] ->Fill(1, evtwt);
  else if (chi2_result.second > 0)
    h1_["chi2_mass"] ->Fill(chi2_result.second, evtwt);
  */
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

 if (filterSignal_){h1_["signalEvts"] = fs->make<TH1D>("signalEvts", "signalEvts", 2, 0.5, 2.5) ;}
  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 16, 0.5, 16.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1, "Trig.+l^{+}l^{-}") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2, "75 #lt M(l^{+}l^{-}) #lt 105") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3, "N(AK4) #geq 3") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4, "H_{T} #geq 150") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5, "Z Pt> 150") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6, "leading jet pt > 100") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7, "2nd jet pt > 50") ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8, "N(b jet) #geq 1") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9, "S_{T} #geq 1000") ; 
 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10, "Z-boosted") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(11, "Z-nonboosted") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(12, "at least one Z") ;
  
  // h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10, "top candidate-D") ;                                                                                                              
  //h1_["cutflow"] -> GetXaxis() -> SetBinLabel(11, "top candidate-BC") ;                                                                                                          
  // h1_["cutflow"] -> GetXaxis() -> SetBinLabel(12, "candidate-A") ;                                                                                                                    
  // h1_["cutflow"] -> GetXaxis() -> SetBinLabel(13, "at least one top ") ;    
  //h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9, "N(AK8) #geq 1") ; //this is not the event cut

  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(13, "top candidate-D") ;                                                                                                                                
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(14, "top candidate-BC") ;                                                                                                                                
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(15, "candidate-A") ;                                                                                                                                    
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(16, "at least one top ") ;    


  TFileDirectory pre = fs->mkdir ("pre");
  TFileDirectory sig = fs->mkdir ("sig");
  TFileDirectory cnt = fs->mkdir ("cnt");
  TFileDirectory cat = fs->mkdir ("cat"); 
  TFileDirectory *bookDir[3]; bookDir[0] = &pre; bookDir[1] = &cnt; bookDir[2] = &sig; bookDir[3]= &cat;
  std::vector<string> suffix = {"_pre", "_cnt", ""};

  //btag == 0 control region
  h1_["nob_pt_z"+lep+lep] = fs->make<TH1D>(("nob_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["b_pt_z"+lep+lep] = fs->make<TH1D>(("b_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["nob_st"] = fs->make<TH1D>("nob_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["b_st"] = fs->make<TH1D>("b_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["pt_zlight_pre"] = fs->make<TH1D>("pt_zlight_pre", "p_{T} (Z + q_{light}) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zb_pre"] = fs->make<TH1D>("pt_zb_pre", "p_{T} (Z + b) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zc_pre"] = fs->make<TH1D>("pt_zc_pre", "p_{T} (Z + c) [GeV]", 100, 0., 2000.) ;
  h1_["nob_ht"]= cnt.make<TH1D>("nob_ht", "HT no b tags", 100, 0., 3000.);



for (int i=0; i<3; i++){
     h1_[("npv_noweight"+suffix[i]).c_str()] = bookDir[i]->make<TH1D>( ("npv_noweight"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
     h1_[("npv"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("npv"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
     h1_[("nak4"+suffix[i]).c_str()] =  bookDir[i]->make<TH1D>( ("nak4"+suffix[i]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
     h1_[("ht"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("ht"+suffix[i]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 3000.) ;
     h1_[("st"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("st"+suffix[i]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
     h1_[("met"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met"+suffix[i]).c_str(), "MET [GeV]", 100, 0., 1000.);
     h1_[("metPhi"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("metPhi"+suffix[i]).c_str(), "#Phi(MET)", 20, -5., 5.);

     //jets
     for(int j=1; j<4; ++j){
        string jetPtName = Form("ptak4jet%d", j)+suffix[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
        h1_[jetPtName.c_str()] = bookDir[i]->make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
        string jetEtaName = Form("etaak4jet%d", j)+suffix[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
        h1_[jetEtaName.c_str()] = bookDir[i]->make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-2.4 ,2.4) ;
        string jetCVSName = Form("cvsak4jet%d", j)+suffix[i]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j); 
        h1_[jetCVSName.c_str()] = bookDir[i]->make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
     }
     string jet1METPhiName = "phi_jet1MET"+suffix[i];
     h1_[jet1METPhiName.c_str()] = bookDir[i]->make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;

     //leptons
     string mass_Z = "mass_z"+lep+lep+suffix[i];
     h1_[mass_Z.c_str()] = bookDir[i]->make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
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

 //additional plots
//***** ht,st and met resolution plots in control, signal regions                                                                                                    
 h1_["ht_cnt1"]   = cnt.make<TH1D>("ht-cnt1", ";H_{T} (AK4 jets) [GeV]", 100, 0., 800.) ;
 h1_["st_cnt1"]   = cnt.make<TH1D>("st-cnt1",";S_{T} [GeV]", 100, 0., 800.) ;
 h1_["met_cnt1"]  = cnt.make<TH1D>("met-cnt1", "MET [GeV]", 100, 0., 400.);

 h1_["htsig1"]   = sig.make<TH1D>("htsig1", ";H_{T} (AK4 jets) [GeV]", 100, 0., 2500.) ;
 h1_["stsig1"]   = sig.make<TH1D>("stsig1",";S_{T} [GeV]", 100, 0., 3000.) ;
 h1_["metsig1"]  = sig.make<TH1D>("metsig1", "MET [GeV]", 100, 0., 500.);




  h1_["nbjets"] = sig.make<TH1D>("nbjets", ";N(b jets);;" , 11, -0.5,10.5) ; 
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

  // plots to study T mass recontruction
  if (optimizeReco_){

    // h1_["genZ"] = sig.make<TH1D>("genZ", ";M (Gen Z Boson) [GeV];;", 20, 0., 200.);
    // h1_["genTMass"] = sig.make<TH1D>("genTMass", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     // h1_["ZJetMass"] = sig.make<TH1D>("ZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
    // h1_["gentMass"] = sig.make<TH1D>("gentMass", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     h1_["hadTJetMass"] = sig.make<TH1D>("TJetMass", ";JetM (Hadronic T quark) [GeV];;", 100, 0., 1500.);
     h1_["lepTJetMass"] = sig.make<TH1D>("TJetMasslep", ";M (T Jet Leptonic);;", 100, 0., 1500.);
     h1_["genH"] = sig.make<TH1D>("genH", ";M (Gen H Boson) [GeV];;", 20, 0., 200.); 
     h1_["HJetMass"] = sig.make<TH1D>("HJetMass", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);

     h1_["hadgenTMass"] = sig.make<TH1D>("hadgenTMass", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["hadgentMass"] = sig.make<TH1D>("hadgentMass", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     h1_["lepgenTMass"] = sig.make<TH1D>("lepgenTMass", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["lepgentMass"] = sig.make<TH1D>("lepgentMass", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
    
     //W tags (category B) 
     h1_["hadTJetMass_B"] = sig.make<TH1D>("TJetMass_B", ";JetM (Hadronic T quark) [GeV];;", 100, 0., 1500.);
     h1_["lepTJetMass_B"] = sig.make<TH1D>("TJetMasslep_B", ";M (T Jet Leptonic);;", 100, 0., 1500.);
     h1_["genH_B"] = sig.make<TH1D>("genH_B", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
     h1_["HJetMass_B"] = sig.make<TH1D>("HJetMass_B", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);

     h1_["hadgenTMass_B"] = sig.make<TH1D>("hadgenTMass_B", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["hadgentMass_B"] = sig.make<TH1D>("hadgentMass_B", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     h1_["lepgenTMass_B"] = sig.make<TH1D>("lepgenTMass_B", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["lepgentMass_B"] = sig.make<TH1D>("lepgentMass_B", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     
     //category C
     h1_["hadTJetMass_C"] = sig.make<TH1D>("TJetMass_C", ";JetM (Hadronic T quark) [GeV];;", 100, 0., 1500.);
     h1_["lepTJetMass_C"] = sig.make<TH1D>("TJetMasslep_C", ";M (T Jet Leptonic);;", 100, 0., 1500.);
     h1_["genH_C"] = sig.make<TH1D>("genH_C", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
     h1_["HJetMass_C"] = sig.make<TH1D>("HJetMass_C", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);

     h1_["hadgenTMass_C"] = sig.make<TH1D>("hadgenTMass_C", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["hadgentMass_C"] = sig.make<TH1D>("hadgentMass_C", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     h1_["lepgenTMass_C"] = sig.make<TH1D>("lepgenTMass_C", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["lepgentMass_C"] = sig.make<TH1D>("lepgentMass_C", ";M(Gen t quark) [GeV];;", 100, 0., 1200);

     //category D
     h1_["hadTJetMass_D"] = sig.make<TH1D>("TJetMass_D", ";JetM (Hadronic T quark) [GeV];;", 100, 0., 1500.);
     h1_["lepTJetMass_D"] = sig.make<TH1D>("TJetMasslep_D", ";M (T Jet Leptonic);;", 100, 0., 1500.);
     h1_["genH_D"] = sig.make<TH1D>("genH_D", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
     h1_["HJetMass_D"] = sig.make<TH1D>("HJetMass_D", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);

     h1_["hadgenTMass_D"] = sig.make<TH1D>("hadgenTMass_D", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["hadgentMass_D"] = sig.make<TH1D>("hadgentMass_D", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     h1_["lepgenTMass_D"] = sig.make<TH1D>("lepgenTMass_D", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
     h1_["lepgentMass_D"] = sig.make<TH1D>("lepgentMass_D", ";M(Gen t quark) [GeV];;", 100, 0., 1200);


     if (bosonMass_ == 91.2){
        //higgs
       //nonboosted
       h1_["hadgenZ_nonb"] = sig.make<TH1D>("hadgenZ_nonb", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass_nonb"] = sig.make<TH1D>("hadZJetMass_nonb", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt_nonb"] = sig.make<TH1D>("hadZJetPt_nonb", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       
       //boosted
       h1_["hadgenZ_b"] = sig.make<TH1D>("hadgenZ_b", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass_b"] = sig.make<TH1D>("hadZJetMass_b", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt_b"] = sig.make<TH1D>("hadZJetPt_b", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);


       //  h1_["genZ"] = sig.make<TH1D>("genZ", ";M (Gen Z Boson) [GeV];;", 20, 0., 200.);
       // h1_["ZJetMass"] = sig.make<TH1D>("ZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       // h1_["ZJetPt"] = sig.make<TH1D>("ZJetPt", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       h1_["hadTJetPt"] = sig.make<TH1D>("hadTJetPt", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["lepTJetPt"] = sig.make<TH1D>("lepTJetPT", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);
       
       h1_["hadgenZ"] = sig.make<TH1D>("hadgenZ", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass"] = sig.make<TH1D>("hadZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt"] = sig.make<TH1D>("hadZJetPt", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       h1_["lepZ"] = sig.make<TH1D>("lepZ", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
       h1_["lepZPt"] = sig.make<TH1D>("lepZPt", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);

       //W tags (Cateogory B)
       h1_["hadTJetPt_B"] = sig.make<TH1D>("hadTJetPt_B", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["lepTJetPt_B"] = sig.make<TH1D>("lepTJetPt_B", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["hadgenZ_B"] = sig.make<TH1D>("hadgenZ_B", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass_B"] = sig.make<TH1D>("hadZJetMass_B", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt_B"] = sig.make<TH1D>("hadZJetPt_B", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       h1_["lepZ_B"] = sig.make<TH1D>("lepZ_B", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
       h1_["lepZPt_B"] = sig.make<TH1D>("lepZPt_B", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);
 
       //category C
       h1_["hadTJetPt_C"] = sig.make<TH1D>("hadTJetPt_C", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["lepTJetPt_C"] = sig.make<TH1D>("lepTJetPt_C", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["hadgenZ_C"] = sig.make<TH1D>("hadgenZ_C", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass_C"] = sig.make<TH1D>("hadZJetMass_C", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt_C"] = sig.make<TH1D>("hadZJetPt_C", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       h1_["lepZ_C"] = sig.make<TH1D>("lepZ_C", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
       h1_["lepZPt_C"] = sig.make<TH1D>("lepZPt_C", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);

       //category D
       h1_["hadTJetPt_D"] = sig.make<TH1D>("hadTJetPt_D", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["lepTJetPt_D"] = sig.make<TH1D>("lepTJetPt_D", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);
       h1_["hadgenZ_D"] = sig.make<TH1D>("hadgenZ_D", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass_D"] = sig.make<TH1D>("hadZJetMass_D", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt_D"] = sig.make<TH1D>("hadZJetPt_D", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       h1_["lepZ_D"] = sig.make<TH1D>("lepZ_D", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
       h1_["lepZPt_D"] = sig.make<TH1D>("lepZPt_D", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);

       //category E(z boson with jet masses equal)

       h1_["hadgenZ_E"] = sig.make<TH1D>("hadgenZ_E", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadZJetMass_E"] = sig.make<TH1D>("hadZJetMass_E", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
       h1_["hadZJetPt_E"] = sig.make<TH1D>("hadZJetPt_E", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       h1_["hadTJetPt_E"] = sig.make<TH1D>("hadTJetPt_E", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);


     }
     else{
       //higgs
       //nonboosted
       h1_["hadgenH_nonb"] = sig.make<TH1D>("hadgenH_nonb", ";M (Gen H Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadHJetMass_nonb"] = sig.make<TH1D>("hadHJetMass_nonb", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
       h1_["hadHJetPt_nonb"] = sig.make<TH1D>("hadHJetPt_nonb", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
       
       //boosted
       h1_["hadgenH_b"] = sig.make<TH1D>("hadgenH_b", ";M (Gen H Boson -hadronic) [GeV];;", 20, 0., 200.);
       h1_["hadHJetMass_b"] = sig.make<TH1D>("hadHJetMass_b", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
       h1_["hadHJetPt_b"] = sig.make<TH1D>("hadHJetPt_b", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);



       //	 h1_["genH"] = sig.make<TH1D>("genH", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
       // h1_["HJetMass"] = sig.make<TH1D>("HJetMass", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
       // h1_["HJetPt"] = sig.make<TH1D>("HJetPt", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
	 h1_["hadTJetPt"] = sig.make<TH1D>("hadTJetPt", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
	 h1_["lepTJetPt"] = sig.make<TH1D>("lepTJetPT", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);


	 h1_["hadgenH"] = sig.make<TH1D>("hadgenH", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass"] = sig.make<TH1D>("hadHJetMass", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt"] = sig.make<TH1D>("hadHJetPt", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
	 h1_["lepH"] = sig.make<TH1D>("lepH", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt"] = sig.make<TH1D>("lepHPt", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);

    

	 //W tags (category B)
	 h1_["hadTJetPt_B"] = sig.make<TH1D>("hadTJetPt_B", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
         h1_["lepTJetPt_B"] = sig.make<TH1D>("lepTJetPT_B", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);


         h1_["hadgenH_B"] = sig.make<TH1D>("hadgenH_B", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass_B"] = sig.make<TH1D>("hadHJetMass_B", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt_B"] = sig.make<TH1D>("hadHJetPt_B", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
         h1_["lepH_B"] = sig.make<TH1D>("lepH_B", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt_B"] = sig.make<TH1D>("lepHPt_B", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);

	 //Category C
	 h1_["hadTJetPt_C"] = sig.make<TH1D>("hadTJetPt_C", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
         h1_["lepTJetPt_C"] = sig.make<TH1D>("lepTJetPT_C", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);


         h1_["hadgenH_C"] = sig.make<TH1D>("hadgenH_C", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass_C"] = sig.make<TH1D>("hadHJetMass_C", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt_C"] = sig.make<TH1D>("hadHJetPt_C", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
         h1_["lepH_C"] = sig.make<TH1D>("lepH_C", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt_C"] = sig.make<TH1D>("lepHPt_C", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);


	 //category D
	 h1_["hadTJetPt_D"] = sig.make<TH1D>("hadTJetPt_D", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
         h1_["lepTJetPt_D"] = sig.make<TH1D>("lepTJetPT_D", ";Pt (leptonic T Prime) [GeV];;", 50, 0., 1000.);


         h1_["hadgenH_D"] = sig.make<TH1D>("hadgenH_D", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass_D"] = sig.make<TH1D>("hadHJetMass_D", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt_D"] = sig.make<TH1D>("hadHJetPt_D", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
         h1_["lepH_D"] = sig.make<TH1D>("lepH_D", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt_D"] = sig.make<TH1D>("lepHPt_D", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);



 }
     


       
     h1_["hadtJetMass"] = sig.make<TH1D>("tJetMass-hadronic", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass"] = sig.make<TH1D>("TJetMasslep", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt"] = sig.make<TH1D>("tJetPt-hadronic", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt"] = sig.make<TH1D>("TJetPtlep", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);
     h2_["genhadtJetMasshad"] = sig.make<TH2D>("genthadthad", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep"] = sig.make<TH2D>("genthadtlep", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);

     //W tags (category B)

     h1_["hadtJetMass_B"] = sig.make<TH1D>("tJetMass-hadronic_B", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass_B"] = sig.make<TH1D>("TJetMasslep_B", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt_B"] = sig.make<TH1D>("tJetPt-hadronic_B", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt_B"] = sig.make<TH1D>("TJetPtlep_B", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);
     h2_["genhadtJetMasshad_B"] = sig.make<TH2D>("genthadthad_B", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep_B"] = sig.make<TH2D>("genthadtlep_B", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);   

     h1_["lepWJetMass_B"]= sig.make<TH1D>("lepWJetMass_B", ";M (leptonic W jet) [GeV];;", 20, 0., 200.);
     h1_["lepWJetPt_B"]= sig.make<TH1D>("lepWJetPt_B", ";Pt (leptonic W jet) [GeV];;", 100, 0., 1200.);
     h1_["hadWJetMass_B"]= sig.make<TH1D>("hadWJetMass_B", ";M (Hadronic W jet) [GeV];;", 20, 0., 200.);
     h1_["hadWJetPt_B"]= sig.make<TH1D>("hadWJetPt_B", ";Pt (hadronic W jet) [GeV];;", 100, 0., 1200.);

     //category C

     h1_["hadtJetMass_C"] = sig.make<TH1D>("tJetMass-hadronic_C", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass_C"] = sig.make<TH1D>("TJetMasslep_C", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt_C"] = sig.make<TH1D>("tJetPt-hadronic_C", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt_C"] = sig.make<TH1D>("TJetPtlep_C", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);

     h1_["hadmergedJetMass_C"] = sig.make<TH1D>("hadmergedJetMass_C", ";merged jet mass (Hadronic) [GeV];;", 100, 0., 1000.);
     h1_["hadmergedJetPt_C"]   = sig.make<TH1D>("hadmergedJetPt_C", ";merged jet Pt (Hadronic) [GeV];;", 100, 0., 1000.);
     h1_["lepmergedJetMass_C"] = sig.make<TH1D>("lepmergedJetMass_C", ";merged jet mass (leptonic) [GeV];;", 100, 0., 1000.);
     h1_["lepmergedJetPt_C"]   = sig.make<TH1D>("lepmergedJetPt_C", ";merged jet Pt (Leptonic) [GeV];;", 100, 0., 1000.);

     h2_["genhadtJetMasshad_C"] = sig.make<TH2D>("genthadthad_C", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep_C"] = sig.make<TH2D>("genthadtlep_C", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);

     //category D

     h1_["hadtJetMass_D"] = sig.make<TH1D>("tJetMass-hadronic_D", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass_D"] = sig.make<TH1D>("TJetMasslep_D", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt_D"] = sig.make<TH1D>("tJetPt-hadronic_D", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt_D"] = sig.make<TH1D>("TJetPtlep_D", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);
     h2_["genhadtJetMasshad_D"] = sig.make<TH2D>("genthadthad_D", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep_D"] = sig.make<TH2D>("genthadtlep_D", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
 
  }
  //ZH candidates
  h1_["ZH_mass_b"]  = sig.make<TH1D>("ZHmass-boosted", ";M(H-boosted) [GeV];;", 100, 0., 400);
  h1_["ZH_Pt_b"]  = sig.make<TH1D>("ZHPt-boosted", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
  h1_["nZHcandidatejets_b"] = sig.make<TH1D>("nZHcandidate-boosted", ";N(H jets-boosted);;" , 21, -0.5, 20.5);




  h1_["ZH_mass_nb"]  = sig.make<TH1D>("ZHmassnb", ";M(H) [GeV];;", 100, 0., 400);
  h1_["ZH_Pt_nb"]  = sig.make<TH1D>("ZHPtnb", ";Pt(H) [GeV];;", 100, 0., 1200);
  h1_["nZHcandidatejets_nb"] = sig.make<TH1D>("nZHcandidatesnb", ";N(H jets);;" , 21, -0.5, 20.5);


  h1_["nZHcandidatejets"] = sig.make<TH1D>("nZHcandidates-tot", ";N(H jets);;" , 21, -0.5, 20.5);
  h1_["nZHcandidatejets1"] = sig.make<TH1D>("nZHcandidates1-tot", ";N(H jets);;" , 21, -0.5, 20.5);




  //H candidates                                                                                                                             
  h1_["H_mass_b"]  = sig.make<TH1D>("Hmass-boosted", ";M(H-boosted) [GeV];;", 100, 0., 400);
  h1_["H_Pt_b"]  = sig.make<TH1D>("HPt-boosted", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
  h1_["nHcandidatejets_b"] = sig.make<TH1D>("nHcandidate-boosted", ";N(H jets-boosted);;" , 21, -0.5, 20.5);




  h1_["H_mass_nb"]  = sig.make<TH1D>("Hmassnb", ";M(H) [GeV];;", 100, 0., 400);
  h1_["H_Pt_nb"]  = sig.make<TH1D>("HPtnb", ";Pt(H) [GeV];;", 100, 0., 1200);
  h1_["nHcandidatejets_nb"] = sig.make<TH1D>("nHcandidatesnb", ";N(H jets);;" , 21, -0.5, 20.5);


  h1_["nHcandidatejets"] = sig.make<TH1D>("nHcandidates-tot", ";N(H jets);;" , 21, -0.5, 20.5);
  h1_["nHcandidatejets1"] = sig.make<TH1D>("nHcandidates1-tot", ";N(H jets);;" , 21, -0.5, 20.5);




  // Z candidates
 

  h1_["Z_mass_a"]  = sig.make<TH1D>("Zmass-boosted", ";M(Z-boosted) [GeV];;", 100, 0., 400);
  h1_["Z_Pt_a"]  = sig.make<TH1D>("ZPt-boosted", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
  h1_["nzcandidatejets_a"] = sig.make<TH1D>("nzcandidate-boosted", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);


  h1_["Z_mass_b"]  = sig.make<TH1D>("Zmass", ";M(Z) [GeV];;", 100, 0., 400);
  h1_["Z_Pt_b"]  = sig.make<TH1D>("ZPt", ";Pt(Z) [GeV];;", 100, 0., 1200);
  h1_["nzcandidatejets_b"] = sig.make<TH1D>("nzcandidates", ";N(Z jets);;" , 21, -0.5, 20.5);

  
  h1_["nzcandidatejets_tot"] = sig.make<TH1D>("nzcandidates-tot", ";N(Z jets);;" , 21, -0.5, 20.5);
  h1_["nzcandidatejets1_tot"] = sig.make<TH1D>("nzcandidates1-tot", ";N(Z jets);;" , 21, -0.5, 20.5);  
// cat A
  h1_["top_mass_a"]  = sig.make<TH1D>("topmas-A", ";M( t quark) [GeV];;", 100, 0., 400);
  h1_["top_Pt_a"]  = sig.make<TH1D>("topPt-A", ";Pt( t quark) [GeV];;", 100, 0., 1200);
  h1_["ntopcandidatejets_a"] = sig.make<TH1D>("ntopcandidate-A", ";N(top jets);;" , 21, -0.5, 20.5);

  h1_["top_mass_bc"]  = sig.make<TH1D>("topmass-Bc", ";M( t quark) [GeV];;", 100, 0., 400);
  h1_["top_Pt_bc"]  = sig.make<TH1D>("topPt-BC", ";Pt( t quark) [GeV];;", 100, 0., 1200);
  h1_["ntopcandidatejets_bc"] = sig.make<TH1D>("ntopcandidate-BC", ";N(top jets);;" , 21, -0.5, 20.5);

  // cat D
  h1_["top_mass_d"]  = sig.make<TH1D>("topmass-D", ";M( t quark) [GeV];;", 100, 0., 400);
  h1_["top_Pt_d"]  = sig.make<TH1D>("topPt-D", ";Pt( t quark) [GeV];;", 100, 0., 1200);
  h1_["ntopcandidatejets_d"] = sig.make<TH1D>("ntopcandidate-D", ";N(top jets);;" , 21, -0.5, 20.5);

  //W and light jet(BC)
  h1_["W_mass_bc"]  = sig.make<TH1D>("Wmass-BC", ";M( W boson) [GeV];;", 100, 0., 400);
  h1_["nWcandidatejets_bc"] = sig.make<TH1D>("nWcandidate-BC", ";N(W candidate jets);;" , 21, -0.5, 20.5);

  h1_["lightjet_mass_bc"]  = sig.make<TH1D>("lightjetmass-BC", ";M( light jet) [GeV];;", 100, 0., 400);
  h1_["nlightjetcandidatejets_bc"] = sig.make<TH1D>("nlightjetcandidate", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
  //total top ( A+ BC+D)
  h1_["ntopcandidatejets"] = sig.make<TH1D>("ntopcandidate-tot", ";N(top jets);;" , 21, -0.5, 20.5);
  h1_["ntopcandidatejets1"] = sig.make<TH1D>("ntopcandidate1-tot", ";N(top jets);;" , 21, -0.5, 20.5);


  h1_["ZJetMasslep"] = sig.make<TH1D>("ZJetMasslep", ";M (Z Jet Leptonic);;", 20, 0., 200.);
  h1_["chi2_chi"] = sig.make<TH1D>("chi2_chi", ";#chi^{2};;", 100, 0., 500.);
  h1_["sqrtChi2"] = sig.make<TH1D>("sqrtChi2", ";#chi;;", 100, 0., 500.);
  h1_["chi2_mass"] = sig.make<TH1D>("chi_mass", ";M_{#chi^{2}}(B);;", 60, 200., 2000.);
  h1_["chi_mass_cnt"] = cnt.make<TH1D>("chi_mass_cnt", ";M_{#chi^{2}}(B);;", 60, 200., 2000.);
  h1_["chi2_chi_cnt"] = cnt.make<TH1D>("chi2_chi_cnt", ";#chi^{2};;", 100, 0., 500.);
  h1_["3jets_cnt"] = cnt.make<TH1D>("3jets_cnt", "events w/ <4 jets", 2, .5, 2.5);
  h1_["3jets"] = sig.make<TH1D>("3jets", "events w/ <4 jets", 2, .5, 2.5);

  // Tprime Tprime Categories                                                                                                                                                               
  //for Z and H categories seperetely along with b
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

  //for ZH(combined category)
  h1_["cutflow3"] = cat.make<TH1D>("cutflow3", "cut flow", 15, 0.5, 15.5) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(1, "no T,Z,b ") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(4, "T1ZH1 ") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(5, "T0ZH1") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(6, "T1ZH0 ") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(7, "T0ZH0") ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(8, "T1ZH1b1 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(9, "T1ZH1b2 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(10, "T0ZH1b1 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(11, "T0ZH0b2 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(12, "T1ZH0b1 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(13, "T1ZH0b2 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(14, "T0ZH0b1 " ) ;
  h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(15, "T0ZH0b2 " ) ;
 



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

  //ZH category (expanded mass window from Z)

  h1_["st_sigT1ZH1"] =cat.make<TH1D>("ST_sigT1ZH1", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT0ZH1"] =cat.make<TH1D>("ST_sigT0ZH1", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT1ZH0"] =cat.make<TH1D>("ST_sigT1ZH0", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT0ZH0"] =cat.make<TH1D>("ST_sigT0ZH0", ";S_{T} [Gev];;" , 100,0.,3000.);

  h1_["st_sigT1ZH1b1"] =cat.make<TH1D>("ST_sigT1ZH1b1", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT1ZH1b2"] =cat.make<TH1D>("ST_sigT1ZH1b2", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT0ZH1b1"] =cat.make<TH1D>("ST_sigT0ZH1b1", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT0ZH1b2"] =cat.make<TH1D>("ST_sigT0ZH1b2", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT1ZH0b1"] =cat.make<TH1D>("ST_sigT1ZH0b1", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT1ZH0b2"] =cat.make<TH1D>("ST_sigT1ZH0b2", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT0ZH0b1"] =cat.make<TH1D>("ST_sigT0ZH0b1", ";S_{T} [Gev];;" , 100,0.,3000.);
  h1_["st_sigT0ZH0b2"] =cat.make<TH1D>("ST_sigT0ZH0b2", ";S_{T} [Gev];;" , 100,0.,3000.);




  h1_["nZll"] = cat.make<TH1D>("nZll", ";N(Zll);;" , 8, -0.5,8.5) ;
  h1_["nak4_A"] = cat.make<TH1D>("nAK4-CatA", ";N(AK4);;" , 21, -0.5, 20.5) ;
  h1_["nak4_B"] = cat.make<TH1D>("nAK4-CatB", ";N(AK4);;" , 21, -0.5, 20.5) ;
  h1_["nak4_C"] = cat.make<TH1D>("nAK4-CatC", ";N(AK4);;" , 21, -0.5, 20.5);


  h1_["ptak4jet1_sig"] =cat.make<TH1D>("ak4Jet1Pt", ";P_{T} [Gev];;" , 50,0.,1000.);
  h1_["ptak4jet2_sig"] =cat.make<TH1D>("ak4Jet2Pt", ";P_{T} [Gev];;" , 50,0.,1000.);
  h1_["ptak4jet3_sig"] =cat.make<TH1D>("ak4Jet3Pt", ";P_{T} [Gev];;" , 50,0.,1000.);
  h1_["ptak4jet4_sig"] =cat.make<TH1D>("ak4Jet4Pt", ";P_{T} [Gev];;" , 50,0.,1000.);
  h1_["ptak4jet5_sig"] =cat.make<TH1D>("ak4Jet5Pt", ";P_{T} [Gev];;" , 50,0.,1000.);
  h1_["ptak4jet6_sig"] =cat.make<TH1D>("ak4Jet6Pt", ";P_{T} [Gev];;" , 50,0.,1000.);

  h1_["nak4jet_sig"] = cat.make<TH1D>("nAK4sig", ";N(AK4);;" , 21, -0.5, 20.5) ;

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
  }

}

double OS2LAna::GetDYNLOCorr(const double dileppt) {

  std::unique_ptr<TFile >file_DYNLOCorr = std::unique_ptr<TFile>(new TFile(fname_DYNLOCorr_.c_str())) ;
  std::unique_ptr<TF1>fun_DYNLOCorr = std::unique_ptr<TF1>(dynamic_cast<TF1*>(file_DYNLOCorr->Get(funname_DYNLOCorr_.c_str()))) ; 
  double EWKNLOkfact = fun_DYNLOCorr->Eval(dileppt);
  return EWKNLOkfact; 

}

double OS2LAna::ZptCorr(vlq::Candidate zll, double p0, double p1){
    double pt= zll.getPt();
    double scale = p1*pt + p0;
    return(scale);
}

double OS2LAna::htCorr(double ht, double p0, double p1){
  return(p1*ht + p0);
}

void OS2LAna::endJob() {

  return ; 
}

DEFINE_FWK_MODULE(OS2LAna);
