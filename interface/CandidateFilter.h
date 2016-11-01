#ifndef ANALYSIS_VLQANA_CANDIDATEFILTER_HH
#define ANALYSIS_VLQANA_CANDIDATEFILTER_HH
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

class CandidateFilter {
  public:
    CandidateFilter (edm::ParameterSet const& iConfig) : 
      massMin_(iConfig.getParameter<double> ("massMin")), 
      massMax_(iConfig.getParameter<double> ("massMax")), 
      ptMin_(iConfig.getParameter<double> ("ptMin")), 
      ptMax_(iConfig.getParameter<double> ("ptMax")) { }
    ~CandidateFilter () {}  

    void operator () (const vlq::CandidateCollection cands, vlq::CandidateCollection& filteredcands) {
      filteredcands.clear() ; 
      for ( auto cand : cands ) {
        double mass = cand.getMass() ;
        double pt = cand.getPt() ; 
        if ( mass > massMin_ && mass < massMax_ && pt > ptMin_ && pt < ptMax_ ) filteredcands.push_back(cand) ; 
      } 
      return ; 
    }

    void operator () (const vlq::CandidateCollection cands, vlq::CandidateCollection& filteredcands, vlq::CandidateCollection& lep1, vlq::CandidateCollection& lep2, vlq::CandidateCollection& filterlep1 ,vlq::CandidateCollection& filterlep2) {
      filteredcands.clear() ;
      for ( auto cand : cands ) {
        double mass = cand.getMass() ;
        double pt = cand.getPt() ;
        if ( mass > massMin_ && mass < massMax_ && pt > ptMin_ && pt < ptMax_ ) filteredcands.push_back(cand) ;
      }
      filterlep1.clear();
      filterlep2.clear();
      TLorentzVector p4l1,p4l2;
      for (unsigned i=0; i<lep1.size();i++){
	p4l1 = lep1.at(i).getP4();
	p4l2 = lep2.at(i).getP4();
	
	if ((p4l1+p4l2).Mag()>massMin_ && (p4l1+p4l2).Mag() < massMax_ && (p4l1+p4l2).Pt() > ptMin_ && (p4l1+p4l2).Pt() < ptMax_){ 
	  
	  
	  filterlep1.push_back(lep1.at(i));
	  filterlep2.push_back(lep2.at(i));
	}

      }
      return;
    }

 private:
    double massMin_ ;
    double massMax_ ; 
    double ptMin_ ; 
    double ptMax_ ; 
}; 
#endif 
