#!/bin/python

import subprocess

dir = 'cat'

if dir == 'pre':
    options = [
       
        ['cutflow'],
        ['mass_zelel_pre'],
        ['mass_Zelel_pre'],
        ['ht_pre'],
        ['mass_zmumu_pre'],
        ['dr_mumu_pre'],
        ['pt_zmumu_pre'],
        ['pt_mu1_pre'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],
        ['Iso04_mu_pre'],
        ['Dz_mu_pre'],
        ['D0_mu_pre'],
        ['mass_zmumu_cnt'],
        ['dr_mumu_cnt'],
        ['pt_zmumu_cnt'],
        ['pt_mu1_cnt'],
        ['eta_mu1_cnt'],
        ['pt_mu2_cnt'],
        ['eta_mu2_cnt'],             
       
        ['Eta_EB_el_pre'],
        ['Eta_EE_el_pre'],
        ['scEta_EB_el_pre'],
        ['scEta_EE_el_pre'],
        ['Iso03_EB_el_pre'],
        ['Iso03_EE_el_pre'],
        ['dEtaIn_EB_el_pre'],
        ['dEtaIn_EE_el_pre'],
        ['dPhiIn_EB_el_pre'],
        ['dPhiIn_EE_el_pre'],
        ['Dz_EB_el_pre'],
        ['Dz_EE_el_pre'],
        ['D0_EB_el_pre'],
        ['D0_EE_el_pre'],
        ['le1E'],
        ['le2E'],
        ['zllE'],
        ['el1E'],
        ['el2E'],
        ['nle1'],
        ['nle2'],
        ['nzll1'],
        ['nel1'],
        ['sumle'],
             
        ['pt_el1'], 
        ['npv_noweight_pre'],
        ['npv_pre'],
        ['nak4_pre'],
        ['ht_pre'],
        ['st_pre'],
        ['met_pre'],
        ['metPhi_pre'],
        ['ptak4jet1_pre'],
        ['etaak4jet1_pre'],
        ['cvsak4jet1_pre'],
        ['ptak4jet2_pre'],
        ['etaak4jet2_pre'],
        ['cvsak4jet2_pre'],
        ['ptak4jet3_pre'],
        ['etaak4jet3_pre'],
        ['cvsak4jet3_pre'],
        ['phi_jet1MET_pre'],
        ['mass_zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],
       
        ]

elif dir == 'cnt':
    options = [
        # ['genZ'],
        # ['genTMass'],
        # ['ZJetMass'],
        # ['hadTJetMass'],
        # ['lepTJetMass'],
        
        # ['chi_mass','4'],
        # ['ht_pre'],           
        # ['nob_pt_zelel'],           
        
        #    ['ptak4jet1_cnt'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],

         ['Iso04_mu_pre'],
         ['Dz_mu_pre'],
        ['D0_mu_pre'],
        ['st_pre'],
        ['st'],
        ['mass_Zelel_pre'],
        ['nob_ht'],
        ['1b_ht'],
        ['ht_pre'],

        ['mass_zelel_pre'],
        ['mass_Zelel_pre'],
        ['le1E'],
        ['le2E'],
        ['zllE'],

        
        ['cutflow'],
        ['cutflow1'],
        ['cutflow2'],
        ['cutflow3'],
        ['cutflow4'],
        


        ['mass_zelel_pre'],
        
        ['mass_Zelel_pre'],
        ['ht_pre'],
        ['nbjets'],
        ['nob_ht'],
        ['1b_ht'],
        ['nob_st'],
        ['nbjets_met_sig'],
        ['nbjets_met_cnt'],
        ['nbjets_met_0btagcnt'],
        ['nbjets_met_1btagcnt'],
        ['ht_met_0btagcnt'],
        ['1b_ht'],
        ['ht_met_1btagcnt'],
        ['lowmet_ht'],
        ['st_met_0btagcnt'],
        ['1b_st'],
        ['st_met_1btagcnt'],
        ['lowmet_st'],
        ['btagSF_0btag'],
        ['btagSF_1btag'],
        ['nbjets_cnt'],
        ['mass_zmumu_pre'],
        ['dr_mumu_pre'],
        ['pt_zmumu_pre'],
        ['pt_mu1_pre'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],
        ['Iso04_mu_pre'],
        ['Dz_mu_pre'],
        ['D0_mu_pre'],
        ['mass_zmumu_cnt'],
        ['dr_mumu_cnt'],
        ['pt_zmumu_cnt'],
        ['pt_mu1_cnt'],
        ['eta_mu1_cnt'],
        ['pt_mu2_cnt'],
        ['eta_mu2_cnt'],             
        
        ['nob_ht'],
        # ['cutflow'],
           #  ['cutflow1'],
        
        # ['ST_sig'],       
        # ['ST_sig1b'],
        # ['ST_sig2b'],
        # ['ST_sigT1Z1'],
        # ['ST_sigT0Z1'],
        # ['ST_sigT1Z0'],
        # ['ST_sigT0Z0'],
        
        # ['ST_sigT1ZH1'],
        # ['ST_sigT0ZH1'],
        # ['ST_sigT1ZH0'],
        # ['ST_sigT0ZH0'],
        
        #  ['ST_sigT1Z1H1'],
        #  ['ST_sigT1Z1H0'],
        #  ['ST_sigT0Z1H1'],
        #  ['ST_sigT0Z1H0'],
        #  ['ST_sigT1Z0H1'],
        #  ['ST_sigT1Z0H0'],
        #  ['ST_sigT0Z0H1'],
        #  ['ST_sigT0Z0H0'],
        
        # ['ST_sigT1Z1H1b1'],
        # ['ST_sigT1Z1H1b2'],
        # ['ST_sigT1Z1H0b1'],
        # ['ST_sigT1Z1H0b2'],
        # ['ST_sigT0Z1H1b1'],
        # ['ST_sigT0Z1H1b2'],
        # ['ST_sigT0Z1H0b1'],
        # ['ST_sigT0Z1H0b2'],
        # ['ST_sigT1Z0H1b1'],
        # ['ST_sigT1Z0H1b2'],
        # ['ST_sigT1Z0H0b1'],
        # ['ST_sigT1Z0H0b2'],
        # ['ST_sigT0Z0H1b1'],
        # ['ST_sigT0Z0H1b2'],
        # ['ST_sigT0Z0H0b1'],
        # ['ST_sigT0Z0H0b2'],
        
        
        # ['ST_sigT1Z1b1'],
        # ['ST_sigT1Z1b2'],
            # ['ST_sigT0Z1b1'],
        # ['ST_sigT0Z1b2'],
        # ['ST_sigT1Z0b1'],
        # ['ST_sigT1Z0b2'],
        # ['ST_sigT0Z0b1'],
        # ['ST_sigT0Z0b2'],
        
        
        # ['ST_1A'],
        #  ['ST_1B'],
        #  ['ST_2A'],
        #  ['ST_2B'],
        #  ['ST_3A'],
        #  ['ST_3B'],
        #  ['ST_4A'],
        #  ['ST_4B'],
        #  ['cutflow3'],
        # ['cutflow1'],
        #  ['ak4Jet1Pt'],            
        # ['ak4Jet2Pt'],
        # ['ak4Jet3Pt'],
        # ['ak4Jet4Pt'],
        #['ak4Jet5Pt'],
        #['ak4Jet6Pt'],
        #['nak4sig'],
        
        
        # ['cutflow'],
        #['b_st'],
        #['nZll'],
        #['nAK4-CatA'],
        #['nAK4-CatB'],
        #['nAK4-CatC'],
          ['ZHmass-boosted'],                                                                 
                                     
         ['ZHPt-boosted'],                                                                                                                                             
         ['nZHcandidate-boosted'],                                                                                                                                     
         ['ZHmassnb'],                                                                                                                                                 
         ['ZHPtnb'],                                                                                                                                                   
         ['nZHcandidatesnb'],                                                                                                                                          
        ['nZHcandidates-tot'],  
        ['nZHcandidates1-tot'],

        
        ['Hmass-boosted-cnt'],
        ['HPt-boosted-cnt'],
        ['nHcandidate-boosted-cnt'],
        ['Hmassnb-cnt'],
        ['HPtnb-cnt'],
        ['nHcandidatesnb-cnt'],
        ['nHcandidates-tot-cnt'],                                                                                                                                             
        ['nHcandidates1-tot-cnt'],
        
        
        
        
        # ['Zmass-boosted'],
        # ['ZPt-boosted'],
        # ['nzcandidate-boosted'],
        
        #  ['Zmass'],
        # ['ZPt'],
        # ['nzcandidates'],
        #  ['nzcandidates-tot'],
        #  ['nzcandidates1-tot'], 
        
        
        
        
        #  ['topmass-D'],
        # ['topPt-D'],
        # ['ntopcandidate-D'],

        
        # ['Wmass-BC'],
        # ['nWcandidate-BC'],
        # ['lightjetmass-BC'],
        # ['nlightjetcandidate'],
        
        #['topmas-A'],
        #['topPt-A'],
        #['ntopcandidate-A'],
        
        #['topmass-Bc'],
        #['topPt-BC'],
           # ['ntopcandidate-BC'],
        
        # ['ntopcandidate-tot'],
        # ['ntopcandidate1-tot'],
        
        
        # ['Eta_EB_el_pre'],
        # ['Eta_EE_el_pre'],
        
        ['scEta_EB_el_pre'],
        ['scEta_EE_el_pre'],
        ['Iso03_EB_el_pre'],
        ['Iso03_EE_el_pre'],
        ['dEtaIn_EB_el_pre'],
        ['dEtaIn_EE_el_pre'],
        ['dPhiIn_EB_el_pre'],
        ['dPhiIn_EE_el_pre'],
        ['Dz_EB_el_pre'],
        ['Dz_EE_el_pre'],
        ['D0_EB_el_pre'],
        ['D0_EE_el_pre'],
        ['Full5x5siee_EB_el_pre'],
        ['Full5x5siee_EE_el_pre'],
        ['HoE_EB_el_pre'],
        ['HoE_EE_el_pre'],
        ['ooEmooP_EB_el_pre'],
        ['ooEmooP_EE_el_pre'],
        ['missHits_EB_el_pre'],
        ['missHits_EE_el_pre'],
        ['conveto_EB_El_pre'],
        ['conveto_EE_el_pre'],


        ['le1E'],
        ['le2E'],
        ['zllE'],
        ['el1E'],
        ['el2E'],
        ['le1E_pre'],
        ['le2E_pre'],
        ['zllE_pre'],
        ['el1E_pre'],
        ['el2E_pre'],



        ['nle1'],
        ['nle2'],
        ['nzll1'],
        ['nel1'],
        ['sumle'],
             
        
        # ['nAK4-CatA'],
        # ['nAK4-CatB'],
        # ['nAK4-CatC'],
        # ['nAK4-CatD'],
        # ['nAK4-CatE'],
        # ['nAK4-CatF'],
        # ['nAK4-CatG'],
        # ['nAK4-CatH'],
        # ['nAK4-CatI'],
        
        # ['ST_CatB'],
        # ['ST_CatC'],
            # ['ST_CatD'],
        # ['ST_CatE'],
        # ['ST_CatF'],
        # ['ST_CatG'],
        # ['ST_CatH'],
        # ['ST_CatI'],
        # ['ST_CatJ'],
        # ['ST_CatK'],
        ['nak8'],
        
        
        
        ['cutflow'],
        ['pt_el1'], 
        ['npv_noweight_pre'],
        ['npv_pre'],
        ['nak4_pre'],
        ['ht_pre'],
        ['st_pre'],
        ['met_pre'],
        ['metPhi_pre'],
        ['ptak4jet1_pre'],
        ['etaak4jet1_pre'],
        ['cvsak4jet1_pre'],
        ['ptak4jet2_pre'],
        ['etaak4jet2_pre'],
        ['cvsak4jet2_pre'],
        ['ptak4jet3_pre'],
        ['etaak4jet3_pre'],
        ['cvsak4jet3_pre'],
        ['phi_jet1MET_pre'],
        ['mass_zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],
        
        ['npv_noweight_cnt'],
        ['npv_cnt'],
        ['nak4_cnt'],
        ['ht_cnt'],
        ['st_cnt'],
        ['met_cnt'],
        ['metPhi_cnt'],
        ['ptak4jet1_cnt'],
        ['etaak4jet1_cnt'],
        ['cvsak4jet1_cnt'],
        ['ptak4jet2_cnt'],
        ['etaak4jet2_cnt'],
        ['cvsak4jet2_cnt'],
        ['ptak4jet3_cnt'],
        ['etaak4jet3_cnt'],
        ['cvsak4jet3_cnt'],
        ['phi_jet1MET_cnt'],
        ['mass_zelel_cnt'],
        ['dr_elel_cnt'],
        ['pt_zelel_cnt'],
        ['pt_el1_cnt'],
        ['eta_el1_cnt'],
        ['pt_el2_cnt'],
        ['eta_el2_cnt'],            
        ['ht-cnt1'],
        ['st-cnt1'],
        ['met-cnt1'],

        ['htsig1'],
        ['stsig1'],
        ['metsig1'],
        ['npv_noweight'],
        ['npv'],
        ['nak4'],
        ['ht'],
        ['st'],
        ['met'],
        ['metPhi'],
        ['ptak4jet1'],
        ['etaak4jet1'],
        ['cvsak4jet1'],
        ['ptak4jet2'],
        ['etaak4jet2'],
        ['cvsak4jet2'],
        ['ptak4jet3'],
        ['etaak4jet3'],
        ['cvsak4jet3'],
        ['phi_jet1MET'],
        ['mass_zelel'],
        ['dr_elel'],
        ['pt_zelel'],
        ['pt_el'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],
        ['nbjets'],
        ['ptbjetleading'],
        ['etabjetleading'],
        ['nak8'],
        ['nwjet'],
        ['nhjet'],
        ['ntjet'],
        ['ptak8leading'],
        ['etaak8leading'],
        ['mak8leading'],
        ['prunedmak8leading'],
        ['trimmedmak8leading'],
        ['softdropmak8leading'],
        ['ptak82nd'],
        ['etamak82nd'],
        ['mak82nd'],
        ['purnedmak82nd'],
        ['trimmedmak82nd'],
        ['softdropmak82nd'],
        ['ptTprime'],
        ['yTprime'],
        ['mTprime'],
        ['ptBprime'],
        ['yBprime'],
        ['mBprime'],
        ['ZJetMasslep'],
        ['chi2_chi'],
        ['sqrtChi2'],
        ['chi_mass'],

        ]

elif dir == 'sig':
    options = [
        # ['genZ'],
        # ['genTMass'],
        # ['ZJetMass'],
        # ['hadTJetMass'],
        # ['lepTJetMass'],
        
        # ['chi_mass','4'],
        # ['ht_pre'],           
        # ['nob_pt_zelel'],           
        
        #    ['ptak4jet1_cnt'],
        ['le1E'],
        ['le2E'],
        ['zllE'],

        
        ['cutflow'],
        ['cutflow1'],
        ['cutflow2'],
        ['cutflow3'],
        ['cutflow4'],
        ['mass_zelel_pre'],
        
        ['mass_Zelel_pre'],
        ['ht_pre'],
        ['nob_ht'],
        ['1b_ht'],
        ['nob_st'],
        ['nbjets_sig'],
        ['nbjets_cnt'],
        ['nbjets_0btagcnt'],
        ['nbjets_1btagcnt'],
        ['ht_0btagcnt'],
        ['1b_ht'],
        ['ht_1btagcnt'],
        ['lowmet_ht'],
        ['st_0btagcnt'],
        ['1b_st'],
        ['st_1btagcnt'],
        ['lowmet_st'],
        ['mass_zmumu_pre'],
        ['dr_mumu_pre'],
        ['pt_zmumu_pre'],
        ['pt_mu1_pre'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],
        ['Iso04_mu_pre'],
        ['Dz_mu_pre'],
        ['D0_mu_pre'],
        ['mass_zmumu_cnt'],
        ['dr_mumu_cnt'],
        ['pt_zmumu_cnt'],
        ['pt_mu1_cnt'],
        ['eta_mu1_cnt'],
        ['pt_mu2_cnt'],
        ['eta_mu2_cnt'],             
        
        ['nob_ht'],
        # ['cutflow'],
           #  ['cutflow1'],
        
        # ['ST_sig'],       
        # ['ST_sig1b'],
        # ['ST_sig2b'],
        # ['ST_sigT1Z1'],
        # ['ST_sigT0Z1'],
        # ['ST_sigT1Z0'],
        # ['ST_sigT0Z0'],
        
        # ['ST_sigT1ZH1'],
        # ['ST_sigT0ZH1'],
        # ['ST_sigT1ZH0'],
        # ['ST_sigT0ZH0'],
        
        #  ['ST_sigT1Z1H1'],
        #  ['ST_sigT1Z1H0'],
        #  ['ST_sigT0Z1H1'],
        #  ['ST_sigT0Z1H0'],
        #  ['ST_sigT1Z0H1'],
        #  ['ST_sigT1Z0H0'],
        #  ['ST_sigT0Z0H1'],
        #  ['ST_sigT0Z0H0'],
        
        # ['ST_sigT1Z1H1b1'],
        # ['ST_sigT1Z1H1b2'],
        # ['ST_sigT1Z1H0b1'],
        # ['ST_sigT1Z1H0b2'],
        # ['ST_sigT0Z1H1b1'],
        # ['ST_sigT0Z1H1b2'],
        # ['ST_sigT0Z1H0b1'],
        # ['ST_sigT0Z1H0b2'],
        # ['ST_sigT1Z0H1b1'],
        # ['ST_sigT1Z0H1b2'],
        # ['ST_sigT1Z0H0b1'],
        # ['ST_sigT1Z0H0b2'],
        # ['ST_sigT0Z0H1b1'],
        # ['ST_sigT0Z0H1b2'],
        # ['ST_sigT0Z0H0b1'],
        # ['ST_sigT0Z0H0b2'],
        
        
        # ['ST_sigT1Z1b1'],
        # ['ST_sigT1Z1b2'],
            # ['ST_sigT0Z1b1'],
        # ['ST_sigT0Z1b2'],
        # ['ST_sigT1Z0b1'],
        # ['ST_sigT1Z0b2'],
        # ['ST_sigT0Z0b1'],
        # ['ST_sigT0Z0b2'],
        
        
        # ['ST_1A'],
        #  ['ST_1B'],
        #  ['ST_2A'],
        #  ['ST_2B'],
        #  ['ST_3A'],
        #  ['ST_3B'],
        #  ['ST_4A'],
        #  ['ST_4B'],
        #  ['cutflow3'],
        # ['cutflow1'],
        #  ['ak4Jet1Pt'],            
        # ['ak4Jet2Pt'],
        # ['ak4Jet3Pt'],
        # ['ak4Jet4Pt'],
        #['ak4Jet5Pt'],
        #['ak4Jet6Pt'],
        #['nak4sig'],
        
        
        # ['cutflow'],
        #['b_st'],
        #['nZll'],
        #['nAK4-CatA'],
        #['nAK4-CatB'],
        #['nAK4-CatC'],
        #  ['ZHmass-boosted'],                                                                                                                                           
        # ['ZHPt-boosted'],                                                                                                                                             
        # ['nZHcandidate-boosted'],                                                                                                                                     
        # ['ZHmassnb'],                                                                                                                                                 
        # ['ZHPtnb'],                                                                                                                                                   
        # ['nZHcandidatesnb'],                                                                                                                                          
        #['nZHcandidates-tot'],  
        # ['nZHcandidates1-tot'],

        
        ['Hmass-boosted'],
        ['HPt-boosted'],
        ['nHcandidate-boosted'],
        ['Hmassnb'],
        ['HPtnb'],
        ['nHcandidatesnb'],
        ['nHcandidates-tot'],                                                                                                                                             
        ['nHcandidates1-tot'],
        
        
        
        
        # ['Zmass-boosted'],
        # ['ZPt-boosted'],
        # ['nzcandidate-boosted'],
        
        #  ['Zmass'],
        # ['ZPt'],
        # ['nzcandidates'],
        #  ['nzcandidates-tot'],
        #  ['nzcandidates1-tot'], 
        
        
        
        
        #  ['topmass-D'],
        # ['topPt-D'],
        # ['ntopcandidate-D'],

        
        # ['Wmass-BC'],
        # ['nWcandidate-BC'],
        # ['lightjetmass-BC'],
        # ['nlightjetcandidate'],
        
        #['topmas-A'],
        #['topPt-A'],
        #['ntopcandidate-A'],
        
        #['topmass-Bc'],
        #['topPt-BC'],
           # ['ntopcandidate-BC'],
        
        # ['ntopcandidate-tot'],
        # ['ntopcandidate1-tot'],
        
        
        # ['Eta_EB_el_pre'],
        # ['Eta_EE_el_pre'],
        
        ['scEta_EB_el_pre'],
        ['scEta_EE_el_pre'],
        ['Iso03_EB_el_pre'],
        ['Iso03_EE_el_pre'],
        ['dEtaIn_EB_el_pre'],
        ['dEtaIn_EE_el_pre'],
        ['dPhiIn_EB_el_pre'],
        ['dPhiIn_EE_el_pre'],
        ['Dz_EB_el_pre'],
        ['Dz_EE_el_pre'],
        ['D0_EB_el_pre'],
        ['D0_EE_el_pre'],
        ['le1E'],
        ['le2E'],
        ['zllE'],
        ['el1E'],
        ['el2E'],
        ['nle1'],
        ['nle2'],
        ['nzll1'],
        ['nel1'],
        ['sumle'],
             
        
        # ['nAK4-CatA'],
        # ['nAK4-CatB'],
        # ['nAK4-CatC'],
        # ['nAK4-CatD'],
        # ['nAK4-CatE'],
        # ['nAK4-CatF'],
        # ['nAK4-CatG'],
        # ['nAK4-CatH'],
        # ['nAK4-CatI'],
        
        # ['ST_CatB'],
        # ['ST_CatC'],
            # ['ST_CatD'],
        # ['ST_CatE'],
        # ['ST_CatF'],
        # ['ST_CatG'],
        # ['ST_CatH'],
        # ['ST_CatI'],
        # ['ST_CatJ'],
        # ['ST_CatK'],
        ['nak8'],
        
        
        
        ['cutflow'],
        ['pt_el1'], 
        ['npv_noweight_pre'],
        ['npv_pre'],
        ['nak4_pre'],
        ['ht_pre'],
        ['st_pre'],
        ['met_pre'],
        ['metPhi_pre'],
        ['ptak4jet1_pre'],
        ['etaak4jet1_pre'],
        ['cvsak4jet1_pre'],
        ['ptak4jet2_pre'],
        ['etaak4jet2_pre'],
        ['cvsak4jet2_pre'],
        ['ptak4jet3_pre'],
        ['etaak4jet3_pre'],
        ['cvsak4jet3_pre'],
        ['phi_jet1MET_pre'],
        ['mass_zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],
        
        ['npv_noweight_cnt'],
        ['npv_cnt'],
        ['nak4_cnt'],
        ['ht_cnt'],
        ['st_cnt'],
        ['met_cnt'],
        ['metPhi_cnt'],
        ['ptak4jet1_cnt'],
        ['etaak4jet1_cnt'],
        ['cvsak4jet1_cnt'],
        ['ptak4jet2_cnt'],
        ['etaak4jet2_cnt'],
        ['cvsak4jet2_cnt'],
        ['ptak4jet3_cnt'],
        ['etaak4jet3_cnt'],
        ['cvsak4jet3_cnt'],
        ['phi_jet1MET_cnt'],
        ['mass_zelel_cnt'],
        ['dr_elel_cnt'],
        ['pt_zelel_cnt'],
        ['pt_el1_cnt'],
        ['eta_el1_cnt'],
        ['pt_el2_cnt'],
        ['eta_el2_cnt'],            
        ['ht-cnt1'],
        ['st-cnt1'],
        ['met-cnt1'],

        ['htsig1'],
        ['stsig1'],
        ['metsig1'],
        ['npv_noweight'],
        ['npv'],
        ['nak4'],
        ['ht'],
        ['st'],
        ['met'],
        ['metPhi'],
        ['ptak4jet1'],
        ['etaak4jet1'],
        ['cvsak4jet1'],
        ['ptak4jet2'],
        ['etaak4jet2'],
        ['cvsak4jet2'],
        ['ptak4jet3'],
        ['etaak4jet3'],
        ['cvsak4jet3'],
        ['phi_jet1MET'],
        ['mass_zelel'],
        ['dr_elel'],
        ['pt_zelel'],
        ['pt_el'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],
        ['nbjets'],
        ['ptbjetleading'],
        ['etabjetleading'],
        ['nak8'],
        ['nwjet'],
        ['nhjet'],
        ['ntjet'],
        ['ptak8leading'],
        ['etaak8leading'],
        ['mak8leading'],
        ['prunedmak8leading'],
        ['trimmedmak8leading'],
        ['softdropmak8leading'],
        ['ptak82nd'],
        ['etamak82nd'],
        ['mak82nd'],
        ['purnedmak82nd'],
        ['trimmedmak82nd'],
        ['softdropmak82nd'],
        ['ptTprime'],
        ['yTprime'],
        ['mTprime'],
        ['ptBprime'],
        ['yBprime'],
        ['mBprime'],
        ['ZJetMasslep'],
        ['chi2_chi'],
        ['sqrtChi2'],
        ['chi_mass'],

        ]
elif dir == 'cat':
    options = [

          ['cutflow4'],
          ['cutflow1'],
          ['cutflow2'],
          ['cutflow3'],
          ['cutflow4'],

        # ['ST_sig'],       
        # ['ST_sig1b'],
        # ['ST_sig2b'],
        # ['ST_sigT1Z1'],
        # ['ST_sigT0Z1'],
        # ['ST_sigT1Z0'],
        # ['ST_sigT0Z0'],
        

        #  ['ST_sigT1Z1H1'],
        #  ['ST_sigT1Z1H0'],
        #  ['ST_sigT0Z1H1'],
        #  ['ST_sigT0Z1H0'],
        #  ['ST_sigT1Z0H1'],
        #  ['ST_sigT1Z0H0'],
        #  ['ST_sigT0Z0H1'],
        #  ['ST_sigT0Z0H0'],
        
        # ['ST_sigT1Z1H1b1'],
        # ['ST_sigT1Z1H1b2'],
        # ['ST_sigT1Z1H0b1'],
        # ['ST_sigT1Z1H0b2'],
        # ['ST_sigT0Z1H1b1'],
        # ['ST_sigT0Z1H1b2'],
        # ['ST_sigT0Z1H0b1'],
        # ['ST_sigT0Z1H0b2'],
        # ['ST_sigT1Z0H1b1'],
        # ['ST_sigT1Z0H1b2'],
        # ['ST_sigT1Z0H0b1'],
        # ['ST_sigT1Z0H0b2'],
        # ['ST_sigT0Z0H1b1'],
        # ['ST_sigT0Z0H1b2'],
        # ['ST_sigT0Z0H0b1'],
        # ['ST_sigT0Z0H0b2'],
        
        
        # ['ST_sigT1Z1b1'],
        # ['ST_sigT1Z1b2'],
        # ['ST_sigT0Z1b1'],
        # ['ST_sigT0Z1b2'],
        # ['ST_sigT1Z0b1'],
        # ['ST_sigT1Z0b2'],
        # ['ST_sigT0Z0b1'],
        # ['ST_sigT0Z0b2'],
        #  ['ZHmass-boosted'],                                                                                                                                           
        # ['ZHPt-boosted'],                                                                                                                                             
        # ['nZHcandidate-boosted'],                                                                                                                                     
        # ['ZHmassnb'],                                                                                                                                                 
        # ['ZHPtnb'],                                                                                                                                                   
        # ['nZHcandidatesnb'],                                                                                                                                          
        #['nZHcandidates-tot'],  
        # ['nZHcandidates1-tot'],

        
        ['Hmass-boosted-cnt'],
        ['HPt-boosted-cnt'],
        ['nHcandidate-boosted-cnt'],
        ['Hmassnb-cnt'],
        ['HPtnb-cnt'],
        ['nHcandidatesnb-cnt'],
        ['nHcandidates-tot-cnt'],                                                                                                                                             
        ['nHcandidates1-tot-cnt'],
        
        
        
        
        ['Zmass-boosted-cnt'],
        ['ZPt-boosted-cnt'],
        ['nzcandidate-boosted-cnt'],
        
        ['Zmass-cnt'],
        ['ZPt-cnt'],
        ['nzcandidates-cnt'],
        ['nzcandidates-tot-cnt'],
        ['nzcandidates1-tot-cnt'], 
   
        ['topmass-D-cnt'],
        ['topPt-D-cnt'],
        ['ntopcandidate-D-cnt'],
        
        
        ['Wmass-BC-cnt'],
        ['nWcandidate-BC-cnt'],
        ['lightjetmass-BC-cnt'],
        ['nlightjetcandidate-cnt'],
        
        ['topmas-A-cnt'],
        ['topPt-A-cnt'],
        ['ntopcandidate-A-cnt'],
        
        ['topmass-Bc-cnt'],
        ['topPt-BC-cnt'],
        ['ntopcandidate-BC-cnt'],
        
        ['ntopcandidate-tot-cnt'],
        ['ntopcandidate1-tot-cnt'],
        
        
        ['Hmass-boosted-sig'],
        ['HPt-boosted-sig'],
        ['nHcandidate-boosted-sig'],
        ['Hmassnb-sig'],
        ['HPtnb-sig'],
        ['nHcandidatesnb-sig'],
        ['nHcandidates-tot-sig'],                                                                                                                                       

        ['nHcandidates1-tot-sig'],




        ['Zmass-boosted-sig'],
        ['ZPt-boosted-sig'],
        ['nzcandidate-boosted-sig'],

        ['Zmass-sig'],
        ['ZPt-sig'],
        ['nzcandidates-sig'],
        ['nzcandidates-tot-sig'],
        ['nzcandidates1-tot-sig'],

        ['topmass-D-sig'],
        ['topPt-D-sig'],
        ['ntopcandidate-D-sig'],


        ['Wmass-BC-sig'],
        ['nWcandidate-BC-sig'],
        ['lightjetmass-BC-sig'],
        ['nlightjetcandidate-sig'],

        ['topmas-A-sig'],
        ['topPt-A-sig'],
        ['ntopcandidate-A-sig'],

        ['topmass-Bc-sig'],
        ['topPt-BC-sig'],
        ['ntopcandidate-BC-sig'],

        ['ntopcandidate-tot-sig'],
        ['ntopcandidate1-tot-sig'],

             

        ]


command = 'python plot.py --var={0:s}'

for option in options :
    s = command.format(
        option[0]
        )

    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo %s"%s,""]                                                                      , shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( [s, ""]                                                                               , shell=True)
