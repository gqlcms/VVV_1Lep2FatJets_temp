#ifndef Process_1Lep2fatjet_OldLeptonID_h
#define Process_1Lep2fatjet_OldLeptonID_h


















void Process_1Lep2fatjet_Lepton(){
    // cout << "num of muon"  <<nt.Muon_p4().size()<<endl; 
    int OldID_nLooseMu = 0;
    int OldID_nLooseEle = 0;
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu)
    {
        int temp_id = nt.Muon_highPtId()[imu]; // have to switch to int
        // cout<<"nt.Muon_highPtId()[imu]:"<<temp_id<<" id origin:  "<<nt.Muon_highPtId()[imu]<<endl;
        if (not (temp_id==2 )) continue;
        if (not (nt.Muon_p4()[imu].pt()>20.0 )) continue;
        if (not (fabs(nt.Muon_p4()[imu].eta())<2.4 )) continue;
        if (not (nt.Muon_tkRelIso()[imu]<0.1 )) continue;
        OldID_nLooseMu++;
        // std::cout<<"find one highPt muon"<<std::endl;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_LooseMuon_highPtId_idxs",imu);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_LooseMuon_highPtId_idxs",imu);
    }
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu)
    {
        int temp_id = nt.Muon_highPtId()[imu]; // have to switch to int
        // cout<<"nt.Muon_highPtId()[imu]:"<<temp_id<<" id origin:  "<<nt.Muon_highPtId()[imu]<<endl;
        if (not (temp_id==2 )) continue;
        if (not (nt.Muon_p4()[imu].pt()>55.0 )) continue;
        if (not (fabs(nt.Muon_p4()[imu].eta())<2.4 )) continue;
        if (not (nt.Muon_tkRelIso()[imu]<0.1 )) continue;
        // std::cout<<"find one highPt muon"<<std::endl;
        //ana.tx.pushbackToBranch<int>("1Lep2fatjet_CoodMuon_highPtId_idxs",imu);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_GoodMuon_highPtId_idxs",imu);
    }
    for (unsigned int iel = 0; iel < nt.Electron_p4().size(); ++iel)
    {
        if (not (nt.Electron_cutBased_HEEP()[iel])) continue;
        if (not (nt.Electron_p4()[iel].pt()>35.0 )) continue;
        if (not (fabs(nt.Electron_p4()[iel].eta())<2.5 )) continue;
        if (fabs(nt.Electron_p4().at(iel).eta())<1.556&&fabs(nt.Electron_p4().at(iel).eta())>1.442) continue;
        OldID_nLooseEle++;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs",iel);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs",iel);
    }
    for (unsigned int iel = 0; iel < nt.Electron_p4().size(); ++iel)
    {
        if (not (nt.Electron_cutBased_HEEP()[iel])) continue;
        if (not (nt.Electron_p4()[iel].pt()>55.0 )) continue;
        if (not (fabs(nt.Electron_p4()[iel].eta())<2.5 )) continue;
        if (fabs(nt.Electron_p4().at(iel).eta())<1.556&&fabs(nt.Electron_p4().at(iel).eta())>1.442) continue;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs",iel);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs",iel);
    }
    ana.tx.setBranch<int>	("Process_1Lep2fatjet_nlooseMu_Old",OldID_nLooseMu);
    ana.tx.setBranch<int>	("Process_1Lep2fatjet_nlooseEle_oldID",OldID_nLooseEle);
    

}



void Process_1Lep2fatjet_Fatjet(){
    //float fjSFvlc(1.), fjSFvlu(1.), fjSFvld(1.), fjSFlc(1.), fjSFlu(1.), fjSFld(1.), fjSFmc(1.), fjSFmu(1.), fjSFmd(1.), fjSFtc(1.), fjSFtu(1.), fjSFtd(1.);
    for (unsigned int usenumber1 = 0; usenumber1 < nt.FatJet_p4().size(); ++usenumber1)
    {
      float fjWPvloose = 0.274; //https://twiki.cern.ch/twiki/bin/view/CMS/DeepAK8Tagging2018WPsSFs
      float fjWPloose  = 0.506;
      float fjWPmedium = 0.731;
      float fjWPtight  = 0.828;
      if(nt.year() == 2017){
        fjWPvloose = 0.258;
        fjWPloose  = 0.506;
        fjWPmedium = 0.739;
        fjWPtight  = 0.838;
      }
      if(nt.year() == 2018){
        fjWPvloose = 0.245;
        fjWPloose  = 0.479;
        fjWPmedium = 0.704;
        fjWPtight  = 0.806;
      }
      cout<<fjWPvloose<<fjWPloose<<fjWPmedium<<fjWPtight<<endl;

        // TODO: What is POG recommendation? do we use nt.FatJet_jetId()?
        // Figure this out
        // For now, accept anything above 250 GeV (TODO: is 250 GeV also ok?)
        float minFatJetPt = 0;
        if (not nt.isData()){
            minFatJetPt = 180;  // relax the pt to 180; apply JER later;
        }
        else{
            minFatJetPt = 200;
        }
        if (not (nt.FatJet_p4()[usenumber1].pt() > minFatJetPt)) 
            continue;
        if (not (abs(nt.FatJet_p4()[usenumber1].eta()) < 2.4))
            continue;

        int jetpuid = nt.FatJet_jetId()[usenumber1]&2;
        if (not (jetpuid ==2))
            continue;


        // Because every muon and electron shows up in PF FatJet collections
        // Need to check against leptons
        bool is_overlapping_with_a_lepton = false;

        // Overlap check against leptons (electrons)
        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.FatJet_p4()[usenumber1], nt.Electron_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"ele clean jet pt"<<nt.FatJet_p4()[usenumber1].pt()<<"ele pt "<<nt.Electron_p4()[ilep_idx].pt()<<endl;
                    break;
                }
        }

        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodMuon_highPtId_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodMuon_highPtId_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.FatJet_p4()[usenumber1], nt.Muon_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"mu clean jet pt"<<nt.FatJet_p4()[usenumber1].pt()<<"  pt "<<nt.Muon_p4()[ilep_idx].pt()<<endl;
                    break;
                }
        }

        if (is_overlapping_with_a_lepton)
            continue;

 
        // For now, accept anything that reaches this point

    ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_fatjet_idxs_new", usenumber1);
    ana.tx.pushbackToBranch<LorentzVector>("t_1Lep2fatjet_fatjet_p4_new", nt.FatJet_p4()[usenumber1]);
    }

    

    // leading pt jet
    int usenumber3 = -1; double pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].eta())<2.4 && inum<4) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt(); 
            usenumber3 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]; 
            continue;
        }
    }
    if (usenumber3>-1) {
        if (not nt.isData()){
    TLorentzVector JERcorrectedP4,JERuncorrectedP4;
    JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[usenumber3].pt(),nt.FatJet_p4()[usenumber3].eta(),nt.FatJet_p4()[usenumber3].phi(),nt.FatJet_p4()[usenumber3].energy());
    JERcorrectedP4 = Process_1Lep2fatjet_JER(usenumber3) * JERuncorrectedP4;

    ana.tx.setBranch<float>("jetAK8puppi_ptJECJER",JERcorrectedP4.Pt());
}
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC", nt.FatJet_p4()[usenumber3].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta", nt.FatJet_p4()[usenumber3].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi", nt.FatJet_p4()[usenumber3].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid", nt.FatJet_jetId()[usenumber3]&2);
    
    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr", nt.FatJet_msoftdrop()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_mass", nt.FatJet_mass()[usenumber3]);
    // cout<<"sd_sub"<<Process_1Lep2fatjet_Fatjet_SD_sub(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3])<<endl;
    // ana.tx.setBranch<float>("jetAK8puppi_sd_sub", Process_1Lep2fatjet_Fatjet_SD_sub(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3]));
    ana.tx.setBranch<float>("jetAK8puppi_sd", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3]));
    ana.tx.setBranch<float>("Jet_rawFactor", nt.Jet_rawFactor()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW", nt.FatJet_deepTagMD_WvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW", nt.FatJet_deepTag_WvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ", nt.FatJet_deepTagMD_ZvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ", nt.FatJet_deepTag_ZvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop", nt.FatJet_deepTagMD_TvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop", nt.FatJet_deepTag_TvsQCD()[usenumber3]);
    // ana.tx.pushbackToBranch<float>("Common_fatjet_deepMD_bb", nt.FatJet_deepTagMD_bbvsLight()[usenumber3]);
    }
    
    // second pt jet
    int usenumber2 = -1; 
    pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum] != usenumber3) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt(); 
            usenumber2 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]; 
            continue;
        }
    }
    if (usenumber2>-1) {

    
    if (not nt.isData()){
        TLorentzVector JERcorrectedP4,JERuncorrectedP4;
        JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[usenumber2].pt(),nt.FatJet_p4()[usenumber2].eta(),nt.FatJet_p4()[usenumber2].phi(),nt.FatJet_p4()[usenumber2].energy());
        JERcorrectedP4 = Process_1Lep2fatjet_JER(usenumber2) * JERuncorrectedP4;

        ana.tx.setBranch<float>("jetAK8puppi_ptJECJER_2",JERcorrectedP4.Pt());
    }
    
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_2", nt.FatJet_p4()[usenumber2].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_2", nt.FatJet_p4()[usenumber2].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_2", nt.FatJet_p4()[usenumber2].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_2", nt.FatJet_jetId()[usenumber2]&2);

    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_2", nt.FatJet_msoftdrop()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_sd_2", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber2],nt.FatJet_subJetIdx2()[usenumber2]));
    ana.tx.setBranch<float>("Jet_rawFactor_2", nt.Jet_rawFactor()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_2", nt.FatJet_deepTagMD_WvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_2", nt.FatJet_deepTag_WvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_2", nt.FatJet_deepTagMD_ZvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_2", nt.FatJet_deepTag_ZvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_2", nt.FatJet_deepTagMD_TvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_2", nt.FatJet_deepTag_TvsQCD()[usenumber2]);
    }

    // third pt jet
    int usenumber1 = -1; 
    pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum] != usenumber3 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum] != usenumber2) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt(); 
            usenumber1 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]; 
            continue;
        }
    }
    if (usenumber1>-1) {

    if (not nt.isData()){
        TLorentzVector JERcorrectedP4,JERuncorrectedP4;
        JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[usenumber1].pt(),nt.FatJet_p4()[usenumber1].eta(),nt.FatJet_p4()[usenumber1].phi(),nt.FatJet_p4()[usenumber1].energy());
        JERcorrectedP4 = Process_1Lep2fatjet_JER(usenumber1) * JERuncorrectedP4;

        ana.tx.setBranch<float>("jetAK8puppi_ptJECJER_3",JERcorrectedP4.Pt());
    }

    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_3", nt.FatJet_p4()[usenumber1].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_3", nt.FatJet_p4()[usenumber1].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_3", nt.FatJet_p4()[usenumber1].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_3", nt.FatJet_jetId()[usenumber1]&2);

    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_3", nt.FatJet_msoftdrop()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_sd_3", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber1],nt.FatJet_subJetIdx2()[usenumber1]));
    ana.tx.setBranch<float>("Jet_rawFactor_3", nt.Jet_rawFactor()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_3", nt.FatJet_deepTagMD_WvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_3", nt.FatJet_deepTag_WvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_3", nt.FatJet_deepTagMD_ZvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_3", nt.FatJet_deepTag_ZvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_3", nt.FatJet_deepTagMD_TvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_3", nt.FatJet_deepTag_TvsQCD()[usenumber1]);
    }

    int numofAK8 = 0;
    if(ana.tx.getBranchLazy<float>("jetAK8puppi_ptJECJER")>200){
        numofAK8++;
    }
    if(ana.tx.getBranchLazy<float>("jetAK8puppi_ptJECJER_2")>200){
        numofAK8++;
    }
    if(ana.tx.getBranchLazy<float>("jetAK8puppi_ptJECJER_3")>0){
        numofAK8++;
    }
    ana.tx.setBranch<int>("t_1Lep2fatjet_Nj8", numofAK8);

}

void Process_1Lep2fatjet_jet(){
for (unsigned int ijet = 0; ijet < nt.Jet_p4().size(); ++ijet)
    {
      
        //  do we use nt.Jet_jetId()?
        //   Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto

        if ( not (nt.Jet_p4()[ijet].pt() > 20.) )
                continue;
        int _Jet_jetId = nt.Jet_jetId()[ijet];
        if (_Jet_jetId &2 !=2)
            continue;
        if (not (abs(nt.Jet_p4()[ijet].eta()) < 5.0))
            continue;

        // Because every muon and electron shows up in PF FatJet collections
        // Need to check against leptons
        bool is_overlapping_with_a_lepton = false;

        // Overlap check against leptons (electrons)
        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.Jet_p4()[ijet], nt.Electron_p4()[ilep_idx]) < 0.4)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"ele clean jet pt"<<nt.Jet_p4()[ijet].pt()<<"ele pt "<<nt.Electron_p4()[ilep_idx].pt()<<"ele ID"<<nt.Electron_mvaFall17V2Iso_WP90()[ilep_idx]<<endl;
                    break;
                }
        }

        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodMuon_highPtId_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_GoodMuon_highPtId_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.Jet_p4()[ijet], nt.Muon_p4()[ilep_idx]) < 0.4)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"mu clean jet pt"<<nt.Jet_p4()[ijet].pt()<<"  pt "<<nt.Muon_p4()[ilep_idx].pt()<<"Muon ID"<<nt.Muon_mediumId()[ilep_idx]<<"Mu Iso"<<nt.Muon_pfRelIso04_all()[ilep_idx]<<endl;
                    break;
                }
        }

        if (is_overlapping_with_a_lepton)
            continue;

 
        // For now, accept anything that reaches this point

    ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_jet_idxs", ijet);
    ana.tx.pushbackToBranch<LorentzVector>("Common_jet_p4_new", nt.Jet_p4()[ijet]);
    if (not nt.isData()){
    ana.tx.pushbackToBranch<int>("ak4jet_hf", nt.Jet_hadronFlavour()[ijet]);
    ana.tx.pushbackToBranch<int>("ak4jet_pf", nt.Jet_partonFlavour()[ijet]);
    }
    ana.tx.pushbackToBranch<float>("ak4jet_pt", nt.Jet_p4()[ijet].pt());
    ana.tx.pushbackToBranch<float>("ak4jet_eta", nt.Jet_p4()[ijet].eta());
    ana.tx.pushbackToBranch<float>("ak4jet_phi", nt.Jet_p4()[ijet].phi());
    ana.tx.pushbackToBranch<float>("ak4jet_e", nt.Jet_p4()[ijet].energy());
    ana.tx.pushbackToBranch<float>("Jet_btagDeepB", nt.Jet_btagDeepB()[ijet]);
    ana.tx.pushbackToBranch<float>("Jet_btagDeepC", nt.Jet_btagDeepC()[ijet]);


    if (not nt.isData()){
    TLorentzVector JERcorrectedP4,JERuncorrectedP4;
    JERuncorrectedP4.SetPtEtaPhiE(nt.Jet_p4()[ijet].pt(),nt.Jet_p4()[ijet].eta(),nt.Jet_p4()[ijet].phi(),nt.Jet_p4()[ijet].energy());
    JERcorrectedP4 = Process_1Lep2fatjet_JER_AK4(ijet) * JERuncorrectedP4;
    ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_JER", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()));
    ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_unJER", LorentzVector(JERuncorrectedP4.Pt(),JERuncorrectedP4.Eta(),JERuncorrectedP4.Phi(),JERuncorrectedP4.M()));
    ana.tx.pushbackToBranch<float>("ak4jet_pt_JER",JERcorrectedP4.Pt() );
    ana.tx.pushbackToBranch<float>("ak4jet_e_JER",JERcorrectedP4.Energy() );
    }

    }
}

void Process_1Lep2fatjet_MET_Lepton_leptonicW(){



    // lepton
    
    float ptlep1, etalep1, philep1, energylep1;
    int lep;
    float RelIso04=-99.;

    
    if ( ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs").size() == 1 ){
        ptlep1 = nt.Electron_p4()[ ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].pt();
        etalep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].eta();
        philep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].phi();
        energylep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].energy();
        lep = nt.Electron_pdgId()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]];
    }
    if ( ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs").size() == 1 ){
        ptlep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].pt();
        etalep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].eta();
        philep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].phi();
        energylep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].energy();
        lep = nt.Muon_pdgId()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]];
        RelIso04 = nt.Muon_pfRelIso04_all()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]];
    }
    ana.tx.setBranch<float>("ptlep1", ptlep1 );
    ana.tx.setBranch<float>("etalep1", etalep1 );
    ana.tx.setBranch<float>("philep1", philep1 );
    ana.tx.setBranch<float>("energylep1", energylep1 );
    ana.tx.setBranch<int>("lep", lep );
    ana.tx.setBranch<float>("RelIso04", RelIso04 );



    // leptonic w
    TLorentzVector  glepton,neutrino,neutrinoP4,WLeptonic;
    glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
    int leptontype = 1;
    // double MET_et_JER = ana.tx.getBranchLazy<float>("MET_et");
    double MET_et_JER = 0;

    if (not nt.isData()){
    MET_et_JER = ana.tx.getBranchLazy<float>("MET_et_JER");}
    else{
        MET_et_JER = ana.tx.getBranchLazy<float>("MET_et");
    }
    double MET_phi = 0;
    MET_phi = ana.tx.getBranchLazy<float>("MET_phi");
    neutrino = getNeutrinoP4(MET_et_JER , MET_phi, glepton, leptontype);
    // neutrino = getNeutrinoP4(MET_et, MET_phi, glepton, leptontype);
    neutrinoP4.SetPtEtaPhiE(neutrino.Pt(),neutrino.Eta(),neutrino.Phi(),neutrino.Energy());
    WLeptonic = glepton+neutrinoP4;
    float ptVlepJEC, yVlepJEC, phiVlepJEC,massVlepJEC,mtVlepJEC;
    ptVlepJEC    = WLeptonic.Pt();
    yVlepJEC     = WLeptonic.Eta();
    phiVlepJEC   = WLeptonic.Phi();
    massVlepJEC  = WLeptonic.M();
    mtVlepJEC    = WLeptonic.Mt();
    ana.tx.setBranch<float>("ptVlepJEC", ptVlepJEC );
    ana.tx.setBranch<float>("yVlepJEC", yVlepJEC );
    ana.tx.setBranch<float>("phiVlepJEC", phiVlepJEC );
    ana.tx.setBranch<float>("massVlepJEC", massVlepJEC );
    ana.tx.setBranch<float>("mtVlepJEC", mtVlepJEC );

    // neutrino
    float ptlep2, etalep2, philep2, energylep2;
    ptlep2 = neutrinoP4.Pt();
    etalep2 = neutrinoP4.Eta();
    philep2 = neutrinoP4.Phi();
    energylep2 = neutrinoP4.Energy();
    ana.tx.setBranch<float>("ptlep2", ptlep2 );
    ana.tx.setBranch<float>("etalep2", etalep2 );
    ana.tx.setBranch<float>("philep2", philep2 );
    ana.tx.setBranch<float>("energylep2", energylep2 );
}




#endif