#ifndef Process_1Lep2fatjet_NewLeptonID_h
#define Process_1Lep2fatjet_NewLeptonID_h








void Process_1Lep2fatjet_Fatjet_NewID(){
    //float fjSFvlc(1.), fjSFvlu(1.), fjSFvld(1.), fjSFlc(1.), fjSFlu(1.), fjSFld(1.), fjSFmc(1.), fjSFmu(1.), fjSFmd(1.), fjSFtc(1.), fjSFtu(1.), fjSFtd(1.);
    for (unsigned int usenumber1 = 0; usenumber1 < nt.FatJet_p4().size(); ++usenumber1)
    {
      

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
        // cout<< "pt" <<nt.FatJet_p4()[usenumber1].pt()<<"eta"<<abs(nt.FatJet_p4()[usenumber1].eta())<<"ID"<<nt.FatJet_jetId()[usenumber1]<<endl;
        // cout<< "pt" <<nt.FatJet_p4()[usenumber1].pt()<<"eta"<<abs(nt.FatJet_p4()[usenumber1].eta())<<"ID"<<(nt.FatJet_jetId()[usenumber1]&2)<<endl;
        // cout<< "pt" <<(nt.FatJet_p4()[usenumber1].pt() > minFatJetPt)<<"eta"<<abs(nt.FatJet_p4()[usenumber1].eta())<<"ID"<<(nt.FatJet_jetId()[usenumber1]&2==2)<<endl;
        int jetpuid = nt.FatJet_jetId()[usenumber1]&2;
        // cout<< "pt" <<(nt.FatJet_p4()[usenumber1].pt() > minFatJetPt)<<"eta"<<(abs(nt.FatJet_p4()[usenumber1].eta())<2.4)<<"ID"<<(jetpuid==2)<<endl;
        if (not (nt.FatJet_p4()[usenumber1].pt() > minFatJetPt)) 
            continue;
        if (not (abs(nt.FatJet_p4()[usenumber1].eta()) < 2.4))
            continue;
        // if (not (nt.FatJet_jetId()[usenumber1]&2 ==2)) # met problem for aQGC sample
        if (not (jetpuid ==2))
            continue;


        // Because every muon and electron shows up in PF FatJet collections
        // Need to check against leptons
       bool is_overlapping_with_a_lepton = false;

        // Overlap check against leptons (electrons)
        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs")[ilep];
            // If electron
            if (abs(ana.tx.getBranchLazy<vector<int>>("Common_lep_pdgid")[ilep]) == 11)
            {
                if (RooUtil::Calc::DeltaR(nt.FatJet_p4()[usenumber1], nt.Electron_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    break;
                }
            }
            // else muon
            else
            {
                if (RooUtil::Calc::DeltaR(nt.FatJet_p4()[usenumber1], nt.Muon_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    break;
                }
            }
        }

        // cout<<"is_overlapping_with_a_lepton"<<endl;
        if (is_overlapping_with_a_lepton)
            continue;

 
        // For now, accept anything that reaches this point

    ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_fatjet_idxs_new_NewID", usenumber1);
    ana.tx.pushbackToBranch<LorentzVector>("t_1Lep2fatjet_fatjet_p4_new_NewID", nt.FatJet_p4()[usenumber1]);
    }

    // leading pt jet
    int usenumber3 = -1; double pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].eta())<2.4 && inum<4) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].pt(); 
            usenumber3 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]; 
            continue;
        }
    }
    cout <<"usenumber3"<<usenumber3<<endl;
    cout<<"njets"<<ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID").size()<<endl;
    if (usenumber3>-1) {
        if (not nt.isData()){
    TLorentzVector JERcorrectedP4,JERuncorrectedP4;
    JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[usenumber3].pt(),nt.FatJet_p4()[usenumber3].eta(),nt.FatJet_p4()[usenumber3].phi(),nt.FatJet_p4()[usenumber3].energy());
    JERcorrectedP4 = Process_1Lep2fatjet_JER(usenumber3) * JERuncorrectedP4;

    ana.tx.setBranch<float>("jetAK8puppi_ptJECJER_NewID",JERcorrectedP4.Pt());
}
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_NewID", nt.FatJet_p4()[usenumber3].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_NewID", nt.FatJet_p4()[usenumber3].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_NewID", nt.FatJet_p4()[usenumber3].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_NewID", nt.FatJet_jetId()[usenumber3]&2);
    
    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_NewID", nt.FatJet_msoftdrop()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_mass_NewID", nt.FatJet_mass()[usenumber3]);
    // cout<<"sd_sub"<<Process_1Lep2fatjet_Fatjet_SD_sub(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3])<<endl;
    // ana.tx.setBranch<float>("jetAK8puppi_sd_sub", Process_1Lep2fatjet_Fatjet_SD_sub(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3]));
    ana.tx.setBranch<float>("jetAK8puppi_sd_NewID", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3]));
    ana.tx.setBranch<float>("Jet_rawFactor_NewID", nt.Jet_rawFactor()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_NewID", nt.FatJet_deepTagMD_WvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_NewID", nt.FatJet_deepTag_WvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_NewID", nt.FatJet_deepTagMD_ZvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_NewID", nt.FatJet_deepTag_ZvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_NewID", nt.FatJet_deepTagMD_TvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_NewID", nt.FatJet_deepTag_TvsQCD()[usenumber3]);
    // ana.tx.pushbackToBranch<float>("Common_fatjet_deepMD_bb", nt.FatJet_deepTagMD_bbvsLight()[usenumber3]);
    }
    
    // second pt jet
    int usenumber2 = -1; 
    pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum] != usenumber3) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].pt(); 
            usenumber2 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]; 
            continue;
        }
    }
    if (usenumber2>-1) {

    
    if (not nt.isData()){
        TLorentzVector JERcorrectedP4,JERuncorrectedP4;
        JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[usenumber2].pt(),nt.FatJet_p4()[usenumber2].eta(),nt.FatJet_p4()[usenumber2].phi(),nt.FatJet_p4()[usenumber2].energy());
        JERcorrectedP4 = Process_1Lep2fatjet_JER(usenumber2) * JERuncorrectedP4;

        ana.tx.setBranch<float>("jetAK8puppi_ptJECJER_2_NewID",JERcorrectedP4.Pt());
    }
    
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_2_NewID", nt.FatJet_p4()[usenumber2].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_2_NewID", nt.FatJet_p4()[usenumber2].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_2_NewID", nt.FatJet_p4()[usenumber2].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_2_NewID", nt.FatJet_jetId()[usenumber2]&2);

    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_2_NewID", nt.FatJet_msoftdrop()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_sd_2_NewID", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber2],nt.FatJet_subJetIdx2()[usenumber2]));
    ana.tx.setBranch<float>("Jet_rawFactor_2_NewID", nt.Jet_rawFactor()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_2_NewID", nt.FatJet_deepTagMD_WvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_2_NewID", nt.FatJet_deepTag_WvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_2_NewID", nt.FatJet_deepTagMD_ZvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_2_NewID", nt.FatJet_deepTag_ZvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_2_NewID", nt.FatJet_deepTagMD_TvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_2_NewID", nt.FatJet_deepTag_TvsQCD()[usenumber2]);
    }

    // third pt jet
    int usenumber1 = -1; 
    pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum] != usenumber3 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum] != usenumber2) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]].pt(); 
            usenumber1 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new_NewID")[inum]; 
            continue;
        }
    }
    if (usenumber1>-1) {

    if (not nt.isData()){
        TLorentzVector JERcorrectedP4,JERuncorrectedP4;
        JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[usenumber1].pt(),nt.FatJet_p4()[usenumber1].eta(),nt.FatJet_p4()[usenumber1].phi(),nt.FatJet_p4()[usenumber1].energy());
        JERcorrectedP4 = Process_1Lep2fatjet_JER(usenumber1) * JERuncorrectedP4;

        ana.tx.setBranch<float>("jetAK8puppi_ptJECJER_3_NewID",JERcorrectedP4.Pt());
    }

    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_3_NewID", nt.FatJet_p4()[usenumber1].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_3_NewID", nt.FatJet_p4()[usenumber1].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_3_NewID", nt.FatJet_p4()[usenumber1].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_3_NewID", nt.FatJet_jetId()[usenumber1]&2);

    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_3_NewID", nt.FatJet_msoftdrop()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_sd_3_NewID", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber1],nt.FatJet_subJetIdx2()[usenumber1]));
    ana.tx.setBranch<float>("Jet_rawFactor_3_NewID", nt.Jet_rawFactor()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_3_NewID", nt.FatJet_deepTagMD_WvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_3_NewID", nt.FatJet_deepTag_WvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_3_NewID", nt.FatJet_deepTagMD_ZvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_3_NewID", nt.FatJet_deepTag_ZvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_3_NewID", nt.FatJet_deepTagMD_TvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_3_NewID", nt.FatJet_deepTag_TvsQCD()[usenumber1]);
    }

    int numofAK8 = 0;
    if(ana.tx.getBranchLazy<float>("jetAK8puppi_ptJECJER_NewID")>200){
        numofAK8++;
    }
    if(ana.tx.getBranchLazy<float>("jetAK8puppi_ptJECJER_2_NewID")>200){
        numofAK8++;
    }
    if(ana.tx.getBranchLazy<float>("jetAK8puppi_ptJECJER_3_NewID")>0){
        numofAK8++;
    }
    ana.tx.setBranch<int>("t_1Lep2fatjet_Nj8_NewID", numofAK8);

}

void Process_1Lep2fatjet_jet_NewID(){
for (unsigned int ijet = 0; ijet < nt.Jet_p4().size(); ++ijet)
    {
      
        //  do we use nt.Jet_jetId()?
        //   Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto

        if ( not (nt.Jet_p4()[ijet].pt() > 20.) )
            continue;
        int _Jet_jetId = nt.Jet_jetId()[ijet];
        // cout<<"_Jet_jetId"<<_Jet_jetId<<(_Jet_jetId &2)<<endl;
        if (_Jet_jetId &2 !=2)
            continue;
        if (not (abs(nt.Jet_p4()[ijet].eta()) < 5.0))
            continue;

        // Because every muon and electron shows up in PF FatJet collections
        // Need to check against leptons
        bool is_overlapping_with_a_lepton = false;

        // Overlap check against leptons (electrons)
        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs")[ilep];
            // If electron
            if (abs(ana.tx.getBranchLazy<vector<int>>("Common_lep_pdgid")[ilep]) == 11)
            {
                if (RooUtil::Calc::DeltaR(nt.Jet_p4()[ijet], nt.Electron_p4()[ilep_idx]) < 0.4)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    break;
                }
            }
            // else muon
            else
            {

                if (RooUtil::Calc::DeltaR(nt.Jet_p4()[ijet], nt.Muon_p4()[ilep_idx]) < 0.4)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    break;
                }
            }
        }

        if (is_overlapping_with_a_lepton)
            continue;
 
        // For now, accept anything that reaches this point

    ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_jet_idxs_NewID", ijet);
    ana.tx.pushbackToBranch<LorentzVector>("Common_jet_p4_new_NewID", nt.Jet_p4()[ijet]);
    if (not nt.isData()){
    ana.tx.pushbackToBranch<int>("ak4jet_hf_NewID", nt.Jet_hadronFlavour()[ijet]);
    ana.tx.pushbackToBranch<int>("ak4jet_pf_NewID", nt.Jet_partonFlavour()[ijet]);
    }
    ana.tx.pushbackToBranch<float>("ak4jet_pt_NewID", nt.Jet_p4()[ijet].pt());
    ana.tx.pushbackToBranch<float>("ak4jet_eta_NewID", nt.Jet_p4()[ijet].eta());
    ana.tx.pushbackToBranch<float>("ak4jet_phi_NewID", nt.Jet_p4()[ijet].phi());
    ana.tx.pushbackToBranch<float>("ak4jet_e_NewID", nt.Jet_p4()[ijet].energy());
    ana.tx.pushbackToBranch<float>("Jet_btagDeepB_NewID", nt.Jet_btagDeepB()[ijet]);
    ana.tx.pushbackToBranch<float>("Jet_btagDeepC_NewID", nt.Jet_btagDeepC()[ijet]);

    if (not nt.isData()){
    
    TLorentzVector JERcorrectedP4,JERuncorrectedP4;
    JERuncorrectedP4.SetPtEtaPhiE(nt.Jet_p4()[ijet].pt(),nt.Jet_p4()[ijet].eta(),nt.Jet_p4()[ijet].phi(),nt.Jet_p4()[ijet].energy());
    JERcorrectedP4 = Process_1Lep2fatjet_JER_AK4(ijet) * JERuncorrectedP4;
    ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_JER_NewID", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()));
    ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_unJER_NewID", LorentzVector(JERuncorrectedP4.Pt(),JERuncorrectedP4.Eta(),JERuncorrectedP4.Phi(),JERuncorrectedP4.M()));
    ana.tx.pushbackToBranch<float>("ak4jet_pt_JER_NewID",JERcorrectedP4.Pt() );
    ana.tx.pushbackToBranch<float>("ak4jet_e_JER_NewID",JERcorrectedP4.Energy() );
    }

    }
}

void Process_1Lep2fatjet_MET_Lepton_leptonicW_NewID_bk(){



    // lepton

    int nloose=0;
    int nlooseEle=0;
    int nlooseMu=0;
    int nreallooseEle=0;
    int nreallooseMu=0;
    //loop over all electrons
    for(unsigned int iel=0;iel<nt.Electron_p4().size();iel++){
	if(!(nt.Electron_mvaFall17V2Iso_WP90().at(iel))) continue;
	if(!(nt.Electron_p4().at(iel).pt()>10)) continue;
	if(!(fabs(nt.Electron_p4().at(iel).eta())<2.5)) continue;
	if(fabs(nt.Electron_p4().at(iel).eta())<1.556&&fabs(nt.Electron_p4().at(iel).eta())>1.442) continue;
	if(!(nt.Electron_ip3d().at(iel)<0.01)) continue;
	if(!(fabs(nt.Electron_dxy().at(iel))<0.05)) continue;
	if(!(fabs(nt.Electron_dz().at(iel))<0.1)) continue;
    nreallooseEle++;
	if(!(nt.Electron_p4().at(iel).pt()>25)) continue;
	nloose++;
	nlooseEle++;
	ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_idx", iel);
    ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_pdgid", nt.Electron_pdgId().at(iel));
    }
    //loop over all muons
    for(unsigned int imu=0;imu<nt.Muon_p4().size();imu++){
	if(!(nt.Muon_mediumId().at(imu))) continue;
	if(!(nt.Muon_p4().at(imu).pt()>10)) continue;
	if(!(fabs(nt.Muon_p4().at(imu).eta())<2.4)) continue;
	if(!(nt.Muon_ip3d().at(imu)<0.015)) continue;
	if(!(nt.Muon_sip3d().at(imu)<4)) continue;
	if(!(fabs(nt.Muon_dxy().at(imu))<0.05)) continue;
	if(!(fabs(nt.Muon_dz().at(imu))<0.1)) continue;
	if(not (nt.Muon_pfRelIso03_all().at(imu)<0.4)) continue;
    nreallooseMu++;
	if(!(nt.Muon_p4().at(imu).pt()>25)) continue;
    nlooseMu++;
    nloose++;
	ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_idx", imu);
    ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_pdgid", nt.Muon_pdgId().at(imu));
   }
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nloose",nloose);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nlooseEle",nlooseEle);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nlooseMu",nlooseMu);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nreallooseMu",nreallooseMu);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nreallooseEle",nreallooseEle);
    
    float ptlep1 =-99., etalep1, philep1, energylep1;
    int lep;
    float RelIso03=-99;
    float Electron_mvaFall17V2Iso_WP80=-99;
    
    if ( ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx").size() == 1 ){
        int ilep_idx = ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0];
            // If electron
            if (abs(ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_pdgid")[0]) == 11)
            {
                ptlep1 = nt.Electron_p4()[ ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].pt();
                etalep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].eta();
                philep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].phi();
                energylep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].energy();
                lep = nt.Electron_pdgId()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];
                Electron_mvaFall17V2Iso_WP80 = nt.Electron_mvaFall17V2Iso_WP80()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];
            }
            if (abs(ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_pdgid")[0]) == 13)
            {
                ptlep1 = nt.Muon_p4()[ ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].pt();
                etalep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].eta();
                philep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].phi();
                energylep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].energy();
                lep = nt.Muon_pdgId()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];
                RelIso03 = nt.Muon_pfRelIso03_all()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];

            }
            
    }
    ana.tx.setBranch<float>("ptlep1_NewID", ptlep1 );
    ana.tx.setBranch<float>("etalep1_NewID", etalep1 );
    ana.tx.setBranch<float>("philep1_NewID", philep1 );
    ana.tx.setBranch<float>("energylep1_NewID", energylep1 );
    ana.tx.setBranch<int>("lep_NewID", lep );
    ana.tx.setBranch<float>("RelIso03_NewID", RelIso03 );
    ana.tx.setBranch<float>("mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80 );


    


 
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
    ana.tx.setBranch<float>("ptVlepJEC_NewID", ptVlepJEC );
    ana.tx.setBranch<float>("yVlepJEC_NewID", yVlepJEC );
    ana.tx.setBranch<float>("phiVlepJEC_NewID", phiVlepJEC );
    ana.tx.setBranch<float>("massVlepJEC_NewID", massVlepJEC );
    ana.tx.setBranch<float>("mtVlepJEC_NewID", mtVlepJEC );

    // neutrino
    float ptlep2, etalep2, philep2, energylep2;
    ptlep2 = neutrinoP4.Pt();
    etalep2 = neutrinoP4.Eta();
    philep2 = neutrinoP4.Phi();
    energylep2 = neutrinoP4.Energy();
    ana.tx.setBranch<float>("ptlep2_NewID", ptlep2 );
    ana.tx.setBranch<float>("etalep2_NewID", etalep2 );
    ana.tx.setBranch<float>("philep2_NewID", philep2 );
    ana.tx.setBranch<float>("energylep2_NewID", energylep2 );
}

void Process_1Lep2fatjet_MET_Lepton_leptonicW_NewID(){



    // lepton

    int nloose=0;
    int nlooseEle=0;
    int nlooseMu=0;
    int nreallooseEle=0;
    int nreallooseMu=0;
    //loop over all electrons
    for(unsigned int iel=0;iel<nt.Electron_p4().size();iel++){
	if(!(nt.Electron_mvaFall17V2Iso_WP90().at(iel))) continue;
	if(!(nt.Electron_p4().at(iel).pt()>10)) continue;
	if(!(fabs(nt.Electron_p4().at(iel).eta())<2.5)) continue;
	if(fabs(nt.Electron_p4().at(iel).eta())<1.556&&fabs(nt.Electron_p4().at(iel).eta())>1.442) continue;
	if(!(nt.Electron_ip3d().at(iel)<0.01)) continue;
	if(!(fabs(nt.Electron_dxy().at(iel))<0.05)) continue;
	if(!(fabs(nt.Electron_dz().at(iel))<0.1)) continue;
    nreallooseEle++;
	if(!(nt.Electron_p4().at(iel).pt()>25)) continue;
	nloose++;
	nlooseEle++;
	ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_idx", iel);
    ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_pdgid", nt.Electron_pdgId().at(iel));
    }
    //loop over all muons
    for(unsigned int imu=0;imu<nt.Muon_p4().size();imu++){
	if(!(nt.Muon_mediumId().at(imu))) continue;
	if(!(nt.Muon_p4().at(imu).pt()>10)) continue;
	if(!(fabs(nt.Muon_p4().at(imu).eta())<2.4)) continue;
	if(!(nt.Muon_ip3d().at(imu)<0.015)) continue;
	if(!(nt.Muon_sip3d().at(imu)<4)) continue;
	if(!(fabs(nt.Muon_dxy().at(imu))<0.05)) continue;
	if(!(fabs(nt.Muon_dz().at(imu))<0.1)) continue;
	if(not (nt.Muon_pfRelIso03_all().at(imu)<0.4)) continue;
    nreallooseMu++;
	if(!(nt.Muon_p4().at(imu).pt()>25)) continue;
    nlooseMu++;
    nloose++;
	ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_idx", imu);
    ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_lep_pdgid", nt.Muon_pdgId().at(imu));
   }
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nloose",nloose);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nlooseEle",nlooseEle);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nlooseMu",nlooseMu);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nreallooseMu",nreallooseMu);
   ana.tx.setBranch<int>	("Process_1Lep2fatjet_nreallooseEle",nreallooseEle);
    
    float ptlep1 =-99., etalep1, philep1, energylep1;
    int lep;
    float RelIso03=-99;
    float Electron_mvaFall17V2Iso_WP80=-99;
    
    if ( ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx").size() == 1 ){
        int ilep_idx = ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0];
            // If electron
            if (abs(ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_pdgid")[0]) == 11)
            {
                ptlep1 = nt.Electron_p4()[ ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].pt();
                etalep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].eta();
                philep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].phi();
                energylep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].energy();
                lep = nt.Electron_pdgId()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];
                Electron_mvaFall17V2Iso_WP80 = nt.Electron_mvaFall17V2Iso_WP80()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];
            }
            if (abs(ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_pdgid")[0]) == 13)
            {
                ptlep1 = nt.Muon_p4()[ ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].pt();
                etalep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].eta();
                philep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].phi();
                energylep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]].energy();
                lep = nt.Muon_pdgId()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];
                RelIso03 = nt.Muon_pfRelIso03_all()[ana.tx.getBranchLazy<vector<int>>("Process_1Lep2fatjet_lep_idx")[0]];

            }
            
    }
    ana.tx.setBranch<float>("ptlep1_NewID", ptlep1 );
    ana.tx.setBranch<float>("etalep1_NewID", etalep1 );
    ana.tx.setBranch<float>("philep1_NewID", philep1 );
    ana.tx.setBranch<float>("energylep1_NewID", energylep1 );
    ana.tx.setBranch<int>("lep_NewID", lep );
    ana.tx.setBranch<float>("RelIso03_NewID", RelIso03 );
    ana.tx.setBranch<float>("mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80 );


    


 
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
    ana.tx.setBranch<float>("ptVlepJEC_NewID", ptVlepJEC );
    ana.tx.setBranch<float>("yVlepJEC_NewID", yVlepJEC );
    ana.tx.setBranch<float>("phiVlepJEC_NewID", phiVlepJEC );
    ana.tx.setBranch<float>("massVlepJEC_NewID", massVlepJEC );
    ana.tx.setBranch<float>("mtVlepJEC_NewID", mtVlepJEC );

    // neutrino
    float ptlep2, etalep2, philep2, energylep2;
    ptlep2 = neutrinoP4.Pt();
    etalep2 = neutrinoP4.Eta();
    philep2 = neutrinoP4.Phi();
    energylep2 = neutrinoP4.Energy();
    ana.tx.setBranch<float>("ptlep2_NewID", ptlep2 );
    ana.tx.setBranch<float>("etalep2_NewID", etalep2 );
    ana.tx.setBranch<float>("philep2_NewID", philep2 );
    ana.tx.setBranch<float>("energylep2_NewID", energylep2 );
}


#endif