#ifndef Process_1Lep2fatjet_UpdateJER_h
#define Process_1Lep2fatjet_UpdateJER_h

void Process_1Lep2fatjet_UpdateJER_AK4(){
    if (not nt.isData()){
        for (unsigned int ijet = 0; ijet < nt.Jet_p4().size(); ++ijet){
            
            int ijetid = ijet;
            ana.res.loadVariable("JetEta", nt.Jet_p4()[ijetid].eta());
            ana.res.loadVariable("Rho", nt.fixedGridRhoFastjetAll());
            ana.res.loadVariable("JetPt", nt.Jet_p4()[ijetid].pt());
            vector<double> GenJetPt;
            for (unsigned int genJetIdx = 0; genJetIdx < nt.GenJet_p4().size(); ++genJetIdx){
                GenJetPt.push_back(nt.GenJet_p4()[genJetIdx].pt());
            }
            auto smearing = ana.res.smear(nt.Jet_p4()[ijetid], nt.GenJet_p4(), GenJetPt, 0, 0.4); // need to clean the leptonic jets?
            auto matching = ana.res.match();

            ana.tx.pushbackToBranch<float>("Common_Jet_pt_JER", smearing[0] );
            ana.tx.pushbackToBranch<float>("Common_Jet_pt_unJER", nt.Jet_p4()[ijetid].pt() );
            ana.tx.pushbackToBranch<float>("Common_JER_matching", matching[0] );

            TLorentzVector JERcorrectedP4,JERuncorrectedP4;
            JERuncorrectedP4.SetPtEtaPhiE(nt.Jet_p4()[ijetid].pt(),nt.Jet_p4()[ijetid].eta(),nt.Jet_p4()[ijetid].phi(),nt.Jet_p4()[ijetid].energy());
            JERcorrectedP4 = (smearing[0]/nt.Jet_p4()[ijetid].pt()) * JERuncorrectedP4;
            ana.tx.pushbackToBranch<LorentzVector>("Common_Jet_JER_p4", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()) );



        }

        for (unsigned int ijet = 0; ijet < nt.FatJet_p4().size(); ++ijet){
            
            int ijetid = ijet;
            ana.resfatJet.loadVariable("JetEta", nt.FatJet_p4()[ijetid].eta());
            ana.resfatJet.loadVariable("Rho", nt.fixedGridRhoFastjetAll());
            ana.resfatJet.loadVariable("JetPt", nt.FatJet_p4()[ijetid].pt());
            vector<double> GenJetPt;
            for (unsigned int genJetIdx = 0; genJetIdx < nt.GenJetAK8_p4().size(); ++genJetIdx){
                GenJetPt.push_back(nt.GenJetAK8_p4()[genJetIdx].pt());
            }
            auto smearing = ana.resfatJet.smear(nt.FatJet_p4()[ijetid], nt.GenJetAK8_p4(), GenJetPt, 0, 0.8); // need to clean the leptonic jets?
            auto matching = ana.resfatJet.match();
            auto matchingpt = ana.resfatJet.matchpt();


            ana.tx.pushbackToBranch<float>("Common_fatJet_pt_JER", smearing[0] );
            ana.tx.pushbackToBranch<float>("Common_fatJet_pt_unJER", nt.FatJet_p4()[ijetid].pt() );
            ana.tx.pushbackToBranch<float>("Common_fatJER_matching", matching[0] );
            ana.tx.pushbackToBranch<int>("Lep1fatJet2_MatchGenJet_num", matching.size() );
            ana.tx.pushbackToBranch<float>("Lep1fatJet2_MatchGenJet_pt", matchingpt[0] );
            ana.tx.pushbackToBranch<float>("Lep1fatJet2_JERSF", ana.resfatJet.getScaleFactor(0) );

            TLorentzVector JERcorrectedP4,JERuncorrectedP4;
            JERuncorrectedP4.SetPtEtaPhiE(nt.FatJet_p4()[ijetid].pt(),nt.FatJet_p4()[ijetid].eta(),nt.FatJet_p4()[ijetid].phi(),nt.FatJet_p4()[ijetid].energy());
            JERcorrectedP4 = (smearing[0]/nt.FatJet_p4()[ijetid].pt()) * JERuncorrectedP4;
            ana.tx.pushbackToBranch<LorentzVector>("Common_fatJet_JER_p4", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()) );
        }

        float corrEx_MET_JER = 0;
        float corrEy_MET_JER = 0;
        for (unsigned int ijet = 0; ijet < ana.tx.getBranchLazy<vector<LorentzVector>>("Common_Jet_JER_p4").size(); ++ijet){
            corrEx_MET_JER += ana.tx.getBranchLazy<vector<LorentzVector>>("Common_Jet_JER_p4")[ijet].px() - ana.tx.getBranchLazy<vector<LorentzVector>>("Common_Jet_JER_p4")[ijet].px();
            corrEy_MET_JER += ana.tx.getBranchLazy<vector<LorentzVector>>("Common_Jet_JER_p4")[ijet].py() - ana.tx.getBranchLazy<vector<LorentzVector>>("Common_Jet_JER_p4")[ijet].py();
        }
        TVector2 rawMET_;
        rawMET_.SetMagPhi (nt.MET_pt(), nt.MET_phi() );
        float rawPx = rawMET_.Px();
        float rawPy = rawMET_.Py();
        float pxcorr = rawPx+corrEx_MET_JER;
        float pycorr = rawPy+corrEy_MET_JER;
        float et     = std::hypot(pxcorr,pycorr);
        TLorentzVector corrmet;
        corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
        ana.tx.setBranch<float>("Common_MET_pt_JER",et );
        ana.tx.setBranch<float>("Common_MET_phi_JER",corrmet.Phi() );
    }
}

TLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){
        float MW_=80.385;

        double leppt = lep.Pt();
        double lepphi = lep.Phi();
        double lepeta = lep.Eta();
        double lepenergy = lep.Energy();

        double metpt = MetPt;
        double metphi = MetPhi;

        double  px = metpt*cos(metphi);
        double  py = metpt*sin(metphi);
        double  pz = 0;
        double  pxl= leppt*cos(lepphi);
        double  pyl= leppt*sin(lepphi);
        double  pzl= leppt*sinh(lepeta);
        double  El = lepenergy;
        double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
        double  b = 2.*pzl;
        double  A = b*b -4.*El*El;
        double  B = 2.*a*b;
        double  C = a*a-4.*(px*px+py*py)*El*El;

        ///////////////////////////pz for fnal
        double M_mu =  0;

        //if(lepType==1)M_mu=0.105658367;//mu
        //if(lepType==0)M_mu=0.00051099891;//electron

        int type=2; // use the small abs real root

        a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
        A = 4.0*(El*El - pzl*pzl);
        B = -4.0*a*pzl;
        C = 4.0*El*El*(px*px + py*py) - a*a;

        double tmproot = B*B - 4.0*A*C;

        if (tmproot<0) {
            //std::cout << "Complex root detected, taking real part..." << std::endl;
            pz = - B/(2*A); // take real part of complex roots
        }
        else {
            double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
            double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
            //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;

            if (type == 0 ) {
                // two real roots, pick the one closest to pz of muon
                if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
                else { pz = tmpsol1; }
                // if pz is > 300 pick the most central root
                if ( abs(pz) > 300. ) {
                    if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                    else { pz = tmpsol2; }
                }
            }
            if (type == 1 ) {
                // two real roots, pick the one closest to pz of muon
                if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
                else {pz = tmpsol1; }
            }
            if (type == 2 ) {
                // pick the most central root.
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
            /*if (type == 3 ) {
             // pick the largest value of the cosine
             TVector3 p3w, p3mu;
             p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
             p3mu.SetXYZ(pxl, pyl, pzl );

             double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
             p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
             double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;

             double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
             double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);

             if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
             else { pz = tmpsol2;otherSol_ = tmpsol1; }

             }*///end of type3

        }//endl of if real root

        //dont correct pt neutrino
        TLorentzVector outP4;
        outP4.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz));
        return outP4;

    }

void Process_1Lep2fatjet_Geninfo(){
    // Gen level info
    
    if (not nt.isData())
    {
        // Z boson info, last copy
        for(size_t ik=0; ik<nt.nGenPart();ik++)
        {
            if (abs(nt.GenPart_pdgId()[ik]) == 23 )
            {
                if (not (nt.GenPart_statusFlags()[ik]&(1<<13))) continue; // isLastCopy
                int indexgenzl; indexgenzl = ik;
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_ptgenzl",nt.GenPart_pt()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_etagenzl",nt.GenPart_eta()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_phigenzl",nt.GenPart_phi()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_energygenzl",nt.GenPart_p4()[ik].energy());
                vector<float> daughter_index;
                for (size_t id=0; id<nt.nGenPart();id++){
                    if (nt.GenPart_genPartIdxMother()[id] == indexgenzl){
                        daughter_index.push_back(id);
                    }
                }
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q1_pt",nt.GenPart_pt()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q1_eta",nt.GenPart_eta()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q1_phi",nt.GenPart_phi()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q1_e",nt.GenPart_p4()[daughter_index[0]].energy());
                ana.tx.pushbackToBranch<int>("Lep1fatJet2_genz_q1_pdg",nt.GenPart_pdgId()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q2_pt",nt.GenPart_pt()[daughter_index[1]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q2_eta",nt.GenPart_eta()[daughter_index[1]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q2_phi",nt.GenPart_phi()[daughter_index[1]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genz_q2_e",nt.GenPart_p4()[daughter_index[1]].energy());
                ana.tx.pushbackToBranch<int>("Lep1fatJet2_genz_q2_pdg",nt.GenPart_pdgId()[daughter_index[1]]);
            }
        }

        // W boson info, first copy
        for(size_t ik=0; ik<nt.nGenPart();ik++)
        {
            if (abs(nt.GenPart_pdgId()[ik]) == 24 )
            {
                if (not (nt.GenPart_statusFlags()[ik]&(1<<13))) continue; // isLastCopy
                int indexgenwl;  indexgenwl = ik;
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_ptgenwl",nt.GenPart_pt()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_etagenwl",nt.GenPart_eta()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_phigenwl",nt.GenPart_phi()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_energygenwl",nt.GenPart_p4()[ik].energy());
                vector<float> daughter_index;
                for (size_t id=0; id<nt.nGenPart();id++){
                
                    if (nt.GenPart_genPartIdxMother()[id] == indexgenwl){
                        daughter_index.push_back(id);
                    }
                }
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q1_pt",nt.GenPart_p4()[daughter_index[0]].pt());
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q1_eta",nt.GenPart_eta()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q1_phi",nt.GenPart_phi()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q1_e",nt.GenPart_p4()[daughter_index[0]].energy());
                ana.tx.pushbackToBranch<int>("Lep1fatJet2_genw_q1_pdg",nt.GenPart_pdgId()[daughter_index[0]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q2_pt",nt.GenPart_pt()[daughter_index[1]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q2_eta",nt.GenPart_eta()[daughter_index[1]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q2_phi",nt.GenPart_phi()[daughter_index[1]]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_genw_q2_e",nt.GenPart_p4()[daughter_index[1]].energy());
                ana.tx.pushbackToBranch<int>("Lep1fatJet2_genw_q2_pdg",nt.GenPart_pdgId()[daughter_index[1]]);
            }
        }
        // top info
        for(size_t ik=0; ik<nt.nGenPart();ik++)
        {
            if (abs(nt.GenPart_pdgId()[ik]) == 6 )
            {
                if (not (nt.GenPart_statusFlags()[ik]&(1<<13))) continue; // isLastCopy
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_ptgentopl",nt.GenPart_pt()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_etagentopl",nt.GenPart_eta()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_phigentopl",nt.GenPart_phi()[ik]);
                ana.tx.pushbackToBranch<float>("Lep1fatJet2_energygentopl",nt.GenPart_p4()[ik].energy());
            }
        }



    }
}

void Process_1Lep2fatjet_Weight_HLT(){
    // for L1 prefiring 
    if ( (not ana.is_EFT_sample) && (nt.year() == 2016)){
    ana.tx.setBranch<float>("Lep1fatJet2_L1PreFiringWeight_Nom", nt.L1PreFiringWeight_Nom());
    }

    // pile-up reweight
    if ( not nt.isData()){
    ana.tx.setBranch<float>("Lep1fatJet2_Pileup_nTrueInt", nt.Pileup_nTrueInt());
    ana.tx.setBranch<float>("Lep1fatJet2_Pileup_nPU", nt.Pileup_nPU());
    }


    // HLT
    ana.tx.setBranch<int>("Lep1fatJet2_HLT_Mu50", nt.HLT_Mu50() );
    ana.tx.setBranch<int>("Lep1fatJet2_HLT_IsoMu24", nt.HLT_IsoMu24() );
    try { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_OldMu100", nt.HLT_OldMu100()); }                catch (std::runtime_error) { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_OldMu100", 0); }
    try { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_TkMu100", nt.HLT_TkMu100()); }                catch (std::runtime_error) { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_TkMu100", 0); }
    ana.tx.setBranch<int>("Lep1fatJet2_HLT_IsoMu27", nt.HLT_IsoMu27() );
    
    ana.tx.setBranch<int>("Lep1fatJet2_HLT_Ele27", nt.HLT_Ele27_WPTight_Gsf() );
    ana.tx.setBranch<int>("Lep1fatJet2_HLT_Ele35_WPTight_Gsf", nt.HLT_Ele35_WPTight_Gsf() );
    ana.tx.setBranch<int>("Lep1fatJet2_HLT_Photon200", nt.HLT_Photon200() );
    try { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_Ele32_WPTight_Gsf_L1DoubleEG", nt.HLT_Ele32_WPTight_Gsf_L1DoubleEG()); }                catch (std::runtime_error) { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_Ele32_WPTight_Gsf_L1DoubleEG", 0); }
    try { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_Ele32_WPTight_Gsf", nt.HLT_Ele32_WPTight_Gsf()); }                catch (std::runtime_error) { ana.tx.setBranch<int>        ("Lep1fatJet2_HLT_Ele32_WPTight_Gsf", 0); }


    // vertex
    // ana.tx.setBranch<int>("Process_1Lep2fatjet_PV_npvsGood", nt.PV_npvsGood() );
}

void Process_1Lep2fatjet_Jet(){
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs").size(); ++inum){
        ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_btagDeepB", nt.Jet_btagDeepB()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]]);
        ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_btagDeepC", nt.Jet_btagDeepC()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]]);
        int _Jet_jetId = nt.Jet_jetId()[inum]; _Jet_jetId = ((_Jet_jetId &2) == 2);
        ana.tx.pushbackToBranch<int>("Lep1fatJet2_Jet_jetId", _Jet_jetId);
        if (not nt.isData()){
            ana.tx.pushbackToBranch<int>("Lep1fatJet2_Jet_hadronFlavour", nt.Jet_hadronFlavour()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]]);
            ana.tx.pushbackToBranch<int>("Lep1fatJet2_Jet_partonFlavour", nt.Jet_partonFlavour()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]]);
            ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_pt_JER", ana.tx.getBranchLazy<vector<float>>("Common_Jet_pt_JER")[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]]);
        }
        else{
            ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_pt_JER", nt.Jet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]].pt());
        }
        ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_eta", nt.Jet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]].eta());
        ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_phi", nt.Jet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]].phi());
        ana.tx.pushbackToBranch<float>("Lep1fatJet2_Jet_e", nt.Jet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_jet_idxs")[inum]].energy());
    }
}

void Process_1Lep2fatjet_fatJet(){
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size(); ++inum){
        int _Jet_jetId = nt.FatJet_jetId()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]];         _Jet_jetId = ((_Jet_jetId &2) == 2);
        ana.tx.pushbackToBranch<int>("Lep1fatJet2_FatJet_jetId", _Jet_jetId);
        if (not nt.isData()){
            ana.tx.pushbackToBranch<float>("Lep1fatJet2_FatJet_pt_JER", ana.tx.getBranchLazy<vector<float>>("Common_fatJet_pt_JER")[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]]);
        }
    }
}

void Process_1Lep2fatjet_MET_Lepton_leptonicW(){

    float ptlep1 =-99., etalep1 =-99., philep1 =-99., energylep1 =-99.;
    
    int lep;
    if ( ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs").size() == 1 ){
        ptlep1 = ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[0].pt();
        etalep1 = ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[0].eta();
        philep1 = ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[0].phi();
        energylep1 = ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[0].energy();
        lep = ana.tx.getBranchLazy<vector<int>>("Common_lep_pdgid")[0];
        ana.tx.setBranch<float>("Lep1fatJet2_LeptonPt", ptlep1 );
        ana.tx.setBranch<float>("Lep1fatJet2_LeptonEta", etalep1 );
        ana.tx.setBranch<float>("Lep1fatJet2_LeptonPhi", philep1 );
        ana.tx.setBranch<float>("Lep1fatJet2_LeptonE", energylep1 );
        ana.tx.setBranch<int>("Lep1fatJet2_LeptonPDGID", lep );
        if (abs(lep) ==13){
        ana.tx.setBranch<float>("Lep1fatJet2_Muon_pfRelIso04_all", nt.Muon_pfRelIso04_all()[ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs")[0]] );
        }
    }
    
    
    TLorentzVector  glepton,neutrino,neutrinoP4,WLeptonic;
    glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
    int leptontype = 1; double MET_et_JER = 0, MET_phi = 0;

    if (not nt.isData()){
        MET_et_JER = ana.tx.getBranchLazy<float>("Common_MET_pt_JER");
        MET_phi = ana.tx.getBranchLazy<float>("Common_MET_phi_JER");
    }
    else{
        MET_et_JER = nt.MET_pt();
        MET_phi = nt.MET_phi();
    }
    
    neutrino = getNeutrinoP4(MET_et_JER , MET_phi, glepton, leptontype);
    neutrinoP4.SetPtEtaPhiE(neutrino.Pt(),neutrino.Eta(),neutrino.Phi(),neutrino.Energy());
    
    ana.tx.setBranch<float>("Lep1fatJet2_NeutrinoPt", neutrinoP4.Pt() );
    ana.tx.setBranch<float>("Lep1fatJet2_NeutrinoEta", neutrinoP4.Eta() );
    ana.tx.setBranch<float>("Lep1fatJet2_Neutrinophi", neutrinoP4.Phi() );
    ana.tx.setBranch<float>("Lep1fatJet2_NeutrinoE", neutrinoP4.Energy() );

    WLeptonic = glepton+neutrinoP4;
    ana.tx.setBranch<float>("Lep1fatJet2_LeptonicWPt", WLeptonic.Pt() );
    ana.tx.setBranch<float>("Lep1fatJet2_LeptonicWEta", WLeptonic.Eta() );
    ana.tx.setBranch<float>("Lep1fatJet2_LeptonicWPhi", WLeptonic.Phi() );
    ana.tx.setBranch<float>("Lep1fatJet2_LeptonicWM", WLeptonic.M() );
    ana.tx.setBranch<float>("Lep1fatJet2_LeptonicWMt", WLeptonic.Mt() );


}

void Process_1Lep2fatjet_StorefatJet(){
    if (nt.isData()){
        int usenumber3 = -1; double pt_larger=0;
        for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size(); ++inum){
            if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].eta())<2.4 && inum<4) {
                pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].pt(); 
                usenumber3 = ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]; 
                continue;
            }
        }
        if (usenumber3>-1) {
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJECJER", nt.FatJet_p4()[usenumber3].pt());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJEC", nt.FatJet_p4()[usenumber3].pt());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_eta", nt.FatJet_p4()[usenumber3].eta());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_phi", nt.FatJet_p4()[usenumber3].phi());
            ana.tx.setBranch<int>("Lep1fatJet2_jetAK8puppi_tightid", nt.FatJet_jetId()[usenumber3]&2);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_sd", nt.FatJet_msoftdrop()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jet_rawFactor", nt.Jet_rawFactor()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrW", nt.FatJet_deepTagMD_WvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnW", nt.FatJet_deepTag_WvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrZ", nt.FatJet_deepTagMD_ZvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnZ", nt.FatJet_deepTag_ZvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrTop", nt.FatJet_deepTagMD_TvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnTop", nt.FatJet_deepTag_TvsQCD()[usenumber3]);
        }
        int usenumber2 = -1; pt_larger=0;
        for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size(); ++inum){
            if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum] != usenumber3) {
                pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].pt(); 
                usenumber2 = ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]; 
                continue;
            }
        }
        if (usenumber2>-1) {
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJECJER_2",nt.FatJet_p4()[usenumber2].pt());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJEC_2", nt.FatJet_p4()[usenumber2].pt());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_eta_2", nt.FatJet_p4()[usenumber2].eta());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_phi_2", nt.FatJet_p4()[usenumber2].phi());
            ana.tx.setBranch<int>("Lep1fatJet2_jetAK8puppi_tightid_2", nt.FatJet_jetId()[usenumber2]&2);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_sd_2", nt.FatJet_msoftdrop()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jet_rawFactor_2", nt.Jet_rawFactor()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrW_2", nt.FatJet_deepTagMD_WvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnW_2", nt.FatJet_deepTag_WvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrZ_2", nt.FatJet_deepTagMD_ZvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnZ_2", nt.FatJet_deepTag_ZvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrTop_2", nt.FatJet_deepTagMD_TvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnTop_2", nt.FatJet_deepTag_TvsQCD()[usenumber2]);
        }
    }
    else {
        int usenumber3 = -1, usenumber3Jer = -1; double pt_larger=0;
        for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size(); ++inum){
            if(ana.tx.getBranchLazy<vector<float>>("Lep1fatJet2_FatJet_pt_JER")[inum] > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].eta())<2.4 && inum<4) {
                pt_larger = ana.tx.getBranchLazy<vector<float>>("Lep1fatJet2_FatJet_pt_JER")[inum]; 
                usenumber3 = ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]; 
                usenumber3Jer = inum; 
                continue;
            }
        }
        if (usenumber3>-1) {
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJECJER", ana.tx.getBranchLazy<vector<float>>("Lep1fatJet2_FatJet_pt_JER")[usenumber3Jer]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJEC", nt.FatJet_p4()[usenumber3].pt());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_eta", nt.FatJet_p4()[usenumber3].eta());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_phi", nt.FatJet_p4()[usenumber3].phi());
            ana.tx.setBranch<int>("Lep1fatJet2_jetAK8puppi_tightid", nt.FatJet_jetId()[usenumber3]&2);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_sd", nt.FatJet_msoftdrop()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jet_rawFactor", nt.Jet_rawFactor()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrW", nt.FatJet_deepTagMD_WvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnW", nt.FatJet_deepTag_WvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrZ", nt.FatJet_deepTagMD_ZvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnZ", nt.FatJet_deepTag_ZvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrTop", nt.FatJet_deepTagMD_TvsQCD()[usenumber3]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnTop", nt.FatJet_deepTag_TvsQCD()[usenumber3]);
        }
        int usenumber2 = -1, usenumber2Jer = -1; pt_larger=0;
        for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size(); ++inum){
            if(ana.tx.getBranchLazy<vector<float>>("Lep1fatJet2_FatJet_pt_JER")[inum] > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum] != usenumber3) {
                pt_larger = ana.tx.getBranchLazy<vector<float>>("Lep1fatJet2_FatJet_pt_JER")[inum]; 
                usenumber2 = ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs")[inum]; 
                usenumber2Jer = inum; 
                continue;
            }
        }
        if (usenumber2>-1) {
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJECJER_2",ana.tx.getBranchLazy<vector<float>>("Lep1fatJet2_FatJet_pt_JER")[usenumber2Jer]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_ptJEC_2", nt.FatJet_p4()[usenumber2].pt());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_eta_2", nt.FatJet_p4()[usenumber2].eta());
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_phi_2", nt.FatJet_p4()[usenumber2].phi());
            ana.tx.setBranch<int>("Lep1fatJet2_jetAK8puppi_tightid_2", nt.FatJet_jetId()[usenumber2]&2);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_sd_2", nt.FatJet_msoftdrop()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jet_rawFactor_2", nt.Jet_rawFactor()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrW_2", nt.FatJet_deepTagMD_WvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnW_2", nt.FatJet_deepTag_WvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrZ_2", nt.FatJet_deepTagMD_ZvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnZ_2", nt.FatJet_deepTag_ZvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrTop_2", nt.FatJet_deepTagMD_TvsQCD()[usenumber2]);
            ana.tx.setBranch<float>("Lep1fatJet2_jetAK8puppi_dnnTop_2", nt.FatJet_deepTag_TvsQCD()[usenumber2]);
        }
    }
}

void Process_1Lep2fatjet_StoreMET(){
    ana.tx.setBranch<float>("Lep1fatJet2_MET_pt" , nt.MET_pt()  );
    ana.tx.setBranch<float>("Lep1fatJet2_MET_phi", nt.MET_phi() );
}

void Process_1Lep2fatjet_Selection(){
    if (ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs").size() == 1 && ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size() == 2){
        ana.tx.setBranch<int>("Lep1fatJet2_Preselection", 1);
    }
    else{
        ana.tx.setBranch<int>("Lep1fatJet2_Preselection", 0);
    }
    if (ana.tx.getBranchLazy<int>("Lep1fatJet2_Preselection") == 1){
    // if (true) {
        Process_1Lep2fatjet_MET_Lepton_leptonicW();
        Process_1Lep2fatjet_Geninfo();
        if(true){
            Process_1Lep2fatjet_StorefatJet();
            Process_1Lep2fatjet_StoreMET();
        }
        if (true){
            ana.tx.fill();
        }
    }
}


#endif