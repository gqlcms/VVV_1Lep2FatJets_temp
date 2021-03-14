#ifndef Process_1Lep2fatjet_common_h
#define Process_1Lep2fatjet_common_h



// JER; FatJet
double Process_1Lep2fatjet_get_JER_corr(float JERSF, bool isMC, Int_t FatJet_index, double conSize, float PtResolution){
        double JER_corrFactor = 1.;
        if(isMC) {
                bool isGenMatched = 0;
                int genJetAK8Idx = nt.FatJet_genJetAK8Idx()[FatJet_index];
                int un_int_zero(0);
                int numofgenjet = nt.FatJet_genJetAK8Idx().size();
                float zero(0.);

                if (genJetAK8Idx>=un_int_zero && genJetAK8Idx < numofgenjet) {
                        TLorentzVector jetp4, genjetp4;
                        jetp4.SetPtEtaPhiE(nt.FatJet_p4()[FatJet_index].pt(), nt.FatJet_p4()[FatJet_index].eta(), nt.FatJet_p4()[FatJet_index].phi(), nt.FatJet_p4()[FatJet_index].energy());
                        genjetp4.SetPtEtaPhiE(nt.GenJetAK8_p4()[genJetAK8Idx].pt(), nt.GenJetAK8_p4()[genJetAK8Idx].eta(), nt.GenJetAK8_p4()[genJetAK8Idx].phi(), nt.GenJetAK8_p4()[genJetAK8Idx].energy());
                        float dR = jetp4.DeltaR(genjetp4);
                        float dPt = nt.FatJet_p4()[FatJet_index].pt()-nt.GenJetAK8_p4()[genJetAK8Idx].pt();
                        if ((dR<conSize/2.0)&&(std::abs(dPt)<(3*PtResolution*nt.FatJet_p4()[FatJet_index].pt()))) {
                                isGenMatched = 1;
                                JER_corrFactor = std::max(zero, 1 + (JERSF     - 1) * dPt / (nt.FatJet_p4()[FatJet_index].pt()));
                        }
                }
                if (!isGenMatched && JERSF>1) {
                        double sigma = std::sqrt(JERSF * JERSF - 1) * PtResolution;
                        TRandom3 rnd_;
                        // rnd_.SetSeed(3);
                        JER_corrFactor = 1 + rnd_.Gaus(0, sigma);
                }
        }
        return JER_corrFactor;
}

float Process_1Lep2fatjet_JER(Int_t FatJet_index){


// res_sf = JME::JetResolutionScaleFactor("/eos/user/q/qiguo/NanoAOD/updateJER/JRDatabase/textFiles/Summer16_25nsV1_MC/Summer16_25nsV1_MC_SF_AK8PFPuppi.txt");

// cout<<"produce JME::JetResolution"<<endl;
JME::JetParameters jetParam;
jetParam.setJetPt(nt.FatJet_p4()[FatJet_index].pt()).setJetEta(nt.FatJet_p4()[FatJet_index].eta()).setRho(nt.fixedGridRhoFastjetAll());

float PtResolution = ana.resolution.getResolution(jetParam);
float JERSF_temp        = ana.res_sf.getScaleFactor(jetParam);

float first_get_JER_corr = Process_1Lep2fatjet_get_JER_corr(JERSF_temp, true, FatJet_index, 0.8, PtResolution);
// cout<<"first_get_JER_corr"<<first_get_JER_corr<<endl;

return first_get_JER_corr;

}

// JER; Jet

double Process_1Lep2fatjet_get_JER_corr_AK4(float JERSF, bool isMC, Int_t Jet_index, double conSize, float PtResolution){
        double JER_corrFactor = 1.;
        if(isMC) {
                bool isGenMatched = 0;
                int genJetAK8Idx = nt.Jet_genJetIdx()[Jet_index];
                int numofgenjet = nt.Jet_genJetIdx().size();
                int un_int_zero(0);
                float zero(0.);

                if (genJetAK8Idx>=un_int_zero && genJetAK8Idx < numofgenjet) {
                        TLorentzVector jetp4, genjetp4;
                        jetp4.SetPtEtaPhiE(nt.Jet_p4()[Jet_index].pt(), nt.Jet_p4()[Jet_index].eta(), nt.Jet_p4()[Jet_index].phi(), nt.Jet_p4()[Jet_index].energy());
                        genjetp4.SetPtEtaPhiE(nt.GenJet_p4()[genJetAK8Idx].pt(), nt.GenJet_p4()[genJetAK8Idx].eta(), nt.GenJet_p4()[genJetAK8Idx].phi(), nt.GenJet_p4()[genJetAK8Idx].energy());
                        float dR = jetp4.DeltaR(genjetp4);
                        float dPt = nt.Jet_p4()[Jet_index].pt()-nt.GenJet_p4()[genJetAK8Idx].pt();
                        if ((dR<conSize/2.0)&&(std::abs(dPt)<(3*PtResolution*nt.Jet_p4()[Jet_index].pt()))) {
                                isGenMatched = 1;
                                JER_corrFactor = std::max(zero, 1 + (JERSF     - 1) * dPt / (nt.Jet_p4()[Jet_index].pt()));
                        }
                }
                if (!isGenMatched && JERSF>1) {
                        double sigma = std::sqrt(JERSF * JERSF - 1) * PtResolution;
                        TRandom3 rnd_;
                        // rnd_.SetSeed(3);
                        JER_corrFactor = 1 + rnd_.Gaus(0, sigma);
                }
        }
        return JER_corrFactor;
}

double Process_1Lep2fatjet_get_JER_corr_addtional_AK4(float JERSF, bool isMC, Int_t Jet_index, double conSize, float PtResolution){
        double JER_corrFactor = 1.;
        double mindr = 0.2;
        int genJetAK8Idx = -1;
        ana.tx.setBranch<int>("Process_1Lep2fatjet_nGenJet",nt.nGenJet() );
        
        LorentzVector rawJetP4= LorentzVector(ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[Jet_index].pt(),ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[Jet_index].eta(), ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[Jet_index].phi(), 0.);
        if(isMC) {
                bool isGenMatched = 0;
                float min_DR_tmp = 2*mindr;
                for(unsigned int igenjet = 0; igenjet < nt.GenJet_p4().size(); ++igenjet){
                    float min_DR_i = RooUtil::Calc::DeltaR(rawJetP4,nt.GenJet_p4()[igenjet]);
                    if(min_DR_i<min_DR_tmp && min_DR_i<mindr){
                        min_DR_tmp = mindr;
                        genJetAK8Idx = igenjet;
                    }
                }

                ana.tx.pushbackToBranch<int>("Process_1Lep2fatjet_addak4_GenJetid", genJetAK8Idx);
                if(genJetAK8Idx>=0){
                ana.tx.pushbackToBranch<float>("Process_1Lep2fatjet_addak4_GenJetpt", nt.GenJet_pt()[genJetAK8Idx]);
                }

                int un_int_zero(0);
                int numofgenjet = nt.Jet_genJetIdx().size();
                float zero(0.);

                if (genJetAK8Idx>=un_int_zero && genJetAK8Idx < numofgenjet) {
                        TLorentzVector jetp4, genjetp4;
                        jetp4.SetPtEtaPhiE(rawJetP4.pt(), rawJetP4.eta(), rawJetP4.phi(), rawJetP4.energy());
                        genjetp4.SetPtEtaPhiE(nt.GenJet_p4()[genJetAK8Idx].pt(), nt.GenJet_p4()[genJetAK8Idx].eta(), nt.GenJet_p4()[genJetAK8Idx].phi(), nt.GenJet_p4()[genJetAK8Idx].energy());
                        float dR = jetp4.DeltaR(genjetp4);
                        float dPt = jetp4.Pt()-genjetp4.Pt();
                        if ((dR<conSize/2.0)&&(std::abs(dPt)<(3*PtResolution*jetp4.Pt()))) {
                                isGenMatched = 1;
                                JER_corrFactor = std::max(zero, 1 + (JERSF     - 1) * dPt / float(jetp4.Pt()));
                        }
                }
                if (!isGenMatched && JERSF>1) {
                        double sigma = std::sqrt(JERSF * JERSF - 1) * PtResolution;
                        TRandom3 rnd_;
                        // rnd_.SetSeed(3);
                        JER_corrFactor = 1 + rnd_.Gaus(0, sigma);
                }
        }
        return JER_corrFactor;
}

float Process_1Lep2fatjet_JER_AK4(Int_t Jet_index){


// cout<<"produce JME::JetResolution"<<endl;
JME::JetParameters jetParam;
jetParam.setJetPt(nt.Jet_p4()[Jet_index].pt()).setJetEta(nt.Jet_p4()[Jet_index].eta()).setRho(nt.fixedGridRhoFastjetAll());

float PtResolution = ana.resolution_AK4.getResolution(jetParam);
float JERSF_temp        = ana.res_sf_AK4.getScaleFactor(jetParam);

float first_get_JER_corr = Process_1Lep2fatjet_get_JER_corr_AK4(JERSF_temp, true, Jet_index, 0.4, PtResolution);
// cout<<"first_get_JER_corr"<<first_get_JER_corr<<endl;

return first_get_JER_corr;

}

float Process_1Lep2fatjet_JER_AK4_addtional_AK4(Int_t Jet_index){


// cout<<"produce JME::JetResolution"<<endl;
JME::JetParameters jetParam;

LorentzVector rawJetP4= LorentzVector(ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[Jet_index].pt(),ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[Jet_index].eta(), ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[Jet_index].phi(), 0.);

jetParam.setJetPt(rawJetP4.pt()).setJetEta(rawJetP4.eta()).setRho(nt.fixedGridRhoFastjetAll());

float PtResolution = ana.resolution_AK4.getResolution(jetParam);
float JERSF_temp        = ana.res_sf_AK4.getScaleFactor(jetParam);

float first_get_JER_corr = Process_1Lep2fatjet_get_JER_corr_addtional_AK4(JERSF_temp, true, Jet_index, 0.4, PtResolution);
// cout<<"first_get_JER_corr"<<first_get_JER_corr<<endl;

return first_get_JER_corr;

}


// JER; MET
void Process_1Lep2fatjet_allAK4_JERCorr_ForMET(){

    for (unsigned int ijet = 0; ijet < nt.Jet_p4().size(); ++ijet){
        TLorentzVector JERcorrectedP4,JERuncorrectedP4;
        int ijetid = ijet;
        JERuncorrectedP4.SetPtEtaPhiE(nt.Jet_p4()[ijetid].pt(),nt.Jet_p4()[ijetid].eta(),nt.Jet_p4()[ijetid].phi(),nt.Jet_p4()[ijetid].energy());
        JERcorrectedP4 = Process_1Lep2fatjet_JER_AK4(ijetid) * JERuncorrectedP4;
        ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_JER_ForMET", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()));
        ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_unJER_ForMET", LorentzVector(JERuncorrectedP4.Pt(),JERuncorrectedP4.Eta(),JERuncorrectedP4.Phi(),JERuncorrectedP4.M()));
    
    }
    for (unsigned int ijet = 0; ijet < nt.CorrT1METJet_rawPt().size(); ++ijet){
        TLorentzVector JERcorrectedP4,JERuncorrectedP4;
        int ijetid = ijet;
        LorentzVector rawJetP4= LorentzVector(ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[ijetid].pt(),ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[ijetid].eta(), ana.tx.getBranchLazy<vector<LorentzVector>>("addtional_ak4jet_p4_unJER_ForMET")[ijetid].phi(), 0.);
        JERuncorrectedP4.SetPtEtaPhiE(rawJetP4.pt(),rawJetP4.eta(),rawJetP4.phi(),rawJetP4.energy());
        JERcorrectedP4 = Process_1Lep2fatjet_JER_AK4_addtional_AK4(ijetid) * JERuncorrectedP4;
        ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_JER_ForMET", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()));
        ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_unJER_ForMET", LorentzVector(JERuncorrectedP4.Pt(),JERuncorrectedP4.Eta(),JERuncorrectedP4.Phi(),JERuncorrectedP4.M()));

    }
}

// void Process_1Lep2fatjet_allAK4_JERCorr_ForMET(){

//     for (unsigned int ijet = 0; ijet < nt.Jet_p4().size(); ++ijet){
//         TLorentzVector JERcorrectedP4,JERuncorrectedP4;
//         int ijetid = ijet;
//         JERuncorrectedP4.SetPtEtaPhiE(nt.Jet_p4()[ijetid].pt(),nt.Jet_p4()[ijetid].eta(),nt.Jet_p4()[ijetid].phi(),nt.Jet_p4()[ijetid].energy());
//         JERcorrectedP4 = Process_1Lep2fatjet_JER_AK4(ijetid) * JERuncorrectedP4;
//         ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_JER_ForMET", LorentzVector(JERcorrectedP4.Pt(),JERcorrectedP4.Eta(),JERcorrectedP4.Phi(),JERcorrectedP4.M()));
//         ana.tx.pushbackToBranch<LorentzVector>("ak4jet_p4_unJER_ForMET", LorentzVector(JERuncorrectedP4.Pt(),JERuncorrectedP4.Eta(),JERuncorrectedP4.Phi(),JERuncorrectedP4.M()));
    
//     }
// }


void Process_1Lep2fatjet_METCorr(){
    float corrEx_MET_JER = 0;
    float corrEy_MET_JER = 0;
    for (unsigned int ijet = 0; ijet < ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_JER_ForMET").size(); ++ijet){
        corrEx_MET_JER += ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_unJER_ForMET")[ijet].px() - ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_JER_ForMET")[ijet].px();
        corrEy_MET_JER += ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_unJER_ForMET")[ijet].py() - ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_JER_ForMET")[ijet].py();

       ana.tx.pushbackToBranch<float>("ak4jet_p4_JER_ForMET_correx",(ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_unJER_ForMET")[ijet].px() - ana.tx.getBranchLazy<vector<LorentzVector>>("ak4jet_p4_JER_ForMET")[ijet].px()));
    }
    TVector2 rawMET_;
    rawMET_.SetMagPhi (nt.MET_pt(), nt.MET_phi() );
    float rawPx = rawMET_.Px();
    float rawPy = rawMET_.Py();
    float pxcorr = rawPx+corrEx_MET_JER;
    float pycorr = rawPy+corrEy_MET_JER;
    float et     = std::hypot(pxcorr,pycorr);
    ana.tx.setBranch<float>("MET_et_JER",et );
    ana.tx.setBranch<int>("Process_1Lep2fatjet_nCorrT1METJet",nt.nCorrT1METJet() );


}

void Process_1Lep2fatjet_corrJEC_ForAdditionalJet(){
    double jetCorrEtaMax           = 9.9;
    float jetCorrFactor = 1;

    for (unsigned int ijet = 0; ijet < nt.CorrT1METJet_rawPt().size(); ++ijet){
        // float rawJetPt = nt.CorrT1METJet_rawPt()[ijet];
        // float rawJeteta = nt.CorrT1METJet_rawPt()[ijet];
        // float rawJetPhi = nt.CorrT1METJet_rawPt()[ijet];

        

        LorentzVector rawJetP4= LorentzVector(nt.CorrT1METJet_rawPt()[ijet], nt.CorrT1METJet_eta()[ijet], nt.CorrT1METJet_phi()[ijet], 0);
        
        float jetArea = nt.CorrT1METJet_area()[ijet];
        float jetrho = nt.fixedGridRhoFastjetAll();
        int nVtx = nt.PV_npvs();

        jetCorrFactor = 1;
        if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
                ana.jecAK4_->setJetEta( rawJetP4.Eta() );
                ana.jecAK4_->setJetPt ( rawJetP4.Pt() );
                ana.jecAK4_->setJetE  ( rawJetP4.energy() );
                ana.jecAK4_->setJetPhi( rawJetP4.Phi()    );
                ana.jecAK4_->setJetA  ( jetArea );
                ana.jecAK4_->setRho   ( jetrho );
                ana.jecAK4_->setNPV   ( nVtx );
                jetCorrFactor = ana.jecAK4_->getCorrection();
            }
        ana.tx.pushbackToBranch<LorentzVector>("addtional_ak4jet_p4_unJEC_ForMET", LorentzVector(rawJetP4.Pt(),rawJetP4.Eta(),rawJetP4.Phi(),rawJetP4.M()));

        LorentzVector JECJetP4= LorentzVector(nt.CorrT1METJet_rawPt()[ijet]*jetCorrFactor, nt.CorrT1METJet_eta()[ijet], nt.CorrT1METJet_phi()[ijet], 0.);
        ana.tx.pushbackToBranch<LorentzVector>("addtional_ak4jet_p4_unJER_ForMET", LorentzVector(JECJetP4.Pt(),JECJetP4.Eta(),JECJetP4.Phi(),JECJetP4.M()));



    
    }
}


float Process_1Lep2fatjet_Fatjet_SD_sub(Int_t FatJet_subJetIdx1,Int_t FatJet_subJetIdx2){
    TLorentzVector subjet1,subjet2,sum_p4;
    subjet1.SetPtEtaPhiE(nt.SubJet_p4()[FatJet_subJetIdx1].pt(),nt.SubJet_p4()[FatJet_subJetIdx1].eta(),nt.SubJet_p4()[FatJet_subJetIdx1].phi(),nt.SubJet_p4()[FatJet_subJetIdx1].energy());
    subjet2.SetPtEtaPhiE(nt.SubJet_p4()[FatJet_subJetIdx2].pt(),nt.SubJet_p4()[FatJet_subJetIdx2].eta(),nt.SubJet_p4()[FatJet_subJetIdx2].phi(),nt.SubJet_p4()[FatJet_subJetIdx2].energy());
    sum_p4 = subjet1+subjet2;
    return sum_p4.M();
}

float Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(Int_t FatJet_subJetIdx1,Int_t FatJet_subJetIdx2){
    TLorentzVector subjet1,subjet2,sum_p4;
    if(FatJet_subJetIdx1>=0&&FatJet_subJetIdx2>=0){
    float pt1 = nt.SubJet_p4()[FatJet_subJetIdx1].pt()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx1]);
    float pt2 = nt.SubJet_p4()[FatJet_subJetIdx2].pt()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx2]);
    float mass1 = nt.SubJet_p4()[FatJet_subJetIdx1].mass()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx1]);
    float mass2 = nt.SubJet_p4()[FatJet_subJetIdx2].mass()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx2]);
    subjet1.SetPtEtaPhiM(pt1,nt.SubJet_p4()[FatJet_subJetIdx1].eta(),nt.SubJet_p4()[FatJet_subJetIdx1].phi(),mass1);
    subjet2.SetPtEtaPhiM(pt2,nt.SubJet_p4()[FatJet_subJetIdx2].eta(),nt.SubJet_p4()[FatJet_subJetIdx2].phi(),mass2);
    sum_p4 = subjet1+subjet2;
    return sum_p4.M();
    }
    else{
        return -99998.0;
    }
}



void Process_1Lep2fatjet_test(){
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu){
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_pt_test",nt.Muon_p4()[imu].pt());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_eta_test",nt.Muon_p4()[imu].eta());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_iso_test",nt.Muon_pfRelIso03_all()[imu]);
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_iso2_test",nt.Muon_pfRelIso04_all()[imu]);
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_iso3_test",nt.Muon_tkRelIso  ()[imu]);
        
    }
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu)
    {
        int temp_id = nt.Muon_highPtId()[imu]; // have to switch to int
        // cout<<"nt.Muon_highPtId()[imu]:"<<temp_id<<" id origin:  "<<nt.Muon_highPtId()[imu]<<endl;
        if (not (temp_id==2 )) continue;
        if (not (nt.Muon_p4()[imu].pt()>50.0 )) continue;
        // std::cout<<"find one highPt muon"<<std::endl;
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_highptMuon_pt_test",nt.Muon_p4()[imu].pt());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_highptMuon_eta_test",nt.Muon_p4()[imu].eta());
    }
    for (unsigned int usenumber1 = 0; usenumber1 < nt.FatJet_p4().size(); ++usenumber1)
    {
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_fatjet_pt_test",nt.FatJet_p4()[usenumber1].pt());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_fatjet_eta_test",nt.FatJet_p4()[usenumber1].eta());
    }
}

void Process_1Lep2fatjet_MET(){
    // Met
    double MET_et,MET_phi;
    MET_et = nt.MET_pt();
    MET_phi = nt.MET_phi();
    ana.tx.setBranch<float>("MET_et", MET_et );
    ana.tx.setBranch<float>("MET_phi", MET_phi );
    if (not nt.isData()){
    Process_1Lep2fatjet_corrJEC_ForAdditionalJet();
    Process_1Lep2fatjet_allAK4_JERCorr_ForMET();
    Process_1Lep2fatjet_METCorr();
}

}


#endif