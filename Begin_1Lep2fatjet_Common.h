#ifndef Begin_1Lep2fatjet_Common_h
#define Begin_1Lep2fatjet_Common_h

void Begin_1Lep2fatjet_test(){
    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_Muon_idxs_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_pt_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_eta_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_iso_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_iso2_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_iso3_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_highptMuon_pt_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_highptMuon_eta_test");

    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_fatjet_pt_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_fatjet_eta_test");



}

#endif


