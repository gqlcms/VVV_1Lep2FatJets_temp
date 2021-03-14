#include "Begin_1Lep2fatjet.h"


void Begin_1Lep2fatjet()
{
    //==============================================
    // Begin_1Lep2fatjet:
    // This function gets called prior to the event looping.
    // This is where one declares variables, histograms, event selections for the category 1Lep2fatjet.
    //==============================================


    // Create variables used in this category.
    // Please follow the convention of <category>_<varname> structure.


    // ana.tx.createBranch<>("");

    // Lepton
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonPt" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonEta" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonPhi" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonE" );
    ana.tx.createBranch<int>("Lep1fatJet2_LeptonPDGID" );
    ana.tx.createBranch<float>("Lep1fatJet2_Muon_pfRelIso04_all" );

    // MET
    ana.tx.createBranch<float>("Lep1fatJet2_MET_pt" );
    ana.tx.createBranch<float>("Lep1fatJet2_MET_phi" );
    ana.tx.createBranch<float>("Common_MET_pt_JER" );
    ana.tx.createBranch<float>("Common_MET_phi_JER" );

    // neutrino
    ana.tx.createBranch<float>("Lep1fatJet2_NeutrinoPt" );
    ana.tx.createBranch<float>("Lep1fatJet2_NeutrinoEta" );
    ana.tx.createBranch<float>("Lep1fatJet2_Neutrinophi" );
    ana.tx.createBranch<float>("Lep1fatJet2_NeutrinoE" );

    // Leptonic W
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonicWPt" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonicWEta" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonicWPhi" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonicWM" );
    ana.tx.createBranch<float>("Lep1fatJet2_LeptonicWMt" );

    // Jet
    ana.tx.createBranch<vector<float>>("Common_Jet_pt_JER");
    ana.tx.createBranch<vector<LorentzVector>>("Common_Jet_JER_p4");
    ana.tx.createBranch<vector<float>>("Common_Jet_pt_unJER");
    ana.tx.createBranch<vector<float>>("Common_JER_matching");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_Jet_btagDeepB");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_Jet_btagDeepC");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_Jet_jetId");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_Jet_hadronFlavour");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_Jet_partonFlavour");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_Jet_pt_JER");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_Jet_eta");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_Jet_phi");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_Jet_e");
   
    // fatJet
    ana.tx.createBranch<vector<LorentzVector>>("Common_fatJet_JER_p4");
    ana.tx.createBranch<vector<float>>("Common_fatJet_pt_JER");
    ana.tx.createBranch<vector<float>>("Common_fatJet_pt_unJER");
    ana.tx.createBranch<vector<float>>("Common_fatJER_matching");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_FatJet_jetId");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_FatJet_pt_JER");

    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_ptJECJER");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_ptJEC");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_eta");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_phi");
    ana.tx.createBranch<int>("Lep1fatJet2_jetAK8puppi_tightid");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_sd");
    ana.tx.createBranch<float>("Lep1fatJet2_jet_rawFactor");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrW");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnW");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrZ");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnZ");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrTop");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnTop");
    
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_ptJECJER_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_ptJEC_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_eta_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_phi_2");
    ana.tx.createBranch<int>("Lep1fatJet2_jetAK8puppi_tightid_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_sd_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jet_rawFactor_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrW_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnW_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrZ_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnZ_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnDecorrTop_2");
    ana.tx.createBranch<float>("Lep1fatJet2_jetAK8puppi_dnnTop_2");

    // Gen
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_ptgenzl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_etagenzl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_phigenzl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_energygenzl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q1_pt");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q1_eta");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q1_phi");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q1_e");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_genz_q1_pdg");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q2_pt");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q2_eta");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q2_phi");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genz_q2_e");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_genz_q2_pdg");

    ana.tx.createBranch<vector<float>>("Lep1fatJet2_ptgenwl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_etagenwl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_phigenwl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_energygenwl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q1_pt");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q1_eta");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q1_phi");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q1_e");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_genw_q1_pdg");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q2_pt");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q2_eta");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q2_phi");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_genw_q2_e");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_genw_q2_pdg");

    ana.tx.createBranch<vector<float>>("Lep1fatJet2_ptgentopl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_etagentopl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_phigentopl");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_energygentopl");

    // selection
    ana.tx.createBranch<int>("Lep1fatJet2_Preselection");

    // weight, HLT
    ana.tx.createBranch<float>("Lep1fatJet2_L1PreFiringWeight_Nom");
    ana.tx.createBranch<float>("Lep1fatJet2_Pileup_nTrueInt");
    ana.tx.createBranch<float>("Lep1fatJet2_Pileup_nPU");

    ana.tx.createBranch<int>("Lep1fatJet2_HLT_Mu50");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_IsoMu24");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_OldMu100");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_TkMu100");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_IsoMu27");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_Ele27");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_Ele32_WPTight_Gsf_L1DoubleEG");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_Ele35_WPTight_Gsf");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_Photon200");
    ana.tx.createBranch<int>("Lep1fatJet2_HLT_Ele32_WPTight_Gsf");

    // check Ntuple
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_CleanJet_pt");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_MatchGenJet_pt");
    ana.tx.createBranch<vector<int>>("Lep1fatJet2_MatchGenJet_num");
    ana.tx.createBranch<vector<float>>("Lep1fatJet2_JERSF");



    ana.res.loadResolutionFile("Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt");
    ana.res.loadScaleFactorFile("Summer16_25nsV1_MC_SF_AK4PFchs.txt");

    ana.resfatJet.loadResolutionFile("Summer16_25nsV1_MC_PtResolution_AK8PFPuppi.txt");
    ana.resfatJet.loadScaleFactorFile("Summer16_25nsV1_MC_SF_AK8PFPuppi.txt");

    
    

    ana.cutflow.getCut("CommonCut");
    // ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has1Lepton_new", [&](){
    //     if (not ((ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs").size()+ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs").size()) == 1))
    //        return false;
    //     return true;

    // }, UNITY);
    // ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has2FatJet_new", [&]() { return (ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs_new").size() >= 2);}, UNITY);




    // Define selections
    // CommonCut will contain selections that should be common to all categories, starting from this cut, add cuts for this category of the analysis.
    ana.cutflow.getCut("CommonCut");
    // ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_GenHT", [&]() { return ana.tx.getBranch<float>("Common_genHT"); }, UNITY);
//     ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has1Lepton", [&]()
//             {
//                 if (not (ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs").size() == 1))
//                     return false;
//                 ana.tx.setBranch<LorentzVector>("1Lep2fatjet_lep_p4", ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[0]);
//                 ana.tx.setBranch<int>("1Lep2fatjet_lep_pdgid", ana.tx.getBranchLazy<vector<int>>("Common_lep_pdgid")[0]);
//                 return true;
//             }, UNITY);
//     ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has2FatJet", [&]() { return (ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size() >= 2);}, UNITY);
//     ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_nbVeto", [&]() { return (ana.tx.getBranch<int>("Common_nb_tight") == 0);}, UNITY);
//     ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_FatJetMassPresel", [&]() { return (nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] > 20.) and (nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] > 20.);}, UNITY);
//     ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_DeepWProd", [&]() { return (nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]) > 0.5;}, UNITY);
//     ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Preselection", UNITY, UNITY);

//     // Create histograms used in this category.
//     // Please follow the convention of h_<category>_<varname> structure.
//     // N.B. Using nbins of size 180 or 360 can provide flexibility as it can be rebinned easily, as 180, 360 are highly composite numbers.
//     RooUtil::Histograms hists_1Lep2fatjet_new;
        
//     hists_1Lep2fatjet_new.addHistogram("h_1Lep2fatjet_new_event_id", 100, 0, 100, [&]() { return ana.tx.getBranch<unsigned long long>("Common_evt"); } );
//     hists_1Lep2fatjet_new.addHistogram("h_1Lep2fatjet_new_fatjet_msoftdrop_0", 180, 0, 250, [&]() { return nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs_new")[0]]; } );

//     RooUtil::Histograms hists_1Lep2fatjet;
// //    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_lep_pt", 180, 0, 500, [&]() { return ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4").pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_event_id", 100, 0, 100, [&]() { return ana.tx.getBranch<unsigned long long>("Common_evt"); } );
// //    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_lep_pdgid", 4, 10, 14, [&]() { return abs(ana.tx.getBranch<int>("1Lep2fatjet_lep_pdgid")); } );
// //    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_lep_eta", 180, -3, 3, [&]() { return ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4").eta(); } );
// //    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt", 180, 0, 250, [&]() { return RooUtil::Calc::mT(ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4"), ana.tx.getBranch<LorentzVector>("Common_met_p4")); } );
// //    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_pt_lv", 180, 0, 1000, [&]() { return (ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_met_pt", 180, 0, 500, [&]() { return ana.tx.getBranch<LorentzVector>("Common_met_p4").pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_met_phi", 180, -3.1416, 3.1416, [&]() { return ana.tx.getBranch<LorentzVector>("Common_met_p4").phi(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nb_loose", 8, 0, 8, [&]() { return ana.tx.getBranch<int>("Common_nb_loose"); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nb_medium", 8, 0, 8, [&]() { return ana.tx.getBranch<int>("Common_nb_medium"); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nb_tight", 8, 0, 8, [&]() { return ana.tx.getBranch<int>("Common_nb_tight"); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_njet", 15, 0, 15, [&]() { return ana.tx.getBranch<vector<int>>("Common_jet_idxs").size(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nfatjet", 5, 0, 5, [&]() { return ana.tx.getBranch<vector<int>>("Common_fatjet_idxs").size(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_difatjet_mass", 180, 0, 3000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1]).mass(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_difatjet_pt", 180, 0, 1000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1]).pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_eta_0", 180, -5, 5, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].eta(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_eta_1", 180, -5, 5, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].eta(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_pt_0", 180, 0, 1000, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_pt_1", 180, 0, 750, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_msoftdrop_0", 180, 0, 250, [&]() { return nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_msoftdrop_1", 180, 0, 250, [&]() { return nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_genjetmass_0", 180, 0, 250, [&]() { return nt.GenJetAK8_mass()[nt.FatJet_genJetAK8Idx()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_genjetmass_1", 180, 0, 250, [&]() { return nt.GenJetAK8_mass()[nt.FatJet_genJetAK8Idx()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_0", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_1", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_sum", 180, -1, 2, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] + nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_prod", 180, -1, 2, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_0", 180, -1, 1, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_1", 180, -1, 1, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_sum", 180, -1, 2, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] + nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_prod", 180, -1, 2, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_prod_zoom", 180, 0.8, 1, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepT_0", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepT_1", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepT_0", 180, -1, 1, [&]() { return nt.FatJet_deepTag_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepT_1", 180, -1, 1, [&]() { return nt.FatJet_deepTag_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepbbvsL_0", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_bbvsLight()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepbbvsL_1", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_bbvsLight()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_lsf3_0", 180, -1, 1, [&]() { return nt.FatJet_lsf3()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_lsf3_1", 180, -1, 1, [&]() { return nt.FatJet_lsf3()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nBHadrons_0", 5, 0, 5, [&]() { return nt.FatJet_nBHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nBHadrons_1", 5, 0, 5, [&]() { return nt.FatJet_nBHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nCHadrons_0", 5, 0, 5, [&]() { return nt.FatJet_nCHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nCHadrons_1", 5, 0, 5, [&]() { return nt.FatJet_nCHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_0", 180, 0, 1, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_1", 180, 0, 1, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_sum", 180, 0, 2, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] + nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_prod", 180, 0, 1, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_0_0", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_1_0", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_0_1", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_1_1", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_0_0", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt(); else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_1_0", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt(); else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_0_1", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt(); else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_1_1", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt(); else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_ratio_0", 180, -20, 10, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0 and nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]] / nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; else return -999.f; } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_ratio_1", 180, -20, 10, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0 and nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]] / nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; else return -999.f; } );
// /*    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_htSum", 180, 0, 4000, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt() + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt() + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4").pt() + ana.tx.getBranch<LorentzVector>("Common_met_p4").pt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mass_lvJJ", 180, 0, 4000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt_lvJJ", 180, 0, 4000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_pt_lvJJ", 180, 0, 4000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mass_lvJ_0", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mass_lvJ_1", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt_lvJ_0", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mt(); } );
//     hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt_lvJ_1", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mt(); } );
// */

//     // Now book cutflow histogram (could be commented out if user does not want.)
//     // N.B. Cutflow histogramming can be CPU consuming.
//     ana.cutflow.bookCutflows();

//     // Book histograms to cuts that user wants for this category.
//     ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Has2FatJet");
// //    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Has2FatJet_new");
// //    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet_new, "Cut_1Lep2fatjet_Has1Lepton_new");
//     ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet_new, "Cut_1Lep2fatjet_Has2FatJet_new");
// //    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Has1Lepton_new");
//     ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_FatJetMassPresel");
//     ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Preselection");
    
    // ana.tx.fill();
}
