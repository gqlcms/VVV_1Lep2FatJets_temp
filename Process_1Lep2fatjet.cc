#include "Process_1Lep2fatjet.h"
// #include "Process_1Lep2fatjet_common.h"
// #include "Process_1Lep2fatjet_OldLeptonID.h"
// #include "Process_1Lep2fatjet_NewLeptonID.h"
#include "Process_1Lep2fatjet_UpdateJER.h"
#include "TMath.h"


void Process_1Lep2fatjet()
{
    //==============================================
    // Process_1Lep2fatjet:
    // This function gets called during the event looping.
    // This is where one sets the variables used for the category 1Lep2fatjet.
    //==============================================

    Process_1Lep2fatjet_UpdateJER_AK4();
    Process_1Lep2fatjet_Weight_HLT();
    Process_1Lep2fatjet_Jet();
    Process_1Lep2fatjet_fatJet();
    Process_1Lep2fatjet_Selection();

    // Set variables used in this category.
    // If histograms are booked with these variables the histograms will be filled automatically.
    // Please follow the convention of <category>_<varname> structure.

    // Example of reading from Nano
    // std::vector<LorentzVector> electron_p4s = nt.Electron_p4(); // nt is a global variable that accesses NanoAOD
    // std::vector<float> electron_mvaTTH = nt.Electron_mvaTTH(); // electron ttH MVA scores from NanoAOD
    // Semi-complete list of NanoAOD for 102X can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html
    // Also consult here: https://github.com/cmstas/NanoTools/blob/d641a6d6c1aa9ecc8094a1af73d5e1bd7d6502ab/NanoCORE/Nano.h#L4875 (if new variables are added they may show up in master)
}
