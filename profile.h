    
//const Int_t nev = -1; // number of events to read, -1 - until CTRL+C
const Int_t nev = 1000000; // Only nev events to read
const Int_t fRunId = 1;

// *********************************** //
// PLEASE CHANGE THE EXPERIMENT NUMBER //
// *********************************** //
const Int_t expId = 467;               // select experiment: 444 or 467
// *********************************** //


std::ifstream RunList("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/RunSummary.csv", std::ios::in);
TString dir_rawfile = "/u/taniuchi/s467/lmd_stitched/";
TString dir_output = "./rootfiles/rootfiletmp/fragment_Nov2021/";
TString dir = gSystem->Getenv("VMCWORKDIR");
TString sofiacaldir = dir + "/sofia/macros/s467_ryotani/parameters/";
TString ucesb_dir = getenv("UCESB_DIR");
TString upexps_dir = ucesb_dir + "/../upexps/";
TString ucesb_path = upexps_dir + "/202002_s467/202002_s467 --allow-errors --input-buffer=100Mi";

// store data or not ------------------------------------
Bool_t fCal_level_califa = true;  // set true if there exists a file with the calibration parameters
Bool_t NOTstoremappeddata = false; // if true, don't store mapped data in the root file
Bool_t NOTstorecaldata = false;    // if true, don't store cal data in the root file
Bool_t NOTstorehitdata = false;    // if true, don't store hit data in the root file
    
// Setup: Selection of detectors ------------------------
Bool_t fFrs = true;      // FRS for production of exotic beams (just scintillators)
Bool_t fFrsTpcs = false; // Tpcs at FRS (S2) for scintillator calibration in position
Bool_t fFrsMws = false;  // MWs at FRS (S8) for beam position
Bool_t fFrsSci = true;   // Start: Plastic scintillators at FRS
Bool_t fMwpc0 = true;    // MWPC0 for tracking at entrance of Cave-C
Bool_t fMusic = true;    // R3B-Music: Ionization chamber for charge-Z
Bool_t fSci = true;      // Start: Plastic scintillator for ToF
//Bool_t fAms = false;     // AMS tracking detectors
Bool_t fCalifa = false;  // Califa calorimeter
Bool_t fMwpc1 = true;    // MWPC1 for tracking of fragments in front of target
Bool_t fMwpc2 = true;    // MWPC2 for tracking of fragments before GLAD
Bool_t fTwim = true;     // Twim: Ionization chamber for charge-Z of fragments
Bool_t fMwpc3 = true;    // MWPC3 for tracking of fragments behind GLAD
Bool_t fTofW = true;     // ToF-Wall for time-of-flight of fragments behind GLAD
Bool_t fScalers = false;  // SIS3820 scalers at Cave C
//Bool_t fNeuland = true;  // NeuLAND for neutrons behind GLAD
//Bool_t fTracking = true; // Tracking of fragments inside GLAD

Bool_t fSkip_tpat0 = true;

// Calibration files ------------------------------------
// Parameters for CALIFA mapping
TString califamapdir = dir + "/macros/r3b/unpack/s467/califa/parameters/";
TString califamapfilename = califamapdir + "CALIFA_mapping.par";
// Parameters for CALIFA calibration in keV
TString califadir = dir + "/macros/r3b/unpack/s467/califa/parameters/";
TString califacalfilename = califadir + "Califa_Cal8Feb2020.root";
