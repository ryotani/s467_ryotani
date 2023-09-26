const Int_t nev = -1; // number of events to read, -1 - until CTRL+C
//const Int_t nev = 10000; // Only nev events to read

const Int_t fRunId = 1;
//
// Select tpat condition : 1=min bias, 2=califaOR(p), 3=neuland, 4=califa&neuland,
// 8=califa off, 9=neuland off.
const Int_t min_tpat = 1;
const Int_t max_tpat = 4;

// *********************************** //
// PLEASE CHANGE THE EXPERIMENT NUMBER //
// *********************************** //
const Int_t expId = 467;               // select experiment: 444 or 467
// *********************************** //

TString dir_ana = "./";
TString dir = gSystem->Getenv("VMCWORKDIR");
//
TString s_runlist = dir_ana + "RunSummary.csv";
std::ifstream RunList(s_runlist, std::ios::in);
TString lmdfilename = "/u/taniuchi/s467/lmd_stitched/mainRUNNUM_*.lmd";
TString filename = lmdfilename; // For backward compatibility
TString dir_output = dir_ana + "rootfiles/Feb2023_FSnov22/";
TString outfile_template = "s467_filltree_neuland_SettingSETID_RUNNUM_DAYMONTH.root"; // SETID, RUNNUM, DAY, and MONTH will be replaced accordingly.
TString sofiacaldir = dir_ana + "parameters/";
//
TString ucesb_dir = getenv("UCESB_DIR");
//TString upexps_dir = "/u/land/latar/upexps";
TString upexps_dir = ucesb_dir + "/../upexps/";
//TString ucesb_path = upexps_dir + "/202002_s467_nov22/202002_s467 --allow-errors --input-buffer=100Mi";
TString ucesb_path = upexps_dir + "/202002_s467/202002_s467 --allow-errors --input-buffer=100Mi";
TString ntuple_options = "RAW"; //"RAW,time-stitch=1000" // For no stitched data

// store data or not ------------------------------------
Bool_t fCal_level_califa = true;  // set true if there exists a file with the calibration parameters
Bool_t NOTstoremappeddata = false; // if true, don't store mapped data in the root file
Bool_t NOTstorecaldata = false;    // if true, don't store cal data in the root file
Bool_t NOTstorehitdata = false;    // if true, don't store hit data in the root file
    
// Setup: Selection of detectors ------------------------
Bool_t fFrsTpcs = false; // Tpcs at FRS (S2) for scintillator calibration in position
Bool_t fFrsSci = false;   // Start: Plastic scintillators at FRS
Bool_t fMwpc0 = true;    // MWPC0 for tracking at entrance of Cave-C
Bool_t fMusic = true;    // R3B-Music: Ionization chamber for charge-Z
Bool_t fSci = true;      // Start: Plastic scintillator for ToF
Bool_t fCalifa = true;  // Califa calorimeter
Bool_t fMwpc1 = true;    // MWPC1 for tracking of fragments in front of target
Bool_t fMwpc2 = true;    // MWPC2 for tracking of fragments before GLAD
Bool_t fTwim = true;     // Twim: Ionization chamber for charge-Z of fragments
Bool_t fMwpc3 = true;    // MWPC3 for tracking of fragments behind GLAD
Bool_t fTofW = true;     // ToF-Wall for time-of-flight of fragments behind GLAD
Bool_t fScalers = false;  // SIS3820 scalers at Cave C
Bool_t fNeuland = true;  // NeuLAND for neutrons behind GLAD


// Calibration files ------------------------------------
TString califamapdir = dir + "/macros/r3b/unpack/s467/califa/parameters/";
TString califamapfilename = califamapdir + "CALIFA_mapping_Feb82020.par";// not used
// Parameters for CALIFA calibration in keV
TString califadir = dir_ana + "parameters/";
TString califacalfilename = califadir + "Califa_CalibParam25012022_60Co_Incomplete.root";
