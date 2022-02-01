using namespace std;

//The number of events you want to process
Long64_t nEvents=-1; //Set to -1 to run over the whole file


//Conditions for the cuts. Here you can adjust the macro to your region of interest
const Double_t EMINPROTON=10000;        //Min Hit energy in keV for a CALIFA hit
//Spread is currently not used
//const Double_t THETASPREAD = 20.;       //Spread around the correlation of the two protons, having a theta of 80 +/- SPREAD degree
//const Double_t PHISPREAD = 180.;        //Spread around the correlation of the two protons, having a phi of 180 +/- SPREAD degree
const Double_t BETA = 0.7652;           //Beta value to calculate the the lab to rest frame energy, depending on the hit crystal. At the moment
                                        //this is just rough, with one theta angle per crystal, in the middle of the crystal axis
/// Ryo: Beta should be the one at the center of the mass of two protons?
/*
const Double_t R3BMusicCutMin = 23.2;       //Min charge for R3BMusic (before target)
const Double_t R3BMusicCutMax = 24.5;       //Max charge for R3BMusic (before target)
const Double_t twimMusicCutMin = 16.4;      //Min charge for TwimMusic (after target)
const Double_t twimMusicCutMax = 17.2;      //Max charge for TwimMusic (after target)
*/

//Input File
const char* RUNLISTNAME="./RunSummary.csv";
const Int_t FRS=13; // calib:0, ToFCalib:6-8, 40Ca:9, 39Ca:10, 38Ca:11,12, 50Ca:13, ToFWcalib:14
//
//const int posmin=1325, posmax=1424; const TString targetname = "empty";
const int posmin=539, posmax=539; const TString targetname = "ch2-24mm";
//const int posmin=362, posmax=362; const TString targetname = "carbon";
//const int posmin=893, posmax=893; const TString targetname = "PP";
//
const char* FILENAME="s467_filltree_Setting13_RUN_25Jan.root"; // RUN will be replaced with runnumber later
const char* INPUTDIR="rootfiles/rootfile_land/";

//Output File
const char* OUTFILENAME="Setting13_25Jan";
const char* OUTPUTDIR="rootfiles/rootfile_land/";

// Defining the result file and creating some histos
TFile* resultFile;
TH1F* histCalifaHighEnTheta= new TH1F("histCalifaHighEnTheta", "The theta for 2xHighEn & phi~180 Ar->Cl; Theta; ",180,0,180);
TH2F* histCalifaTheta1vsTheta2 = new TH2F("histCalifaTheta1vsTheta2","Theta1 vs Theta2 for two high energy hits Ar->Cl; Theta1 ; Theta2",180,0,180,180,0,180);
TH2F* histCalifaPhi1vsPhi2 = new TH2F("histCalifaPhi1vsPhi2","Phi1 vs Phi2 for two high energy hits Ar->Cl; Phi1 ; Phi2",360,-180,180,360,-180,180);


// ===== REQUIRED FOR ANALYSIS LOOP ===============================
// Processing the input file
//TFile* eventFile; 
//TTree* eventTree;
TChain* eventTree;
TTree* OutTree;
TString InputFileName, OutputFileName;
//char inputFile[300];    //Holds the name of the input file.
char outputFile[300];   //Holds the name of the output file.

//For cleaning up variables
const Double_t DEFAULTDOUBLE = NAN;
const Int_t DEFAULTINTEGER=0;
// Integers for the loops through the data
Int_t nCalifaHits=DEFAULTINTEGER;
//Int_t nFragmentHits=DEFAULTINTEGER;
//Int_t nR3BMusicHits=DEFAULTINTEGER;
Int_t nFrsHits=DEFAULTINTEGER, nFragmentHits=DEFAULTINTEGER;

//R3BEventHeaderPropagator* EventHeader;
//FairEventHeader* EventHeader;
//Int_t ftpat;//=DEFAULTINTEGER;
// Trying to obtain the tpat info but not working right now.

R3BSofFrsData* FrsHit;
R3BSofTrackingData* FragmentHit; 
R3BCalifaHitData* califaHit;        //Data on hit level. It also includes a cluster algorithm

// TClonesArrays to hold the data structures
TClonesArray* fEventHeader;
TClonesArray* fCalifaHitData;
//TClonesArray* fR3BMusicHitData;
TClonesArray* fSofFrsData;
TClonesArray* fSofTrackingData;

// Variables to hold the values during the loop
// CALIFA
Double_t energy[2]={DEFAULTDOUBLE};    //Energy of proton 1 in the rest frame
Double_t theta[2]={DEFAULTDOUBLE};     //Theta of proton 1
Double_t phi[2]={DEFAULTDOUBLE};       //Phi of proton 1
Int_t mult=DEFAULTINTEGER;
const int MAXMULT =10;
Int_t mult_dist[MAXMULT]={DEFAULTINTEGER};

Double_t thetaP2P=DEFAULTDOUBLE;            //The theta between proton 1 and 2 (only true, if phi is ~180 degree)
Double_t phiP2P=DEFAULTDOUBLE;              //The phi between proton 1 and 2
Double_t sumEnergy=DEFAULTDOUBLE;           //The sum energy in the rest frame of the two proton hits

// MUSICs
Double_t r3bMusicCharge=DEFAULTDOUBLE;      //Charge of the ion before the target
Double_t TwimCharge=DEFAULTDOUBLE;     //Charge of the ion after the target

//Addition
Int_t multiEventFlag=0;                     //If 1, we have more then two proton hits and the event is discarded
Int_t multiEventCounter=DEFAULTINTEGER;     //For curiosity, count the amount of events, which have more then two proton hits (this events will not be used).
Int_t musicsChargeMatch=DEFAULTINTEGER;     //Count, how often we have the right charges
Int_t califaP2PCandidate=DEFAULTINTEGER;    //Count, how often a good p2p event was found in CALIFA

//==================================================


// Forward decleration of functions
void analysisLoop(int max_events=-1);
void findP2PCandidate();
void writeHistos();
void SetTreeBranches();
void resetVariables();
void loadtrees();
Double_t dopplerCorrection(Double_t theta, Double_t BETA);


