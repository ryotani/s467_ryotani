//#define test
#ifdef test
const int NUMTARGET=2, NUMPADDLE=11, MAXEVENT=100000;
#else
const int NUMTARGET=3, NUMPADDLE=28, MAXEVENT=-1;
#endif
/*
const int MINZ=17, MAXZ=23, MAXA=55, MININCLZ=18, MAXINCLZ=23, BINAOQ=500, BINZ=500;
const double MINAOQ=1.9, MAXAOQ=2.75;
TString frsname = "50Ca";
double zoff =0;
const bool frs13 = true;
*/
const int MINZ=17, MAXZ=23, MAXA=55, MININCLZ=18, MAXINCLZ=23, BINAOQ=500, BINZ=500;
const double MINAOQ=1.65, MAXAOQ=2.5;
TString frsname = "38Ca";
const bool frs13 = false;
double zoff=0.3; // tentative.

const Double_t d_targ[3] = {0, 1.23E+23, 9.85E+22}; // Num / cm^2
const Double_t to_mb = 1.e27;
const Double_t significance = 3.;
//
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TTimer.h"
#include <fstream>
#include <iostream>
using namespace std;
//
TString dir_rootfile = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/fragment/output/" ;
TString filename = "tofw_calib_FRS_TARG_Apr.root"; // TARG will be replaced
#ifdef test
TString outfilename = "incl_out_test_indivfit_30May_FRS_TARG.root";
#else
TString outfilename = "incl_out_indivfit_30May_FRS_TARG.root";
#endif
TString outcsvname = "incl_list_FRS_TARG.csv";
TString targetname[3] = {"empty","carbon","ch2"};
TString reactionname[2] = {"-1n", "-1p"};
TString targetelement[3] = {"carbon", "CH2", "proton"};
TFile *fin[NUMTARGET], *fin2[NUMTARGET], *fout;
TTree *tree[NUMTARGET], *tree2[NUMTARGET];
ofstream fcsv;
//
// Def variables
Long64_t nentry[NUMTARGET]={0};
Long64_t nentry2[NUMTARGET]={0};
Double_t pid_fraggate[NUMTARGET][NUMPADDLE+1][MAXZ][MAXA][4]={NAN}; // 0:Z, 1:sigmaZ, 2:aq, 3:sigmaaq
Double_t trans[NUMTARGET][MAXZ][MAXA]={NAN}, sigma[MAXZ][MAXA][3][2]={NAN}; // 0: C, 1: CH2, 2: H, 0: -1n, 1: -1p
Double_t Etrans[NUMTARGET][MAXZ][MAXA]={NAN}, Esigma[MAXZ][MAXA][3][2]={NAN}; // 0: C, 1: CH2, 2: H, 0: -1n, 1: -1p
TGraphErrors *g_incl[MAXZ][3][2]; // 0: C, 1: CH2, 2: H, 0: -1n, 1: -1p
//
// Def histograms
TH1D* h_music[NUMTARGET][NUMPADDLE+1], *h_aoq[NUMTARGET][NUMPADDLE+1][MAXA];
TH1I* h_aoq_gated[NUMTARGET][NUMPADDLE+1][MAXZ][MAXA][MAXZ];
TH2I* h_pid[NUMTARGET][NUMPADDLE+1];
Int_t num_pid[NUMTARGET][NUMPADDLE+1][MAXZ][MAXA]={0};
//TH2S* h_gated_pid[NUMTARGET][NUMPADDLE][MAXZ][MAXA];
TF1* f_music[NUMTARGET][NUMPADDLE+1], *f_aoq[NUMTARGET][NUMPADDLE+1][MAXZ][MAXA][MAXZ];
//
// Def functions
void calc_incl(int index = -1);
int loadfile();
void fill_histos(int i), fit_music(int i), fill_aoq_gated(int i), fit_aoq(int i);
void fill_histos_38Ca(int i), fill_aoq_gated_38Ca(int i);
Int_t fit_music_paddle(int i, int p), fit_aoq_paddle(int i, int p, int z), fit_aoq_gated_paddle(int i, int p, int zfrs, int afrs, int zfrag);
Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag=false, Int_t i=0, Int_t targ=-1);
TF1* generate_func(TString name, int numpeaks, double minpos, double interval, double par_limit_range, double height_limit=1e10);
void delete_histos(int i=0), calcincl();
//
ProcInfo_t pinfo;
void get_mem_usage();
/////////////////////////////////////////////
void calc_incl(int index){
//  ROOT::EnableImplicitMT(10);
//  ROOT::EnableThreadSafety();
  get_mem_usage();
  if(index <0){
    outfilename.ReplaceAll("TARG","all");
    outcsvname.ReplaceAll("TARG","all");
  }else if(index<3){
    outfilename.ReplaceAll("TARG",targetname[index]);
    outcsvname.ReplaceAll("TARG",targetname[index]);
  }else{
    return;
  }
  outfilename.ReplaceAll("FRS",frsname);
  outcsvname.ReplaceAll("FRS",frsname);
  //
  if(loadfile()!=0) return;
  //
  for(int i=0; i<NUMTARGET; i++){
    if(index>=0 && index != i) continue;
    fill_histos(i);
    if(!frs13){
      fill_histos_38Ca(i);
    }
    fit_music(i);
    fit_aoq(i);
    fill_aoq_gated(i);
    if(!frs13){
      fill_aoq_gated_38Ca(i);
    }
    //
    ROOT::DisableImplicitMT();
    for(int p =0; p<NUMPADDLE+1; p++){
      cout<<"Targ: "<< targetname[i] <<", Fit for Paddle"<<p+1;
      for(int ZA=MININCLZ*MAXA; ZA<MAXINCLZ*MAXA; ZA++){//FRS: Z: ZA/MAXA, A: ZA%MAXA
	for(int k =MINZ; k<=ZA/MAXA; k++){//Frag
	  fit_aoq_gated_paddle(i,p,ZA/MAXA,ZA%MAXA,k);
	  //if(0!=fit_aoq_gated_paddle(i,p,ZA/MAXA,ZA%MAXA,k))
	    //cout<<"Not fitted: "<<i<<", "<<p<<", "<<ZA/MAXA<<", "<<ZA%MAXA<<", "<<k<<endl;
	}
      }
    }
  }
  if(index>=0) return;
  cout<<"csv out"<<endl;
  calcincl();
}
/////////////////////////////////////////////
void fill_histos(int i){
  Float_t FRSAoQ, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, ROLU_X, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
  UChar_t Tofw_Paddle;
  //
  tree[i]->SetBranchAddress("FRS_AoQ",&FRSAoQ);
  tree[i]->SetBranchAddress("FRS_Brho",&FRSBrho);
  tree[i]->SetBranchAddress("FRS_Z", &MusicZ);
  tree[i]->SetBranchAddress("Frag_AoQ", &FragAoQ_corr);
  tree[i]->SetBranchAddress("Frag_AoQ_original", &FragAoQ);
  tree[i]->SetBranchAddress("Frag_Brho",&FragBrho);
  tree[i]->SetBranchAddress("Frag_Z", &FragZ);
  tree[i]->SetBranchAddress("Frag_Tof", &FragTof);
  tree[i]->SetBranchAddress("MusicE", &MusicE);
  tree[i]->SetBranchAddress("TwimE", &TwimE);
  tree[i]->SetBranchAddress("TwimTheta", &TwimTheta);
  tree[i]->SetBranchAddress("ROLU_X", &ROLU_X);
  tree[i]->SetBranchAddress("Mw1_X", &Mw1_X);
  tree[i]->SetBranchAddress("Mw2_X", &Mw2_X);
  tree[i]->SetBranchAddress("Mw3_X", &Mw3_X);
  tree[i]->SetBranchAddress("Mw1_Y", &Mw1_Y);
  tree[i]->SetBranchAddress("Mw2_Y", &Mw2_Y);
  tree[i]->SetBranchAddress("Mw3_Y", &Mw3_Y);
  tree[i]->SetBranchAddress("Tofw_Y", &Tofw_Y);
  tree[i]->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  //
  for(int p =0; p<NUMPADDLE+1; p++){
    cout<<"Preparing histos for paddle"<<p+1<<" ";
    TString dummy = Form("h_gated_pid_%i_",p+1) + targetname[i];
    h_pid[i][p] = new TH2I(dummy,dummy,BINAOQ, MINAOQ, MAXAOQ, BINZ, MINZ-0.5, MAXZ+0.5);
    /*
    for(int j=0; j<MAXZ*MAXA; j++){
      dummy = Form("h_gated_pid_paddle%i_Z%i_A%i_Z%s", p+1, j/MAXA, j%MAXA, targetname[i].Data());
      h_gated_pid[i][p][j/MAXA][j%MAXA] = new TH2S(dummy, dummy, BINAOQ, MINAOQ, MAXAOQ, BINZ, MINZ-0.5, MAXZ+0.5);
      }*/
    get_mem_usage();
  }
  //
  ROOT::DisableImplicitMT();
  int neve=0;
  for(Long64_t n=0; n<nentry[i]; n++){
    if(++neve%10000==0)
      cout<<"\r"<<neve<<" enrties done in "<<nentry[i]<<flush;
    //
    tree[i]->GetEntry(n);
    //
    if(ROLU_X < 0) continue;
    int p = Tofw_Paddle-1;
    if(p >= NUMPADDLE || p < 0) continue;
    h_pid[i][NUMPADDLE]->Fill(FragAoQ_corr,FragZ);
    h_pid[i][p]->Fill(FragAoQ_corr,FragZ);
    //
    /* // h_gated_pid is not used
    double tmpzet = MusicZ, tmpaoq = FRSAoQ;
    double sigma=pid(tmpzet, tmpaoq);
    if(isnan(sigma) || sigma > significance) continue; // Get FRS PID
    int tmpA = tmpaoq*tmpzet + 0.1;
    int tmpZ = tmpzet + 0.1;//0.1 is to make sure the conversion from double to int to be correct
    if(tmpZ<MAXZ && tmpA<MAXA)
      h_gated_pid[i][p][tmpZ][tmpA]->Fill(FragAoQ_corr,FragZ);
    */
  }
  cout<<endl;
}
//
void fill_histos_38Ca(int i){
  Float_t FRSAoQ, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, ROLU_X, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
  UChar_t Tofw_Paddle;
  //
  tree2[i]->SetBranchAddress("FRS_AoQ",&FRSAoQ);
  tree2[i]->SetBranchAddress("FRS_Brho",&FRSBrho);
  tree2[i]->SetBranchAddress("FRS_Z", &MusicZ);
  tree2[i]->SetBranchAddress("Frag_AoQ", &FragAoQ_corr);
  tree2[i]->SetBranchAddress("Frag_AoQ_original", &FragAoQ);
  tree2[i]->SetBranchAddress("Frag_Brho",&FragBrho);
  tree2[i]->SetBranchAddress("Frag_Z", &FragZ);
  tree2[i]->SetBranchAddress("Frag_Tof", &FragTof);
  tree2[i]->SetBranchAddress("MusicE", &MusicE);
  tree2[i]->SetBranchAddress("TwimE", &TwimE);
  tree2[i]->SetBranchAddress("TwimTheta", &TwimTheta);
  tree2[i]->SetBranchAddress("ROLU_X", &ROLU_X);
  tree2[i]->SetBranchAddress("Mw1_X", &Mw1_X);
  tree2[i]->SetBranchAddress("Mw2_X", &Mw2_X);
  tree2[i]->SetBranchAddress("Mw3_X", &Mw3_X);
  tree2[i]->SetBranchAddress("Mw1_Y", &Mw1_Y);
  tree2[i]->SetBranchAddress("Mw2_Y", &Mw2_Y);
  tree2[i]->SetBranchAddress("Mw3_Y", &Mw3_Y);
  tree2[i]->SetBranchAddress("Tofw_Y", &Tofw_Y);
  tree2[i]->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  //
  int neve=0;
  for(Long64_t n=0; n<nentry2[i]; n++){
    if(++neve%10000==0)
      cout<<"\r"<<neve<<" enrties done in "<<nentry2[i]<<flush;
    //
    tree2[i]->GetEntry(n);
    //
    if(ROLU_X < 0) continue;
    int p = Tofw_Paddle-1;
    if(p >= NUMPADDLE || p < 0) continue;
    h_pid[i][NUMPADDLE]->Fill(FragAoQ_corr,FragZ+zoff);// offset only for FRS122 settings
    h_pid[i][p]->Fill(FragAoQ_corr,FragZ+zoff);
  }
  cout<<endl;
}
//
void delete_histos(int i){
  cout<<"delete_histos("<<i<<")"<<endl;
  get_mem_usage();
  for(int p =0; p<NUMPADDLE; p++){
    //h_pid[i][p]->Clear();
    delete h_pid[i][p];
    for(int j=0; j<MAXZ*MAXA; j++){
      //h_gated_pid[i][p][j/MAXA][j%MAXA]->Clear();
      //delete h_gated_pid[i][p][j/MAXA][j%MAXA];
    }
    cout<<"Paddle"<<p+1<<endl;
  }
  cout<<"Deleted filled histos."<<endl;
  get_mem_usage();
}
//
void fit_music(int i){
  //
  fout->cd();
  for(int p =0; p<NUMPADDLE+1; p++)
    {
      cout<<"Paddle"<< p+1 << ", max bin: " << h_pid[i][p]->GetMaximum()<<endl;
      if(h_pid[i][p]->GetMaximum()>30000){//limit for TH2S
	cout<<"Too many events for TH2S"<<endl;
      }
      h_pid[i][p]->Write();
      /*
      for(int j=0; j<MAXZ*MAXA; j++){
	h_gated_pid[i][p][j/MAXA][j%MAXA]->Write();
      }*/
    }
  //TTimer *timer = new TTimer();
  //timer->Start();
  //for(int p =0; p<NUMPADDLE; p++)
  //fit_music_paddle(i, p);
  // Multi threading
  //timer->Stop();
  //cout<<"Time: "<<timer->GetTime().AsString()<<endl;
  //cout<<"Going into TTaskGroup"<<endl;
  //ROOT::EnableImplicitMT(4);
  ROOT::EnableImplicitMT(10);
  //ROOT::EnableThreadSafety();
  //ROOT::Experimental::TTaskGroup tg;
  //int i_paddles[NUMPADDLE];
  for(int p =0; p<NUMPADDLE+1; p++){
    //i_paddles[]
    //tg.Run([&, p, i] { cout<<"Paddle:"<<p+1<<" is set"<<endl; sleep(2); fit_music_paddle(i, p); } );
    //tg.Run([&, p] { cout<<"test"<<p<<endl; } );
    //cout<<"Paddle:"<<p+1<<" is set"<<endl;
    fit_music_paddle(i, p);
  }
  //ROOT::TThreadExecutor pool; auto return_music = pool.Foreach([&] (int i_paddle) { fit_music_paddle(i, i_paddle); }, ROOT::TSeqI(NUMPADDLE));
  //ROOT::TThreadExecutor pool; pool.Foreach([&, i] (int i_paddle) { fit_music_paddle(i, i_paddle); }, ROOT::TSeqI(NUMPADDLE));
  //tg.Wait(); cout<<"TTask end"<<endl;
  //timer->Stop();
  //cout<<"Time: "<<timer->GetTime().AsString()<<endl;
  //
  for(int p =0; p<NUMPADDLE+1; p++)
    h_music[i][p]->Write();
  ROOT::DisableImplicitMT();
}
Int_t fit_music_paddle(int i, int p){
  cout<<"fit_music_paddle for paddle"<<p+1<<endl;
  h_music[i][p] = h_pid[i][p]->ProjectionY();
  int numpeaks = MAXZ-MINZ;
  cout<<"Entries for paddle"<<p+1<<": "<<h_music[i][p]->GetEntries()<<endl;
  if(h_music[i][p]->GetEntries()<1){
    cout<<"Skippling paddle"<<p+1<<endl;
    return -1;
  }
  TString f_name = "f_music_" + targetname[i] + "_" + TString::Itoa(p+1,10);
  cout<<f_name<<endl;
  f_music[i][p] = generate_func(f_name, numpeaks, (double)MINZ, 1, 0.2, 0); // all peaks are fixed to 0 first
  f_music[i][p] ->FixParameter(0,0.);
  for(int j=0; j<numpeaks; j++){
    f_music[i][p] ->SetParLimits(3*j+1,0., 0.2 * h_music[i][p]->GetMaximum());
  }
  //
  TSpectrum *sp = new TSpectrum(numpeaks);
  auto S_numpeaks = sp->Search(h_music[i][p]);
  auto posX = sp->GetPositionX();
  auto posY = sp->GetPositionY();
  for(int j=S_numpeaks-1; 0<=j; j--){
    int bin_z = (int)(posX[j]+0.5) - MINZ;
    f_music[i][p] ->FixParameter(3*bin_z +1, posY[j]);
    f_music[i][p] ->FixParameter(3*bin_z +2, posX[j]);
  }
  /*
    Int_t Mx, My, Mz;
    h_music[i][p]->GetMaximumBin(Mx, My,Mz);
    Double_t MCentre = h_music[i][p]->GetXaxis()->GetBinCenter(Mx);
    Int_t IMax = (Int_t)MCentre - MINZ;
    if(IMax<numpeaks){
    f_music[i][p] ->FixParameter(3*IMax+1, h_music[i][p]->GetMaximum());
    f_music[i][p] ->FixParameter(3*IMax+2, MCentre);
    }*/
  //
  for(int j=S_numpeaks-1; 0<=j; j--){
    int bin_z = (int)(posX[j]+0.5) - MINZ;
    if(numpeaks <= bin_z || bin_z < 0) continue;
    f_music[i][p] ->SetParLimits(3*bin_z +1, 0.5 * posY[j], 1.5* posY[j]);
    h_music[i][p] ->Fit(f_music[i][p],"Q 0","",(double)(bin_z+MINZ)-0.3,(double)(bin_z+MINZ)+0.3);
    f_music[i][p] ->SetParLimits(3*bin_z +2, posX[j]-0.1, posX[j]+0.1);
    h_music[i][p] ->Fit(f_music[i][p],"Q 0","goff",(double)(bin_z+MINZ)-0.3,(double)(bin_z+MINZ)+0.3);
    //
    f_music[i][p] ->FixParameter(3*bin_z +1,f_music[i][p] ->GetParameter(3*bin_z +1));
    f_music[i][p] ->FixParameter(3*bin_z +2,f_music[i][p] ->GetParameter(3*bin_z +2));
    f_music[i][p] ->FixParameter(3*bin_z +3,f_music[i][p] ->GetParameter(3*bin_z +3));
    cout<<"Frag_Z Peak"<<j<<": at "<<f_music[i][p] ->GetParameter(3*bin_z +2)<<endl;

  }
  for(int j=0; j<numpeaks; j++){
    if(f_music[i][p] ->GetParameter(3*j+1)>0) continue;
    f_music[i][p] ->SetParLimits(3*j+1,0., 0.2 * h_music[i][p]->GetMaximum());
    h_music[i][p] ->Fit(f_music[i][p],"Q 0","",MINZ+j-0.3,MINZ+j+0.3);
    f_music[i][p] ->FixParameter(3*j+1,f_music[i][p] ->GetParameter(3*j+1));
  }
  //
  /*
    for(int j=0; j<numpeaks; j++){
    //cout<<j<<"-th peak: "<<params[3*j+1]<<", ";
    f_music[i][p] ->SetParameter(3*j+1,posY[j]);
    f_music[i][p] ->SetParameter(3*j+2,posX[j]);
    f_music[i][p] ->SetParLimits(3*j+2,posX[j]-0.3,posX[j]+0.3);
    f_music[i][p] ->SetParameter(3*j+3,0.2);
    }*/
  //h_music[i][p] ->Fit(f_music[i][p],"LL","",MINZ-0.5, MAXZ+0.5);
  //
  //f_music[i][p]->SetParLimits(3*IMax+2, MCentre-0.3, MCentre+0.3);
  //f_music[i][p]->SetParLimits(0,0.,0.01*h_music[i][p] ->GetMaximum());
  h_music[i][p] ->Fit(f_music[i][p],"LL Q 0","",MINZ-0.3, MAXZ+0.3); // Fit rest of peaks
  //
  auto params = f_music[i][p] ->GetParameters();
  for(int j=0; j<numpeaks; j++){
    //cout<<j<<"-th peak: "<<params[3*j+1]<<", ";
    if(params[3*j+1]<1.){
      f_music[i][p] ->SetParameter(3*j+1,0.);
      f_music[i][p] ->SetParameter(3*j+2,0.);
      f_music[i][p] ->SetParameter(3*j+3,0.1);
      f_music[i][p] ->FixParameter(3*j+1,0.);
      f_music[i][p] ->FixParameter(3*j+2,0.);
      f_music[i][p] ->FixParameter(3*j+3,0.1);
    }else{
      f_music[i][p] ->SetParLimits(3*j+1,0.75*params[3*j+1],1.5*params[3*j+1]);
      f_music[i][p] ->SetParameter(3*j+1,params[3*j+1]);
      f_music[i][p] ->SetParLimits(3*j+2,params[3*j+2]-0.05,params[3*j+2]+0.05);
      f_music[i][p] ->SetParameter(3*j+2,params[3*j+2]);
      f_music[i][p] ->SetParLimits(3*j+3,0.75*params[3*j+3],1.5*params[3*j+3]);
      f_music[i][p] ->SetParameter(3*j+3,params[3*j+3]);
    }
  }
  h_music[i][p] ->Fit(f_music[i][p],"LL Q 0","",MINZ-0.3, MAXZ+0.3);
  params = f_music[i][p] ->GetParameters();
  //
  cout<<"Fit result for paddle"<<p+1<<". Chi-sq:"<<f_music[i][p]->GetChisquare()<<", NDF:"<<f_music[i][p]->GetNDF()<<endl<<endl;
  //
  for(int j=0; j<numpeaks; j++){
    if(params[3*j+1]>1.){
      for(int i_aq=0; i_aq<MAXA; i_aq++){
	pid_fraggate[i][p][j+MINZ][i_aq][0] = params[3*j+2];
	pid_fraggate[i][p][j+MINZ][i_aq][1] = params[3*j+3];
      }
    }
  }
  return 0;
}

void fit_aoq(int i){
  ROOT::EnableImplicitMT(10);
  for(int p =0; p<NUMPADDLE+1; p++){
    for(int z=MINZ; z<MAXZ; z++){
      if(0 != fit_aoq_paddle(i, p, z))
	continue;
      if(h_aoq[i][p][z]->GetEntries()==0)
	continue;
      f_music[i][p] ->Draw("same");
      f_music[i][p] ->Write();
      h_aoq[i][p][z]->Write();
    }
  }

  ROOT::DisableImplicitMT();
}
Int_t fit_aoq_paddle(int i, int p, int z){
  cout<<"fit_aoq_paddle for paddle"<<p+1<<" and Z"<<z<<endl;
  int k = z; // just because I copied code from different part..
  if(isnan(pid_fraggate[i][p][z][0][0]) || isnan(pid_fraggate[i][p][z][0][1]) || pid_fraggate[i][p][z][0][0]<MINZ)
    return 1;
  //
  TString dummy = Form("%s_z%i",h_pid[i][p]->GetTitle(), z);
  int min_bin = h_pid[i][p]->GetYaxis()->FindBin(pid_fraggate[i][p][z][0][0] - significance * pid_fraggate[i][p][z][0][1]);
  int max_bin = h_pid[i][p]->GetYaxis()->FindBin(pid_fraggate[i][p][z][0][0] + significance * pid_fraggate[i][p][z][0][1]);
  //int min_bin = h_pid[i][p]->FindBin(pid_fraggate[i][p][z][0][0] - significance * pid_fraggate[i][p][z][0][1]);
  //int max_bin = h_pid[i][p]->FindBin(pid_fraggate[i][p][z][0][0] + significance * pid_fraggate[i][p][z][0][1]);
  cout<<"zgate"<<pid_fraggate[i][p][z][0][0] - significance * pid_fraggate[i][p][z][0][1]<<" to "<<pid_fraggate[i][p][z][0][0] + significance * pid_fraggate[i][p][z][0][1]<<endl;
  cout<<"min_bin"<<min_bin<<", max_bin"<<max_bin<<endl;
  auto h_tmp = h_pid[i][p]->ProjectionX(dummy, min_bin, max_bin);
  h_tmp ->SetTitle(dummy);
  //
  int rangeA[2] = {(int)(MINAOQ*(double)z+0.3), (int)min(MAXAOQ*(double)z-0.3,(double)(MAXA))};
  int numpeaks = rangeA[1]-rangeA[0]+1;
  if(numpeaks<1)
    return 1;
  dummy = Form("f_aoq_%s_paddle%i_z%i",targetname[i].Data(),p,z);
  auto f_tmp = generate_func(dummy, numpeaks, (double)rangeA[0]/(double)z, 1./(double)z, 0.1/(double)z, 0);
  f_tmp->FixParameter(0,0.);
  if(h_tmp->GetEntries()<1) return 1;

  cout<<"Paddle"<<p+1<<", Z"<<z <<endl;
  int num_peaks_guess=4;
  TSpectrum *sp = new TSpectrum(num_peaks_guess);
  auto S_numpeaks = sp->Search(h_tmp);
  auto posX = sp->GetPositionX();
  auto posY = sp->GetPositionY();
  for(int j=S_numpeaks-1; 0<=j; j--){
    int bin_z = (int)(posX[j]+0.5) - MINZ;
    f_tmp ->FixParameter(3*bin_z +1, posY[j]);
    f_tmp ->FixParameter(3*bin_z +2, posX[j]);
  }
  for(int j=0; j<S_numpeaks; j++){
    //for(int j=S_numpeaks-1; 0<=j; j--){
    int bin_aoq = (int)(posX[j] * (double)k + 0.5) - rangeA[0];
    if(numpeaks <= bin_aoq ) continue;
    if(f_tmp->GetParameter(3*bin_aoq +1) > 0 ) continue;
    f_tmp ->SetParLimits(3*bin_aoq +1, 0.5 * posY[j], 1.5* posY[j]);
    h_tmp ->Fit(f_tmp,"Q","",((double)(bin_aoq+ rangeA[0])-0.3)/(double)k,((double)(bin_aoq+ rangeA[0])+0.3)/(double)k);
  }
  // divide fits
  for(int j=0; j<S_numpeaks; j++){
    int bin_aoq = (int)(posX[j] * (double)k + 0.5) - rangeA[0];
    f_tmp ->SetParLimits(3*bin_aoq +2, posX[j]-0.05, posX[j]+0.05);
    h_tmp ->Fit(f_tmp,"Q","",((double)(bin_aoq+ rangeA[0])-0.3)/(double)k,((double)(bin_aoq+ rangeA[0])+0.3)/(double)k);
    //
    //f_tmp ->FixParameter(3*bin_aoq +1,f_tmp ->GetParameter(3*bin_aoq +1));
    f_tmp ->FixParameter(3*bin_aoq +2,f_tmp ->GetParameter(3*bin_aoq +2));
    f_tmp ->FixParameter(3*bin_aoq +3,f_tmp ->GetParameter(3*bin_aoq +3));
    cout<<"Frag Mass Peak"<<j<<" for z="<<k<<": at "<<f_tmp ->GetParameter(3*bin_aoq +2)<<endl;
  }
  for(int j=numpeaks-1; j>=0; j--){
    if(f_tmp->GetParameter(3*j+1)>0) continue;
    f_tmp->SetParameter(3*j+1,0.2* h_tmp->GetMaximum());
    f_tmp->SetParLimits(3*j+1,0., 2.* h_tmp->GetMaximum());
    h_tmp->Fit(f_tmp,"QL","",((double)(rangeA[0]+j)-0.3)/(double)k,((double)(rangeA[0]+j)+0.3)/(double)k);
    f_tmp->FixParameter(3*j+1,f_tmp->GetParameter(3*j+1));
  }
  auto params = f_tmp->GetParameters();
  for(int j=0; j<numpeaks; j++){
    if(params[3*j+1]<0.1){
      f_tmp ->SetParameter(3*j+1,0.);
      f_tmp ->SetParameter(3*j+2,0.);
      f_tmp ->SetParameter(3*j+3,0.1);
      f_tmp ->FixParameter(3*j+1,0.);
      f_tmp ->FixParameter(3*j+2,0.);
      f_tmp ->FixParameter(3*j+3,0.1);
    }else{
      f_tmp ->SetParLimits(3*j+1,0.*params[3*j+1],1.3*params[3*j+1]);
      f_tmp ->SetParameter(3*j+1,params[3*j+1]);
      f_tmp ->SetParLimits(3*j+2,params[3*j+2] - 0.1/(double)k, params[3*j+2] + 0.1/(double)k);
      f_tmp ->SetParameter(3*j+2,params[3*j+2]);
      //f_tmp ->SetParLimits(3*j+3,0.75*params[3*j+3],3.5*params[3*j+3]);
      f_tmp ->SetParLimits(3*j+3,0.75*params[3*j+3], significance*params[3*j+3]);
      f_tmp ->SetParameter(3*j+3,params[3*j+3]);
    }
  }
  h_tmp->Fit(f_tmp,"Q LL","",((double)rangeA[0]-0.3)/(double)k, ((double)rangeA[1]+0.7)/(double)k);
  //
  cout<<"AoQ fit result for target "<<targetname[i]<<", paddle"<<p+1<<", Z_frag"<<k<<". Chi-sq:"<<f_tmp->GetChisquare()<<", NDF:"<<f_tmp->GetNDF()<<endl<<endl;
  //
  for(int j=0; j<numpeaks; j++){
    if(params[3*j+1]<0.1)continue;
    int i_A = rangeA[0] + j;
    pid_fraggate[i][p][k][i_A][2] = params[3*j+2];
    pid_fraggate[i][p][k][i_A][3] = params[3*j+3];
    for(int i_param=0; i_param<4; i_param++) cout<<pid_fraggate[i][p][k][i_A][i_param]<<" ";
    cout<<"A:"<<pid_fraggate[i][p][k][i_A][0]*pid_fraggate[i][p][k][i_A][2]<<endl;
  }
  //
  f_aoq[i][p][z][0][0] = f_tmp;
  h_aoq[i][p][z] = h_tmp;
  f_tmp->Draw("same");
  f_tmp->Write();
  h_aoq[i][p][z] ->Write();
  return 0;
}

//
void fill_aoq_gated(int i){
  Float_t FRSAoQ, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, ROLU_X, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
  UChar_t Tofw_Paddle;
  //
  tree[i]->SetBranchAddress("FRS_AoQ",&FRSAoQ);
  tree[i]->SetBranchAddress("FRS_Brho",&FRSBrho);
  tree[i]->SetBranchAddress("FRS_Z", &MusicZ);
  tree[i]->SetBranchAddress("Frag_AoQ", &FragAoQ_corr);
  tree[i]->SetBranchAddress("Frag_AoQ_original", &FragAoQ);
  tree[i]->SetBranchAddress("Frag_Brho",&FragBrho);
  tree[i]->SetBranchAddress("Frag_Z", &FragZ);
  tree[i]->SetBranchAddress("Frag_Tof", &FragTof);
  tree[i]->SetBranchAddress("MusicE", &MusicE);
  tree[i]->SetBranchAddress("TwimE", &TwimE);
  tree[i]->SetBranchAddress("TwimTheta", &TwimTheta);
  tree[i]->SetBranchAddress("ROLU_X", &ROLU_X);
  tree[i]->SetBranchAddress("Mw1_X", &Mw1_X);
  tree[i]->SetBranchAddress("Mw2_X", &Mw2_X);
  tree[i]->SetBranchAddress("Mw3_X", &Mw3_X);
  tree[i]->SetBranchAddress("Mw1_Y", &Mw1_Y);
  tree[i]->SetBranchAddress("Mw2_Y", &Mw2_Y);
  tree[i]->SetBranchAddress("Mw3_Y", &Mw3_Y);
  tree[i]->SetBranchAddress("Tofw_Y", &Tofw_Y);
  tree[i]->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  //
  for(int p =0; p<NUMPADDLE+1; p++){
    cout<<"Preparing histos for aoq fit"<<p+1<<" ";
    for(int j=MININCLZ*MAXA; j<MAXINCLZ*MAXA; j++){
      for(int k =MINZ; k<MAXZ; k++){
	if(isnan(pid_fraggate[i][p][k][0][0])) continue; // Check the Twim gates.
	TString dummy = Form("h_aoq_gated_paddle%i_Z%i_A%i_Z%i_%s", p+1, j/MAXA, j%MAXA, k, targetname[i].Data());
	if(p==NUMPADDLE){
	  dummy = Form("h_aoq_gated_all_Z%i_A%i_Z%i_%s", j/MAXA, j%MAXA, k, targetname[i].Data());
	}
	h_aoq_gated[i][p][j/MAXA][j%MAXA][k] = new TH1I(dummy,dummy,BINAOQ, MINAOQ, MAXAOQ);
      }
    }
    get_mem_usage();
  }
  //
  int neve=0;
  for(Long64_t n=0; n<nentry[i]; n++){
    if(++neve%10000==0)
      cout<<"\r"<<neve<<" enrties done in "<<nentry[i]<<flush;
    //
    tree[i]->GetEntry(n);
    //incoming cut
    if(ROLU_X < 0) continue;
    double tmpzet = MusicZ, tmpaoq = FRSAoQ;
    double sigma=pid(tmpzet, tmpaoq);
    if(isnan(sigma) || sigma > significance) continue; // Get FRS PID
    int tmpA = tmpaoq*tmpzet + 0.1;
    int tmpZ = tmpzet + 0.1;//0.1 is to make sure the conversion from double to int to be correct
    //
    num_pid[i][NUMPADDLE][tmpZ][tmpA]++;
    //
    int p = Tofw_Paddle-1;
    if(p >= NUMPADDLE || p < 0) continue;
    num_pid[i][p][tmpZ][tmpA]++;
    //
    if(isnan(MusicZ*FRSAoQ*FragZ*FragAoQ_corr)) continue;
    int tmpZ_frag = FragZ + zoff + 0.5, tmpA_frag = ((FragZ+zoff)*FragAoQ_corr) + 0.5;
    if(tmpZ<MININCLZ || tmpZ>=MAXINCLZ || tmpA>=MAXA || tmpZ_frag>=MAXZ || tmpA_frag>=MAXA) continue;
    if(isnan(pid_fraggate[i][p][tmpZ_frag][0][0])) continue;
    if(TMath::Abs(FragZ+zoff-pid_fraggate[i][p][tmpZ_frag][0][0]) > significance * pid_fraggate[i][p][tmpZ_frag][0][1]) continue; // For Z, 3 sigma cut applied.
    h_aoq_gated[i][NUMPADDLE][tmpZ][tmpA][tmpZ_frag]->Fill(FragAoQ_corr);
    h_aoq_gated[i][p][tmpZ][tmpA][tmpZ_frag]->Fill(FragAoQ_corr);
  }
}
//
void fill_aoq_gated_38Ca(int i){
  Float_t FRSAoQ, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, ROLU_X, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
  UChar_t Tofw_Paddle;
  //
  tree2[i]->SetBranchAddress("FRS_AoQ",&FRSAoQ);
  tree2[i]->SetBranchAddress("FRS_Brho",&FRSBrho);
  tree2[i]->SetBranchAddress("FRS_Z", &MusicZ);
  tree2[i]->SetBranchAddress("Frag_AoQ", &FragAoQ_corr);
  tree2[i]->SetBranchAddress("Frag_AoQ_original", &FragAoQ);
  tree2[i]->SetBranchAddress("Frag_Brho",&FragBrho);
  tree2[i]->SetBranchAddress("Frag_Z", &FragZ);
  tree2[i]->SetBranchAddress("Frag_Tof", &FragTof);
  tree2[i]->SetBranchAddress("MusicE", &MusicE);
  tree2[i]->SetBranchAddress("TwimE", &TwimE);
  tree2[i]->SetBranchAddress("TwimTheta", &TwimTheta);
  tree2[i]->SetBranchAddress("ROLU_X", &ROLU_X);
  tree2[i]->SetBranchAddress("Mw1_X", &Mw1_X);
  tree2[i]->SetBranchAddress("Mw2_X", &Mw2_X);
  tree2[i]->SetBranchAddress("Mw3_X", &Mw3_X);
  tree2[i]->SetBranchAddress("Mw1_Y", &Mw1_Y);
  tree2[i]->SetBranchAddress("Mw2_Y", &Mw2_Y);
  tree2[i]->SetBranchAddress("Mw3_Y", &Mw3_Y);
  tree2[i]->SetBranchAddress("Tofw_Y", &Tofw_Y);
  tree2[i]->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  //
  int neve=0;
  for(Long64_t n=0; n<nentry2[i]; n++){
    if(++neve%10000==0)
      cout<<"\r"<<neve<<" enrties done in "<<nentry2[i]<<flush;
    //
    tree2[i]->GetEntry(n);
    //incoming cut
    if(ROLU_X < 0) continue;
    double tmpzet = MusicZ + zoff, tmpaoq = FRSAoQ;
    double sigma=pid(tmpzet, tmpaoq);
    if(isnan(sigma) || sigma > significance) continue; // Get FRS PID
    int tmpA = tmpaoq*tmpzet + 0.1;
    int tmpZ = tmpzet + 0.1;//0.1 is to make sure the conversion from double to int to be correct
    //
    num_pid[i][NUMPADDLE][tmpZ][tmpA]++;
    //
    int p = Tofw_Paddle-1;
    if(p >= NUMPADDLE || p < 0) continue;
    num_pid[i][p][tmpZ][tmpA]++;
    //
    if(isnan(MusicZ*FRSAoQ*FragZ*FragAoQ_corr)) continue;
    int tmpZ_frag = FragZ + zoff + 0.5, tmpA_frag = ((FragZ+zoff)*FragAoQ_corr) + 0.5;
    if(tmpZ<MININCLZ || tmpZ>=MAXINCLZ || tmpA>=MAXA || tmpZ_frag>=MAXZ || tmpA_frag>=MAXA) continue;
    if(isnan(pid_fraggate[i][p][tmpZ_frag][0][0])) continue;
    if(TMath::Abs(FragZ+zoff-pid_fraggate[i][p][tmpZ_frag][0][0]) > significance * pid_fraggate[i][p][tmpZ_frag][0][1]) continue; // For Z, 3 sigma cut applied.
    h_aoq_gated[i][NUMPADDLE][tmpZ][tmpA][tmpZ_frag]->Fill(FragAoQ_corr);
    h_aoq_gated[i][p][tmpZ][tmpA][tmpZ_frag]->Fill(FragAoQ_corr);
  }
}
//
int fit_aoq_gated_paddle(int i, int p, int zfrs, int afrs, int zfrag){
  //
  //ROOT::EnableImplicitMT();
  fout->cd();
  //cout<<"Fitting AoQ for Paddle"<<p+1<<endl;
  int ZA = zfrs*afrs;
  int k = zfrag;
  int rangeA[2] = {(int)(MINAOQ*(double)k+0.3), (int)min(MAXAOQ*(double)k-0.3, (double)(afrs - zfrs + k))}; // Max range is determined for -1p reaction
  int numpeaks = rangeA[1]-rangeA[0]+1;
  if(numpeaks<1) return 1;
  auto h_tmp = h_aoq_gated[i][p][zfrs][afrs][k];
  if(h_tmp->GetEntries()<1) return 1;
  //
  //cout<<h_tmp->GetEntries()<<endl;
  cout<<endl;
  TString dummy = Form("f_aoq_%s_%i_%i_%i_%i",targetname[i].Data(),p, zfrs, afrs,k);
  if(p==NUMPADDLE) dummy =  Form("f_aoq_all_%s_%i_%i_%i",targetname[i].Data(), zfrs, afrs,k);
  auto f_tmp = generate_func(dummy, numpeaks,
			     (double)rangeA[0]/(double)k, 1./(double)k, 0.1/(double)k, 0);
  f_aoq[i][p][zfrs][afrs][k] = f_tmp;
  f_tmp->FixParameter(0,0.);
  cout<<"Paddle"<<p+1<<", Z"<<zfrs<<", A"<<afrs<<", Zfrag"<<k<<endl;
  //

  /*
  int num_peaks_guess=1;
  if(i>0){
    num_peaks_guess+=2;
    if(k<zfrs) num_peaks_guess+=2;
  }
  TSpectrum *sp = new TSpectrum(num_peaks_guess);
  auto S_numpeaks = sp->Search(h_tmp);
  auto posX = sp->GetPositionX();
  auto posY = sp->GetPositionY();
  for(int j=S_numpeaks-1; 0<=j; j--){
    int bin_z = (int)(posX[j]+0.5) - MINZ;
    f_tmp ->FixParameter(3*bin_z +1, posY[j]);
    f_tmp ->FixParameter(3*bin_z +2, posX[j]);
  }
  / *
    Int_t Mx, My, Mz;
    h_tmp->GetMaximumBin(Mx, My,Mz);
    Double_t MCentre = h_tmp->GetXaxis()->GetBinCenter(Mx);
    Int_t AMax = (Int_t)(MCentre*(Double_t)k) - rangeA[0];
    if(AMax<numpeaks){
    f_tmp ->FixParameter(3*AMax+1, h_tmp->GetMaximum());
    f_tmp ->SetParameter(3*AMax+2, MCentre);
    f_tmp ->SetParLimits(3*AMax+2, MCentre-0.1/(double)k, MCentre+0.1/(double)k);
    }
  * /
  //
  for(int j=0; j<S_numpeaks; j++){
    //for(int j=S_numpeaks-1; 0<=j; j--){
    int bin_aoq = (int)(posX[j] * (double)k + 0.5) - rangeA[0];
    if(numpeaks <= bin_aoq ) continue;
    if(f_tmp->GetParameter(3*bin_aoq +1) > 0 ) continue;
    f_tmp ->SetParLimits(3*bin_aoq +1, 0.5 * posY[j], 1.5* posY[j]);
    h_tmp ->Fit(f_tmp,"Q","",((double)(bin_aoq+ rangeA[0])-0.1)/(double)k,((double)(bin_aoq+ rangeA[0])+0.1)/(double)k);
    f_tmp ->SetParLimits(3*bin_aoq +2, posX[j]-0.1, posX[j]-0.1);
    h_tmp ->Fit(f_tmp,"Q","",((double)(bin_aoq+ rangeA[0])-0.1)/(double)k,((double)(bin_aoq+ rangeA[0])+0.1)/(double)k);
    //
    f_tmp ->FixParameter(3*bin_aoq +1,f_tmp ->GetParameter(3*bin_aoq +1));
    f_tmp ->FixParameter(3*bin_aoq +2,f_tmp ->GetParameter(3*bin_aoq +2));
    f_tmp ->FixParameter(3*bin_aoq +3,f_tmp ->GetParameter(3*bin_aoq +3));
    cout<<"Frag Mass Peak"<<j<<" for z="<<k<<": at "<<f_tmp ->GetParameter(3*bin_aoq +2)<<endl;
  }
  */
  for(int j=numpeaks-1; j>=0; j--){
    //if(f_tmp->GetParameter(3*j+1)>0) continue;
    int afrag = rangeA[0] + j;
    if(isnan(pid_fraggate[i][p][k][afrag][2]) || isnan(pid_fraggate[i][p][k][afrag][3]) ||
       pid_fraggate[i][p][k][afrag][0] < MINZ || pid_fraggate[i][p][k][afrag][2]*pid_fraggate[i][p][k][afrag][0] < MINAOQ*MINZ) continue;
    f_tmp->SetParLimits(3*j+1,0., 2.* h_tmp->GetMaximum());
    f_tmp->FixParameter(3*j+2,pid_fraggate[i][p][k][afrag][2]);
    f_tmp->FixParameter(3*j+3,pid_fraggate[i][p][k][afrag][3]);
  }
  h_tmp->Fit(f_tmp,"QLL","",(double)(rangeA[0]-1)/(double)k, (double)(rangeA[1]+1)/(double)k);
  //
  //f_tmp->SetParLimits(3*AMax+2, MCentre-0.3/(double)k, MCentre+0.3/(double)k);
  //f_tmp->SetParLimits(0,0.,0.01*h_tmp->GetMaximum());
  auto params = f_tmp->GetParameters();
  for(int j=0; j<numpeaks; j++){
    if(params[3*j+1]<0.1){
      f_tmp ->SetParameter(3*j+1,0.);
      f_tmp ->SetParameter(3*j+2,0.);
      f_tmp ->SetParameter(3*j+3,0.1);
      f_tmp ->FixParameter(3*j+1,0.);
      f_tmp ->FixParameter(3*j+2,0.);
      f_tmp ->FixParameter(3*j+3,0.1);
    }else{
      //f_tmp ->SetParLimits(3*j+1,0.1*params[3*j+1],1.3*params[3*j+1]);
      f_tmp ->SetParameter(3*j+1,params[3*j+1]);
      //f_tmp ->SetParLimits(3*j+2,params[3*j+2] - 0.1/(double)k, params[3*j+2] + 0.1/(double)k);
      //f_tmp ->SetParameter(3*j+2,params[3*j+2]);
      //f_tmp ->SetParLimits(3*j+3,0.75*params[3*j+3],3.5*params[3*j+3]);
      //f_tmp ->SetParLimits(3*j+3,0.75*params[3*j+3], significance*params[3*j+3]);
      //f_tmp ->SetParameter(3*j+3,params[3*j+3]);
    }
  }
  h_tmp->Fit(f_tmp,"Q LL","",((double)rangeA[0]-0.3)/(double)k, ((double)rangeA[1]+0.7)/(double)k);
  //
  cout<<"AoQ fit result for target "<<targetname[i]<<", paddle"<<p+1<<", Z"<<zfrs<<", A"<<afrs<< ", Z_frag"<<k<<". Chi-sq:"<<f_tmp->GetChisquare()<<", NDF:"<<f_tmp->GetNDF()<<endl<<endl;
  //
  for(int j=0; j<numpeaks; j++){
    if(params[3*j+1]<0.1)continue;
    int i_A = rangeA[0] + j;
    //for(int i_aq=0; i_aq<MAXA; i_aq++){
    //pid_fraggate[i][p][k][i_A][2] = params[3*j+2];
    //pid_fraggate[i][p][k][i_A][3] = params[3*j+3];
    for(int i_param=0; i_param<4; i_param++) cout<<pid_fraggate[i][p][k][i_A][i_param]<<" ";
    cout<<"A:"<<pid_fraggate[i][p][k][i_A][0]*pid_fraggate[i][p][k][i_A][2]<<endl;
    //}
  }
  f_aoq[i][p][zfrs][afrs][k]->Write();
  h_tmp->Write();
  return 0;
}
//
//
//////////
//TF1* generate_func(char*name, int numpeaks, double minpos, double interval, double par_limit_range, double height_limit){
TF1* generate_func(TString name, int numpeaks, double minpos, double interval, double par_limit_range, double height_limit){
  if(numpeaks<1||interval<0) return NULL;
  TString formula="pol0(0)";
  for(int i=0; i<numpeaks; i++){
    formula.Append(Form("+gaus(%i)",3*i+1));
  }
  TF1* f_tmp = new TF1(name, formula);
  for(int i=0; i<numpeaks; i++){
    double peakpos = minpos + (double)i*interval;
    f_tmp->SetParameter(3*i+1, 0);
    f_tmp->SetParameter(3*i+2, peakpos);
    f_tmp->SetParameter(3*i+3, 0.2*par_limit_range);
    //
    if(height_limit>0){
      f_tmp->SetParLimits(3*i+1, 0, height_limit);
    }else{
      f_tmp->FixParameter(3*i+1, 0);
      //f_tmp->FixParameter(3*i+3, 0.2*par_limit_range);
    }
    f_tmp->SetParLimits(3*i+2, peakpos-par_limit_range, peakpos+par_limit_range);
    f_tmp->SetParLimits(3*i+3, 0.01*par_limit_range, par_limit_range);
  }
  //cout<<formula<<" func generated"<<endl;
  //f_tmp->SetVectorized(true);
  return f_tmp;
}

void calcincl(){
  TString temp = dir_rootfile + outcsvname;
  cout<<"CSV: "<<temp<<endl;
  fcsv.open(temp, ofstream::out );
  ///
  fcsv<<"FRSZ, FRSA";
  cout<<"FRSZ, FRSA";
  double width_bin = (MAXAOQ - MINAOQ)/(double)BINAOQ;
  for(int i =0; i<NUMTARGET; i++){
    cout<<", "<<targetname[i]<<" counts in FRS, "<<targetname[i]<<" transmission, +/-, error";
    fcsv<<", "<<targetname[i]<<" counts in FRS, "<<targetname[i]<<" transmission, +/-, error";
    /*
    cout<<", "<<targetname[i]<<"counts in FRS, "<<targetname[i]<<" unreacted, "<<targetname[i]<<" unreacted err, ";
    cout<<targetname[i]<<" -1p, "<<targetname[i]<<" -1p err, "<<targetname[i]<<" -1n, "<<targetname[i]<<" -1n err";
    fcsv<<", "<<targetname[i]<<"counts in FRS, "<<targetname[i]<<" unreacted, "<<targetname[i]<<" unreacted err, ";
    fcsv<<targetname[i]<<" -1p, "<<targetname[i]<<" -1p err, "<<targetname[i]<<" -1n, "<<targetname[i]<<" -1n err";
    */
  }
  cout<<", sigma-1n C, +/-, error, sigma-1n CH2, +/-, error, sigma-1n H, +/-, error, sigma-1p C, +/-, error, sigma-1p CH2, +/-, error, sigma-1p H +/-, error,";
  fcsv<<", sigma-1n C, +/-, error, sigma-1n CH2, +/-, error, sigma-1n H, +/-, error, sigma-1p C, +/-, error, sigma-1p CH2, +/-, error, sigma-1p H +/-, error,";
  for(int FZ=MINZ; FZ<MAXZ; FZ++){
    ///
    // Setup TGraphErrors
    for(int i = 0; i<6; i++){// i%2: reaction, i/2: target
      g_incl[FZ][i/2][i%2] = new TGraphErrors();
      temp = "Systematic inclusive cross sections of " + reactionname[i%2] + " channel for Z=" + TString::Itoa(FZ,10) + " isotopic chain for " + targetelement[i/2] + " target";
      g_incl[FZ][i/2][i%2] -> SetName(Form("g_incl_%i_%i_%i",FZ,i/2,i%2));
      g_incl[FZ][i/2][i%2] -> SetTitle(temp);
    }
    for(int FA=0; FA<MAXA; FA++){
      if(num_pid[0][NUMPADDLE][FZ][FA]<1 || FZ<MINZ+1 || FZ>MAXINCLZ) continue;
      fcsv<<endl<<FZ<<", "<<FA<<", ";
      cout<<endl<<FZ<<", "<<FA<<", ";
      double Iunreacted[NUMTARGET], Eunreacted[NUMTARGET], I1n[NUMTARGET], E1n[NUMTARGET], I1p[NUMTARGET], E1p[NUMTARGET];
      for(int i =0; i<NUMTARGET; i++){
	//counts in frs
	fcsv<<num_pid[i][NUMPADDLE][FZ][FA]<<", ";
	cout<<num_pid[i][NUMPADDLE][FZ][FA]<<", ";
	//unreacted
	int k = FZ;
	int ZA = FZ * MAXA + FA;
	int rangeA[2] = {(int)(MINAOQ*(double)k+0.3), (int)min(MAXAOQ*(double)k-0.3,(double)(ZA%MAXA - ZA/MAXA + k))};
	int numpeaks = rangeA[1]-rangeA[0]+1;
	int i_unreacted = FA - rangeA[0];
	//
	Iunreacted[i] = 0.;
	double E_tmp = 0.;
	Eunreacted[i] = 0.;
	if((f_aoq[i][NUMPADDLE][FZ][FA][k]) && i_unreacted > 0){
	  for(int p =0; p<NUMPADDLE; p++){
	    if(!f_aoq[i][p][FZ][FA][k]) continue;
	    Iunreacted[i] += TMath::Sqrt(2*TMath::Pi()) / width_bin
	      * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_unreacted + 1)
	      * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_unreacted + 3);
	    E_tmp += 2*TMath::Pi() / (width_bin * width_bin)
	      * (TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_unreacted + 1) 
			      * f_aoq[i][p][FZ][FA][k]->GetParError(3*i_unreacted + 3),2)
		 + TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParError(3*i_unreacted + 1)
				* f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_unreacted + 3),2));
	    /*Eunreacted[i] = TMath::Sqrt(2*TMath::Pi()) / width_bin
	     * TMath::Sqrt(TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_unreacted + 1) 
	     * f_aoq[i][p][FZ][FA][k]->GetParError(3*i_unreacted + 3),2)
	     +TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParError(3*i_unreacted + 1)
	     * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_unreacted + 3),2));*/
	  }
	  /*
	  if(Iunreacted[i] == 0){
	    Iunreacted[i] = TMath::Sqrt(2*TMath::Pi()) / width_bin
	      * f_aoq[i][NUMPADDLE][FZ][FA][k]->GetParameter(3*i_unreacted + 1)
	      * f_aoq[i][NUMPADDLE][FZ][FA][k]->GetParameter(3*i_unreacted + 3);
	    E_tmp = 2*TMath::Pi() / (width_bin * width_bin)
	      * (TMath::Power(f_aoq[i][NUMPADDLE][FZ][FA][k]->GetParameter(3*i_unreacted + 1) 
			      * f_aoq[i][NUMPADDLE][FZ][FA][k]->GetParError(3*i_unreacted + 3),2)
		 + TMath::Power(f_aoq[i][NUMPADDLE][FZ][FA][k]->GetParError(3*i_unreacted + 1)
				* f_aoq[i][NUMPADDLE][FZ][FA][k]->GetParameter(3*i_unreacted + 3),2));
	  }*/
	  Eunreacted[i] = TMath::Sqrt(E_tmp);
	  //
	  if(i==0){
	    trans[i][FZ][FA] = Iunreacted[i] / (double)num_pid[i][NUMPADDLE][FZ][FA];
	    Etrans[i][FZ][FA] = trans[i][FZ][FA] *
	      TMath::Sqrt(TMath::Power(Eunreacted[i]/Iunreacted[i],2.)
			  + TMath::Power(1./(double)num_pid[i][NUMPADDLE][FZ][FA],2.));
	  }else{
	    trans[i][FZ][FA] = 0.5 * Iunreacted[i] / (double)num_pid[i][NUMPADDLE][FZ][FA] + 0.5 * trans[0][FZ][FA];
	    Etrans[i][FZ][FA] = 0.5 *
	      TMath::Sqrt(TMath::Power(trans[i][FZ][FA]*Eunreacted[i]/Iunreacted[i],2.)
			  + TMath::Power(trans[i][FZ][FA]/(double)num_pid[i][NUMPADDLE][FZ][FA],2.)
			  + Etrans[0][FZ][FA]*Etrans[0][FZ][FA]);
	  }
	}
	// transmission for each target
	cout<<trans[i][FZ][FA]<<", +/-, "<<Etrans[i][FZ][FA]<<", ";
	fcsv<<trans[i][FZ][FA]<<", +/-, "<<Etrans[i][FZ][FA]<<", ";
	// crorss section
	// -1n
	I1n[i] = 0.;
	E_tmp = 0.;
	E1n[i] = 0.;
	//
	int i_1n = FA - rangeA[0] - 1;
	if(i_1n >= 0){
	  for(int p =0; p<NUMPADDLE; p++){
	    if(!f_aoq[i][p][FZ][FA][k]) continue;
	    I1n[i] += TMath::Sqrt(2*TMath::Pi()) / width_bin
	      * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1n + 1)
	      * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1n + 3);
	    E_tmp += (2*TMath::Pi()) / (width_bin*width_bin)
	      * (TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1n + 1) 
			      * f_aoq[i][p][FZ][FA][k]->GetParError(3*i_1n + 3),2)
		 + TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParError(3*i_1n + 1)
				* f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1n + 3),2));
	  }
	  E1n[i] = TMath::Sqrt(E_tmp);
	  //
	  if(i>0){
	    double frac_targ = I1n[i]/(double)num_pid[i][NUMPADDLE][FZ][FA]/trans[i][FZ][FA];
	    double frac_empty = I1n[0]/(double)num_pid[0][NUMPADDLE][FZ][FA]/trans[0][FZ][FA];
	    if(frac_targ >0 && frac_empty > 0){
	      sigma[FZ][FA][i-1][0] = (frac_targ - frac_empty)/d_targ[i]*to_mb;
	      Esigma[FZ][FA][i-1][0]= to_mb/d_targ[i] *
		TMath::Sqrt(TMath::Power(frac_targ*E1n[i]/I1n[i],2) + frac_targ*frac_targ/(double)num_pid[i][NUMPADDLE][FZ][FA]
			    + TMath::Power(frac_targ*Etrans[i][FZ][FA]/trans[i][FZ][FA],2)
			    + TMath::Power(frac_empty*E1n[0]/I1n[0],2) + frac_empty*frac_empty/(double)num_pid[0][NUMPADDLE][FZ][FA]
			    + TMath::Power(frac_empty*Etrans[0][FZ][FA]/trans[0][FZ][FA],2));
	    }else if(frac_empty <= 0){
	      sigma[FZ][FA][i-1][0] = (frac_targ)/d_targ[i]*to_mb;
	      Esigma[FZ][FA][i-1][0]= to_mb/d_targ[i] *
		TMath::Sqrt(TMath::Power(frac_targ*E1n[i]/I1n[i],2) + frac_targ*frac_targ/(double)num_pid[i][NUMPADDLE][FZ][FA]
			    + TMath::Power(frac_targ*Etrans[i][FZ][FA]/trans[i][FZ][FA],2));
	    }
	    if(sigma[FZ][FA][i-1][0] > 0){
	      int N = g_incl[FZ][i-1][0]->GetN();
	      g_incl[FZ][i-1][0]->SetPoint(N, FA, sigma[FZ][FA][i-1][0]);
	      g_incl[FZ][i-1][0]->SetPointError(N, 0, Esigma[FZ][FA][i-1][0]);
	    }else{
	      if(TMath::Abs(sigma[FZ][FA][i-1][0])>TMath::Abs(Esigma[FZ][FA][i-1][0])) Esigma[FZ][FA][i-1][0] = TMath::Abs(sigma[FZ][FA][i-1][0]);
	      sigma[FZ][FA][i-1][0] = 0.;
	    }
	  }
	  if(i==2){//for proton
	    if(sigma[FZ][FA][1][0] >= 0 && sigma[FZ][FA][0][0] >= 0){
	      sigma[FZ][FA][2][0] = 0.5 * (sigma[FZ][FA][1][0] - sigma[FZ][FA][0][0]);
	      Esigma[FZ][FA][2][0] = 0.5 * TMath::Sqrt(Esigma[FZ][FA][0][0]*Esigma[FZ][FA][0][0] + Esigma[FZ][FA][1][0]*Esigma[FZ][FA][1][0]);
	      if(sigma[FZ][FA][2][0]<0){
		if(TMath::Abs(sigma[FZ][FA][2][0]) < TMath::Abs(Esigma[FZ][FA][2][0])) Esigma[FZ][FA][2][0] = TMath::Abs(sigma[FZ][FA][2][0]);
		sigma[FZ][FA][2][0] = 0.;
	      }
	      int N = g_incl[FZ][2][0]->GetN();
	      g_incl[FZ][2][0]->SetPoint(N, FA, sigma[FZ][FA][2][0]);
	      g_incl[FZ][2][0]->SetPointError(N, 0, Esigma[FZ][FA][2][0]);
	    }
	  }
	}
	//
	// -1p
	k = FZ - 1;
	//
	I1p[i] = 0.;
	E_tmp = 0.;
	E1p[i] = 0.;
	//
	rangeA[0] = (int)(MINAOQ*(double)k+0.3);
	int i_1p = FA - rangeA[0] - 1;
	if((f_aoq[i][NUMPADDLE][FZ][FA][k]) && i_1p >= 0){
	  for(int p =0; p<NUMPADDLE; p++){
	    if(!f_aoq[i][p][FZ][FA][k]) continue;
	    I1p[i] += TMath::Sqrt(2*TMath::Pi()) / width_bin
	      * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1p + 1)
	      * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1p + 3);
	    E_tmp += 2*TMath::Pi() / (width_bin * width_bin)
	      * (TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1p + 1) 
			      * f_aoq[i][p][FZ][FA][k]->GetParError(3*i_1p + 3),2)
		 +TMath::Power(f_aoq[i][p][FZ][FA][k]->GetParError(3*i_1p + 1)
			       * f_aoq[i][p][FZ][FA][k]->GetParameter(3*i_1p + 3),2));
	  }
	  E1p[i] = TMath::Sqrt(E_tmp);
	  //
	  if(i>0){
	    double frac_targ = I1p[i]/(double)num_pid[i][NUMPADDLE][FZ][FA]/trans[i][FZ][FA];
	    double frac_empty = I1p[0]/(double)num_pid[0][NUMPADDLE][FZ][FA]/trans[0][FZ][FA];
	    if(frac_targ >0 && frac_empty > 0){
	      sigma[FZ][FA][i-1][1] = (frac_targ - frac_empty)/d_targ[i]*to_mb;
	      Esigma[FZ][FA][i-1][1]= to_mb/d_targ[i] *
		TMath::Sqrt(TMath::Power(frac_targ*E1p[i]/I1p[i],2) + frac_targ*frac_targ/(double)num_pid[i][NUMPADDLE][FZ][FA]
			    + TMath::Power(frac_targ*Etrans[i][FZ][FA]/trans[i][FZ][FA],2)
			    + TMath::Power(frac_empty*E1p[0]/I1p[0],2) + frac_empty*frac_empty/(double)num_pid[0][NUMPADDLE][FZ][FA]
			    + TMath::Power(frac_empty*Etrans[0][FZ][FA]/trans[0][FZ][FA],2));
	    }else if(frac_empty <= 0){ // Tentatively assume empty contribution is 0 for errors
	      sigma[FZ][FA][i-1][1] = (frac_targ)/d_targ[i]*to_mb;
	      Esigma[FZ][FA][i-1][1]= to_mb/d_targ[i] *
		TMath::Sqrt(TMath::Power(frac_targ*E1p[i]/I1p[i],2) + frac_targ*frac_targ/(double)num_pid[i][NUMPADDLE][FZ][FA]
			    + TMath::Power(frac_targ*Etrans[i][FZ][FA]/trans[i][FZ][FA],2));
	    }
	    if(sigma[FZ][FA][i-1][1] >= 0){
	      int N = g_incl[FZ][i-1][1]->GetN();
	      g_incl[FZ][i-1][1]->SetPoint(N, FA, sigma[FZ][FA][i-1][1]);
	      g_incl[FZ][i-1][1]->SetPointError(N, 0, Esigma[FZ][FA][i-1][1]);
	    }else{
	      if(TMath::Abs(sigma[FZ][FA][i-1][1])>TMath::Abs(Esigma[FZ][FA][i-1][1])) Esigma[FZ][FA][i-1][1] = TMath::Abs(sigma[FZ][FA][i-1][1]);
	      sigma[FZ][FA][i-1][1] = 0.;
	    }
	  }
	  if(i==2){//for proton
	    if(sigma[FZ][FA][1][1] >= 0 && sigma[FZ][FA][0][1] >= 0){
	      sigma[FZ][FA][2][1] = 0.5 * (sigma[FZ][FA][1][1] - sigma[FZ][FA][0][1]);
	      Esigma[FZ][FA][2][1] = 0.5 * TMath::Sqrt(Esigma[FZ][FA][0][1]*Esigma[FZ][FA][0][1] + Esigma[FZ][FA][1][1]*Esigma[FZ][FA][1][1]);
	      if(sigma[FZ][FA][2][1]<0){
		if(TMath::Abs(sigma[FZ][FA][2][1]) < TMath::Abs(Esigma[FZ][FA][2][1])) Esigma[FZ][FA][2][1] = TMath::Abs(sigma[FZ][FA][2][1]);
		sigma[FZ][FA][2][1] = 0.;
	      }
	      int N = g_incl[FZ][2][1]->GetN();
	      g_incl[FZ][2][1]->SetPoint(N, FA, sigma[FZ][FA][2][1]);
	      g_incl[FZ][2][1]->SetPointError(N, 0, Esigma[FZ][FA][2][1]);
	    }
	  }
	}
	////
	//cout <<Iunreacted[i]<<", "<<Eunreacted[i]<<", "<<I1p[i]<<", "<<E1p[i]<<", "<<I1n[i]<<", "<<E1n[i]<<", ";
	//
      }
      for(int i =0;i<NUMTARGET*2;i++){
	cout<<sigma[FZ][FA][i%NUMTARGET][i/NUMTARGET]<<", +/-, "<<Esigma[FZ][FA][i%NUMTARGET][i/NUMTARGET]<<", ";
	fcsv<<sigma[FZ][FA][i%NUMTARGET][i/NUMTARGET]<<", +/-, "<<Esigma[FZ][FA][i%NUMTARGET][i/NUMTARGET]<<", ";
      }
      /*
	for(int p=0; p<NUMPADDLE; p++){
	fcsv<<num_pid[i][p][FZ][FA]<<", ";
	cout<<num_pid[i][p][FZ][FA]<<", ";
	}*/
    }
    for(int i =0;i<NUMTARGET*2;i++)
      g_incl[FZ][i/2][i%2]->Write();// i%2: reaction, i/2: target
  }// for FZ
  fcsv<<endl;
  cout<<endl;
  fcsv.close();
}

int loadfile(){
  for(int i=0; i<NUMTARGET; i++){
    TString temp = dir_rootfile + filename;
    temp.ReplaceAll("FRS",frsname);
    temp.ReplaceAll("TARG",targetname[i]);
    cout<<"Load file: "<<temp<<endl;
    fin[i] = new TFile(temp,"READ");
    if(fin[i]->IsOpen()){
      cout<<"Successfully open."<<endl;
    }else{
      cerr<<"Failed to open file."<<endl;
      return 1;
    }
    //
    tree[i] = (TTree*)fin[i]->Get("FragTree");
    tree[i] ->SetTitle(targetname[i]);
    //tree[i]->Print();
    if(MAXEVENT<1){
      nentry[i] = tree[i]->GetEntries();
    }else{
      nentry[i] = MAXEVENT;
    }
    cout<<nentry[i]<<"entries."<<endl;
  }
  //for 38Ca FRS122
  if(!frs13){
    for(int i=0; i<NUMTARGET; i++){
      TString temp = dir_rootfile + filename;
      temp.ReplaceAll("FRS","38Ca_122");
      temp.ReplaceAll("TARG",targetname[i]);
      cout<<"Load file: "<<temp<<endl;
      fin2[i] = new TFile(temp,"READ");
      if(fin2[i]->IsOpen()){
	cout<<"Successfully open."<<endl;
      }else{
	cerr<<"Failed to open file."<<endl;
	return 1;
      }
      //
      tree2[i] = (TTree*)fin2[i]->Get("FragTree");
      tree2[i] ->SetTitle(targetname[i]);
      //tree[i]->Print();
      if(MAXEVENT<1){
	nentry2[i] = tree2[i]->GetEntries();
      }else{
	nentry2[i] = MAXEVENT;
      }
      cout<<nentry2[i]<<"entries."<<endl;
    }
  }
  //
  TString temp = dir_rootfile + outfilename;
  fout = new TFile(temp, "RECREATE");
  return 0;
}



Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag, Int_t i, Int_t targ){
  Int_t tmpzet = (Int_t) (zet + 0.5);
  Int_t tmpmass = (Int_t) ((Double_t)tmpzet * aoq + 0.5);
  Double_t tmpaoq = (Double_t)tmpmass/((Double_t)tmpzet);
  //
  Double_t rangezet = 0.;
  Double_t rangeaoq = 0.;
  Double_t value = 0.;
  if(!isfrag){
    rangezet = 1.11e-1;
    rangeaoq = 1.34e-3;
    value = pow(((Double_t)tmpzet - zet)/rangezet, 2) + pow((tmpaoq - aoq)/rangeaoq, 2);
  }else if(0 <= i && i < NUMPADDLE && targ>=0 && targ<NUMTARGET){
    if(tmpmass<0||tmpmass>=MAXA||tmpzet<MINZ||tmpzet>=MAXZ) return NAN;
    //pid_fraggate[NUMTARGET][NUMPADDLE][MAXZ][MAXA][4]={NAN}; // 0:lowZ, 1:highZ, 2:lowaq, 3:highaq
    rangezet = pid_fraggate[targ][i][tmpmass][tmpzet][1];
    rangeaoq = pid_fraggate[targ][i][tmpmass][tmpzet][3];
    value = pow((pid_fraggate[targ][i][tmpmass][tmpzet][0] - zet)/rangezet, 2) + pow((pid_fraggate[targ][i][tmpmass][tmpzet][2] - aoq)/rangeaoq, 2);
  }else{
    return NAN;
  }
  //
  if(!isnan(value)){
    zet = (Double_t) tmpzet;
    aoq = (Double_t) tmpaoq;
  }else{
    zet = NAN;
    aoq = NAN;
  }
  return TMath::Sqrt(value);
}

void get_mem_usage(){
  gSystem->GetProcInfo(&pinfo);
  cout<<" Memory usage (MB) " << pinfo.fMemResident /1024. << ", Virtual" << pinfo.fMemVirtual/1024.<<endl;
}
