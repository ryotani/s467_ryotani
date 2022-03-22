const int NUMTARGET=3, NUMPADDLE=28, MAXEVENT=1000000;
const int MINZ=10, MAXZ=30, MAXA=60, BINAOQ=500, BINZ=500;
const double MINAOQ=1.9, MAXAOQ=2.75;
//
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
//#include "TSpectrum.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
using namespace std;
//
TString dir_rootfile = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/fragment/output/" ;
TString filename = "tofw_beta_offset_paddle_50Ca_TARG_nofrsgate_Feb.root"; // TARG will be replaced
TString outfilename = "incl_out_test_indivfit.root";
TString targetname[3] = {"empty","carbon","ch2"};
TFile *fin[NUMTARGET], *fout;
TTree *tree[NUMTARGET];
//
// Def variables
Long64_t nentry[NUMTARGET]={0};
Double_t pid_fraggate[NUMTARGET][NUMPADDLE][MAXZ][MAXA][4]={NAN}; // 0:Z, 1:sigmaZ, 2:aq, 3:sigmaaq
//
// Def histograms
TH1D* h_music[NUMTARGET][NUMPADDLE];
TH1I* h_aoq_gated[NUMTARGET][NUMPADDLE][MAXZ][MAXA][MAXZ];
TH2I* h_pid[NUMTARGET][NUMPADDLE];
TH2S* h_gated_pid[NUMTARGET][NUMPADDLE][MAXZ][MAXA];
TF1* f_music[NUMTARGET][NUMPADDLE], *f_aoq[NUMTARGET][NUMPADDLE][MAXZ][MAXA][MAXZ];
//
// Def functions
void calc_incl();
int loadfile();
void fill_histos(int i), fit_music(int i), fill_aoq_gated(int i);
Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag=false, Int_t i=0, Int_t targ=-1);
TF1* generate_func(char* name, int numpeaks, double minpos, double interval, double par_limit_range, double height_limit=1e10);
void delete_histos(int i=0);
//
ProcInfo_t pinfo;
void get_mem_usage();
/////////////////////////////////////////////
void calc_incl(){
  ROOT::EnableImplicitMT();
  get_mem_usage();
  if(loadfile()!=0) return;
  //
  for(int i=0; i<NUMTARGET; i++){
    fill_histos(i);
    fit_music(i);
    //if(NUMTARGET>1) delete_histos(i); // To reduce memory usage
    fill_aoq_gated(i);
  }
}
/////////////////////////////////////////////
void fill_histos(int i){
  Float_t FRSAoQ, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
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
  tree[i]->SetBranchAddress("Mw1_X", &Mw1_X);
  tree[i]->SetBranchAddress("Mw2_X", &Mw2_X);
  tree[i]->SetBranchAddress("Mw3_X", &Mw3_X);
  tree[i]->SetBranchAddress("Mw1_Y", &Mw1_Y);
  tree[i]->SetBranchAddress("Mw2_Y", &Mw2_Y);
  tree[i]->SetBranchAddress("Mw3_Y", &Mw3_Y);
  tree[i]->SetBranchAddress("Tofw_Y", &Tofw_Y);
  tree[i]->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  //
  for(int p =0; p<NUMPADDLE; p++){
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
  int neve=0;
  for(Long64_t n=0; n<nentry[i]; n++){
    if(++neve%10000==0)
      cout<<"\r"<<neve<<" enrties done in "<<nentry[i]<<flush;
    //
    tree[i]->GetEntry(n);
    //
    int p = Tofw_Paddle-1;
    if(p >= NUMPADDLE || p < 0) continue;
    h_pid[i][p]->Fill(FragAoQ,FragZ);
    //
    /* // h_gated_pid is not used
    double tmpzet = MusicZ, tmpaoq = FRSAoQ;
    double sigma=pid(tmpzet, tmpaoq);
    if(isnan(sigma) || sigma > 3.) continue; // Get FRS PID
    int tmpA = tmpaoq*tmpzet + 0.1;
    int tmpZ = tmpzet + 0.1;//0.1 is to make sure the conversion from double to int to be correct
    if(tmpZ<MAXZ && tmpA<MAXA)
      h_gated_pid[i][p][tmpZ][tmpA]->Fill(FragAoQ_corr,FragZ);
    */
  }
  cout<<endl;
  //
  fout->cd();
  for(int p =0; p<NUMPADDLE; p++)
    {
      cout<<"Paddle"<< p+1 << ", max bin: " << h_pid[i][p]->GetMaximum()<<endl;
      if(h_pid[i][p]->GetMaximum()>30000){//limit for TH2S
	cout<<"Too many events for TH2S"<<endl;
      }
      /*
      h_pid[i][p]->Write();
      for(int j=0; j<MAXZ*MAXA; j++){
	h_gated_pid[i][p][j/MAXA][j%MAXA]->Write();
      }
      */
    }
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
      delete h_gated_pid[i][p][j/MAXA][j%MAXA];
    }
    cout<<"Paddle"<<p+1<<endl;
  }
  cout<<"Deleted filled histos."<<endl;
  get_mem_usage();
}
//
void fit_music(int i){
  //h_music[i] = new TH1F(dummy,dummy,MAXZ-5,10.5,MAXZ+5.5);
  //h_music[i]->Fill(FragZ);
  for(int p =0; p<NUMPADDLE; p++)
    {
      h_music[i][p] = h_pid[i][p]->ProjectionY();
      int numpeaks = MAXZ-MINZ;
      cout<<"Entries for paddle"<<p+1<<": "<<h_music[i][p]->GetEntries()<<endl;
      if(h_music[i][p]->GetEntries()<1){
	cout<<"Skippling paddle"<<p+1<<endl;
	continue;
      }
      /*
      TSpectrum *sp = new TSpectrum(numpeaks);
      //numpeaks = sp->Search(h_music[i][p]);
      auto posX = sp->GetPositionX();
      auto posY = sp->GetPositionY();
      */
      f_music[i][p] = generate_func(Form("f_music_%s_%i",targetname[i].Data(),p), numpeaks, MINZ, 1, 0.2, 0);
      f_music[i][p] ->FixParameter(0,0.);
      //
      Int_t Mx, My, Mz;
      h_music[i][p]->GetMaximumBin(Mx, My,Mz);
      Double_t MCentre = h_music[i][p]->GetXaxis()->GetBinCenter(Mx);
      Int_t IMax = (Int_t)MCentre - MINZ;
      if(IMax<numpeaks){
	f_music[i][p] ->FixParameter(3*IMax+1, h_music[i][p]->GetMaximum());
	f_music[i][p] ->FixParameter(3*IMax+2, MCentre);
      }
      //
      for(int j=0; j<numpeaks; j++){
	f_music[i][p] ->SetParLimits(3*j+1,0., 1.2 * h_music[i][p]->GetMaximum());
	h_music[i][p] ->Fit(f_music[i][p],"Q","",MINZ+j-0.3,MINZ+j+0.3);
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
      f_music[i][p]->SetParLimits(3*IMax+2, MCentre-0.3, MCentre+0.3);
      f_music[i][p]->SetParLimits(0,0.,0.01*h_music[i][p] ->GetMaximum());
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
	  f_music[i][p] ->SetParLimits(3*j+1,0.2*params[3*j+1],10.*params[3*j+1]);
	  f_music[i][p] ->SetParameter(3*j+1,params[3*j+1]);
	}
      }
      h_music[i][p] ->Fit(f_music[i][p],"LL","",MINZ-0.3, MAXZ+0.3);
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
      h_music[i][p]->Write();
    }
}
//
void fill_aoq_gated(int i){
    Float_t FRSAoQ, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
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
  tree[i]->SetBranchAddress("Mw1_X", &Mw1_X);
  tree[i]->SetBranchAddress("Mw2_X", &Mw2_X);
  tree[i]->SetBranchAddress("Mw3_X", &Mw3_X);
  tree[i]->SetBranchAddress("Mw1_Y", &Mw1_Y);
  tree[i]->SetBranchAddress("Mw2_Y", &Mw2_Y);
  tree[i]->SetBranchAddress("Mw3_Y", &Mw3_Y);
  tree[i]->SetBranchAddress("Tofw_Y", &Tofw_Y);
  tree[i]->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  //
  for(int p =0; p<NUMPADDLE; p++){
    cout<<"Preparing histos for aoq fit"<<p+1<<" ";
    for(int j=MINZ*MAXA; j<MAXZ*MAXA; j++){
      for(int k =MINZ; k<MAXZ; k++){
	if(isnan(pid_fraggate[i][p][k][0][0])) continue; // Check the Twim gates.
	TString dummy = Form("h_aoq_gated_paddle%i_Z%i_A%i_Z%i_%s", p+1, j/MAXA, j%MAXA, k, targetname[i].Data());
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
    //
    int p = Tofw_Paddle-1;
    if(p >= NUMPADDLE || p < 0) continue;
    //
    double tmpzet = MusicZ, tmpaoq = FRSAoQ;
    double sigma=pid(tmpzet, tmpaoq);
    if(isnan(sigma) || sigma > 3.) continue; // Get FRS PID
    int tmpA = tmpaoq*tmpzet + 0.1;
    int tmpZ = tmpzet + 0.1;//0.1 is to make sure the conversion from double to int to be correct
    if(isnan(MusicZ*FRSAoQ*FragZ*FragAoQ_corr)) continue;
    int tmpZ_frag = FragZ + 0.5, tmpA_frag = (FragZ*FragAoQ_corr) + 0.5;
    if(tmpZ>=MAXZ || tmpA>=MAXA || tmpZ_frag>=MAXZ || tmpA_frag>=MAXA) continue;
    if(isnan(pid_fraggate[i][p][tmpZ_frag][0][0])) continue;
    if(TMath::Abs(FragZ-pid_fraggate[i][p][tmpZ_frag][0][0]) > 3. * pid_fraggate[i][p][tmpZ_frag][0][1]) continue; // For Z, 3 sigma cut applied.
    h_aoq_gated[i][p][tmpZ][tmpA][tmpZ_frag]->Fill(FragAoQ_corr);
  }
  cout<<endl;
  //
  fout->cd();
  //
  for(int p =0; p<NUMPADDLE; p++){
    cout<<"Fitting AoQ for Paddle"<<p+1<<endl;
    for(int ZA=MINZ*MAXA; ZA<MAXZ*MAXA; ZA++){//FRS
      for(int k =MINZ; k<=ZA/MAXA; k++){//Frag
	int rangeA[2] = {(int)(MINAOQ*(double)k+0.3), (int)min(MAXAOQ*(double)k-0.3,(double)MAXA-1)};
	int numpeaks = rangeA[1]-rangeA[0]+1;
	if(numpeaks<1) continue;
	auto h_tmp = h_aoq_gated[i][p][ZA/MAXA][ZA%MAXA][k];
	//cout<<h_tmp->GetEntries()<<endl;
	if(h_tmp->GetEntries()<numpeaks) continue;
	auto f_tmp = generate_func(Form("f_aoq_%s_%i_%i_%i_%i",targetname[i].Data(),p, ZA/MAXA, ZA%MAXA,k), numpeaks,
			      (double)rangeA[0]/(double)k, 1./(double)k, 0.4/(double)k, 0);
	cout<<"Paddle"<<p+1<<", Z"<<ZA/MAXA<<", A"<<ZA%MAXA<<", Zfrag"<<k<<endl;
	f_aoq[i][p][ZA/MAXA][ZA%MAXA][k] = f_tmp;
	f_tmp->FixParameter(0,0.);
	Int_t Mx, My, Mz;
	h_tmp->GetMaximumBin(Mx, My,Mz);
	Double_t MCentre = h_tmp->GetXaxis()->GetBinCenter(Mx);
	Int_t AMax = (Int_t)(MCentre*(Double_t)k) - rangeA[0];
	if(AMax<numpeaks){
	  f_tmp ->FixParameter(3*AMax+1, h_tmp->GetMaximum());
	  f_tmp ->SetParameter(3*AMax+2, MCentre);
	  f_tmp ->SetParLimits(3*AMax+2, MCentre-0.1/(double)k, MCentre+0.1/(double)k);
	}
	//
	for(int j=numpeaks-1; j>=0; j--){
	  f_tmp->SetParLimits(3*j+1,0., 2.* h_tmp->GetMaximum());
	  h_tmp->Fit(f_tmp,"QL","",((double)(rangeA[0]+j)-0.3)/(double)k,((double)(rangeA[0]+j)+0.3)/(double)k);
	  f_tmp->FixParameter(3*j+1,f_tmp->GetParameter(3*j+1));
	}
	//h_tmp->Fit(f_tmp,"QLL","",(double)(rangeA[0]-1)/(double)k, (double)(rangeA[1]+1)/(double)k);
	//
	f_tmp->SetParLimits(3*AMax+2, MCentre-0.3/(double)k, MCentre+0.3/(double)k);
	f_tmp->SetParLimits(0,0.,0.01*h_tmp->GetMaximum());
	auto params = f_tmp->GetParameters();
	for(int j=0; j<numpeaks; j++){
	  if(params[3*j+1]<1.){
	    f_tmp ->SetParameter(3*j+1,0.);
	    f_tmp ->SetParameter(3*j+2,0.);
	    f_tmp ->SetParameter(3*j+3,0.1);
	    f_tmp ->FixParameter(3*j+1,0.);
	    f_tmp ->FixParameter(3*j+2,0.);
	    f_tmp ->FixParameter(3*j+3,0.1);
	  }else{
	    f_tmp ->SetParLimits(3*j+1,0.1*params[3*j+1],2.*params[3*j+1]);
	    f_tmp ->SetParameter(3*j+1,params[3*j+1]);
	  }
	}
	h_tmp->Fit(f_tmp,"LL","",((double)rangeA[0]-0.3)/(double)k, ((double)rangeA[1]+0.3)/(double)k);
	//
	cout<<"AoQ fit result for paddle"<<p+1<<", Z"<<ZA/MAXA<<", A"<<ZA%MAXA<< ", Z_frag"<<k<<". Chi-sq:"<<f_music[i][p]->GetChisquare()<<", NDF:"<<f_music[i][p]->GetNDF()<<endl<<endl;
	//
	for(int j=0; j<numpeaks; j++){
	  if(params[3*j+1]<1.)continue;
	  int i_A = rangeA[0] + j;
	  //for(int i_aq=0; i_aq<MAXA; i_aq++){
	  pid_fraggate[i][p][k][i_A][2] = params[3*j+2];
	  pid_fraggate[i][p][k][i_A][3] = params[3*j+3];
	  for(int i_param=0; i_param<4; i_param++) cout<<pid_fraggate[i][p][k][i_A][i_param]<<" ";
	  cout<<"A:"<<pid_fraggate[i][p][k][i_A][0]*pid_fraggate[i][p][k][i_A][2]<<endl;
	  //}
	}
	h_tmp->Write();
      }
    }
  }
}
//
//
//////////
TF1* generate_func(char*name, int numpeaks, double minpos, double interval, double par_limit_range, double height_limit){
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
    }
    f_tmp->SetParLimits(3*i+2, peakpos-par_limit_range, peakpos+par_limit_range);
    f_tmp->SetParLimits(3*i+3, 0.01*par_limit_range, par_limit_range);
  }
  //cout<<formula<<" func generated"<<endl;
  f_tmp->SetVectorized(true);
  return f_tmp;
}

int loadfile(){
  for(int i=0; i<NUMTARGET; i++){
    TString temp = dir_rootfile + filename;
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
