#define MINPADDLE 0
#define NUMPADDLE 28

#include "fragmentana_twim.h"
int zetdiff = 0;
double min_aoq = 1.8, max_aoq = 2.7;
//
//for trees
Float_t FRSAoQ, FRSBeta, FRSBrho, FRSGamma, MusicE, TwimE, MusicZ, FragZ, MusicTheta, TwimTheta, Mw0_X, Mw1_X, Mw2_X, Mw3_X, Mw0_Y, Mw1_Y, Mw2_Y, Mw3_Y, ROLU_X, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho, FragBeta, FragGamma;
UChar_t Tofw_Paddle;
Long64_t nentry=0;
TString targ, infile;
//
TString outpdf ;
TString outcsv ;
TString outroot;
//
TString sofiacaldir ="./parameters/";
TString FRS;
TString inpar;
TString outpar;
//
bool IsEmpty=false;
const Double_t sigma_frs = 4., sigma = 2., sigma3 = 3.;//Set conditions for PID gating
const Int_t MINA=36, MAXA=55, MINZ=16, MAXZ=23;
Int_t NumGoodPaddles=0;
Int_t GoodPaddles[NUMPADDLE]={0};
TGraph *g_efflength, *g_tofoffset;
TF1* f_efflength, *f_tofoffset;
//
ofstream fcsv;
TFile *fout;
TF1 *fit_prof[NUMPADDLE];
TH2F *h_frsbeta_betacalc[NUMPADDLE], *h_deltabeta_tof[NUMPADDLE], *h_deltabeta_twime[NUMPADDLE], *h_deltabeta_theta[NUMPADDLE], *h_deltabeta_pid[NUMPADDLE];
TH2F *h_gatedpid[NUMPADDLE][100][50], *h_zgated_mw3x_aoq[MAXZ];
TH1F *h_zgatedpid[NUMPADDLE][MAXZ];
Int_t NumGated[NUMPADDLE][100][50][100][50] = {0}, NumGated3[NUMPADDLE][100][50][100][50] = {0}; // Count hits in 2sigma and 3sigma
TH1F *h_noreaction[NUMPADDLE], *h_counts_paddle[100][50][3][3], *h_counts_mw3x[100][50][3][3], *h_counts_mw3x3[100][50][3][3];
TF1 *f_FitGated[100][50][3][3];
TF2 *f_fragfit[NUMPADDLE][100][50];
TF1 *f_aoq[NUMPADDLE], *f_betatof[NUMPADDLE], *f_betatheta[NUMPADDLE];
TTree *tree;

void frag_tofw_calib(int i_target = 0, int FRS_setting = 13);
void define_conditions();
void delta_beta_method();
//void fit_beta_paddle();
void draw_aoq_paddle(Double_t zet);
void initialise();
void filltree();
Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag=false, Int_t i=0);
Int_t initbranch();
void writecsv();
void writeparam();

//////////
void frag_tofw_calib(int i_target, int FRS_setting){
  //ROOT::EnableImplicitMT();
  if (FRS_setting==13){
    FRS="50Ca";
    inpar = sofiacaldir + "FRS13.par";
  }else if(FRS_setting==11 || FRS_setting==12){
    FRS="38Ca";
    inpar = sofiacaldir + "FRS11.par";
  }else if(FRS_setting==122){
    FRS="38Ca_122"; // only the latter half
    inpar = sofiacaldir + "FRS122.par";
  }else{
    cerr<<"Incorrect FRS number"<<endl;
    return;
  }
  //
  switch (i_target){
  case 0:
    targ = "empty";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_empty.root";
    IsEmpty = true;
    if(FRS_setting==13){
      NumGoodPaddles=4;
      for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 12 + i;
    }else{
      NumGoodPaddles=5;
      for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 9 + i;
    }
    break;
  case 1:
    targ = "ch2";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_ch2-24mm.root";
    if(FRS_setting==13){
      /*
	NumGoodPaddles=6;
	for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 14 + i;
      */
      NumGoodPaddles=4;
      for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 14 + i;
    }else{
      NumGoodPaddles=8;
      for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 12 + i;
    }
    break;
  case 2:
    targ = "carbon";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_carbon.root";
    if(FRS_setting==13){
      NumGoodPaddles=5;
      for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 13 + i;
    }else{
      NumGoodPaddles=6;
      for(int i=0; i<NumGoodPaddles; i++) GoodPaddles[i] = 12 + i;
    }
    break;
  case 3:
    targ = "PP";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+ FRS +"_PP.root";
    break;
  default:
    cerr<<"No target info"<<endl;
    return;
  }
  outpdf  = "./fragment/output/tofw_calib_" + FRS + "_"+ targ + Form("_May") + ".pdf";
  outcsv  = "./fragment/output/tofw_calib_" + FRS + "_"+ targ + Form("_May") + ".csv";
  outroot = "./fragment/output/tofw_calib_" + FRS + "_"+ targ + Form("_May") + ".root";
  outpar  = "./fragment/output/frag_tofw_calib_" + FRS + "_" + targ + ".par";
  //
  fcsv.open(outcsv, ofstream::out);
  fout = new TFile(outroot,"RECREATE");
  
  cout<< targ<< " "<<infile<<endl;

  initialise();
  c = new TCanvas();
  c->Divide(4,4);
  c->Print(outpdf + "[");
  //
  define_conditions();
  //
  initbranch();
  delta_beta_method();
  //
  c->Print(outpdf + "]");
  //
  filltree();
  writecsv();
  writeparam();
}

void define_conditions(){
}

void delta_beta_method(){
  // Unreacted cut
  nentry = ch->GetEntries();
  Long64_t neve = 0;
  for(Long64_t n=0; n<nentry; n++){
    ch->GetEntry(n);
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    if( TMath::Abs(FragZ-tmpzet)>0.4) continue;
    if((Mw0_X + 405.9 * MusicTheta < 0)||(Mw0_X + 405.9 * MusicTheta > 15)) continue; // ROLU events are not good always
    h_noreaction[i]->Fill(tmpaoq - FragAoQ);
    Int_t i_fragZ=(0.5+FragZ);
    if(abs(FragZ-(double)i_fragZ)<0.3&&i_fragZ<MAXZ&&i_fragZ>=MINZ){
      h_zgatedpid[i][i_fragZ]->Fill(FragAoQ);
      h_zgated_mw3x_aoq[i_fragZ]->Fill(Mw3_X,FragAoQ);
      //
    }
  }
  for(int i_fragZ =MINZ; i_fragZ<MAXZ; i_fragZ++){
    h_zgated_mw3x_aoq[i_fragZ]->Write();
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
      //c -> cd((i%(NUMPADDLE/2))+1);
      TSpectrum *tspec = new TSpectrum(10);
      Int_t npeaks = tspec->Search(h_zgatedpid[i][i_fragZ]);
      Double_t *peakXspec = tspec->GetPositionX();
      //peakX[i]=2000.;
      for(int j = 0; j<npeaks; j++){
	//peakX[i] = (abs(peakX[i]) < abs(peakXspec[j]))? peakX[i] : peakXspec[j];
	//if(TMath::Abs(peakX[i]) > TMath::Abs(peakXspec[j])) peakX[i] = peakXspec[j];
	cout<<peakXspec[j]<<endl;//" "<<peakX[i]<<endl;
      }
      h_zgatedpid[i][i_fragZ]->Write();
    }
  }
  //
  Double_t peakX[NUMPADDLE]={1000.};
  int i_maxpaddle=0;
  int maxhits=0;
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    if(maxhits < h_noreaction[i]->GetEntries()){
      maxhits = h_noreaction[i]->GetEntries();
      i_maxpaddle = i;
    }
  }
  //
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    TSpectrum *tspec = new TSpectrum(5);
    Int_t npeaks = tspec->Search(h_noreaction[i]);
    Double_t *peakXspec = tspec->GetPositionX();
    peakX[i]=2000.;
    for(int j = 0; j<npeaks; j++){
      if(TMath::Abs(peakX[i]) > TMath::Abs(peakXspec[j])) peakX[i] = peakXspec[j];
      cout<<peakXspec[j]<<" "<<peakX[i]<<endl;
    }
    //peakX[i] = h_noreaction[i]->GetBinCenter(h_noreaction[i]->GetMaximumBin());
    f_aoq[i] = new TF1(Form("f_aoqdiff%i",i+1),"gaus");
    cout<<"peakX"<<i+1<<": "<<peakX[i]<<endl;
    f_aoq[i]->SetParameter(0,h_noreaction[i]->GetMaximum());
    f_aoq[i]->FixParameter(1,peakX[i]);
    f_aoq[i]->SetParameter(2,0.02);
    h_noreaction[i]->Fit(f_aoq[i],"L","",peakX[i]-0.03,peakX[i]+0.03);
    //
    //
    h_noreaction[i]->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  //
  neve = 0;
  for(Long64_t n=0; n<nentry; n++){
  //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    if((Mw0_X + 405.9 * MusicTheta < 0)||(Mw0_X + 405.9 * MusicTheta > 15)) continue;
    //if((Mw0_X + 405.9 * MusicTheta > 0)&&(Mw0_X + 405.9 * MusicTheta < 15)) continue;
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    //cout<<i<<": "<<MusicZ<<" "<<FRSAoQ<<" "<<tmpzet<<" "<<tmpaoq<<endl;
    //if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2)  || TMath::Abs(FragZ-tmpzet)>0.4) continue;
    if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2)) continue; // Z gate not needed?
    Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    h_deltabeta_tof[i]->Fill(FragTof, betacalc);
    h_frsbeta_betacalc[i]->Fill(FRSBeta, betacalc);
    //
    if(++neve%100000==0)
      cout<<"\r"<<i<<" entries done in "<<nentry<<flush;//endl;//
  }
  cout<<endl;
  //Draw
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_deltabeta_tof[i]->Draw("colz");
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_deltabeta_tof[i]->Draw("colz");
    f_betatof[i]->FixParameter(2,0);
    h_deltabeta_tof[i]->Fit(Form("f_betatof%i",i+1),"","",h_deltabeta_tof[i]->GetMean() - 4.* h_deltabeta_tof[i]->GetStdDev(), h_deltabeta_tof[i]->GetMean() + 4.* h_deltabeta_tof[i]->GetStdDev());
    h_deltabeta_tof[i]->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_frsbeta_betacalc[i]->Draw("colz");
    h_frsbeta_betacalc[i]->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //
  g_efflength = new TGraph();
  g_tofoffset = new TGraph();
  g_efflength->SetName("g_efflength");
  g_tofoffset->SetName("g_tofoffset");
  for(int i = 0; i<NumGoodPaddles; i++){
    g_efflength->SetPoint(g_efflength->GetN(),GoodPaddles[i], f_betatof[GoodPaddles[i]-1]->GetParameter(0));
    g_tofoffset->SetPoint(g_tofoffset->GetN(),GoodPaddles[i], f_betatof[GoodPaddles[i]-1]->GetParameter(1));
  }
  f_efflength = new TF1("f_efflength","pol1");
  g_efflength->Fit(f_efflength);
  g_efflength->Write();
  f_tofoffset = new TF1("f_tofoffset","pol1");
  g_tofoffset->Fit(f_tofoffset);
  g_tofoffset->Write();
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    f_betatof[i]->SetParameter(0, f_efflength->Eval((double)i+1.));
    f_betatof[i]->SetParameter(1, f_tofoffset->Eval((double)i+1.));
    //Ryo
  }
  //
  // Draw PID
  //
  for(Long64_t n=0; n<nentry; n++){
  //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    if(n%100000==0)
      cout<<"\r"<<neve<< " in "<<n<<" entries done in "<<nentry<<flush;//endl;//
    if((Mw0_X + 405.9 * MusicTheta < 0)||(Mw0_X + 405.9 * MusicTheta > 15)) continue;
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    //Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    Double_t beta_fit = f_betatof[i]->GetParameter(0)/(FragTof + f_betatof[i]->GetParameter(1)) + f_betatof[i]->GetParameter(2);
    Double_t aoq_fit = (FragBrho * TMath::Sqrt(1.-beta_fit*beta_fit))/(beta_fit * mc_e);
    //cout<<i << " "<<beta_fit<< " "<<aoq_fit<<" "<<FragZ<<endl;
    h_deltabeta_pid[i]->Fill(aoq_fit,FragZ);
    ++neve;
  }
  cout<<endl;
  //Draw
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    //c -> cd((i%(NUMPADDLE/2))+1);
    c->cd(0);
    for (int A = MINA; A < MAXA; A++){
      for (int Z = MINZ; Z < MAXZ; Z++){
	h_deltabeta_pid[i]->Fit(Form("f_fragfit%i_%i_%i",i,A,Z),"R N 0","");
      }
    }
    h_deltabeta_pid[i]->Draw("colz");
    c->Print(outpdf);
    for (int A = MINA; A < MAXA; A++){
      for (int Z = MINZ; Z < MAXZ; Z++){
	if(f_fragfit[i][A][Z]->GetParameter(0)<=0.001*h_deltabeta_pid[i]->GetMaximum()) continue;
	TEllipse *el = new TEllipse(f_fragfit[i][A][Z]->GetParameter(1),
				    f_fragfit[i][A][Z]->GetParameter(3),
				    sigma*f_fragfit[i][A][Z]->GetParameter(2),
				    sigma*f_fragfit[i][A][Z]->GetParameter(4));
	el->SetLineColor(2);
	el->SetLineWidth(1);
	el->SetFillColor(0);
	el->SetFillStyle(0);
	//el->Draw("same");
	el->Draw();
      }
    }
    h_deltabeta_pid[i]->Write();
    c->Print(outpdf);
    //if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  // c->Print(outpdf);
  /*
  // Draw gated pid
  c->Clear();
  c->cd(0);
  c->Divide(4,4);
  for(Long64_t n=0; n<nentry; n++){
    //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    if((Mw0_X + 405.9 * MusicTheta < 0)||(Mw0_X + 405.9 * MusicTheta > 15)) continue;
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    Double_t beta_fit = f_betatof[i]->GetParameter(0)*(1. + 1000.*TwimTheta*f_betatheta[i]->GetParameter(1))/(FragTof + f_betatof[i]->GetParameter(1)) + f_betatof[i]->GetParameter(2);
    Double_t aoq_fit = (FragBrho * TMath::Sqrt(1.-beta_fit*beta_fit))/(beta_fit * mc_e);
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    Int_t tmpA = tmpaoq*tmpzet;
    Int_t tmpZ = tmpzet;
    if(tmpA<MINA||tmpA>=MAXA||tmpZ<MINZ||tmpZ>=MAXZ) continue;
    h_gatedpid[i][tmpA][tmpZ]->Fill(aoq_fit,FragZ);
    //
    NumGated[i][tmpA][tmpZ][0][0]++; // Incoming particles within acceptance
    //
    Double_t tmpfragzet = FragZ, tmpfragaoq = aoq_fit;
    Double_t sigma_tmp = pid(tmpfragzet, tmpfragaoq, true, i);
    if(sigma_tmp < sigma3){
      Int_t tmpfragA = tmpfragaoq*tmpfragzet;
      Int_t tmpfragZ = tmpfragzet;
      Int_t DZ = tmpZ-tmpfragZ;
      Int_t DN = tmpA-tmpfragA - DZ;
      NumGated3[i][tmpA][tmpZ][tmpfragA][tmpfragZ]++;
      if(DN>=0 && DZ>=0 && DN<3 && DZ<3){
	h_counts_mw3x3[tmpA][tmpZ][DN][DZ]->Fill(Mw3_X);
      }
      if(sigma_tmp < sigma){ // stricter gating
	NumGated[i][tmpA][tmpZ][tmpfragA][tmpfragZ]++;
	if(DN>=0 && DZ>=0 && DN<3 && DZ<3){
	  h_counts_mw3x[tmpA][tmpZ][DN][DZ]->Fill(Mw3_X);
	}
      }
    }
    / *
    else{
      NumGated[i][tmpA][tmpZ][0][0]++; // outside of the gate
      }* /
    //
    if(++neve%100000==0)
      cout<<"\r"<<n<<" entries done in "<<nentry<<flush;//endl;//
  }
  for (int Z1 = MINZ; Z1 < MAXZ; Z1++){
    for (int A1 = MINA; A1 < MAXA; A1++){
      for(int i = MINPADDLE; i<NUMPADDLE; i++){
	c -> cd((i%(NUMPADDLE/2))+1);     
	h_gatedpid[i][A1][Z1]->Draw("colz");
	h_gatedpid[i][A1][Z1]->Write();
	for (int A = MINA; A < MAXA; A++){
	  for (int Z = MINZ; Z < MAXZ; Z++){
	    if(f_fragfit[i][A][Z]->GetParameter(0)<=0.001*h_deltabeta_pid[i]->GetMaximum()) continue;
	    TEllipse *el = new TEllipse(f_fragfit[i][A][Z]->GetParameter(1),
					f_fragfit[i][A][Z]->GetParameter(3),
					sigma*f_fragfit[i][A][Z]->GetParameter(2),
					sigma*f_fragfit[i][A][Z]->GetParameter(4));
	    el->SetLineColor(2);
	    el->SetLineWidth(1);
	    el->SetFillColor(0);
	    el->SetFillStyle(0);
	    //el->Draw("same");
	    el->Draw();
	  }
	}
	//cout<<"paddle:"<<i+1<<", A:"<<A<<", Z:"<<Z<<" counts:"<<NumGated[i][A][Z][A][Z]<<endl;
	if(i==NUMPADDLE/2-1) c->Print(outpdf);
	//
	for(int D = 0; D<9; D++){
	  int DN=D%3, DZ=D/3;
	  h_counts_paddle[A1][Z1][D%3][D/3]->SetBinContent(i, NumGated[i][A1][Z1][A1-DN-DZ][Z1-DZ]);
	}
      }
      c->cd(15);
      h_counts_paddle[A1][Z1][0][0]->Draw();
      c->cd(16);
      h_counts_paddle[A1][Z1][1][1]->Draw(); //-1p
      c->Print(outpdf);
    }
  }
  //
  c->Clear();
  c->cd(0);
  c->Divide(3,3);
  for (int Z = MINZ; Z < MAXZ; Z++){
    for (int A = MINA; A < MAXA; A++){
      if(h_counts_paddle[A][Z][0][0]->Integral()==0) continue;
      for(int D = 0; D<9; D++){
	c->cd(D+1);
	h_counts_paddle[A][Z][D%3][D/3]->Draw();
      }
      c->Print(outpdf);
    }
  }
  for (int Z = MINZ; Z < MAXZ; Z++){
    for (int A = MINA; A < MAXA; A++){
      if(h_counts_mw3x[A][Z][0][0]->GetEntries()==0) continue;
      for(int D = 0; D<9; D++){
	c->cd(D+1);
	int DN=D%3, DZ=D/3;
  	f_FitGated[A][Z][DN][DZ] = new TF1(Form("f_FitGated%i%i%i%i",A,Z,DN,DZ), "gausn", -300, 300);
	h_counts_mw3x3[A][Z][DN][DZ]->Draw();
	h_counts_mw3x[A][Z][DN][DZ]->Draw("same");
	h_counts_mw3x[A][Z][DN][DZ]->Fit(f_FitGated[A][Z][DN][DZ], "LL R","");
      }
      c->Print(outpdf);
    }
  }  */
}

Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag, Int_t i){
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
  }else if(MINPADDLE <= i && i < NUMPADDLE){
    if(tmpmass<MINA||tmpmass>=MAXA||tmpzet<MINZ||tmpzet>=MAXZ) return NAN;
    rangezet = f_fragfit[i][tmpmass][tmpzet]->GetParameter(4);
    rangeaoq = f_fragfit[i][tmpmass][tmpzet]->GetParameter(2);
    value = pow((f_fragfit[i][tmpmass][tmpzet]->GetParameter(3) - zet)/rangezet, 2) + pow((f_fragfit[i][tmpmass][tmpzet]->GetParameter(1) - aoq)/rangeaoq, 2);
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
///
//void fit_beta_paddle(Int_i i, Double_t &offset, Double_t &eff_length){
//void fit_beta_paddle(){ // not used
  /* / Fit max hit paddle
  f_aoq[i_maxpaddle] = new TF1(Form("f_aoqdiff%i",i_maxpaddle+1),"gaus");
  peakX[i_maxpaddle] = h_noreaction[i_maxpaddle]->GetBinCenter(h_noreaction[i_maxpaddle]->GetMaximumBin());
  cout<<"peakX"<<i_maxpaddle+1<<": "<<peakX[i_maxpaddle]<<endl;
  f_aoq[i_maxpaddle]->SetParameter(0,h_noreaction[i_maxpaddle]->GetMaximum());
  f_aoq[i_maxpaddle]->FixParameter(1,peakX[i_maxpaddle]);
  f_aoq[i_maxpaddle]->SetParameter(2,0.02);
  h_noreaction[i_maxpaddle]->Fit(f_aoq[i_maxpaddle],"L","",peakX[i_maxpaddle]-0.03,peakX[i_maxpaddle]+0.03);
  h_noreaction[i_maxpaddle]->Write();
  //
  //
  for(Long64_t n=0; n<nentry; n++){
    ch->GetEntry(n);
    Int_t i = Tofw_Paddle - 1;
    if(i != i_maxpaddle) continue;
    if((Mw0_X + 405.9 * MusicTheta < 0)||(Mw0_X + 405.9 * MusicTheta > 15)) continue;
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    //if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2)  || TMath::Abs(FragZ-tmpzet)>0.4) continue;
    if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2)) continue; // Z gate not needed?
    Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    h_deltabeta_tof[i]->Fill(FragTof, betacalc);
    h_frsbeta_betacalc[i]->Fill(FRSBeta, betacalc);
  }
  //
  h_deltabeta_tof[i_maxpaddle]->Draw("colz");
  f_betatof[i_maxpaddle]->FixParameter(2,0);
  h_deltabeta_tof[i_maxpaddle]->Fit(Form("f_betatof%i",i_maxpaddle+1),"","",h_deltabeta_tof[i_maxpaddle]->GetMean() - 2.6* h_deltabeta_tof[i_maxpaddle]->GetStdDev(), h_deltabeta_tof[i_maxpaddle]->GetMean() + 2.6* h_deltabeta_tof[i_maxpaddle]->GetStdDev());
  */
//}
void draw_aoq_paddle(Double_t zet){
  TH2F* h_aoq_paddle = new TH2F(Form("h_aoq_paddle%i",(int)zet),Form("h_aoq_paddle%i",(int)zet),NUMPADDLE,0.5,NUMPADDLE,500,min_aoq,max_aoq);
  nentry = ch->GetEntries();
  for(Long64_t n=0; n<nentry; n++){
    ch->GetEntry(n);
    Int_t j = Tofw_Paddle-1;
    Double_t beta_fit = f_betatof[j]->GetParameter(0)/(FragTof + f_betatof[j]->GetParameter(1)) + f_betatof[j]->GetParameter(2);
    Double_t aoq_fit = (FragBrho * TMath::Sqrt(1.-beta_fit*beta_fit))/(beta_fit * mc_e);
    //FragBeta = beta_fit;
    //FragGamma = 1./TMath::Sqrt(1.-beta_fit*beta_fit);
    h_aoq_paddle->Fill(Tofw_Paddle,aoq_fit);
  }
}
///
void filltree(){
  nentry = ch->GetEntries();
  cout<<"FILLTREE"<<endl;
  for(Long64_t i=0; i<nentry; i++){
    //for(Long64_t i=0; i<10000; i++){
    ch->GetEntry(i);
    //if((Mw0_X + 405.9 * MusicTheta < 0)||(Mw0_X + 405.9 * MusicTheta > 15)) continue;
    //
    ROLU_X = Mw0_X + 405.9 * MusicTheta;
    //
    if(Tofw_Paddle - 1 < MINPADDLE || Tofw_Paddle > NUMPADDLE){
      Tofw_Paddle = -1;
      //FragAoQ = NAN;
      FragAoQ_corr = NAN;
      FragBeta = NAN;
      FragGamma = NAN;
      FragTof = NAN;
    }else{
      //cout<<(int)Tofw_Paddle<<endl;
      //    f_betatof[i] = new TF1(Form("f_betatof%i",i+1), "[0]/(x + [1])+[2]",40,600);
      Int_t j = Tofw_Paddle-1;
      Double_t beta_fit = f_betatof[j]->GetParameter(0)/(FragTof + f_betatof[j]->GetParameter(1)) + f_betatof[j]->GetParameter(2);
      Double_t aoq_fit = (FragBrho * TMath::Sqrt(1.-beta_fit*beta_fit))/(beta_fit * mc_e);
      FragBeta = beta_fit;
      FragGamma = 1./TMath::Sqrt(1.-beta_fit*beta_fit);
      FragAoQ_corr = aoq_fit;
      FragTof += f_betatof[j]->GetParameter(1);
	/*
	FragAoQ + fit_prof[Tofw_Paddle-1]->GetParameter(0)
	+ fit_prof[Tofw_Paddle-1]->GetParameter(1) * TwimTheta * 1000.;
	//+ fit_prof[Tofw_Paddle-1]->GetParameter(2) * pow(TwimTheta*1000.,2.); // pol2
	*/
    }
    tree->Fill();
    if(i%100000==0)
      cout<<"\r"<< i<<" entries done in "<<nentry<<flush;//endl;//
  }
  cout<<endl;
  tree->Write();
}
///
void writecsv(){
  //  fcsv<<"test"<<endl;
  fcsv << "FRS Z, FRS A, Frag Z, Frag A";
  for(int i = MINPADDLE; i<NUMPADDLE; i++) fcsv<<", Paddle "<<i+1;
  fcsv << ", total, total3, total_ratio 2sigma-3sigma"<<endl;
  for(int FZ = MINZ; FZ< MAXZ; FZ++){
    for(int FA= MINA; FA< MAXA; FA++){
      int total=0, total3=0;
      if(h_counts_mw3x[FA][FZ][0][0]->GetEntries()==0) continue;
      //
      fcsv<<FZ<<", "<<FA<<", "<<"0"<<", 0";
      for(int i = MINPADDLE; i<NUMPADDLE; i++){
	fcsv<<", "<<NumGated[i][FA][FZ][0][0];
	total +=NumGated[i][FA][FZ][0][0];
	//total3 +=NumGated3[i][FA][FZ][0][0];
      }// no gate in fragment
      fcsv<<", "<<total<<",0 ,0"<<endl; //total3<<", "<<(double)total/(double)total3<<endl;
      //
      for(int CZ = MINZ; CZ<=FZ; CZ++){
	for(int CA= MINA; CA<=FA; CA++){
	  total = 0;
	  total3= 0;
	  fcsv<<FZ<<", "<<FA<<", "<<CZ<<", "<<CA;
	  for(int i = MINPADDLE; i<NUMPADDLE; i++){
	    fcsv<<", "<<NumGated[i][FA][FZ][CA][CZ];
	    total +=NumGated[i][FA][FZ][CA][CZ];
	    total3 +=NumGated3[i][FA][FZ][CA][CZ];
	  }
	  Double_t total_fit = 0;
	  Int_t DZ = FZ-CZ, DN = (FA-CA)-(FZ-CZ);
	  if(DN<3 && DZ<3 && DN>=0)
	    total_fit = f_FitGated[FA][FZ][DN][DZ]->GetParameter(0)/10.; // bin width =10mm
	  //fcsv<<", "<<total<<", "<<(double)total/(double)total3<<", "<<total_fit<<endl;
	  fcsv<<", "<<total<<", "<<total3<<", "<<(double)total/(double)total3<<endl;
	}
      }
    }
  }
}
///
void writeparam(){
  //FairRunOnline* run = new FairRunOnline(source);
  FairRunOnline* run = new FairRunOnline();
  //FairRun* run = new FairRun();
  run->SetRunId(1);

  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  //FairRuntimeDb* rtdb = FairRuntimeDb::instance() ;
  if (!rtdb)
    {
      return;
    }
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii
  FairParAsciiFileIo* parOut = new FairParAsciiFileIo();
  parIo1->open(inpar, "in");
  //
  rtdb->setFirstInput(parIo1);
  rtdb->print();
  rtdb->getContainer("SofSciRawPosPar");
  rtdb->getContainer("tofwHitPar");
  auto fFragmentPar = (R3BSofFragmentAnaPar*)rtdb->getContainer("soffragmentAnaPar");
  cout<<"Parameter file: "<<inpar<<" is open."<<endl;
  rtdb->initContainers(1);
  cout<<"Old FragmentPar:"<<endl;
  fFragmentPar->printParams();
  //
  for(int i = 0; i<fFragmentPar->GetNumSci(); i++){
    int Paddle = i+1;
    if (fFragmentPar->GetInUse(Paddle) != 1) continue;
    fFragmentPar->SetEffectivLength(f_betatof[i]->GetParameter(0), Paddle);
    fFragmentPar->SetTofWOffset(-1. * f_betatof[i]->GetParameter(1), Paddle); // sign should be inverted
  }
  cout<<endl<<"New FragmentPar:"<<endl;
  fFragmentPar->printParams();
  //
  parOut->open(outpar,"out");
  rtdb->setOutput(parOut);
  rtdb->saveOutput();
}
///
Int_t initbranch(){
  ch->SetBranchAddress("FRSAoQ",&FRSAoQ);
  ch->SetBranchAddress("FRSBrho",&FRSBrho);
  ch->SetBranchAddress("FRSBeta",&FRSBeta);
  ch->SetBranchAddress("FRSGamma",&FRSGamma);
  ch->SetBranchAddress("MusicE", &MusicE);
  ch->SetBranchAddress("MusicZ", &MusicZ);
  ch->SetBranchAddress("MusicTheta", &MusicTheta);
  ch->SetBranchAddress("TwimE", &TwimE);
  ch->SetBranchAddress("FragZ", &FragZ);
  ch->SetBranchAddress("FragBrho",&FragBrho);
  ch->SetBranchAddress("TwimTheta", &TwimTheta);
  ch->SetBranchAddress("Mw0_X", &Mw0_X);
  ch->SetBranchAddress("Mw1_X", &Mw1_X);
  ch->SetBranchAddress("Mw2_X", &Mw2_X);
  ch->SetBranchAddress("Mw3_X", &Mw3_X);
  ch->SetBranchAddress("Mw0_Y", &Mw0_Y);
  ch->SetBranchAddress("Mw1_Y", &Mw1_Y);
  ch->SetBranchAddress("Mw2_Y", &Mw2_Y);
  ch->SetBranchAddress("Mw3_Y", &Mw3_Y);
  ch->SetBranchAddress("Tofw_Y", &Tofw_Y);
  ch->SetBranchAddress("FragTof", &FragTof);
  ch->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  ch->SetBranchAddress("FragAoQ", &FragAoQ);
  //
  tree = new TTree("FragTree", "Fragment tree");
  tree->Branch("FRS_AoQ",&FRSAoQ);
  tree->Branch("FRS_Brho",&FRSBrho);
  tree->Branch("FRSBeta",&FRSBeta);
  tree->Branch("FRSGamma",&FRSGamma);
  tree->Branch("FRS_Z", &MusicZ);
  tree->Branch("Frag_AoQ", &FragAoQ_corr);
  tree->Branch("Frag_AoQ_original", &FragAoQ);
  tree->Branch("Frag_Brho",&FragBrho);
  tree->Branch("Frag_Z", &FragZ);
  tree->Branch("Frag_Tof", &FragTof);
  tree->Branch("Frag_Beta", &FragBeta);
  tree->Branch("Frag_Gamma", &FragGamma);
  tree->Branch("MusicE", &MusicE);
  tree->Branch("MusicTheta", &MusicTheta);
  tree->Branch("TwimE", &TwimE);
  tree->Branch("TwimTheta", &TwimTheta);
  tree->Branch("Mw0_X", &Mw0_X);
  tree->Branch("Mw1_X", &Mw1_X);
  tree->Branch("Mw2_X", &Mw2_X);
  tree->Branch("Mw3_X", &Mw3_X);
  tree->Branch("Mw0_Y", &Mw0_Y);
  tree->Branch("Mw1_Y", &Mw1_Y);
  tree->Branch("Mw2_Y", &Mw2_Y);
  tree->Branch("Mw3_Y", &Mw3_Y);
  tree->Branch("ROLU_X", &ROLU_X);
  tree->Branch("Tofw_Y", &Tofw_Y);
  tree->Branch("Tofw_Paddle", &Tofw_Paddle);
  //
  return 0;
}


void initialise(){
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleFont(42,"T");

  /*
  gStyle->SetLabelSize(25,"XYZ");
  gStyle->SetTitleSize(25,"XYZ");
  */
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");

  gStyle->SetTitleOffset(.9,"X");
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetLabelOffset(0.01,"X");
  gStyle->SetLabelOffset(0.01,"Y");
 

  ch = new TChain("Tree");
  ch -> Add(infile);
  //ch -> Add(infile);
  //ch -> SetProof();
  //
  //frspidgate = Form("abs(MusicZ-%i)<0.4 && abs(FRSAoQ -%f)<%f",zetdiff, (double)mass/(double)zetdiff, 0.22/(double)zetdiff);
  //
  brho = "FRSBrho";
  beta_tofw_mod = "(FragBeta)";
  fragaoqstring = "FragAoQ";
  //cut_mw = "Mw0_X > 0";
  cut_mw = "(Mw0_X + 405.9 * MusicTheta > 0) && (Mw0_X + 405.9 * MusicTheta < 15)"; // ROLU condition; // not used in the later analysis
  Zgate=Form("abs(MusicZ-FragZ-%f)<0.4",(float)zetdiff);
  //
    for(int i=MINPADDLE; i<NUMPADDLE; i++){
    h_noreaction[i] = new TH1F(Form("h_noreaction_%i",i+1),Form("AoQ difference of Paddle%i",i+1),500,-0.6,0.6);
    h_frsbeta_betacalc[i] = new TH2F(Form("h_frsbeta_betacalc%i",i+1), Form("h_frsbeta_betacalc%i",i+1), 500, 0.60, 0.85,500, 0.60, 0.85);
    h_deltabeta_tof[i] = new TH2F(Form("h_deltabeta_tof%i",i+1), Form("#beta vs raw ToF with limieted theta of Paddle %i", i+1), 1000, 40, 60, 500, 0.60, 0.85);
    f_betatof[i] = new TF1(Form("f_betatof%i",i+1), "[0]/(x + [1])+[2]",40,600);
    h_deltabeta_twime[i] = new TH2F(Form("h_deltabeta_twime%i;Twim E / ch; #beta deviation in percent",i+1), Form("#Energy_{Twim} vs #beta * ToF of Paddle %i", i+1), 500, 0, 10000, 500, -2, 2);
    h_deltabeta_theta[i] = new TH2F(Form("h_deltabeta_theta%i",i+1), Form("#Theta_{Twim} vs #beta * ToF of Paddle %i", i+1), 500, -22, 22, 500, -0.02, 0.02);
    //f_betatheta[i] = new TF1(Form("f_betatheta%i",i+1), "pol1", -0.5, 0.5);
    h_deltabeta_pid[i] = new TH2F(Form("h_deltabeta_pid%i",i+1), Form("PID of Paddle  %i", i+1), 500,min_aoq,max_aoq,500,10,25);
    for (int A = MINA; A < MAXA; A++){
      for (int Z = MINZ; Z < MAXZ; Z++){
	h_gatedpid[i][A][Z] = new TH2F(Form("h_gatedpid%i_%i_%i",i,A,Z),Form("PID of Paddle%i gated in FRS, A=%i, Z=%i",i+1,A,Z),
				       500, min_aoq, max_aoq, 500, 10, 25);
	if(A==MINA){//only once
	  h_zgatedpid[i][Z] = new TH1F(Form("h_zgatedpid%i_%i",i+1,Z),Form("PID of Paddle%i gated with frag Z=%i",i+1,Z),
					 500, min_aoq, max_aoq);
	  if(i==MINPADDLE)
	    h_zgated_mw3x_aoq[Z] = new TH2F(Form("h_zgated_mw3x_aoq_%i",Z), Form("h_zgated_mw3x_aoq_%i",Z), 500, -500,500,500, min_aoq, max_aoq);
	}
	double tmpaoq = (double)A/(double)Z;
	f_fragfit[i][A][Z] = new TF2(Form("f_fragfit%i_%i_%i",i,A,Z),"[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",
				     tmpaoq-0.025,tmpaoq+0.025,(double)Z-0.4,(double)Z+0.4);
	f_fragfit[i][A][Z] -> SetParameter(0,0);
	f_fragfit[i][A][Z] -> SetParameter(1,tmpaoq);
      	f_fragfit[i][A][Z] -> SetParameter(2,1./200.);
      	f_fragfit[i][A][Z] -> SetParameter(3,(double)Z);
      	f_fragfit[i][A][Z] -> SetParameter(4,0.2);
	//
      	f_fragfit[i][A][Z] -> SetParLimits(0,0,1e8);
      	f_fragfit[i][A][Z] -> SetParLimits(1,tmpaoq-0.08,tmpaoq+0.008);
      	f_fragfit[i][A][Z] -> SetParLimits(2,0.002,0.01);
      	f_fragfit[i][A][Z] -> SetParLimits(3,(double)Z-0.1,(double)Z+0.1);
	f_fragfit[i][A][Z] -> SetParLimits(4,0.03,0.15);
	//
	if(i>MINPADDLE) continue; // call below only once
	for(int DN=0; DN<3; DN++){ // Delta N
	  for(int DZ=0; DZ<3; DZ++){ // Delta Z
	    h_counts_paddle[A][Z][DN][DZ] = new TH1F(Form("h_counts_paddle_%i_%i_%i_%i",A,Z,DN,DZ),
						     Form("Counts for each paddle from A=%i Z=%i with #DeltaN=%i #DeltaZ=%i reaction",A,Z,DN,DZ),
						     NUMPADDLE, 0.5, NUMPADDLE+0.5);
	    h_counts_mw3x[A][Z][DN][DZ] = new TH1F(Form("h_counts_mw3x_%i_%i_%i_%i",A,Z,DN,DZ),
						   Form("Distribution of hits in X for the reaction from A=%i Z=%i with #DeltaN=%i #DeltaZ=%i reaction",A,Z,DN,DZ),
						   60,-300,300);
	    h_counts_mw3x3[A][Z][DN][DZ] = new TH1F(Form("h_counts_mw3x3_%i_%i_%i_%i",A,Z,DN,DZ),
						   Form("Distribution of hits in X for the reaction from A=%i Z=%i with #DeltaN=%i #DeltaZ=%i reaction with larger pid gate",A,Z,DN,DZ),
						   60,-300,300);
	    h_counts_mw3x3[A][Z][DN][DZ]->SetLineColor(kBlue);
	  }
	}
	h_counts_paddle[A][Z][0][0]->SetLineColor(1);
      	h_counts_paddle[A][Z][1][0]->SetLineColor(2);
      	h_counts_paddle[A][Z][0][1]->SetLineColor(3);
      }
    }
  }
}
