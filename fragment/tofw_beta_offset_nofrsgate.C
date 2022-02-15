#define MINPADDLE 0
#define NUMPADDLE 28

#include "fragmentana_twim.h"
int zetdiff = 0;
double min_aoq = 1.8, max_aoq = 2.7;
TString FRS="50Ca";
//
//for trees
Float_t FRSAoQ, FRSBeta, FRSBrho, MusicE, TwimE, MusicZ, FragZ, TwimTheta, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
UChar_t Tofw_Paddle;
Long64_t nentry=0;
TString targ, infile;
//
/*
TString targ = "empty";
TString infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_empty.root";
*/
//
/*
TString targ = "ch2";
TString infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_ch2-24mm.root";
*/
//
/*
TString targ = "carbon";
TString infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_carbon.root";
*/
//
/*
TString targ = "PP";
TString infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+ FRS +"_PP.root";
*/
TString outpdf ;
TString outcsv ;
TString outroot;
//
bool IsEmpty=false;
const Double_t sigma_frs = 4., sigma = 2., sigma3 = 3.;//Set conditions for PID gating
const Int_t MINA=36, MAXA=55, MINZ=16, MAXZ=22;
//
ofstream fcsv;
TFile *fout;
TF1 *fit_prof[NUMPADDLE];
TH2F *h_frsbeta_betacalc[NUMPADDLE], *h_deltabeta_tof[NUMPADDLE], *h_deltabeta_twime[NUMPADDLE], *h_deltabeta_theta[NUMPADDLE], *h_deltabeta_pid[NUMPADDLE];
TH2F *h_gatedpid[NUMPADDLE][100][50];
Int_t NumGated[NUMPADDLE][100][50][100][50] = {0}, NumGated3[NUMPADDLE][100][50][100][50] = {0}; // Count hits in 2sigma and 3sigma
TH1F *h_counts_paddle[100][50][3][3], *h_counts_mw3x[100][50][3][3], *h_counts_mw3x3[100][50][3][3];
TF1 *f_FitGated[100][50][3][3];
TF2 *f_fragfit[NUMPADDLE][100][50];
TF1 *f_aoq[NUMPADDLE], *f_betatof[NUMPADDLE], *f_betatheta[NUMPADDLE];
TTree *tree;

void tofw_beta_offset_nofrsgate(int i_target=1);
void define_conditions();
void delta_beta_method();
void initialise();
void filltree();
Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag=false, Int_t i=0);
Int_t initbranch();
void writecsv();

//////////
void tofw_beta_offset_nofrsgate(int i_target){
  switch (i_target){
  case 0:
    targ = "empty";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_empty.root";
    break;
  case 1:
    targ = "ch2";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_ch2-24mm.root";
    break;
  case 2:
    targ = "carbon";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+FRS+"_carbon.root";
    break;
  case 3:
    targ = "PP";
    infile = "./rootfiles/rootfile_land/mktree/mktree_fragment_"+ FRS +"_PP.root";
    break;
  default:
    cerr<<"No target info"<<endl;
    return;
  }
  outpdf = "./fragment/output/tofw_beta_offset_paddle_" + FRS + "_"+ targ + Form("_nofrsgate_Feb") + ".pdf";
  outcsv = "./fragment/output/tofw_beta_offset_paddle_" + FRS + "_"+ targ + Form("_nofrsgate_Feb") + ".csv";
  outroot = "./fragment/output/tofw_beta_offset_paddle_" + FRS + "_"+ targ + Form("_nofrsgate_Feb") + ".root";
  //
  fcsv.open(outcsv, ofstream::out);
  fout = new TFile(outroot,"RECREATE");
  //TString recobrho_outpdf = "./fragment/output/tofw_beta_offset_paddlebrho_" + targ + Form("_nofrsgate_zdiff%i",zetdiff) + ".pdf";


  cout<< targ<< " "<<infile<<endl;

  initialise();
  c = new TCanvas();
  c->Divide(4,4);
  c->Print(outpdf + "[");
  //
  define_conditions(); // tentatively commented
  //
  initbranch();
  delta_beta_method();
  //
  c->Print(outpdf + "]");
  //
  filltree();
  writecsv();
}

void define_conditions(){
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h_aoqdiff_twim%i",i),Form("Correlation in AoQ and X-angle of Paddle%i;AoQ_{FRS} - AoQ_{Cave};Twim Theta (mrdad)",i+1),500,-0.6,0.6,500,-55,55);
    ch->Draw(Form("TwimTheta*1000.:(FRSAoQ-FragAoQ)>>h_aoqdiff_twim%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  /*
  // Angle change
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h2_%i",i),Form("Correlation in AoQ and X-angle change of Paddle%i;AoQ_{FRS} - AoQ_{Cave};#Theta_{Twim} - #Theta_{R3BMusic} (mrad)",i+1),500,-0.6,0.6,500,-55,55);
    ch->Draw(Form("(TwimTheta-MusicTheta)*1000.:(FRSAoQ-FragAoQ)>>h2_%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  / */
  // AoQ Zet
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h3_%i",i),Form("Z vs AoQ difference of Paddle%i;AoQ_{FRS} - AoQ_{Cave};Twim Z",i+1),500,-0.6,0.6,500,15,25);
    ch->Draw(Form("FragZ:(FRSAoQ-FragAoQ)>>h3_%i",i),Form("Tofw_Paddle==%i ",i+1),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Theta Tof
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h4_%i",i),Form("Tof and Twim Theta of Paddle%i; ToF in Cave (ns); #Theta_{Twim} (mrad)",i+1),500,40,60,500,-55,55);
    ch->Draw(Form("(TwimTheta)*1000.:(FragTof)>>h4_%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  
  // Unreacted cut
  //
  Double_t peakX[NUMPADDLE]={1000.};
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH1F(Form("h_noreaction_%i",i),Form("AoQ difference of Paddle%i",i+1),500,-0.6,0.6);
    ch->Draw(Form("(FRSAoQ-FragAoQ)>>h_noreaction_%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    TSpectrum *tspec = new TSpectrum(10);
    Int_t npeaks = tspec->Search(h);
    Double_t *peakXspec = tspec->GetPositionX();
    peakX[i]=2000.;
    for(int j = 0; j<npeaks; j++){
      //peakX[i] = (abs(peakX[i]) < abs(peakXspec[j]))? peakX[i] : peakXspec[j];
      if(TMath::Abs(peakX[i]) > TMath::Abs(peakXspec[j])) peakX[i] = peakXspec[j];
      cout<<peakXspec[j]<<" "<<peakX[i]<<endl;
    }      
    //
    f_aoq[i] = new TF1(Form("f_aoqdiff%i",i),"gaus");
    cout<<"peakX"<<i+1<<": "<<peakX[i]<<endl;
    f_aoq[i]->SetParameter(0,100);
    f_aoq[i]->FixParameter(1,peakX[i]);
    f_aoq[i]->SetParameter(2,0.01);
    h->Fit(f_aoq[i],"L","",peakX[i]-0.03,peakX[i]+0.03);
    //
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Theta Tof (gated with aoq)
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h4_gated_%i",i),Form("Tof and Twim Theta with aoq gated of Paddle%i",i+1),500,40,60,500,-55,55);
    TCut cut = Form("abs(FRSAoQ-FragAoQ-%f)<2.*%f && abs(TwimTheta-6e-3)*1000.<14.", f_aoq[i]->GetParameter(1), f_aoq[i]->GetParameter(2)); // cut with +/- 20 mrad
    ch->Draw(Form("(TwimTheta)*1000.:(FragTof)>>h4_gated_%i",i),cut && Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );

}

void delta_beta_method(){
  for(int i=MINPADDLE; i<NUMPADDLE; i++){
    h_frsbeta_betacalc[i] = new TH2F(Form("h_frsbeta_betacalc%i",i+1), Form("h_frsbeta_betacalc%i",i+1), 500, 0.60, 0.85,500, 0.60, 0.85);
    h_deltabeta_tof[i] = new TH2F(Form("h_deltabeta_tof%i",i+1), Form("#beta vs raw ToF with limieted theta of Paddle %i", i+1), 1000, 40, 60, 500, 0.60, 0.85);
    f_betatof[i] = new TF1(Form("f_betatof%i",i+1), "[0]/(x + [1])+[2]",40,600);
    h_deltabeta_twime[i] = new TH2F(Form("h_deltabeta_twime%i;Twim E / ch; #beta deviation in percent",i+1), Form("#Energy_{Twim} vs #beta * ToF of Paddle %i", i+1), 500, 0, 10000, 500, -2, 2);
    h_deltabeta_theta[i] = new TH2F(Form("h_deltabeta_theta%i",i+1), Form("#Theta_{Twim} vs #beta * ToF of Paddle %i", i+1), 500, -22, 22, 500, -0.02, 0.02);
    //f_betatheta[i] = new TF1(Form("f_betatheta%i",i+1), "[0]*x", -14, 14);
    f_betatheta[i] = new TF1(Form("f_betatheta%i",i+1), "pol1", -0.5, 0.5);
    h_deltabeta_pid[i] = new TH2F(Form("h_deltabeta_pid%i",i+1), Form("PID of Paddle  %i", i+1), 500,min_aoq,max_aoq,500,10,25);
    for (int A = MINA; A < MAXA; A++){
      for (int Z = MINZ; Z < MAXZ; Z++){
	h_gatedpid[i][A][Z] = new TH2F(Form("h_gatedpid%i_%i_%i",i,A,Z),Form("PID of Paddle%i gated in FRS, A=%i, Z=%i",i+1,A,Z),
				       500, min_aoq, max_aoq, 500, 10, 25);
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
  //
  nentry = ch->GetEntries();
  cout<<"FILLTREE"<<endl;
  Long64_t neve = 0;
  for(Long64_t n=0; n<nentry; n++){
  //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    //cout<<i<<": "<<MusicZ<<" "<<FRSAoQ<<" "<<tmpzet<<" "<<tmpaoq<<endl;
    if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2)  || TMath::Abs(FragZ-tmpzet)>0.4 ||
       TMath::Abs(TwimTheta-6e-3)*1000. > 14. || TMath::Abs(FragTof-40-0.4*(Double_t)i)>2.) continue;
    //if(TMath::Abs(TwimTheta * 1000.) > 5) continue;
    Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    h_deltabeta_tof[i]->Fill(FragTof, betacalc);
    h_frsbeta_betacalc[i]->Fill(FRSBeta, betacalc);
    //
    if(++neve%100000==0)
      cout<<"\r"<<i<<" enrties done in "<<nentry<<flush;//endl;//
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
    //TProfile *prof = h_deltabeta_tof[i]->ProfileX();
    //prof->Draw("same");
    //prof
    f_betatof[i]->FixParameter(2,0);
    h_deltabeta_tof[i]->Fit(Form("f_betatof%i",i+1),"","",h_deltabeta_tof[i]->GetMean() - 2.6* h_deltabeta_tof[i]->GetStdDev(), h_deltabeta_tof[i]->GetMean() + 2.6* h_deltabeta_tof[i]->GetStdDev());
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_frsbeta_betacalc[i]->Draw("colz");
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //////////////////////////////
  ////deltabeta_TwimE  ////////////////////////////
  //////////////////////////////
  for(Long64_t n=0; n<nentry; n++){
  //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    //if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2) || TMath::Abs(FragZ-tmpzet)>0.4 || TMath::Abs(FragTof-40-0.4*(Double_t)i)>2.) continue; //||
    //if(TMath::Abs(44./18.-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.5*f_aoq[i]->GetParameter(2) || TMath::Abs(FragZ-tmpzet)>0.4 || TMath::Abs(FragZ-18.)>0.4 || TMath::Abs(Mw3_Y)>10.) continue; //||select only one isotope
    //if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.5*f_aoq[i]->GetParameter(2) || TMath::Abs(FragZ-tmpzet)>0.4 || TMath::Abs(FragZ-18.)>0.4 || TMath::Abs(Mw3_Y)>10.) continue; //||select only one isotope
    if(TMath::Abs(2.5-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.5*f_aoq[i]->GetParameter(2) || TMath::Abs(FragZ-tmpzet)>0.4 || ((int)(tmpzet/2))*2!=((int)tmpzet)) continue; //||select only one isotope
    Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    //for(int j=0;j<3;j++) cout<<f_betatof[i]->GetParameter(j)<<" ";
    //cout<<endl;
    h_deltabeta_twime[i]->Fill(TwimE,
			       - 100. * (betacalc) * (FragTof + f_betatof[i]->GetParameter(1)) / f_betatof[i]->GetParameter(0) + 100.);
    //
    if(++neve%100000==0)
      cout<<"\r"<<n<<" enrties done in "<<nentry<<flush;//endl;//
  }
  cout<<endl;
  //Draw
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_deltabeta_twime[i]->Draw("colz");
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_deltabeta_twime[i]->Draw("colz");
    //TProfile *prof = h_deltabeta_twime[i]->ProfileX();
    //prof->Draw("same");
    //prof
    h_deltabeta_twime[i]->Fit(Form("f_betatheta%i",i+1),"","",h_deltabeta_twime[i]->GetMean() - 1.6* h_deltabeta_twime[i]->GetStdDev(), h_deltabeta_twime[i]->GetMean() + 1.6* h_deltabeta_twime[i]->GetStdDev());
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  ////////////////////////////
  ////deltabeta_theta  ////////////////////////////
  ////////////////////////////
  for(Long64_t n=0; n<nentry; n++){
  //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    Double_t tmpzet = MusicZ, tmpaoq = FRSAoQ;
    if(pid(tmpzet, tmpaoq, false)>sigma_frs) continue;
    if(TMath::Abs(FRSAoQ-FragAoQ-f_aoq[i]->GetParameter(1)) > 2.*f_aoq[i]->GetParameter(2) || TMath::Abs(FragZ-tmpzet)>0.4 || TMath::Abs(FragTof-40-0.4*(Double_t)i)>2.) continue; //||
       //TMath::Abs(TwimTheta-6e-3)*1000. > 14. ) continue;
    //if(TMath::Abs(Mw3_Y - Mw1_Y -20.)>5) continue; // Y cut
    Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    //for(int j=0;j<3;j++) cout<<f_betatof[i]->GetParameter(j)<<" ";
    //cout<<endl;
    h_deltabeta_theta[i]->Fill(TwimTheta  *1000.,
    //h_deltabeta_theta[i]->Fill(Mw3_Y - Mw1_Y -20.,
			       //(betacalc - f_betatof[i]->GetParameter(2)) * (FragTof + f_betatof[i]->GetParameter(1)) - f_betatof[i]->GetParameter(0));
			       (betacalc) * (FragTof + f_betatof[i]->GetParameter(1)) / f_betatof[i]->GetParameter(0) - 1.);
    //
    if(++neve%100000==0)
      cout<<"\r"<<n<<" enrties done in "<<nentry<<flush;//endl;//
  }
  cout<<endl;
  //Draw
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_deltabeta_theta[i]->Draw("colz");
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    h_deltabeta_theta[i]->Draw("colz");
    //TProfile *prof = h_deltabeta_theta[i]->ProfileX();
    //prof->Draw("same");
    //prof
    h_deltabeta_theta[i]->Fit(Form("f_betatheta%i",i+1),"","",h_deltabeta_theta[i]->GetMean() - 1.6* h_deltabeta_theta[i]->GetStdDev(), h_deltabeta_theta[i]->GetMean() + 1.6* h_deltabeta_theta[i]->GetStdDev());
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //////////////
  /*
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    TProfile *prof = h_deltabeta_theta[i]->ProfileY();
    prof->Draw();
    prof->Fit("pol1","","",h_deltabeta_theta[i]->GetMean(2) - 1.6* h_deltabeta_theta[i]->GetStdDev(2), h_deltabeta_theta[i]->GetMean(2) + 1.6* h_deltabeta_theta[i]->GetStdDev(2));
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  */
  // Draw PID
  //
  for(Long64_t n=0; n<nentry; n++){
  //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
    Int_t i = Tofw_Paddle - 1;
    if(i < MINPADDLE || i >= NUMPADDLE) continue;
    //Double_t betacalc = TMath::Sqrt((FragBrho*FragBrho) / (pow(tmpaoq*mc_e,2) + pow(FragBrho,2)));
    Double_t beta_fit = f_betatof[i]->GetParameter(0)*(1. + 1000.*TwimTheta*f_betatheta[i]->GetParameter(1))/(FragTof + f_betatof[i]->GetParameter(1)) + f_betatof[i]->GetParameter(2);
    Double_t aoq_fit = (FragBrho * TMath::Sqrt(1.-beta_fit*beta_fit))/(beta_fit * mc_e);
    //cout<<i << " "<<beta_fit<< " "<<aoq_fit<<" "<<FragZ<<endl;
    h_deltabeta_pid[i]->Fill(aoq_fit,FragZ);
    //
    if(++neve%100000==0)
      cout<<"\r"<<n<<" enrties done in "<<nentry<<flush;//endl;//
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
    c->Print(outpdf);
    //if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  // c->Print(outpdf);

  // Draw gated pid
  c->Clear();
  c->cd(0);
  c->Divide(4,4);
  for(Long64_t n=0; n<nentry; n++){
    //for(Long64_t n=0; n<100; n++){
    ch->GetEntry(n);
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
    /*
    else{
      NumGated[i][tmpA][tmpZ][0][0]++; // outside of the gate
      }*/
    //
    if(++neve%100000==0)
      cout<<"\r"<<n<<" enrties done in "<<nentry<<flush;//endl;//
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
  }
}

Double_t pid(Double_t &zet, Double_t &aoq, bool isfrag=false, Int_t i){
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

void filltree(){
  nentry = ch->GetEntries();
  cout<<"FILLTREE"<<endl;
  for(Long64_t i=0; i<nentry; i++){
    //for(Long64_t i=0; i<10000; i++){
    ch->GetEntry(i);
    if(Tofw_Paddle - 1 < MINPADDLE || Tofw_Paddle > NUMPADDLE){
      Tofw_Paddle = -1;
      //FragAoQ = NAN;
      FragAoQ_corr = NAN;
    }else{
      //cout<<(int)Tofw_Paddle<<endl;
      //    f_betatof[i] = new TF1(Form("f_betatof%i",i+1), "[0]/(x + [1])+[2]",40,600);
      Int_t j = Tofw_Paddle-1;
      Double_t beta_fit = f_betatof[j]->GetParameter(0)*(1. + 1000.*TwimTheta*f_betatheta[j]->GetParameter(1))/(FragTof + f_betatof[j]->GetParameter(1)) + f_betatof[j]->GetParameter(2);
      Double_t aoq_fit = (FragBrho * TMath::Sqrt(1.-beta_fit*beta_fit))/(beta_fit * mc_e);
      FragAoQ_corr = aoq_fit;
	/*
	FragAoQ + fit_prof[Tofw_Paddle-1]->GetParameter(0)
	+ fit_prof[Tofw_Paddle-1]->GetParameter(1) * TwimTheta * 1000.;
	//+ fit_prof[Tofw_Paddle-1]->GetParameter(2) * pow(TwimTheta*1000.,2.); // pol2
	*/
    }
    tree->Fill();
    if(i%100000==0)
      cout<<"\r"<< i<<" enrties done in "<<nentry<<flush;//endl;//
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
Int_t initbranch(){
  ch->SetBranchAddress("FRSAoQ",&FRSAoQ);
  ch->SetBranchAddress("FRSBrho",&FRSBrho);
  ch->SetBranchAddress("FRSBeta",&FRSBeta);
  ch->SetBranchAddress("MusicE", &MusicE);
  ch->SetBranchAddress("MusicZ", &MusicZ);
  ch->SetBranchAddress("TwimE", &TwimE);
  ch->SetBranchAddress("FragZ", &FragZ);
  ch->SetBranchAddress("FragBrho",&FragBrho);
  ch->SetBranchAddress("TwimTheta", &TwimTheta);
  ch->SetBranchAddress("Mw1_X", &Mw1_X);
  ch->SetBranchAddress("Mw2_X", &Mw2_X);
  ch->SetBranchAddress("Mw3_X", &Mw3_X);
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
  tree->Branch("FRS_Z", &MusicZ);
  tree->Branch("Frag_AoQ", &FragAoQ_corr);
  tree->Branch("Frag_AoQ_original", &FragAoQ);
  tree->Branch("Frag_Brho",&FragBrho);
  tree->Branch("Frag_Z", &FragZ);
  tree->Branch("Frag_Tof", &FragTof);
  tree->Branch("MusicE", &MusicE);
  tree->Branch("TwimE", &TwimE);
  tree->Branch("TwimTheta", &TwimTheta);
  tree->Branch("Mw1_X", &Mw1_X);
  tree->Branch("Mw2_X", &Mw2_X);
  tree->Branch("Mw3_X", &Mw3_X);
  tree->Branch("Mw1_Y", &Mw1_Y);
  tree->Branch("Mw2_Y", &Mw2_Y);
  tree->Branch("Mw3_Y", &Mw3_Y);
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
  cut_mw = Form("(Mw1_X>%f)&&(Mw1_X<%f)&&", -30.,25.);
  for(int i=0;i<(IsEmpty?4:3);i++) cut_mw += Form("(%s>%f)&&(%s<%f)&&",axis_mw3[i].Data(),range_cut_mw3_low[i],axis_mw3[i].Data(),range_cut_mw3_high[i]);
  Zgate=Form("abs(MusicZ-FragZ-%f)<0.4",(float)zetdiff);
}
