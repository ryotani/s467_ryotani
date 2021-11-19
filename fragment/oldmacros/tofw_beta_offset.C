#define MINPADDLE 0
#define NUMPADDLE 28

#include "fragmentana_twim.h"
int zetfrs = 18, zetdiff = 0;//, mass = 44;
double min_aoq = 1.8, max_aoq = 2.7;
//
/*
TString targ = "empty";
TString infile = "./fragment/output/mktree_fragment_Nov_empty.root";
double beta_offset = 0.0;
*/
//

TString targ = "ch2";
TString infile = "./fragment/output/mktree_fragment_Nov_ch2-24mm.root";
//TString infile = "./rootfiles/rootfiletmp/fragment_Sep2021/s467_filltree_Setting13_0354_28Sep.root";//1hr CH2 run
//double beta_offset = 0.006;
double beta_offset = 0.000;

//
/*
TString targ = "carbon";
TString infile = "./fragment/output/mktree_fragment_Nov_carbon.root";
double beta_offset = 0.000;
*/
//
/*
TString targ = "PP";
TString infile = "./fragment/output/mktree_fragment_Nov_PP.root";
double beta_offset = 0.000;
*/
//
TString outpdf = "./fragment/output/tofw_beta_offset_" + targ + Form("_zfrs%i_zdiff%i",zetfrs,zetdiff) + "_beta_offset_paddle.pdf";
TString outcsv = "./fragment/output/tofw_beta_offset_" + targ + Form("_zfrs%i_zdiff%i",zetfrs,zetdiff) + "_beta_offset_paddle.csv";
TString outroot = "./fragment/output/tofw_beta_offset_" + targ + Form("_zfrs%i_zdiff%i",zetfrs,zetdiff) + "_beta_offset_paddle.root";
TString recobrho_outpdf = "./fragment/output/tofw_beta_offset_brho_" + targ + Form("_zfrs%i_zdiff%i",zetfrs,zetdiff) + "_beta_offset_paddle.pdf";
bool IsEmpty=false;
ofstream fcsv(outcsv, ofstream::out);
TFile *fout = new TFile(outroot,"RECREATE");

void tofw_beta_offset();
void initialise();

//////////
void tofw_beta_offset(){
  initialise();
  TCanvas *c = new TCanvas();
  c->Divide(4,4);
  c->Print(outpdf + "[");
  //

  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h%i",i),Form("h%i",i),500,-0.6,0.6,500,-55,55);
    ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ)>>h%i",i),Form("Tofw_Paddle==%i && %s && abs(MusicZ-%i)<0.4",i+1, Zgate.Data(), zetfrs),"colz");
    h->Write();
    //ch->Scan("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ):Tofw_Paddle","","",10);
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Angle change
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h2_%i",i),Form("h2_%i",i),500,-0.6,0.6,500,-55,55);
    ch->Draw(Form("(TwimTheta-MusicTheta)*1000.:(AoQ_S2_Cave-FragAoQ)>>h2_%i",i),Form("Tofw_Paddle==%i && %s && abs(MusicZ-%i)<0.4",i+1, Zgate.Data(), zetfrs),"colz");
    h->Write();
    //ch->Scan("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ):Tofw_Paddle","","",10);
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // AoQ Zet
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h3_%i",i),Form("h3_%i",i),500,-0.6,0.6,500,15,25);
    ch->Draw(Form("FragZ:(AoQ_S2_Cave-FragAoQ)>>h3_%i",i),Form("Tofw_Paddle==%i && abs(MusicZ-%i)<0.4",i+1, zetfrs),"colz");
    h->Write();
    //ch->Scan("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ):Tofw_Paddle","","",10);
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Theta Tof
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h4_%i",i),Form("h4_%i",i),500,40,60,500,-55,55);
    ch->Draw(Form("(TwimTheta)*1000.:(FragTof)>>h4_%i",i),Form("Tofw_Paddle==%i && %s && abs(MusicZ-%i)<0.4",i+1, Zgate.Data(), zetfrs),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //

  // Unreacted cut
  TF1 *fit[NUMPADDLE];
  Double_t peakX[NUMPADDLE]={1000.};
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH1F(Form("h_noreaction_%i",i),Form("h_noreaction_%i",i),500,-0.6,0.6);
    ch->Draw(Form("(AoQ_S2_Cave-FragAoQ)>>h_noreaction_%i",i),Form("Tofw_Paddle==%i && %s && abs(MusicZ-%i)<0.4",i+1, Zgate.Data(), zetfrs),"colz");
    TSpectrum *tspec = new TSpectrum(10);
    Int_t npeaks = tspec->Search(h);
    Double_t *peakXspec = tspec->GetPositionX();
    peakX[i]=2000.;
    for(int j = 0; j<npeaks; j++){
      //peakX[i] = (abs(peakX[i]) < abs(peakXspec[j]))? peakX[i] : peakXspec[j];
      if(abs(peakX[i]) > abs(peakXspec[j])) peakX[i] = peakXspec[j];
      cout<<peakXspec[j]<<" "<<peakX[i]<<endl;
    }      
    //
    fit[i] = new TF1(Form("f_aoqdiff%i",i),"gaus");
    cout<<"peakX"<<i+1<<": "<<peakX[i]<<endl;
    fit[i]->SetParameter(0,100);
    fit[i]->FixParameter(1,peakX[i]);
    fit[i]->SetParameter(2,0.01);
    h->Fit(fit[i],"LL","",peakX[i]-0.03,peakX[i]+0.03);
    //
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Angle change (gated)
  TF1 *fit_prof[NUMPADDLE];
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h2_gated__%i",i),Form("h2_gated__%i",i),500,-55,55,500,-0.6,0.6);
    TCut cut = Form("abs(AoQ_S2_Cave-FragAoQ-%f)<2.*%f && abs(TwimTheta)*1000.<8.", fit[i]->GetParameter(1), fit[i]->GetParameter(2));
    ch->Draw(Form("(AoQ_S2_Cave-FragAoQ):TwimTheta*1000.>>h2_gated__%i",i),cut && Form("Tofw_Paddle==%i && %s && abs(MusicZ-%i)<0.4",i+1, Zgate.Data(), zetfrs),"colz");
    h->Write();
    TProfile* prof = h->ProfileX();
    prof->Draw("same");
    //fit_prof[i] = new TF1(Form("fit_prof%i",i),"pol2");
    fit_prof[i] = new TF1(Form("fit_prof%i",i),"pol1");
    fit_prof[i]->SetLineWidth(1);
    prof->Fit(fit_prof[i]);
    //
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //
  //Again the AoQDiff
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h2_corrected__%i",i),Form("h2_corrected__%i",i),500,-0.6,0.6,500,-55,55);
    //ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ-%f-TwimTheta*%f-TwimTheta*TwimTheta*%f)>>h2_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000.,fit_prof[i]->GetParameter(2)*1.e6, i),
    ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ-%f-TwimTheta*%f)>>h2_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000., i),
	     Form("Tofw_Paddle==%i && %s && abs(MusicZ-%i)<0.4",i+1, Zgate.Data(), zetfrs),"colz");
    h->Write();
    //
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //
  //Corrected AoQ
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h_aoq_corrected__%i",i),Form("h_aoq_corrected__%i",i),500,1.8,2.8,500,10,22);
    ch->Draw(Form("FragZ:(FragAoQ+%f+TwimTheta*%f)>>h_aoq_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000., i),
	     Form("Tofw_Paddle==%i && abs(MusicZ-%i)<0.4",i+1, zetfrs),"colz");
    h->Write();
    //
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //
  //Corrected AoQ
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h_aoq_all_corrected__%i",i),Form("h_aoq_all_corrected__%i",i),500,1.8,2.8,500,10,22);
    ch->Draw(Form("FragZ:(FragAoQ+%f+TwimTheta*%f)>>h_aoq_all_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000., i),
	     Form("Tofw_Paddle==%i",i+1),"colz");
    h->Write();
    //
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //
  c->Print(outpdf + "]");
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
  //ch -> SetProof();
  //
  //frspidgate = Form("abs(MusicZ-%i)<0.4 && abs(AoQ_S2_Cave -%f)<%f",zetdiff, (double)mass/(double)zetdiff, 0.22/(double)zetdiff);
  //
  brho = "Brho_S2_Cave";
  beta_tofw_mod = "(FragBeta)";
  fragaoqstring = "FragAoQ";
  cut_mw = Form("(Mw1_X>%f)&&(Mw1_X<%f)&&", -30.,25.);
  for(int i=0;i<(IsEmpty?4:3);i++) cut_mw += Form("(%s>%f)&&(%s<%f)&&",axis_mw3[i].Data(),range_cut_mw3_low[i],axis_mw3[i].Data(),range_cut_mw3_high[i]);
  Zgate=Form("abs(MusicZ-FragZ-%f)<0.4",(float)zetdiff);
}
