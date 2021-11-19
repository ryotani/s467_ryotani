#define MINPADDLE 0
#define NUMPADDLE 28

#include "fragmentana_twim.h"
int zetdiff = 0;
double min_aoq = 1.8, max_aoq = 2.7;
//
//for trees
Float_t AoQ_S2_Cave, brho_frs, MusicZ, FragZ, TwimTheta, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ, FragAoQ_corr, FragBrho;
UChar_t Tofw_Paddle;
Long64_t nentry=0;
//
/*
TString targ = "empty";
TString infile = "./fragment/output/mktree_fragment_Nov_empty.root";
*/
//

TString targ = "ch2";
TString infile = "./fragment/output/mktree_fragment_Nov_ch2-24mm.root";

//
/*
TString targ = "carbon";
TString infile = "./fragment/output/mktree_fragment_Nov_carbon.root";
*/
//
/*
TString targ = "PP";
TString infile = "./fragment/output/mktree_fragment_Nov_PP.root";
*/
//
TString outpdf = "./fragment/output/tofw_beta_offset_paddle" + targ + Form("_nofrsgate_zdiff%i",zetdiff) + ".pdf";
TString outcsv = "./fragment/output/tofw_beta_offset_paddle" + targ + Form("_nofrsgate_zdiff%i",zetdiff) + ".csv";
TString outroot = "./fragment/output/tofw_beta_offset_paddle" + targ + Form("_nofrsgate_zdiff%i",zetdiff) + ".root";
//TString recobrho_outpdf = "./fragment/output/tofw_beta_offset_paddlebrho_" + targ + Form("_nofrsgate_zdiff%i",zetdiff) + ".pdf";
bool IsEmpty=false;
ofstream fcsv(outcsv, ofstream::out);
TFile *fout = new TFile(outroot,"RECREATE");
TF1 *fit_prof[NUMPADDLE];
TTree *tree;

void tofw_beta_offset_nofrsgate();
void initialise();
void filltree();
Int_t initbranch();
void writecsv();

//////////
void tofw_beta_offset_nofrsgate(){
  initialise();
  TCanvas *c = new TCanvas();
  c->Divide(4,4);
  c->Print(outpdf + "[");
  //
  
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h_aoqdiff_twim%i",i),Form("Correlation in AoQ and X-angle of Paddle%i;AoQ_{FRS} - AoQ_{Cave};Twim Theta (mrdad)",i+1),500,-0.6,0.6,500,-55,55);
    ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ)>>h_aoqdiff_twim%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Angle change
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h2_%i",i),Form("Correlation in AoQ and X-angle change of Paddle%i;AoQ_{FRS} - AoQ_{Cave};#Theta_{Twim} - #Theta_{R3BMusic} (mrad)",i+1),500,-0.6,0.6,500,-55,55);
    ch->Draw(Form("(TwimTheta-MusicTheta)*1000.:(AoQ_S2_Cave-FragAoQ)>>h2_%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // AoQ Zet
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h3_%i",i),Form("Z vs AoQ difference of Paddle%i;AoQ_{FRS} - AoQ_{Cave};Twim Z",i+1),500,-0.6,0.6,500,15,25);
    ch->Draw(Form("FragZ:(AoQ_S2_Cave-FragAoQ)>>h3_%i",i),Form("Tofw_Paddle==%i ",i+1),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Theta Tof
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h4_%i",i),Form("Tof and Twim Theta of Paddle%i; Tof in Cave (ns); #Theta_{Twim} (mrad)",i+1),500,40,60,500,-55,55);
    ch->Draw(Form("(TwimTheta)*1000.:(FragTof)>>h4_%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
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
    auto *h = new TH1F(Form("h_noreaction_%i",i),Form("AoQ difference of Paddle%i",i+1),500,-0.6,0.6);
    ch->Draw(Form("(AoQ_S2_Cave-FragAoQ)>>h_noreaction_%i",i),Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
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
  // Theta Tof (gated with aoq)
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h4_gated_%i",i),Form("Tof and Twim Theta with aoq gated of Paddle%i",i+1),500,40,60,500,-55,55);
    TCut cut = Form("abs(AoQ_S2_Cave-FragAoQ-%f)<2.*%f && abs(TwimTheta-6e-3)*1000.<14.", fit[i]->GetParameter(1), fit[i]->GetParameter(2)); // cut with +/- 20 mrad
    ch->Draw(Form("(TwimTheta)*1000.:(FragTof)>>h4_gated_%i",i),cut && Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf );
  //
  // Angle change (gated)
    for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h2_gated__%i",i),Form("AoQ difference and Twim theta of Paddle%i; #Theta_{Twim} (mrad); AoQ_{FRS} - AoQ_{Cave}",i+1),500,-55,55,500,-0.6,0.6);
    TCut cut = Form("abs(AoQ_S2_Cave-FragAoQ-%f)<2.*%f && abs(TwimTheta-6e-3)*1000.<14.", fit[i]->GetParameter(1), fit[i]->GetParameter(2)); // cut with +/- 20 mrad
    TCut tofcut = Form("abs(FragTof-40-0.4*%i)<2.",i+1);
    ch->Draw(Form("(AoQ_S2_Cave-FragAoQ):TwimTheta*1000.>>h2_gated__%i",i),cut && tofcut && Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    TProfile* prof = h->ProfileX();
    prof->Draw("same");
    //fit_prof[i] = new TF1(Form("fit_prof%i",i),"pol3");
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
    auto *h = new TH2F(Form("h2_corrected__%i",i),Form("Corrected AoQ difference and Twim theta of Paddle%i; #Theta_{Twim} (mrad); AoQ_{FRS} - Corrected AoQ_{Cave}",i+1),500,-0.6,0.6,500,-55,55);
    //ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ-%f-TwimTheta*%f-TwimTheta*TwimTheta*%f-(TwimTheta**3) *%f)>>h2_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000.,fit_prof[i]->GetParameter(2)*1.e6,fit_prof[i]->GetParameter(3)*1.e9, i),
    //ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ-%f-TwimTheta*%f-TwimTheta*TwimTheta*%f)>>h2_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000.,fit_prof[i]->GetParameter(2)*1.e6, i),
    ch->Draw(Form("TwimTheta*1000.:(AoQ_S2_Cave-FragAoQ-%f-TwimTheta*%f)>>h2_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000., i),
	     Form("Tofw_Paddle==%i && %s ",i+1, Zgate.Data()),"colz");
    h->Write();
    //
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  c->Print(outpdf);
  //
  //Corrected AoQ
  TH2F *htot;
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd((i%(NUMPADDLE/2))+1);
    auto *h = new TH2F(Form("h_aoq_corrected__%i",i),Form("Corrected AoQ and Twim theta of Paddle%i; #Theta_{Twim} (mrad); Corrected AoQ_{Cave}",i+1),500,min_aoq,max_aoq,500,10,25);
    ch->Draw(Form("FragZ:(FragAoQ+%f+TwimTheta*%f)>>h_aoq_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000., i),
    //ch->Draw(Form("FragZ:(FragAoQ+%f+TwimTheta*%f+TwimTheta*TwimTheta*%f+(TwimTheta**3) *%f)>>h_aoq_corrected__%i",fit_prof[i]->GetParameter(0),fit_prof[i]->GetParameter(1)*1000.,fit_prof[i]->GetParameter(2)*1.e6,fit_prof[i]->GetParameter(3)*1.e9, i),
	     Form("Tofw_Paddle==%i ",i+1),"colz");
    if(i==MINPADDLE){
      htot = (TH2F*)h->Clone();
    }else{
      htot->Add(h);
    }
    h->Write();
    //
    if(i==NUMPADDLE/2-1) c->Print(outpdf);
  }
  //
  c->cd(NUMPADDLE+1);
  htot->Draw("colz");
  htot->Write();
  c->Print(outpdf);
  c->Print(outpdf + "]");
  //
  filltree();
  writecsv();
}

void filltree(){
  initbranch();
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
      FragAoQ_corr =
	FragAoQ + fit_prof[Tofw_Paddle-1]->GetParameter(0)
	+ fit_prof[Tofw_Paddle-1]->GetParameter(1) * TwimTheta * 1000.;
	//+ fit_prof[Tofw_Paddle-1]->GetParameter(2) * pow(TwimTheta*1000.,2.); // pol2
    }
    tree->Fill();
    if(i%100000==0)
      cout<<i<<" enrties done in "<<nentry<<flush;//endl;//
  }
  cout<<endl;
  tree->Write();
}
///
void writecsv(){
  fcsv<<"test"<<endl;
  
}
///
Int_t initbranch(){
  ch->SetBranchAddress("AoQ_S2_Cave",&AoQ_S2_Cave);
  ch->SetBranchAddress("Brho_S2_Cave",&brho_frs);
  ch->SetBranchAddress("MusicZ", &MusicZ);
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
  tree->Branch("FRS_AoQ",&AoQ_S2_Cave);
  tree->Branch("FRS_Brho",&brho_frs);
  tree->Branch("FRS_Z", &MusicZ);
  tree->Branch("Frag_AoQ", &FragAoQ_corr);
  tree->Branch("Frag_AoQ_original", &FragAoQ);
  tree->Branch("Frag_Brho",&FragBrho);
  tree->Branch("Frag_Z", &FragZ);
  tree->Branch("Frag_Tof", &FragTof);
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
