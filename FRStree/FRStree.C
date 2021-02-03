#include "TProof.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"
//
TChain *ch = new TChain("FrsChain");
TFile *fout;
TCanvas *c;
TH2F *hPID[4], *hAoQcorr[3], *hMusic, *hMult[3], *hXs2, *hPos[3];
TString dummy, infilename, outfilename;
char *focus[3] = {"S2", "S8", "Cave"};
Double_t AoQmin = 2.2, AoQmax=2.7;

void FRStree(int FRSsetting){
  //TProof *proof = TProof::Open("lite://", "workers=20");
  //TProof::Open("");
  if(FRSsetting != 13){
    AoQmin = 1.75, AoQmax=2.25;
  }
  if(FRSsetting<100){
    infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/tillOct2020/s467_FRSTree_Setting%i_0*_FRSTree.root", FRSsetting);
    //TofW/s467_FRSTree_Setting%i_0*_FragmentTree.root", FRSsetting);
    //infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/TofW/s467_FRSTree_Setting%i_0356_FragmentTree.root", FRSsetting);
    //dummy = infilename + "/Tree";
    dummy = infilename + "/FrsTree";
    ch->Add(dummy);
  }else{
    infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/tillOct2020/s467_FRSTree_Setting%i_0*_FRSTree.root", FRSsetting%100);
    //infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/TofW/s467_FRSTree_Setting%i_0*_FragmentTree.root", FRSsetting%100);
    //dummy = infilename + "/Tree";
    dummy = infilename + "/FrsTree";
    ch->Add(dummy);
    //
    infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/tillOct2020/s467_FRSTree_Setting%i_0*_FRSTree.root", FRSsetting/100);
    //infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/TofW/s467_FRSTree_Setting%i_0*_FragmentTree.root", (int)(FRSsetting/100));
    //dummy = infilename + "/Tree";
    dummy = infilename + "/FrsTree";
    ch->Add(dummy);
  }
  outfilename = Form("./FRStree/Feb2021_Histograms_FRS%i", FRSsetting);
  //
  //ch->SetProof();
  //
  cout << ch->GetListOfFiles()->GetEntries() << "files were imported." <<endl;
  for(int i = 0; i < ch->GetListOfFiles()->GetEntries(); i++){
    cout << "File" << i + 1 << ": " << ch->GetListOfFiles()->At(i)->GetTitle() << endl;
  }
  cout << "Number of entries: " << ch->GetEntries() <<endl;
  //ch->Print();
  //
  //fout = new TFile(Form("./FRSsetting/Histograms_FRS%i.root",FRSseting), "recreate");
  dummy = outfilename + ".root";
  fout = new TFile(dummy, "recreate");
  c = new TCanvas("c", "c", 1000, 1000);
  c -> Divide(2,2);
  dummy = outfilename + ".pdf[";
  c->Print(dummy);
  hPID[0] = new TH2F("hPID28", Form("PID (S2-S8) of FRS setting %i", FRSsetting), 500, AoQmin, AoQmax, 500, 12.5, 27.5);
  hPID[1] = new TH2F("hPID2C", Form("PID (S2-CaveC) of FRS setting %i", FRSsetting), 500, AoQmin, AoQmax, 500, 12.5, 27.5);
  hPID[2] = new TH2F("hPID8C", Form("PID (S8-CaveC) of FRS setting %i", FRSsetting), 500, AoQmin, AoQmax, 500, 12.5, 27.5);
  hPID[3] = new TH2F("hPID", Form("PID of FRS setting %i", FRSsetting), 500, AoQmin, AoQmax, 500, 12.5, 27.5);
  for(int i=0; i<4; i++){
    hPID[i]->GetXaxis()->SetTitle("A/Q");
    hPID[i]->GetYaxis()->SetTitle("R3BMusic Z");
  }
  for(int i=0; i<3; i++){
    hMult[i] = new TH2F(Form("hMult%i",i), Form("Multiplicity of Scintillator %i", i+1), 32, -0.5, 31.5, 32, -0.5, 31.5);
    hMult[i]->GetXaxis()->SetTitle("Mult Wixhausen");
    hMult[i]->GetYaxis()->SetTitle("Mult Messel");
    hPos[i] = new TH2F(Form("hPos%i",i), Form("X position at %s vs A/Q %s (Z=20)", focus[i], i==1?"S2-S8":"S2-Cave"), 500, AoQmin, AoQmax, 500, -20, 20);
    hPos[i]->GetXaxis()->SetTitle("A/Q");
    hPos[i]->GetYaxis()->SetTitle("X position / ns");
  }  
  hAoQcorr[0] = new TH2F("hAoQcorr0", "AoQ correlation: S2-S8 vs S2-Cave", 500, AoQmin, AoQmax, 500, AoQmin, AoQmax);
  hAoQcorr[1] = new TH2F("hAoQcorr1", "AoQ correlation: S2-Cave vs S8-Cave", 500, AoQmin, AoQmax, 500, AoQmin, AoQmax);
  hAoQcorr[2] = new TH2F("hAoQcorr2", "AoQ correlation: S8-Cave vs S2-S8", 500, AoQmin, AoQmax, 500, AoQmin, AoQmax);
  hMusic = new TH2F("hMusic", "MusicE vs MusicZ;MusicZ;MusicE", 500, 0, 100, 500, 0, 10000);
  hXs2 = new TH2F("hXs2", "X-position at S2;A/Q;Xs2 / ns", 500, AoQmin, AoQmax, 500, -100., 100.);
  //////
  c->cd(1);
  gPad->SetLogz();
  ch->Draw("MusicZ:AoQ_S2_S8>>hPID28","MusicZ>0","colz");
  hPID[0]->Write();
  //
  c->cd(2);
  gPad->SetLogz();
  ch->Draw("MusicZ:AoQ_S2_Cave+xs2*1.e-3>>hPID2C","MusicZ>0","colz"); // small fix of dispersion
  hPID[1]->Write();
  //
  c->cd(3);
  gPad->SetLogz();
  ch->Draw("MusicZ:AoQ_S8_Cave>>hPID8C","MusicZ>0","colz");
  hPID[2]->Write();
  //
  c->cd(4);
  gPad->SetLogz();
  //ch->Draw("MusicZ:TheAoQ>>hPID","","colz");
  ch->Draw("MusicZ:AoQ_S2_Cave+xs2*1.e-3>>hPID","MusicZ>0 && abs(Beta_S2_Cave-Beta_S2_S8)<0.003 && abs(Beta_S2_Cave-Beta_S8_Cave)<0.005","colz"); // fix of dispersion + ToF correlation
  hPID[3]->Write();
  //
  c->Write();
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  dummy = outfilename + "_PID.png";
  c->Print(dummy);
  ///
  c->cd(1);
  ch->Draw("AoQ_S2_S8:AoQ_S2_Cave>>hAoQcorr0","","colz");
  //
  c->cd(2);
  ch->Draw("AoQ_S2_Cave:AoQ_S8_Cave>>hAoQcorr1","","colz");
  //
  c->cd(3);
  ch->Draw("AoQ_S8_Cave:AoQ_S2_S8>>hAoQcorr2","","colz");
  //
  c->cd(4);
  //ch->Draw("MusicE:MusicZ>>hMusic","MusicZ>0","colz");
  //hMusic->Write();
  //
  c->Write();
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  dummy = outfilename + "_AoQCorr.png";
  c->Print(dummy);
  ///
  for(int i = 0; i < 3; i++){
    c->cd(1+i);
    //ch->Draw(Form("multMapSci[%i]:multMapSci[%i]>>hMult%i", 3*i + 3, 3*i + 4, i), "MusicZ>0", "colz");
  }
  //
  c->cd(4);
  ch->Draw("xs2:TheAoQ>>hXs2", "MusicZ>0", "colz");
  //
  c->Write();
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  dummy = outfilename + "_mult.png";
  c->Print(dummy);
  /////  
  for(int i = 0; i < 3; i++){
    c->cd(1+i);
    ch->Draw(Form("xpos[%i]:AoQ_S2_%s>>hPos%i", i, i==1?"S8":"Cave", i), "abs(MusicZ-20.)<0.4", "colz");
  }
  //
  c->cd(4);
  //ch->Draw("xs2:TheAoQ>>hXs2", "MusicZ>0", "colz");
  //
  c->Write();
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  dummy = outfilename + "_pos.png";
  c->Print(dummy);
  ////
  dummy = outfilename + ".pdf]";
  c->Print(dummy);
}


void FRStree() {
  FRStree(13);
}
