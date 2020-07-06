#include "TProof.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"
//
TChain *ch = new TChain("FrsChain");
TFile *fout;
TCanvas *c;
TH2F *hPID28, *hPID2C, *hPID8C, *hMusic;
TString dummy, infilename, outfilename;

void FRStree(int FRSsetting){
  //TProof *proof = TProof::Open("lite://", "workers=20");
  TProof::Open("");
  infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles_local/s467_FRSTree_Setting%i_0*.root", FRSsetting);
  dummy = infilename + "/FrsTree";
  outfilename = Form("./FRStree/Histograms_FRS%i", FRSsetting);
  //
  ch->Add(dummy);
  //ch->SetProof();
  //
  cout << ch->GetListOfFiles()->GetEntries() << "files were imported." <<endl;
  for(int i = 0; i < ch->GetListOfFiles()->GetEntries(); i++){
    cout << "File" << i + 1 << ": " << ch->GetListOfFiles()->At(i)->GetTitle() << endl;
  }
  cout << "Number of entries: " << ch->GetEntries() <<endl;
  //ch->Print();
  //
  
  //
  //fout = new TFile(Form("./FRSsetting/Histograms_FRS%i.root",FRSseting), "recreate");
  dummy = outfilename + ".root";
  fout = new TFile(dummy, "recreate");
  c = new TCanvas("c", "c", 1000, 1000);
  c -> Divide(2,2);
  dummy = outfilename + ".pdf[";
  c->Print(dummy);
  hPID28 = new TH2F("hPID28", Form("PID (S2-S8) of FRS setting %i", FRSsetting), 1000, 1.48, 2.82, 1000, 12.5, 27.5);
  hPID2C = new TH2F("hPID2C", Form("PID (S2-CaveC) of FRS setting %i", FRSsetting), 1000, 1.48, 2.82, 1000, 12.5, 27.5);
  hPID8C = new TH2F("hPID8C", Form("PID (S8-CaveC) of FRS setting %i", FRSsetting), 1000, 1.48, 2.82, 1000, 12.5, 27.5);
  hMusic = new TH2F("hMusic", "MusicE vs MusicZ", 1000, 0, 100, 1000, 0, 10000);
  //////
  c->cd(1);
  ch->Draw("MusicZ:AoQ_S2_S8>>hPID28","MusicE>0","colz");
  hPID28->Write();
  dummy = outfilename + ".pdf";
  //c->Write();
  //c->Print(dummy);
  //
  c->cd(2);
  ch->Draw("MusicZ:AoQ_S2_Cave>>hPID2C","MusicE>0","colz");
  hPID2C->Write();
  //
  c->cd(3);
  ch->Draw("MusicZ:AoQ_S8_Cave>>hPID8C","MusicE>0","colz");
  hPID8C->Write();
  //
  c->cd(4);
  ch->Draw("MusicE:MusicZ>>hMusic","MusicE>0","colz");
  hMusic->Write();
  //
  c->Write();
  dummy = outfilename + ".pdf";
  c->Print(dummy);

  ////
  dummy = outfilename + ".pdf]";
  c->Print(dummy);
}


void FRStree() {
  FRStree(11);
}
