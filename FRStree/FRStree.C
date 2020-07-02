#include "TProof.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"
//
TChain *ch = new TChain("FrsChain");
TFile *fout;
TCanvas *c;
TH2F *hPID;
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
  dummy = outfilename + ".pdf[";
  c->Print(dummy);
  hPID = new TH2F("hPID", Form("PID of FRS setting %i", FRSsetting), 1000, 1.48, 2.82, 1000, 12.5, 27.5);
  
  ch->Draw("MusicZ:AoQ_S2_S8>>hPID","","colz");
  //////
  c->Write();
  hPID->Write();
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  //

  ////
  dummy = outfilename + ".pdf]";
  c->Print(dummy);
}


void FRStree() {
  FRStree(11);
}
