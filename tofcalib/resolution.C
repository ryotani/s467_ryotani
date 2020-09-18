#include "TProof.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"
//
TChain *ch = new TChain("FrsChain");
TFile *fout;
TCanvas *c;
TH2F *hTofCorr[3];
TProfile* Prof[3], *Profy[3];
TString dummy, infilename, outfilename;
Double_t avex[3], avey[3], tofres[3], profres[3], profresy[3];
//char *focus[3] = {"S2", "S8", "Cave"};
//char *setname[5] = {"^{40}Ca", "^{39}Ca", "^{38}Ca", "^{38}Ca", "^{50}Ca"};
//Double_t aoqmin[5] = {1.7, 1.7, 1.7, 1.7, 2.2};
//Double_t aoqmax[5] = {2.2, 2.2, 2.2, 2.2, 2.7};
char *tofname[3] = {"S2-Cave", "S2-S8", "S8-Cave"};

void resolution(int runnum){
  gStyle->SetLabelSize(12,"XYZ");
  gStyle->SetTitleSize(12,"XYZ");
  //TProof *proof = TProof::Open("lite://", "workers=20");
  //TProof::Open("");

  infilename = Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/s467_FRSTree_Setting*_%04d_ToFResolution.root", runnum);
  dummy = infilename + "/FrsTree";
  ch->Add(dummy);
  outfilename = Form("./tofcalib/tof_%04d", runnum);
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
  hTofCorr[0] = new TH2F("hTof2C-28","Tof correlation S2-Cave (X) vs S2-S8 (Y)", 1000, 1909, 1913, 1000, 1050, 1054);
  hTofCorr[1] = new TH2F("hTof28-8C","Tof correlation S2-S8 (X) vs S8-Cave (Y)", 1000, 1050, 1054, 1000, 857, 861);
  hTofCorr[2] = new TH2F("hTof8C-2C","Tof correlation S8-Cave (X) vs S2-Cave (Y)", 1000, 857, 861, 1000, 1909, 1913);
  dummy = outfilename + ".pdf[";
  c->Print(dummy);

  c->cd(1);
  gPad->SetLogz();

  ch->Draw(">>elist","multMapSci[3]*multMapSci[4]*multMapSci[11]==1");
  auto *evtlist = (TEventList*)gROOT->FindObject("elist");
  ch->SetEventList(evtlist);
  //
  c->cd(1);
  ch->Draw("Tof_wTref_S2_S8:Tof_wTref_S2_Cave>>hTof2C-28","","colz");
  c->Write();
  Prof[0] = hTofCorr[0] -> ProfileX("Prof_S2S8",1,-1,"s");
  Prof[0] -> GetYaxis() -> SetRangeUser(1050,1054);
  Profy[0] = hTofCorr[0] -> ProfileY("Profy_S2Cave",1,-1,"s");
  Profy[0] -> GetYaxis() -> SetRangeUser(1909,1913);
  
  Double_t stats0[7];
  hTofCorr[0] -> GetStats(stats0);
  avex[0] = stats0[2]/stats0[0];
  avey[0] = stats0[4]/stats0[0];

  //
  c->cd(2);
  ch->Draw("Tof_wTref_S8_Cave:Tof_wTref_S2_S8>>hTof28-8C","","colz");
  c->Write();
  Prof[1] = hTofCorr[1] -> ProfileX("Prof_S8Cave",1,-1,"s");
  Prof[1] -> GetYaxis() -> SetRangeUser(857,861);
  Profy[1] = hTofCorr[1] -> ProfileY("Profy_S2S8",1,-1,"s");
  Profy[1] -> GetYaxis() -> SetRangeUser(1050, 1054);
  Double_t stats1[7];
  hTofCorr[1] -> Draw("colz");
  hTofCorr[1] -> GetStats(stats1);
  avex[1] = stats1[2]/stats1[0];
  avey[1] = stats1[4]/stats1[0];
  //
  c->cd(3);
  ch->Draw("Tof_wTref_S2_Cave:Tof_wTref_S8_Cave>>hTof8C-2C","","colz");
  c->Write();
  Prof[2] = hTofCorr[2] -> ProfileX("Prof_S2Cave",1,-1,"s");
  Prof[2] -> GetYaxis() -> SetRangeUser(1909,1913);
  Profy[2] = hTofCorr[2] -> ProfileY("Profy_S8Cave",1,-1,"s");
  Profy[2] -> GetYaxis() -> SetRangeUser(857,861);
  
  Double_t stats2[7];
  hTofCorr[2] -> GetStats(stats2);
  avex[2] = stats2[2]/stats2[0];
  avey[2] = stats2[4]/stats2[0];
  
  //
  c->cd(4);
  TText *t;
  for(int i=0; i<3; i++){
    t = new TText(0.05,0.2*(4.-(Double_t)i),Form("Average%i: X:%.2f,  Y:%.2f", i, avex[i], avey[i]));
    t->SetTextAlign(12);
    t->SetTextColor(1);
    t->SetTextFont(43);
    t->SetTextSize(20);
    t->Draw();
  }
  //
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  //
  c->Clear();
  c->Divide(2,2);
  for(int i=0; i<3; i++){
    c->cd(i+1);
    Prof[i]->Draw();
    //
    c->cd(4);
    Int_t bin = Prof[i]-> GetXaxis() -> FindBin(avex[i]);
    profres[i] = Prof[i]->GetBinError(bin);
    //cout<<bin<< " "<<tofres[i]<<endl;
    t = new TText(0.05,0.2*(4.-(Double_t)i),Form("Std.Dev. in Pad %i: %.2f ps\n (ToF: %s)", i, profres[i]*1000., tofname[(i+1)%3]));
    t->SetTextAlign(12);
    t->SetTextColor(1);
    t->SetTextFont(43);
    t->SetTextSize(20);
    t->Draw();
    //
    /*tofres[i] = ((pow(profres[(i+2)%3],2.) - pow(avey[(i+2)%3]/avex[(i+2)%3]*profres[(i+1)%3], 2.)
		      + pow(avey[(i+1)%3]*avey[(i+2)%3]/(avex[(i+1)%3]*avex[(i+2)%3])*profres[(i+2)%3], 2.))
		     /(1.+pow(avey[i]*avey[(i+1)%3]*avey[(i+2)%3]/(avex[i]*avex[(i+1)%3]*avex[(i+2)%3]), 2.)));
    t = new TText(0.55,0.2*(4.-(Double_t)i),Form("Std.Dev.%i: %.3f ns", i, tofres[i]));
    t->SetTextAlign(12);
    t->SetTextColor(1);
    t->SetTextFont(43);
    t->SetTextSize(20);
    t->Draw();*/
  }
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  //
  c->Clear();
  c->Divide(2,2);
  for(int i=0; i<3; i++){
    c->cd(i+1);
    Profy[i]->Draw();
    //
    c->cd(4);
    Int_t bin = Profy[i]-> GetXaxis() -> FindBin(avey[i]);
    profresy[i] = Profy[i]->GetBinError(bin);
    t = new TText(0.05,0.2*(4.-(Double_t)i),Form("Std.Dev. in Pad %i: %.2f ps\n (ToF: %s)", i, profresy[i]*1000., tofname[(i)%3]));
    t->SetTextAlign(12);
    t->SetTextColor(1);
    t->SetTextFont(43);
    t->SetTextSize(20);
    t->Draw();
  }  
  dummy = outfilename + ".pdf";
  c->Print(dummy);
  dummy = outfilename + ".pdf]";
  c->Print(dummy);
}


void resolution() {
  resolution(237);
}
