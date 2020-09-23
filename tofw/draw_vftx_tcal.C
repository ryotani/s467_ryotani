/*
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
*/

TH1F *hft[60];
TH1F *hcal[60];
TCanvas *cft[4], *ccal[4];

void draw_vftx_tcal(int runnum){
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle();
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat();
  TString filename = runnum<1000?Form("/u/taniuchi/s467/rootfiles/calibVftx%d.root",  runnum): "/u/taniuchi/s467/rootfiles/calibVftx_tofw.root";
  TString pdfout =   runnum<1000?Form("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/vftxout%04d.pdf", runnum): "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/vftxout_tofw.pdf";
  
  TFile *f = TFile::Open(filename, "READ");
  for(int i = 0; i<4; i++){
    ccal[i] = new TCanvas(Form("ccal%i",i), Form("ccal%i",i));
    cft[i] = new TCanvas(Form("cft%i",i), Form("cft%i",i));
    ccal[i]->Divide(4,4);
    cft[i]->Divide(4,4);
  }
  for(int i = 0 ; i<56; i++){
    hcal[i] = (TH1F*)(f->Get(Form("TimeFineNs_TofW_P%i_Pmt%i_Sig%i", i/2+1, i%2+1, i))->Clone());
    hft[i] = (TH1F*)(f->Get(Form("TimeFineBin_TofW_P%i_Pmt%i_Sig%i", i/2+1, i%2+1, i))->Clone());
    hcal[i]->SetStats(0);
    //hft[i]->SetOptStat(1111);
    //cout<<hcal[i]->GetName()<<" "<<hft[i]->GetName()<<endl;
    ccal[i/14]->cd(i%14+1);
    hcal[i]->Draw();
    cft[i/14]->cd(i%14+1);
    hft[i]->Draw();
  }
  ccal[0]->Print(pdfout+"[");
  //ccal[0]->Print(Form("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/vftxout%04d.png", runnum));
  for(int i = 0; i<4; i++){
    ccal[i]->Print(pdfout);
  }
  for(int i = 0; i<4; i++){
    cft[i]->Print(pdfout);
  }
  ccal[0]->Print(pdfout+"]");
}


void draw_vftx_tcal(){
  draw_vftx_tcal(1000);
}
