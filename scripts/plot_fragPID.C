TString dir = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/fragment/output/";
TString infile[2] = {"incl_califa_Mar2023_38Ca_all.root", "incl_califa_Mar2023_50Ca_all.root"};
TString outfile = "incl_califa_Mar2023.pdf";

TFile *fin[2];
TCanvas *c;

void plot_fragPID(){
  for(int i=0; i<2; i++){
    fin[i] = TFile::Open(infile[i], "READ");
  }
  c = new TCanvas("c","",1000,1000);
  fin[1]->cd();
  (TH1F*)fin[1]->Get("h_aoq_gated_all_Z20_A48_Z19_ch2")->Draw();

}
