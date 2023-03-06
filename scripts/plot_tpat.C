TFile *f[3], *fout;
TCanvas *c;
TLegend *lg;
TTree *tr[3];
TH1F *h_tpat[3];
TH2F *h_califa[3], *h_califa_time1mev[3], *h_califa_time10mev[3];
TString targ[3] = {"empty", "carbon", "CH2"};
int color[3] = {kBlack, kBlue, kRed};
void plot_tpat(){
  f[0] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/Feb2023_FSnov22/s467_filltree_Setting13_0340_12Feb.root", "READ"); // empty
  f[1] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/Feb2023_FSnov22/s467_filltree_Setting13_0339_12Feb.root", "READ"); // carbon
  f[2] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/Feb2023_FSnov22/s467_filltree_Setting13_0342_12Feb.root", "READ"); // CH2
  fout = TFile::Open("plot_tpat.root","RECREATE");
  //
  c = new TCanvas("c","c");
  lg = new TLegend(0.7,0.8,0.9,0.9);
  for(int i=0; i<3; i++){
    cout<<endl<<targ[i]<<endl;
    tr[i] = dynamic_cast<TTree*>(f[i]->Get("evt"));
    //
    TString dummy = "htpat_" + targ[i];
    //TString title = "On-spill tpat for " + targ[i] + "; tpat; counts";
    TString title = "";
    h_tpat[i] = new TH1F(dummy, title, 15, 0.5, 15.5);
    dummy = "h_califa_" + targ[i];
    title = "CALIFA Energy vs Timing; Timestamp difference (CALIFA - Master); Energy / 10 keV";
    h_califa[i] = new TH2F(dummy, title, 500, -1000, 4000, 5e4, 0, 5.e5);
    //
    dummy = "h_califa_time1_" + targ[i];
    title = "CALIFA Timing vs ID for 1 MeV hits for " +targ[i] + "; time difference; id";
    h_califa_time1mev[i] = new TH2F(dummy, title, 500, -1000, 4000, 5000, -0.5, 4999.5);
    dummy = "h_califa_time10_" + targ[i];
    title = "CALIFA Timing vs ID for 10 MeV hits for " +targ[i] + "; time difference; id";
    h_califa_time10mev[i] = new TH2F(dummy, title, 500, -1000, 4000, 5000, -0.5, 4999.5);
    //
    auto *b_eh = new R3BEventHeader();
    auto *ca_clu = new TClonesArray("R3BCalifaClusterData");
    auto *ca_cal = new TClonesArray("R3BCalifaCrystalCalData");
    tr[i]->SetBranchAddress("EventHeader.", &b_eh);
    tr[i]->SetBranchAddress("CalifaCrystalCalData", &ca_cal);
    tr[i]->SetBranchAddress("CalifaClusterData", &ca_clu);
    //
    auto nevt = tr[i]->GetEntries();
    //nevt = 10000;
    int n_verbose = 3;
    for(int j = 0; j < nevt; j++){
      if(j%1000 == 0) cout<<j<<" / "<<nevt<< endl;
      tr[i]->GetEntry(j);
      int tpat = b_eh->GetTpat();
      int n_cal = ca_cal->GetEntries();
      int n_clu = ca_clu->GetEntries();
      if(j<n_verbose) cout << j << " tpat:" << tpat<<",CALIFA crystals"<<n_cal<<",  CALIFA cluster hits: "<<n_clu<<endl;
      h_tpat[i]->Fill(tpat);
      for(int i_cal = 0; i_cal<n_cal; i_cal++){
	auto *hit_cal = dynamic_cast<R3BCalifaCrystalCalData*>(ca_cal->At(i_cal));
	/*if(isnan(hit_cal->GetEnergy())){
	  cout<< "NaN"<<endl;
	  continue;
	  }*/
	Long64_t tdiff = hit_cal ->GetTime() - b_eh->GetTimeStamp();
	auto ene = hit_cal->GetEnergy();
	auto id = hit_cal->GetCrystalId();
	if(j<n_verbose)cout << "CALIFA: "<<id<< ", Tdiff="<<tdiff <<endl;
	h_califa[i]->Fill(tdiff, ene);
	if(abs(ene-1000.)<10.) h_califa_time1mev[i]->Fill(tdiff, id);
	if(abs(ene-1.e4)<100.) h_califa_time10mev[i]->Fill(tdiff, id);
      }
      /*
      if(j<n_verbose) cout<<"Cluster"<<endl;
      for(int i_clu = 0; i_clu<n_clu; i_clu++){
	auto *hit_clu = dynamic_cast<R3BCalifaClusterData*>(ca_clu->At(i_clu));
	if(j<n_verbose)cout << "CALIFA: "<<hit_clu->GetMotherCrystal()<<endl;
      }
      */
    }
  }
  //
  fout->cd();
  for(int i=0; i<3; i++){
    h_tpat[i]->Write();
    h_califa[i]->Write();
    h_califa_time1mev[i]->Write();
    h_califa_time10mev[i]->Write();
    //
  }
  c->cd();
  c->SetLogy();
  for(int i=2; i>=0; i--){
    lg->AddEntry(h_tpat[i], targ[i].Data(),"l");
    h_tpat[i]->SetLineColor(color[i]);
    h_tpat[i]->Draw("same");
  }
  lg->Draw();
  c->Write();
  fout->Close();
}
