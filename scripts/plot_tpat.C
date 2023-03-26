const Long64_t search_range = 40000; // 4us for time correlation?

TFile *f[3], *fout;
TCanvas *c;
TLegend *lg, *lg2[3];
TTree *tr[3];
TH1F *h_tpat[3], *h_califa_hits[3][3];
TH2F *h_califa[3], *h_califa_time1mev[3], *h_califa_time10mev[3], *h_ts_tstart[3];
TString targ[3] = {"empty", "carbon", "CH2"};
TString cluster_type[3] = {"proton", "gamma", "saturated"};
int color[3] = {kBlack, kBlue, kRed};
const Int_t CALIFA_Tmin = 1500, CALIFA_Tmax = 2500;
TString outdir = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/scripts/califa/";
//
void plot_tpat(){
  f[0] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/Feb2023_FSnov22/s467_filltree_Setting13_0340_24Mar.root", "READ"); // empty
  f[1] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/Feb2023_FSnov22/s467_filltree_Setting13_0339_24Mar.root", "READ"); // carbon
  f[2] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/Feb2023_FSnov22/s467_filltree_Setting13_0342_24Mar.root", "READ"); // CH2
  //f[2] = TFile::Open("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiles/CALIFA2023/s467_filltree_Setting13_0338_24Mar.root", "READ"); // CH2
  fout = TFile::Open(outdir + "plot_tpat.root","RECREATE");
  //
  c = new TCanvas("c","c");
  lg = new TLegend(0.7,0.8,0.9,0.9);
  for(int i=0; i<3; i++){
    lg2[i] = new TLegend(0.7,0.8,0.9,0.9);
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
    title = "CALIFA Timing vs ID for 1 MeV hits for " +targ[i] + " target in CAL-level; time difference; id";
    h_califa_time1mev[i] = new TH2F(dummy, title, 500, -1000, 4000, 5000, -0.5, 4999.5);
    dummy = "h_califa_time10_" + targ[i];
    //title = "CALIFA Timing vs ID for 10 MeV hits for " +targ[i] + "; time difference; id";
    title = "CALIFA Timing vs ID for proton hits for " +targ[i] + " target in Cluster-level; time difference; id";
    h_califa_time10mev[i] = new TH2F(dummy, title, 500, -1000, 4000, 5000, -0.5, 4999.5);
    //
    /*
    dummy = "h_ts_tstart" + targ[i];
    title = "Time correlations" + targ[i] + "; TStart / ns; Time-Stamp difference (Califa-Master)";
    h_ts_tstart[i] = new TH2F(dummy, title, 10000, 1000, 4000, 10000, 1000, 4000);
    */
    //
    for(int j=0; j<3; j++){
      dummy = "h_califa_hits_" + targ[i] + "_" + cluster_type[j];
      title = "Califa multiplicity for " + cluster_type[j] + " hits " + targ[i] + " target; Multiplicity; Counts";
      h_califa_hits[i][j] = new TH1F(dummy, title, 10, -0.5, 9.5);
    }
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
      h_tpat[i]->Fill(tpat);
      //
      int n_cal = ca_cal->GetEntries();
      int n_clu = ca_clu->GetEntries();
      if(j<n_verbose) cout << j << " tpat:" << tpat<<",CALIFA crystals"<<n_cal<<",  CALIFA cluster hits: "<<n_clu<<endl;
      for(int i_cal = 0; i_cal<n_cal; i_cal++){
	auto *hit_cal = dynamic_cast<R3BCalifaCrystalCalData*>(ca_cal->At(i_cal));
	Long64_t tdiff = hit_cal ->GetTime() - b_eh->GetTimeStamp();
	auto ene = hit_cal->GetEnergy();
	auto id = hit_cal->GetCrystalId();
	if(abs(ene-1000.)<10.) h_califa_time1mev[i]->Fill(tdiff, id);
      }
      //
      int n_califa_hits[3] = {0,0,0};
      for(int i_clu = 0; i_clu<n_clu; i_clu++){
	auto *hit_clu = dynamic_cast<R3BCalifaClusterData*>(ca_clu->At(i_clu));
	Long64_t tdiff = hit_clu ->GetTime() - b_eh->GetTimeStamp();
	//Double_t ts_califa = (Double_t)remainder(tdiff,search_range);
	//Double_t tcalifa = remainder(ts_califa-b_eh->GetTStart()/1000., (Double_t)search_range);
	//h_ts_tstart[i]->Fill(tcalifa,tdiff);
	//
	auto ene = hit_clu->GetEnergy();
	auto id = hit_clu->GetMotherCrystal();
	auto type = hit_clu->GetClusterType();
	if(type<0 || type>2) continue;
	n_califa_hits[type]++;
	if(j<n_verbose)cout << "CALIFA: "<<id<< ", Tdiff="<<tdiff <<endl;
	h_califa[i]->Fill(tdiff, ene);
	//if(abs(ene-1.e4)<100.) h_califa_time10mev[i]->Fill(tdiff, id);
	if(type==0) h_califa_time10mev[i]->Fill(tdiff, id);
      }
      for(int type=0; type<3; type++)
	h_califa_hits[i][type]->Fill(n_califa_hits[type]);
    }
  }
  //
  fout->cd();
  for(int i=0; i<3; i++){
    h_tpat[i]->Write();
    h_califa[i]->Write();
    h_califa_time1mev[i]->Write();
    h_califa_time10mev[i]->Write();
    //h_ts_tstart[i]->Write();
    //
    for(int j=0; j<3; j++){
      h_califa_hits[i][j]->Write();
    }
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
  c->Print(outdir + "tpat.png");
  //
  for(int i=0; i<3; i++){
    c->SetLogy(0);
    h_califa[i]->Draw("colz");
    c->Print(outdir + h_califa[i]->GetName() + ".png");
    h_califa_time1mev[i]->Draw("colz");
    c->Print(outdir + h_califa_time1mev[i]->GetName() + ".png");
    h_califa_time10mev[i]->Draw("colz");
    c->Print(outdir + h_califa_time10mev[i]->GetName() + ".png");
    //
    c->SetLogy(1);
    for(int j=0; j<3; j++){
      h_califa_hits[i][j]->SetLineColor(color[j]);
      lg2[i]->AddEntry(h_califa_hits[i][j],cluster_type[j],"l");
      h_califa_hits[i][j]->Draw(j==0?"":"same");
    }
    lg2[i]->Draw();
    c->Print(outdir + h_califa_hits[i][0]->GetName() + ".png");
  }
  
  fout->Close();
}
