const char *settingname[3] = {"","",""};
const int run[] = {340,341,342,343,344,345,346,347,348};
const int setting[] = {0,0,0,0,0,0,0,0,0}; //if 0 use
const int NUMRUNS = sizeof(run)/sizeof(run[0]);

TCanvas *c, *c2, *cpid, *cproj, *cprojst;
TPad *p1;
TH1D *h[NUMRUNS], *hproj[NUMRUNS];
TH2D *hpid[NUMRUNS];
THStack hsum("hs","sum");
THStack hsumproj("hsp","R3BMUSIC Z: Run 340-348");
TF1 *f = new TF1("func","[0]*TMath::Gaus(x,[1],[2])");
TF1 *fline = new TF1("line","[0]+x*[1]");
TFile *fin[NUMRUNS];
TFile *fout;
TString suffix = "";
TString fileout = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/www/R3BMUSICCalibOut"+suffix+".root";
TString pdfout = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/r3bmusiccalib/R3BMUSICCalibOut"+suffix+".pdf";
int cutd=0, cutu=0;
TLine *lined, *lineu;

int  gethist(int runid){
  TString filename = Form("/u/taniuchi/s467/rootfiles/data_s467_%04d.root",  run[runid]);
  //
  fin[runid] = TFile::Open(filename, "READ");
  if(!fin[runid]) cerr<<"err";
  h[runid]= (TH1D*)(fin[runid]->Get("fh1_Mus_charge_z")->Clone());
  h[runid]->SetName(Form("Run%4d",run[runid]));
  h[runid]->SetTitle(Form("%s%s",settingname[setting[runid]],h[runid]->GetName()));
  h[runid]->SetFillStyle(0);
  //
  hpid[runid] = (TH2D*)(fin[runid]->Get("fh2v_Aq_vs_q_frs")->Clone());
  if(cutu==0){
    cutd = hpid[runid]->GetXaxis()->FindBin(2.550);
    cutu = hpid[runid]->GetXaxis()->FindBin(2.555);
    lined = new TLine(cutd, 8., cutd, 39.5);
    lineu = new TLine(cutu, 8., cutu, 39.5);
  }
  hproj[runid] = hpid[runid]->ProjectionY(Form("Run%i",run[runid]), cutd, cutu);  
  hproj[runid]->SetFillStyle(0);
  hproj[runid]->SetTitle("MUSIC charge distribution in {}^{50}Ca setting");
  /*
  double maxX=h[runid]->GetBinCenter(h[runid]->GetMaximumBin());
    cerr << maxX <<endl;
    f->SetParLimits(0,0.1*h[runid]->GetMaximum(),10.*h[runid]->GetMaximum());
    f->FixParameter(1,maxX);
    f->SetParLimits(2,0.001,10);
    h[runid]->Fit(f,"","",maxX-1,maxX+1);
    //
    f->SetParLimits(0,0.1*f->GetParameter(0),10.*f->GetParameter(0));
    f->SetParameter(1,f->GetParameter(1));
    f->SetParLimits(1,f->GetParameter(1)-2.,2.+f->GetParameter(1));
    f->SetParLimits(2,0.5*f->GetParameter(2),10.*f->GetParameter(2));
    h[runid]->Fit(f,"","",maxX-1,maxX+1);
    //
    h[runid]->GetXaxis()->SetRangeUser(f->GetParameter(1)-5.,f->GetParameter(1)+5.);
    tof[runid] = f->GetParameter(1);
    toferr[runid] = sqrt(pow(f->GetParameter(2),2.)+pow(f->GetParError(1),2.));
    */
  return 0;
}

void musiccalib(){
  if(NUMRUNS!=sizeof(setting)/sizeof(setting[0])){ cerr <<"Size of runnum and setting are different"<<endl; return;}
  fout= TFile::Open(fileout,"RECREATE");
  //
  for(int i=0; i<NUMRUNS; i++){
    if(setting[i]!=0) continue;
    if(gethist(i)!=0) cerr<<"Read Error in run "<<run[i]<<endl;
  }
  //
  c = new TCanvas("c_music","music",1000,1000);
  int div = (int)(sqrt(NUMRUNS-1)+1.);
  c -> Divide(div,div);
  for(int i=0; i<NUMRUNS; i++){
    if(setting[i]!=0) continue;
    fout->cd();
    h[i]->GetXaxis()->SetRangeUser(15,28);
    h[i]->Write();
    p1 = (TPad*)(c->cd(i+1));
    p1->SetLogy();
    h[i]->Draw();
    hsum.Add(h[i]);
  }
  c->Write();
  c->Print(pdfout+"(");
  c->Print(pdfout);

  c2 = new TCanvas("c_musicsum","musicsum",1000,1000);
  c2->cd();
  c2->SetLogy();
  //(hsum.GetXaxis())->SetRangeUser(15,25);
  hsum.Draw();
  hsum.Write();
  c2->Write();
  c2->Print(pdfout);
  //

  cpid = new TCanvas("pid","pid",1000,1000);
  cpid -> Divide(div,div);
  cproj = new TCanvas("pidproj","PID Proj",1000,1000);
  cproj-> Divide(div,div);
  cprojst = new TCanvas("pidprojstack","PID Proj",1000,1000);
  for(int i=0; i<NUMRUNS; i++){
    if(setting[i]!=0) continue;
    cpid->cd(i+1);
    hpid[i]->Draw("colz");
    lined->Draw("same");
    lineu->Draw("same");
    //
    cproj->cd(i+1);
    p1 = (TPad*)(cproj->cd(i+1));
    p1->SetLogy();
    hproj[i]->GetXaxis()->SetRangeUser(15,28);
    hproj[i]->Draw();
    hsumproj.Add(hproj[i]);
    //
    fout->cd();
    hpid[i]->Write();
    hproj[i]->Write();
  }
  cpid->Write();
  cproj->Write();
  cpid->Print(pdfout);
  cproj->Print(pdfout);
  //
  cprojst->cd();
  cprojst->SetLogy();
  //hsumproj.Draw();
  hproj[2]->Draw();//run342 is 1hr run
  hsumproj.Write();
  cprojst->Write();
  cprojst->Print(pdfout);
  /*
  c2 = new TCanvas("ToFCalib","ToF Calib", 1000, 1000);
  c2->Divide(2,2);  
  for(int j=0; j<NUMTOF; j++){
    double beta[NUMRUNS], tofbetaerr[NUMRUNS], tofbeta[NUMRUNS];
    for(int i=0; i<NUMRUNS; i++){
      beta[i]=veloc28[setting[i]];
      tofbeta[i]=tof[j][i] * beta[i];
      tofbetaerr[i]=toferr[j][i] * beta[i];
    }
    g[j] = new TGraphErrors(NUMRUNS, beta, tofbeta, 0,tofbetaerr);
    g[j]->SetName(tofname[j]);
    g[j]->SetTitle(tofname[j]);
    g[j]->SetMarkerStyle(6); //7 small dot, 8 circle
    //g[j]->SetMarkerSize(1);
    g[j]->GetXaxis()->SetTitle("#beta");
    g[j]->GetYaxis()->SetTitle("#beta #times ToF");
    g[j]->Fit(fline,"","");
    g[j]->Write();
    c2->cd(j+1);
    g[j]->Draw("AP");
  }
  c2->Write();
  c2->Print(pdfout);  
  //
  */
  c->Print(pdfout+")");
  fout->Close();
}
