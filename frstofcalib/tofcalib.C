const int NUMTOF = 3;
const char *tofname[NUMTOF] = {"S2-CaveC", "S8-CaveC", "S2-S8"};
const char *settingname[3] = {"580AMeV no Target", "450AMeV no Target", "580AMeV Be Target"};
const int run[] = {237, 238, 239, 240, 241, 242};//, 245, 246, 248}; // no 247
const int setting[] = {0,0,1,2,2,2};//,0,0,0};
const int NUMRUNS = sizeof(run)/sizeof(run[0]);
const double veloc28[3] = {0.7812423, 0.7494637, 0.7286703}; // s2-s8
const double veloc8c[3] = {0.7796387, 0.7473292, 0.7261188}; // s8-caveC

double tof[NUMTOF][NUMRUNS] = {0.};  double toferr[NUMTOF][NUMRUNS] = {0.};
//double tof[NUMRUNS][NUMTOF] = {0.};  double toferr[NUMRUNS][NUMTOF] = {0.};
TCanvas *c[NUMTOF], *c2;
TH1D *h[NUMRUNS][NUMTOF];
TGraphErrors *g[NUMTOF];
TF1 *f = new TF1("func","[0]*TMath::Gaus(x,[1],[2])");
TF1 *fline = new TF1("line","[0]+x*[1]");
TFile *fin[NUMRUNS];
TFile *fout;
TString suffix = "_m1TRef_237-242";
TString fileout = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/www/ToFCalibOut"+suffix+".root";
TString pdfout = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofcalib/ToFCalibOut"+suffix+".pdf";
TLatex *latex = new TLatex();

int  gettoffit(int runid){
  //TString filename = "/u/taniuchi/s467/rootfiles/data_s467_0" + run[runid] +".root";
  TString filename = Form("/u/taniuchi/s467/rootfiles/data_s467_%04d.root",  run[runid]);
  //
  fin[runid] = TFile::Open(filename, "READ");
  if(!fin[runid]) cerr<<"err";
  //Mult ==1 in reference timing
  h[runid][0]= (TH1D*)(fin[runid]->Get("RawTofNs_m1_wTref_Sci02_to_Sci04")->Clone());
  h[runid][1]= (TH1D*)(fin[runid]->Get("RawTofNs_m1_wTref_Sci03_to_Sci04")->Clone());
  h[runid][2]= (TH1D*)(fin[runid]->Get("RawTofNs_m1_wTref_Sci02_to_Sci03")->Clone());
  //RawTofNs_wTref_Sci should be better
  //h[runid][0]= (TH1D*)(fin[runid]->Get("RawTofNs_wTref_Sci02_to_Sci04")->Clone());
  //h[runid][1]= (TH1D*)(fin[runid]->Get("RawTofNs_wTref_Sci03_to_Sci04")->Clone());
  //h[runid][2]= (TH1D*)(fin[runid]->Get("RawTofNs_wTref_Sci02_to_Sci03")->Clone());
  //RawTofNs without Tref
  //h[runid][0]= (TH1D*)(fin[runid]->Get("RawTofNs_m1_Sci02_to_Sci04")->Clone());
  //h[runid][1]= (TH1D*)(fin[runid]->Get("RawTofNs_m1_Sci03_to_Sci04")->Clone());
  //h[runid][2]= (TH1D*)(fin[runid]->Get("RawTofNs_m1_Sci02_to_Sci03")->Clone());
  for(int i=0; i<NUMTOF; i++){
    h[runid][i]->SetName(Form("%4d_%s",run[runid],tofname[i]));
    h[runid][i]->SetTitle(Form("%s %s",settingname[setting[i]],h[runid][i]->GetName()));
    double maxX=h[runid][i]->GetBinCenter(h[runid][i]->GetMaximumBin());
    cerr << maxX <<endl;
    f->SetParLimits(0,0.1*h[runid][i]->GetMaximum(),10.*h[runid][i]->GetMaximum());
    f->FixParameter(1,maxX);
    f->SetParLimits(2,0.001,10);
    h[runid][i]->Fit(f,"","",maxX-1,maxX+1);
    //
    f->SetParLimits(0,0.1*f->GetParameter(0),10.*f->GetParameter(0));
    f->SetParameter(1,f->GetParameter(1));
    f->SetParLimits(1,f->GetParameter(1)-2.,2.+f->GetParameter(1));
    f->SetParLimits(2,0.5*f->GetParameter(2),10.*f->GetParameter(2));
    h[runid][i]->Fit(f,"","",maxX-1,maxX+1);
    //
    h[runid][i]->GetXaxis()->SetRangeUser(f->GetParameter(1)-5.,f->GetParameter(1)+5.);
    tof[i][runid] = f->GetParameter(1);
    toferr[i][runid] = sqrt(pow(f->GetParameter(2),2.)+pow(f->GetParError(1),2.));
  }
  return 0;
}

void tofcalib(){
  if(NUMRUNS!=sizeof(setting)/sizeof(setting[0])){ cerr <<"Size of runnum and setting are different"<<endl; return;}
  fout= TFile::Open(fileout,"RECREATE");
  //
  for(int i=0; i<NUMRUNS; i++){
    //if(i>0) break;
    if(gettoffit(i)!=0) cerr<<"Read Error in run "<<run[i]<<endl;
  }
  //
  for(int j=0; j<NUMTOF; j++){
    c[j] = new TCanvas(tofname[j],tofname[j],1000,1000);
    //    c[j] -> SetTitle(tofname[j]);
    int div = (int)(sqrt(NUMRUNS)+0.5);
    c[j] -> Divide(div,div);
    for(int i=0; i<NUMRUNS; i++){
      fout->cd();
      h[i][j]->Write();
      c[j]->cd(i+1);
      h[i][j]->Draw();
    }
    c[j]->Write();
    if(j==0) c[j]->Print(pdfout+"[","pdf");
    c[j]->Print(pdfout,"pdf");
  }
  //
  //
  //
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
    //latex->SetTextHeight(0.03);
    latex->DrawLatex(0.75,tofbeta[0],Form("Offset: %f",-1.*fline->GetParameter(1)));
    latex->DrawLatex(0.75,tofbeta[0]-0.01*fline->GetParameter(1),Form("PATH/c: %f",fline->GetParameter(0)));
  }
  c2->Write();
  c2->Print(pdfout,"pdf");  
  //
  c2->Print(pdfout+"]","pdf");
  fout->Close();
}
