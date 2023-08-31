#include "pid_gates.h"
TString dir = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/";
TString str_infile_template = dir + "rootfiles/mktree_2023/evttree_20Jul23_SETTING_TARGET.root";
//TString str_infile_template = dir + "rootfiles/mktree_2023/evttree_8Jun23_SETTING_TARGET.root";
//TString str_infile_template = dir + "rootfiles/mktree_2023/mktree_califa_May23_rolufix_SETTING_TARGET.root";
TString str_runlist = dir + "RunSummary.csv";
//TString str_pid0 = dir + "fragment/pid_gate_Settings12_empty.csv"; // Should be updated
//TString str_pid1 = dir + "fragment/pid_gate_Settings13_carbon.csv"; // Should be updated
TString str_pid_template = dir + "fragment/pid_gate_SettingsSET_TARGET.csv";
TString pdfout = "test_paddle_targradius.pdf";
constexpr int division =10;
vector<int> v_mass0 = {39,40,41};
vector<int> v_mass1 = {47,48,49,50,51};
vector<int> v_mass;
vector<TString> targetname = {"empty","carbon","ch2"};
//
vector<pid_gate> v_pid;
std::map<pid_gate::KeyType, pid_gate> pidMap;
TCanvas *c;
vector<vector<TGraphErrors*>> vv_gtrans(3);
vector<vector<TGraphErrors*>> vv_gtransint(3);
vector<vector<TGraphErrors*>> vv_gsigma(3);
vector<TGraphErrors*> v_gtrans_mass(3);
vector<TGraphErrors*> v_gtotal_mass(2);
//
vector<vector<TH1D*>> vv_h_paddle(3);
vector<vector<TGraphErrors*>> vv_gtrans_pad(3);
vector<vector<TGraphErrors*>> vv_gtransint_pad(3);
vector<vector<TGraphErrors*>> vv_gcentroid_pad(3);
vector<vector<TGraphErrors*>> vv_gsigma_pad(3);
vector<TGraphErrors*> v_gtrans_mass_pad(3);
vector<TGraphErrors*> v_gtotal_mass_pad(2);

void fithisto(TH1D* htmp, TF1 *ftmp,  pid_gate PID);
void transmission(TString outfilename = "Transmission.csv");
int transmission_div_evt(int setting = 13, int targ = 2);
int transmission_paddle(int setting = 13, int targ = 2);
void draw_div_evt(vector<int> &v_mass, vector<vector<TGraphErrors*>> &vv_gtrans_tmp, vector<TGraphErrors*> &v_gtrans_mass_tmp, vector<TGraphErrors*> &v_gtotal_mass_pad, double range_min = 45., double range_max = 105.);
//void draw_paddle(vector<int> &v_mass);

auto sigma_tot(double trans1, double trans2, double ntarg) {
  return TMath::Log(trans1/trans2)/ntarg * 1e27;
}
auto Esigma_tot(double trans1, double trans2, double Etrans1, double Etrans2, double ntarg) {
  double e1=Etrans1/trans1;
  double e2=Etrans2/trans2;
  return TMath::Sqrt(e1*e1+e2*e2)/ntarg * 1e27;
}

void fithisto(TH1D* htmp, TF1 *ftmp, pid_gate PID) {
  ftmp->SetNpx(500);
  double tmpmax = htmp->GetMaximum() * TMath::Sqrt(2 * TMath::Pi()) * PID.getFragaqSigma();
  cout<<"Peak height: "<<tmpmax<<endl;
  ftmp->SetParameter(0,tmpmax);
  ftmp->SetParLimits(0,0,2.*tmpmax);
  //ftmp->SetParameter(1,PID.getFragaqCentre());
  ftmp->FixParameter(1,PID.getFragaqCentre());
  //ftmp->SetParLimits(1,PID.getFragaqCentre()-0.02,PID.getFragaqCentre()+0.02);
  ftmp->SetParLimits(2,0.75*PID.getFragaqSigma(), 1.05*PID.getFragaqSigma());
  //
  //ftmp->FixParameter(2,PID.getFragaqSigma());
  //
  htmp->Fit(ftmp,"LL B","",PID.getFragaqCentre()-1./30, PID.getFragaqCentre()+1./30);
  htmp->GetXaxis()->SetRangeUser(PID.getFragaqCentre()-1./10.,PID.getFragaqCentre()+1./10.);
  //return ftmp;
};

void fithisto_paddle(TH1D* htmp, TF1 *ftmp, pid_gate PID) {
  ftmp->SetNpx(500);
  double tmpmax = htmp->GetMaximum() * TMath::Sqrt(2 * TMath::Pi()) * PID.getFragaqSigma();
  cout<<"Peak height: "<<tmpmax<<endl;
  ftmp->SetParameter(0,tmpmax);
  ftmp->SetParLimits(0,0,2.*tmpmax);
  ftmp->SetParameter(1,PID.getFragaqCentre());
  //ftmp->FixParameter(1,PID.getFragaqCentre());
  ftmp->SetParLimits(1,PID.getFragaqCentre()-0.02,PID.getFragaqCentre()+0.02);
  ftmp->SetParLimits(2,0.5*PID.getFragaqSigma(), 1.15*PID.getFragaqSigma());
  htmp->Fit(ftmp,"LL B","",PID.getFragaqCentre()-1./40, PID.getFragaqCentre()+1./40);
  htmp->GetXaxis()->SetRangeUser(PID.getFragaqCentre()-1./10.,PID.getFragaqCentre()+1./10.);
  //
  // second fit
  ftmp->SetParLimits(2,0.5*PID.getFragaqSigma(), 1.5*PID.getFragaqSigma());
  htmp->Fit(ftmp,"LL B","",ftmp->GetParameter(1)-1./50.,ftmp->GetParameter(1)+1./50.); 
};

void transmission(TString outfilename){
  c = new TCanvas();
  //c->Divide(2,2);
  c->SaveAs(pdfout + "[");
  // Set the number of threads for multi-threading
  ROOT::EnableImplicitMT(10);
  for(auto &gtmp: v_gtrans_mass){
    gtmp = new TGraphErrors();
  }
  for(auto &gtmp: v_gtotal_mass){
    gtmp = new TGraphErrors();
  }
  for(auto &gtmp: v_gtrans_mass_pad){
    gtmp = new TGraphErrors();
  }
  for(auto &gtmp: v_gtotal_mass_pad){
    gtmp = new TGraphErrors();
  }
  //
  v_mass = v_mass0;
  transmission_paddle(12,0);
  transmission_paddle(12,1);
  transmission_paddle(12,2);

  transmission_div_evt(12,0);
  transmission_div_evt(12,1);
  transmission_div_evt(12,2);
  
  v_mass = v_mass1;
  transmission_paddle(13,0);
  transmission_paddle(13,1);
  transmission_paddle(13,2);
  transmission_div_evt(13,0);
  transmission_div_evt(13,1);
  transmission_div_evt(13,2);
  v_mass.insert(v_mass.begin(),v_mass0.begin(),v_mass0.end());
  
  draw_div_evt(v_mass, vv_gtrans, v_gtrans_mass, v_gtotal_mass);
  draw_div_evt(v_mass, vv_gtrans_pad, v_gtrans_mass_pad, v_gtotal_mass_pad, 0.);
  //draw_paddle(v_mass);
  c->SaveAs(pdfout + "]");
}

int transmission_div_evt(int setting, int targ){
  TString set_dummy;
  TString str_infile = str_infile_template;
  TString str_pid = str_pid_template;
  if(setting == 13){
    set_dummy="50Ca";
  }else{
    set_dummy="38Ca";
  }
  str_infile.ReplaceAll("SETTING",set_dummy);
  str_infile.ReplaceAll("TARGET", targetname.at(targ));
  str_pid.ReplaceAll("SET", Form("%i",setting));
  str_pid.ReplaceAll("TARGET", targetname.at(targ));
  std::cout<<"Input file: "<<str_infile<<std::endl<<"PID file: "<<str_pid<<std::endl;
  //
  v_pid = readCSV(str_pid);
  for (const auto& p : v_pid) {
    pidMap[p.getKey()] = p;
  }

  int Z = 20;
  auto *fin = TFile::Open(str_infile, "READ");
  auto *evt = (TTree*)fin->Get("evt");
  //auto *evt = (TTree*)fin->Get("Tree");
  auto nentries = evt->GetEntriesFast();
  cout<<"Entries: "<<nentries;

  auto Nmass = v_mass.size();
  vector<TGraphErrors*> &v_gtrans = vv_gtrans.at(targ);
  vector<TGraphErrors*> &v_gtransint = vv_gtransint.at(targ);
  vector<TGraphErrors*> &v_gsigma = vv_gsigma.at(targ);
  auto mass_it = v_gtrans.size();
  vector<int> nin_tot(Nmass +  mass_it);
  for(auto A:v_mass){
    auto *gtmp = new TGraphErrors();
    v_gtrans.push_back(gtmp);
    auto *gtmpint = new TGraphErrors();
    v_gtransint.push_back(gtmpint);
    auto *gtmp2 = new TGraphErrors();
    v_gsigma.push_back(gtmp2);
  }
  //
  vector<TH1D*> v_htot(Nmass + mass_it);
  for(int i=0; i<division; i++){
    auto nfirst = i* nentries/division;
    auto nplot = (i==division-1)?(nentries-nfirst):nentries/division;
    int i_m= mass_it -1;
    for(auto A:v_mass){
      i_m++;
      vector<double> v_trans; //= &vv_trans.at(i_m);
      vector<double> v_etrans; //= &vv_etrans.at(i_m);
      auto gtmp = v_gtrans.at(i_m);
      auto gtmpint = v_gtransint.at(i_m);
      auto gtmp2 = v_gsigma.at(i_m);
      //c->cd(++i);
      TString hname = Form("h%i_%i", A, i);
      auto htmp = new TH1D(hname,hname, 500,1.8,2.8);
      //evt->Draw("FrsData.fZ:FrsData.fAq>>h","","col",nplot,nfirst);
      TCut condition = Form("PID_FRS_Z==%i && PID_FRS_A==%i", Z, A);
      auto nin = evt->Draw("",condition,"",nplot,nfirst);
      nin_tot.at(i_m) += nin;
      //
      pid_gate::KeyType targetKey = std::make_tuple(Z,A,Z,A);
      // Check if the key exists in the map
      auto it = pidMap.find(targetKey);
      if (it == pidMap.end()) {
	continue;
      }
      // Key exists, so retrieve the pid object
      pid_gate PID = it->second;
      //float frszcentreValue = targetPid.getFrszCentre();
      std::cout << std::endl << "Target: "<< targetname[targ] << " fragZ centre value for frsz=" << PID.getFrsz() << ", frsa=" << PID.getFrsa()
		<< ", fragz=" << PID.getFragz() << ", fraga=" << PID.getFraga() << ": "
		<< PID.getFragaqCentre() <<"+/-" << PID.getFragaqSigma() << std::endl;
      //
      TCut twim_cond = Form("abs((Frag_Z-%f)/%f)<3.", PID.getFragzCentre(), PID.getFragzSigma());
      evt->Draw(Form("Frag_AoQ>>%s",hname.Data()),condition && twim_cond,"",nplot,nfirst);
			/////
      auto draw_trans = [&](TH1D* htmp, TF1 *ftmp, int div, int nin, pid_gate PID) {
			//
			TLatex Tl;
			Tl.SetTextAlign(12);
			Tl.SetTextSize(0.04);
			Tl.DrawLatex(PID.getFragaqCentre(), htmp->GetMaximum(), Form("Entries: %lli - %lli", nfirst, nfirst + nplot));
			auto nout = ftmp->GetParameter(0) * 500.;//multiply by binning
			auto xpos = [&](double sigma) {return htmp->GetXaxis()->FindBin(ftmp->GetParameter(1)+sigma*ftmp->GetParameter(2));};
			auto nintegral = htmp->Integral(xpos(-3.),xpos(3.));
			cout<<"Integral bin range: "<<xpos(-3)<< " "<<xpos(3)<<" Integral:"<<nout<<endl;
			auto Enout = ftmp->GetParError(0) * 500.;
			auto trans = (nout/(double)nin)*100.;
			auto etrans = (Enout/(double)nin)*100;
			auto transint = (double)nintegral /(double)nin * 100.;
			auto etransint = transint / TMath::Sqrt(nintegral + nin);
			v_trans.push_back(trans);
			v_etrans.push_back(etrans);
			Tl.DrawLatex(PID.getFragaqCentre(), 0.75 * htmp->GetMaximum(), Form("Incoming: %i, Outgoing: %0.1f, Trans: %0.1f (%i)",(int)nin, nout, trans, (int)(10.*etrans)));
			gtmp->AddPoint(div, trans);
			gtmp->SetPointError(gtmp->GetN()-1, 0, etrans);
			gtmpint->AddPoint(div, transint);
			gtmpint->SetPointError(gtmpint->GetN()-1, 0, etransint);
			gtmp2->AddPoint(div, ftmp->GetParameter(2)/PID.getFragaqSigma());
			gtmp2->SetPointError(gtmp2->GetN()-1, 0, ftmp->GetParError(2)/PID.getFragaqSigma());
		      };
      auto *ftmp = new TF1(Form("f%i_%i",A,i),"gausn",1.8,2.8);
      fithisto(htmp, ftmp, PID);
      draw_trans(htmp, ftmp, i, nin, PID);
      //c->SaveAs(pdfout);
      delete ftmp;
      //
      if(i==0){ // division == 0;
	v_htot.at(i_m) = (TH1D*)htmp->Clone();
	v_htot.at(i_m)->SetTitle(Form("htotal_%i",A));
      }else{
	v_htot.at(i_m)->Add(htmp,1.);
      }
      //
      if(i == division - 1){
	ftmp = new TF1(Form("f%i_%i",A,division),"gausn",1.8,2.8);
	fithisto(v_htot.at(i_m), ftmp, PID);
	draw_trans(v_htot.at(i_m), ftmp, division, nin_tot.at(i_m), PID);
	c->SaveAs(pdfout);
      }
      delete htmp;
    }
  }
  //
  
  auto lg = TLegend(0.3,0.2);//TLegend(0.7*division, 70, 0.8*division, 85);
  int i_m=-1;
  for(auto A:v_mass){
    i_m++;
    auto gtmp = v_gtrans.at(i_m);
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4);//i_m+1); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(45,105);
    gtmp->Draw(i_m==0?"apl":"pl");
    lg.AddEntry(gtmp, Form("Mass %i", A), "pe");
  }
  lg.Draw();
  c->SaveAs(pdfout);
  //
  auto lgint = TLegend(0.3,0.2);//TLegend(0.7*division, 70, 0.8*division, 85);
  i_m=-1;
  for(auto A:v_mass){
    i_m++;
    auto gtmp = v_gtransint.at(i_m);
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4);//i_m+1); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(45,105);
    gtmp->Draw(i_m==0?"apl":"pl");
    lgint.AddEntry(gtmp, Form("Integral: Mass %i", A), "pe");
  }
  lgint.Draw();
  c->SaveAs(pdfout);
  //
  //sigma
  auto lg2 = TLegend(0.3,0.2);
  i_m=-1;
  for(auto A:v_mass){
    i_m++;
    auto gtmp = v_gsigma.at(i_m);
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4);//i_m+1); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(0.5,1.1);
    gtmp->Draw(i_m==0?"apl":"pl");
    lg2.AddEntry(gtmp, Form("Sigma reduction %i", A), "pe");
  }
  lg2.Draw();
  c->SaveAs(pdfout);
  return 1;
}

int transmission_paddle(int setting, int targ){
  TString set_dummy;
  TString str_infile = str_infile_template;
  TString str_pid = str_pid_template;
  if(setting == 13){
    set_dummy="50Ca";
  }else{
    set_dummy="38Ca";
  }
  str_infile.ReplaceAll("SETTING",set_dummy);
  str_infile.ReplaceAll("TARGET", targetname.at(targ));
  str_pid.ReplaceAll("SET", Form("%i",setting));
  str_pid.ReplaceAll("TARGET", targetname.at(targ));
  std::cout<<"Input file: "<<str_infile<<std::endl<<"PID file: "<<str_pid<<std::endl;
  //
  v_pid = readCSV(str_pid);
  for (const auto& p : v_pid) {
    pidMap[p.getKey()] = p;
  }

  int Z = 20;
  auto *fin = TFile::Open(str_infile, "READ");
  auto *evt = (TTree*)fin->Get("evt");
  //auto *evt = (TTree*)fin->Get("Tree");
  auto nentries = evt->GetEntriesFast();
  cout<<"Entries: "<<nentries;

  auto Nmass = v_mass.size();
  vector<TGraphErrors*> &v_gtrans = vv_gtrans_pad.at(targ);
  vector<TGraphErrors*> &v_gtransint = vv_gtransint_pad.at(targ);
  vector<TGraphErrors*> &v_gcentroid = vv_gcentroid_pad.at(targ);
  vector<TGraphErrors*> &v_gsigma = vv_gsigma_pad.at(targ);
  auto mass_it = v_gtrans.size();
  vector<int> ntot(Nmass +  mass_it);
  vector<double> sqentot(Nmass +  mass_it);
  for(auto A:v_mass){
    auto *gtmp = new TGraphErrors();
    v_gtrans.push_back(gtmp);
    auto *gtmpint = new TGraphErrors();
    v_gtransint.push_back(gtmpint);
    auto *gtmp3 = new TGraphErrors();
    v_gcentroid.push_back(gtmp3);
    auto *gtmp2 = new TGraphErrors();
    v_gsigma.push_back(gtmp2);
  }
  int i_m = mass_it - 1;
  for(auto A:v_mass){
    i_m++;
    vector<double> v_trans;
    vector<double> v_etrans;
    auto gtmp = v_gtrans.at(i_m);
    auto gtmpint = v_gtransint.at(i_m);
    auto gtmp3 = v_gcentroid.at(i_m);
    auto gtmp2 = v_gsigma.at(i_m);
    //
    TString hname = Form("h_paddle_%i",A);
    auto* htmp = new TH2D(hname, hname, 28, 0, 28, 500, 1.8, 2.8);
    TCut condition = Form("PID_FRS_Z==%i && PID_FRS_A==%i", Z, A);
    auto nin = evt->Draw("", condition);
    //
    pid_gate::KeyType targetKey = std::make_tuple(Z,A,Z,A);
    // Check if the key exists in the map
    auto it = pidMap.find(targetKey);
    if (it == pidMap.end()) {
      continue;
    }
    // Key exists, so retrieve the pid object
    pid_gate PID = it->second;
    //float frszcentreValue = targetPid.getFrszCentre();
    std::cout << std::endl << "Target: "<< targetname[targ] << " fragZ centre value for frsz=" << PID.getFrsz() << ", frsa=" << PID.getFrsa()
	      << ", fragz=" << PID.getFragz() << ", fraga=" << PID.getFraga() << ": "
	      << PID.getFragaqCentre() <<"+/-" << PID.getFragaqSigma() << std::endl;
    //
    TCut twim_cond = Form("abs((Frag_Z-%f)/%f)<3.", PID.getFragzCentre(), PID.getFragzSigma());
    evt->Draw(Form("Frag_AoQ:Tofw_Paddle>>%s",hname.Data()),condition && twim_cond,"colz");
    c->SaveAs(pdfout);
    //
    auto draw_trans = [&](TH1D* htmp, TF1 *ftmp, int p, int nin, pid_gate PID) {
			//
			TLatex Tl;
			Tl.SetTextAlign(12);
			Tl.SetTextSize(0.04);
			Tl.DrawLatex(PID.getFragaqCentre(), htmp->GetMaximum(), Form("Paddle: %i",p));
			auto nout = ftmp->GetParameter(0) * 500.;//multiply by binning
			ntot.at(i_m) += nout;
			auto xpos = [&](double sigma) {return htmp->GetXaxis()->FindBin(ftmp->GetParameter(1)+sigma*ftmp->GetParameter(2));};
			auto nintegral = htmp->Integral(xpos(-3.),xpos(3.));
			cout<<"Integral bin range: "<<xpos(-3)<< " "<<xpos(3)<<" Integral:"<<nout<<endl;
			auto Enout = ftmp->GetParError(0) * 500.;
			sqentot.at(i_m) += Enout*Enout;
			auto trans = (nout/(double)nin)*100.;
			auto etrans = (Enout/(double)nin)*100;
			auto transint = (double)nintegral /(double)nin * 100.;
			auto etransint = transint / TMath::Sqrt(nintegral + nin);
			v_trans.push_back(trans);
			v_etrans.push_back(etrans);
			Tl.DrawLatex(PID.getFragaqCentre(), 0.75 * htmp->GetMaximum(), Form("Incoming: %i, Outgoing: %0.1f, Trans: %0.1f (%i)",(int)nin, nout, trans, (int)(10.*etrans)));
			gtmp->AddPoint(p, trans);
			gtmp->SetPointError(gtmp->GetN()-1, 0, etrans);
			gtmpint->AddPoint(p, transint);
			gtmpint->SetPointError(gtmpint->GetN()-1, 0, etransint);
			gtmp3->AddPoint(p, ftmp->GetParameter(1)-PID.getFragaqCentre());
			gtmp3->SetPointError(gtmp3->GetN()-1, 0, ftmp->GetParError(1));
			gtmp2->AddPoint(p, ftmp->GetParameter(2)/PID.getFragaqSigma());
			gtmp2->SetPointError(gtmp2->GetN()-1, 0, ftmp->GetParError(2)/PID.getFragaqSigma());
		      };
    //
    for(int p=0; p<28; p++){
      hname = Form("%s_paddle%i",htmp->GetName(),p);
      auto *hpaddle = htmp->ProjectionY(hname,p,p);
      //if(hpaddle->GetMaximum()<10) continue;
      hpaddle->SetTitle(hname);
      cout<<hname<<endl;
      auto *ftmp = new TF1(Form("f%i_p%i",A,p),"gausn",1.8,2.8);
      fithisto_paddle(hpaddle, ftmp, PID);
      draw_trans(hpaddle, ftmp, p, nin, PID);
      //c->SaveAs(pdfout);
      delete ftmp;
    }
    //
    gtmp->AddPoint(28, (double)ntot.at(i_m)/(double)nin *100.);
    gtmp->SetPointError(gtmp->GetN()-1, 0, sqrt(sqentot.at(i_m))/(double)nin *100.);
    gtmp->SetMarkerStyle(4); // circle
    gtmp->SetMarkerColor(targ+1);
    gtmp->SetMarkerSize(1);
    gtmp->Draw("apl");
    c->SaveAs(pdfout);
  }
  return 1;
}

void draw_div_evt(vector<int> &v_mass, vector<vector<TGraphErrors*>> &vv_gtrans_tmp, vector<TGraphErrors*> &v_gtrans_mass_tmp,vector<TGraphErrors*> &v_gtotal_mass_tmp, double range_min, double range_max){
  int i_m=-1;
  for(auto A:v_mass){
    i_m++;
    auto lg = TLegend(0.25,0.15);//TLegend(0.7*division, 70, 0.8*division, 85);
    for(int targ=0; targ<3; targ++){
      auto v_gtrans = vv_gtrans_tmp.at(targ);
      //cout<<"v_gtrans.size: "<<v_gtrans.size()<<endl;
      auto gtmp = v_gtrans.at(i_m);
      gtmp->SetLineColor(i_m+1);
      gtmp->SetMarkerStyle(4); // circle
      gtmp->SetMarkerColor(targ+1);
      gtmp->SetMarkerSize(1);
      gtmp->GetYaxis()->SetRangeUser(range_min,range_max);
      gtmp->Draw(targ==0?"apl":"pl");
      lg.AddEntry(gtmp, Form("Mass %i %s", A, targetname.at(targ).Data()), "pe");
      //
      Double_t x, y, ex=0., ey;
      gtmp->GetPoint(gtmp->GetN()-1, x, y);
      v_gtrans_mass_tmp.at(targ)->AddPoint(v_mass.at(i_m), y);
      ey = gtmp->GetErrorY(gtmp->GetN()-1);
      v_gtrans_mass_tmp.at(targ)->SetPointError(i_m, ex, ey);
      //
      if(targ<1) continue;
      Double_t CS = sigma_tot(v_gtrans_mass_tmp.at(targ-1)->GetPointY(i_m), y, targ==1?1e23:2e23);
      Double_t eCS = Esigma_tot(v_gtrans_mass_tmp.at(targ-1)->GetPointY(i_m), y, v_gtrans_mass_tmp.at(targ-1)->GetErrorY(i_m), ey, targ==1?1e23:2e23);
      v_gtotal_mass_tmp.at(targ-1)->AddPoint(v_mass.at(i_m), CS);
      v_gtotal_mass_tmp.at(targ-1)->SetPointError(i_m, ex, eCS);
    }
    lg.Draw();
    c->SaveAs(pdfout);
  }
  i_m=0;
  auto lg2 = TLegend(0.25,0.15);
  for(auto &gtmp: v_gtrans_mass_tmp){
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(45,105);
    gtmp->Draw(i_m==0?"apl":"pl");
    lg2.AddEntry(gtmp, Form("Transmission for %s", targetname.at(i_m).Data()), "pe");
    i_m++;
  }
  lg2.Draw();
  c->SaveAs(pdfout);
  //
  vector<double> v_residual(3);
  vector<double> v_sigma_syst(3);
  //vector<double> v_trans, v_etrans, v_systrans, v_dummy;
  TVectorD Tv_mass, v_trans, v_etrans[2], v_systrans, v_dummy;
  auto N = v_mass.size();
  Tv_mass.ResizeTo(N);
  v_trans.ResizeTo(N);
  v_etrans[0].ResizeTo(N);
  v_etrans[1].ResizeTo(N);
  v_systrans.ResizeTo(N);
  v_dummy.ResizeTo(N);
  i_m=0;
  auto lg3 = TLegend(0.25,0.15);
  for(auto &gtmp: v_gtotal_mass_tmp){
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(0,1800);
    gtmp->Draw(i_m==0?"apl":"pl");
    lg3.AddEntry(gtmp, Form("Total reaction cross section (a.u.) for %s", i_m==0?"Carbon":"Proton"), "pe");
    auto *ftmp = new TF1(Form("f%i",i_m), i_m==0?"[0]*pow(pow(x,1./3.)+pow(12,1./3.),2.)":"[0]*pow(pow(x,1./3.)+1.,2.)",38,52);
    gtmp->Fit(ftmp,"RL","");
    for(int n = 0; n<gtmp->GetN(); n++){
      //v_residual.at(i_m) +=
      auto sq_res = pow(gtmp->GetPointY(n) - ftmp->Eval(gtmp->GetPointX(n)) ,2.);
      auto sq_err = pow(gtmp->GetErrorY(n),2.);
      if(sq_err < sq_res) v_residual.at(i_m) += sq_res - sq_err;
    }
    v_sigma_syst.at(i_m) = sqrt(v_residual.at(i_m))/(double)gtmp->GetN();
    //
    i_m++;
  }
  lg3.Draw();
  c->SaveAs(pdfout);
  //
  // Calculate standard deviation for systematic errors
  auto RelSyst = 1.e23 * 1.e-27 * v_sigma_syst.at(0);
  auto* gtmp = v_gtrans_mass_tmp.at(0);
  for(int n = 0; n<gtmp->GetN(); n++){
    Tv_mass[n] = gtmp->GetPointX(n);
    v_trans[n] = (gtmp->GetPointY(n));
    v_etrans[0][n] = (gtmp->GetErrorY(n));
    v_etrans[1][n] = //gtmp->GetPointY(n) * RelSyst;
      (sqrt(pow(gtmp->GetPointY(n) * RelSyst,2.) + pow(gtmp->GetErrorY(n),2)));
    v_dummy[n] = (0);
  }
  auto *g_residual = new TGraphMultiErrors(2, Tv_mass, v_trans, v_dummy, v_dummy, v_etrans, v_etrans);
  g_residual ->SetLineColor(1);
  g_residual ->SetMarkerStyle(4);
  g_residual ->SetMarkerSize(1);
  g_residual ->GetAttFill(1)->SetFillColor(kBlue);
  g_residual ->GetAttFill(1)->SetFillStyle(3144);//3001
  g_residual ->GetYaxis()->SetRangeUser(45, 105);
  g_residual ->Draw("apsl; ;3 s=0.");
  g_residual ->SaveAs("g_residual.root");
  //
  i_m=0;
  auto lg4 = TLegend(0.25,0.15);
  for(auto &gtmp: v_gtrans_mass_tmp){
    //if(i_m==0) continue;
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(45,105);
    gtmp->Draw("pl");
    lg4.AddEntry(gtmp, Form("Transmission for %s", targetname.at(i_m).Data()), "pe");
    i_m++;
  }
  lg4.Draw();
  //v_gtrans_mass_tmp.at(1)->Draw("pl");
  //v_gtrans_mass_tmp.at(2)->Draw("pl");
  c->SaveAs(pdfout);
};

/*
void draw_paddle(vector<int> &v_mass){
  int i_m=-1;
  for(auto A:v_mass){
    i_m++;
    auto lg = TLegend(0.25,0.15);//TLegend(0.7*division, 70, 0.8*division, 85);
    for(int targ=0; targ<3; targ++){
      auto v_gtrans = vv_gtrans_pad.at(targ);
      //cout<<"v_gtrans.size: "<<v_gtrans.size()<<endl;
      auto gtmp = v_gtrans.at(i_m);
      gtmp->SetLineColor(i_m+1);
      gtmp->SetMarkerStyle(4); // circle
      gtmp->SetMarkerColor(targ+1);
      gtmp->SetMarkerSize(1);
      gtmp->GetYaxis()->SetRangeUser(0,105);
      gtmp->Draw(targ==0?"apl":"pl");
      lg.AddEntry(gtmp, Form("Mass %i %s", A, targetname.at(targ).Data()), "pe");
      //
      Double_t x, y, ex=0., ey;
      gtmp->GetPoint(gtmp->GetN()-1, x, y);
      v_gtrans_mass_pad.at(targ)->AddPoint(v_mass.at(i_m), y);
      ey = gtmp->GetErrorY(gtmp->GetN()-1);
      v_gtrans_mass_pad.at(targ)->SetPointError(i_m, ex, ey);
      //
      if(targ<1) continue;
      //Double_t CS = TMath::Exp(-y/v_gtrans_mass_pad.at(targ-1)->GetPointY(i_m));
      //Double_t eCS = CS * TMath::Sqrt(ey*ey/(y*y) + pow(v_gtrans_mass_pad.at(targ-1)->GetErrorY(i_m)/v_gtrans_mass_pad.at(targ-1)->GetPointY(i_m),2));
      Double_t CS = sigma_tot(v_gtrans_mass_pad.at(targ-1)->GetPointY(i_m), y, targ==1?1e23:2e23);
      Double_t eCS = Esigma_tot(v_gtrans_mass_pad.at(targ-1)->GetPointY(i_m), y, v_gtrans_mass_pad.at(targ-1)->GetErrorY(i_m), ey, targ==1?1e23:2e23);
      //Double_t eCS = CS * TMath::Sqrt(ey*ey/y*y);
      v_gtotal_mass_pad.at(targ-1)->AddPoint(v_mass.at(i_m), CS);
      v_gtotal_mass_pad.at(targ-1)->SetPointError(i_m, ex, eCS);
    }
    lg.Draw();
    c->SaveAs(pdfout);
  }
  i_m=0;
  auto lg2 = TLegend(0.25,0.15);
  for(auto &gtmp: v_gtrans_mass_pad){
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(45,105);
    gtmp->Draw(i_m==0?"apl":"pl");
    lg2.AddEntry(gtmp, Form("Transmission for %s", targetname.at(i_m).Data()), "pe");
    i_m++;
  }
  lg2.Draw();
  c->SaveAs(pdfout);
  //
  i_m=0;
  auto lg3 = TLegend(0.25,0.15);
  for(auto &gtmp: v_gtotal_mass_pad){
    gtmp->SetLineColor(i_m+1);
    gtmp->SetMarkerStyle(4); // circle
    gtmp->SetMarkerColor(i_m+1);
    gtmp->SetMarkerSize(1);
    gtmp->GetYaxis()->SetRangeUser(0,1800);
    gtmp->Draw(i_m==0?"apl":"pl");
    lg3.AddEntry(gtmp, Form("Total reaction cross section (a.u.) for %s", i_m==0?"Carbon":"2Protons"), "pe");
    auto *ftmp = new TF1(Form("f%i",i_m), "[0]*pow(x,2./3.)",38,52);
    gtmp->Fit(ftmp,"R","");
    i_m++;
  }
  lg3.Draw();
  c->SaveAs(pdfout);
};
*/
