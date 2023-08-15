#include "pid_gates.h"
TString dir = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/";
TString str_infile_template = dir + "rootfiles/mktree_2023/evttree_20Jul23_SETTING_TARGET.root";
//TString str_infile_template = dir + "rootfiles/mktree_2023/evttree_8Jun23_SETTING_TARGET.root";
//TString str_infile_template = dir + "rootfiles/mktree_2023/mktree_califa_May23_rolufix_SETTING_TARGET.root";
TString str_runlist = dir + "RunSummary.csv";
//TString str_pid0 = dir + "fragment/pid_gate_Settings12_empty.csv"; // Should be updated
//TString str_pid1 = dir + "fragment/pid_gate_Settings13_carbon.csv"; // Should be updated
TString str_pid_template = dir + "fragment/pid_gate_SettingsSET_TARGET.csv";
TString pdfout = "test_3sigma_ll_rolux5.pdf";
constexpr int division =10;
vector<int> v_mass0 = {39,40,41};
vector<int> v_mass1 = {47,48,49,50,51};
vector<int> v_mass;
vector<TString> targetname = {"empty","carbon","ch2"};
//
vector<pid_gate> v_pid;
TCanvas *c;
vector<vector<TGraphErrors*>> vv_gtrans(3);
vector<vector<TGraphErrors*>> vv_gtransint(3);
vector<vector<TGraphErrors*>> vv_gsigma(3);
vector<TGraphErrors*> v_gtrans_mass(3);
vector<TGraphErrors*> v_gtotal_mass(2);

void transmission(TString outfilename = "Transmission.csv");
int transmission_run(int setting = 13, int targ = 2);

auto sigma_tot(double trans1, double trans2, double ntarg) {
  return TMath::Log(trans1/trans2)/ntarg * 1e27;
}
auto Esigma_tot(double trans1, double trans2, double Etrans1, double Etrans2, double ntarg) {
  double e1=Etrans1/trans1;
  double e2=Etrans2/trans2;
  return TMath::Sqrt(e1*e1+e2*e2)/ntarg * 1e27;
}


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
  //
  v_mass = v_mass0;
  transmission_run(12,0);
  transmission_run(12,1);
  transmission_run(12,2);
  v_mass = v_mass1;
  transmission_run(13,0);
  transmission_run(13,1);
  transmission_run(13,2);
  
  //
  v_mass.insert(v_mass.begin(),v_mass0.begin(),v_mass0.end());
  int i_m=-1;
  for(auto A:v_mass){
    i_m++;
    auto lg = TLegend(0.25,0.15);//TLegend(0.7*division, 70, 0.8*division, 85);
    for(int targ=0; targ<3; targ++){
      auto v_gtrans = vv_gtrans.at(targ);
      //cout<<"v_gtrans.size: "<<v_gtrans.size()<<endl;
      auto gtmp = v_gtrans.at(i_m);
      gtmp->SetLineColor(i_m+1);
      gtmp->SetMarkerStyle(4); // circle
      gtmp->SetMarkerColor(targ+1);
      gtmp->SetMarkerSize(1);
      gtmp->GetYaxis()->SetRangeUser(45,105);
      gtmp->Draw(targ==0?"apl":"pl");
      lg.AddEntry(gtmp, Form("Mass %i %s", A, targetname.at(targ).Data()), "pe");
      //
      Double_t x, y, ex=0., ey;
      gtmp->GetPoint(gtmp->GetN()-1, x, y);
      v_gtrans_mass.at(targ)->AddPoint(v_mass.at(i_m), y);
      ey = gtmp->GetErrorY(gtmp->GetN()-1);
      v_gtrans_mass.at(targ)->SetPointError(i_m, ex, ey);
      //
      if(targ<1) continue;
      //Double_t CS = TMath::Exp(-y/v_gtrans_mass.at(targ-1)->GetPointY(i_m));
      //Double_t eCS = CS * TMath::Sqrt(ey*ey/(y*y) + pow(v_gtrans_mass.at(targ-1)->GetErrorY(i_m)/v_gtrans_mass.at(targ-1)->GetPointY(i_m),2));
      Double_t CS = sigma_tot(v_gtrans_mass.at(targ-1)->GetPointY(i_m), y, targ==1?1e23:2e23);
      Double_t eCS = Esigma_tot(v_gtrans_mass.at(targ-1)->GetPointY(i_m), y, v_gtrans_mass.at(targ-1)->GetErrorY(i_m), ey, targ==1?1e23:2e23);
      //Double_t eCS = CS * TMath::Sqrt(ey*ey/y*y);
      v_gtotal_mass.at(targ-1)->AddPoint(v_mass.at(i_m), CS);
      v_gtotal_mass.at(targ-1)->SetPointError(i_m, ex, eCS);
    }
    lg.Draw();
    c->SaveAs(pdfout);
  }
  i_m=0;
  auto lg2 = TLegend(0.25,0.15);
  for(auto &gtmp: v_gtrans_mass){
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
  for(auto &gtmp: v_gtotal_mass){
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
  c->SaveAs(pdfout + "]");
}

int transmission_run(int setting, int targ){
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
  
  TString dummy = TString(str_infile);

  v_pid = readCSV(str_pid);
  std::map<pid_gate::KeyType, pid_gate> pidMap;
  for (const auto& p : v_pid) {
    pidMap[p.getKey()] = p;
  }


  int Z = 20;
  auto *fin = TFile::Open(dummy, "READ");
  auto *evt = (TTree*)fin->Get("evt");
  //auto *evt = (TTree*)fin->Get("Tree");
  auto nentries = evt->GetEntriesFast();
  cout<<"Entries: "<<nentries;
  //auto* histo = new TH2F("in48","in48",400,1.8,2.8,400,12,22);

  auto Nmass = v_mass.size();
  //vector<vector<double>> vv_trans(Nmass);
  //vector<vector<double>> vv_etrans(Nmass);
  //vv_gtrans.at(targ) = v_gtrans;
  //vv_gsigma.at(targ) = v_gsigma;
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
  vector<TH1F*> v_htot(Nmass + mass_it);
  //for(int i=0; i<1; i++){
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
      auto htmp = new TH1F(hname,hname, 500,1.8,2.8);
      //evt->Draw("FrsData.fZ:FrsData.fAq>>h","","col",nplot,nfirst);
      TCut condition = Form("PID_FRS_Z==%i && PID_FRS_A==%i && ROLU_X>5", Z, A);
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
      //
      auto fithisto = [&](TH1F* htmp, int i, int nin) {
			auto *ftmp = new TF1(Form("f%i_%i",A,i),"gausn",1.8,2.8);
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
			gtmp->AddPoint(i, trans);
			gtmp->SetPointError(gtmp->GetN()-1, 0, etrans);
			gtmpint->AddPoint(i, transint);
			gtmpint->SetPointError(gtmpint->GetN()-1, 0, etransint);
			gtmp2->AddPoint(i, ftmp->GetParameter(2)/PID.getFragaqSigma());
			gtmp2->SetPointError(gtmp2->GetN()-1, 0, ftmp->GetParError(2)/PID.getFragaqSigma());
			delete ftmp;
		      };
      fithisto(htmp, i, nin);
      //c->SaveAs(pdfout);
      //
      if(i==0){ // division == 0;
	v_htot.at(i_m) = (TH1F*)htmp->Clone();
	v_htot.at(i_m)->SetTitle(Form("htotal_%i",A));
      }else{
	v_htot.at(i_m)->Add(htmp,1.);
      }
      //
      if(i == division - 1){
	fithisto(v_htot.at(i_m), division, nin_tot.at(i_m));
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

  //vv_gtrans.at(targ) = v_gtrans;
  //vv_gsigma.at(targ) = v_gsigma;
  return 1;
}
