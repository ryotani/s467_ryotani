#define NUMPADDLE 27

#include "fragmentana.h"
// *
TString infile ="s467_FRSTree_Setting13_0340_FragmentTree.root" ;// "./tofw/output/mktree_tofw_frs_empty.root";
TString rootfiledir ="~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfiletmp/TofW/";
//TString outpdf = "./tofw/output/fragment_draw_run.pdf";
//TString brho_outpdf = "./tofw/output/fragment_draw_brho.pdf";
TString recobrho_outpdf = "./tofw/output/fragment_reco_brho_draw_run.pdf";
/* /
TString infile = "./tofw/output/mktree_tofw_frs_ch2-24mm.root";
TString outpdf = "./tofw/output/fragment_draw_ch2.pdf";
TString brho_outpdf = "./tofw/output/fragment_draw_brho_ch2.pdf";
TString recobrho_outpdf = "./tofw/output/fragment_reco_brho_ch2.pdf";
//TString recobrho_outpdf = "./tofw/output/fragment_reco_brho_ch2_49Ca.pdf";
*/

void initialise();
int tofw_calib();
int loadtofpara();
int transfer_mat();
int draw_pidgate(TString conditions);
int draw_transfer_corr(TString conditions);
int draw_toflength_corr(TString conditions);
int brho_corr(Int_t runnum);

int fragment_draw_run(){
  if (NUMPADDLE>28) return 1;
  p->SetRealTimeLog(kFALSE);
  p->SetProgressDialog(kFALSE);
  //
  std::ifstream RunList("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/RunSummary.csv", std::ios::in);
  if(!RunList.is_open()) std::cerr <<"No run summary found\n";
  int runnumcsv[500], targetpos[500], musicgain[500], junk[500];
  int FRSsetting[500]; // calib:0, ToFCalib:6-8, 40Ca:9, 39Ca:10, 38Ca:11,12, 50Ca:13, ToFWcalib:14
  string dummyline;
  char dumchar;
  double brhocsv[500];
  std::getline (RunList, dummyline);
  //std::cout << dummyline << std::endl;
  Int_t i=0;
  while(true){
    RunList>>runnumcsv[i]>>dumchar>>FRSsetting[i]>>dumchar>>brhocsv[i]>>dumchar>>targetpos[i]>>dumchar>>musicgain[i]>>dumchar>>junk[i];
    //std::cout<<runnumcsv[i]<<dumchar<<FRSsetting[i]<<dumchar<<brhocsv[i]<<dumchar<<targetpos[i]<<dumchar<<musicgain[i]<<dumchar<<junk[i]<<std::endl;
    if(i > 400 || !RunList.good()){
      break;
    }
    i++;
  }
  initialise();
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(4,4);
  c -> Print(recobrho_outpdf + "[");
  for(int j=0; j<i; j++){
    if(FRSsetting[j] !=13) continue;
    if(junk[j] == 1) continue;
    //if(targetpos[j]!=1424) continue;
    infile = Form("s467_FRSTree_Setting13_%04d_FragmentTree.root", runnumcsv[j]);
    cout << "Draw: "<<rootfiledir + infile <<endl;
    ch = new TChain("Tree");
    ch -> Add(rootfiledir + infile);
    //ch -> SetProof(kTRUE,kTRUE);
    brho_corr(runnumcsv[j]);
    //
    //p->ClearCache();
    //p->ClearData();
    //p->ClearDataSetCache();
    //ch->SetProof(kFALSE);
    cout<<"Ana end run: "<<runnumcsv[j]<<endl<<endl;
    ch->Reset();
    delete ch;
    //break;
  }
  //
  c -> Print(recobrho_outpdf + "]");
  p->Close();
  //delete ch;
  return 0;
}

int brho_corr(Int_t runnum){
  TString conditions ="";
  TString dummystring ="";
  TString mw3_dummy[3] = {"Mw3_X", "Mw3_Y", "Corrected Mw3_X"};
  brho = "Brho_S2_Cave";
  beta_tofw_mod = "FragBeta";
  for(int i=0; i<3; i++){
    c->cd(1+i);
    conditions = Zgate;
    if(i<2){
      h_brho_mw3[cond][i] = new TH2D(Form("h_brho_mw3%i%i",cond,i),Form("Brho MW3 correlation;%s /mm; Brho in FRS /Tm",mw3_dummy[i].Data()), 500, 8.7, 9.3, 500, -500, 500);
    }else{
      h_brho_mw3[cond][2] = new TH2D(Form("h_brho_brho%i%i",cond, i),Form("Brho-Brho correlation;FRS /Tm; Brho in FRS /Tm"), 500, 8.7, 9.3, 500, 8.7, 9.3);
    }
    if(i==0){
      dummystring = Form("Mw3_X:%s>>h_brho_mw3%i%i",brho.Data(),cond,i);
    }else if(i==1){
      dummystring = Form("Mw3_Y:%s>>h_brho_mw3%i%i",brho.Data(),cond,i);
    }else{
      //dummystring = Form(":>>h_brho_mw3%i%i", Mw3_X_mod.Data(),brho.Data(), cond, i);
      dummystring = Form("Brho_S2_Cave:FragBrho>>h_brho_brho%i%i",cond, i);
      continue;
    }
    ch -> Draw(dummystring,conditions,"col");
    ////
    c->cd(5+i);
    f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"[0]+[1]*x",9.03,9.1);
    f_brho_mw3[cond][i]->SetLineWidth(1);
    prof_brho_mw3[cond][i] = h_brho_mw3[cond][i]->ProfileX();
    prof_brho_mw3[cond][i] ->Draw();
    prof_brho_mw3[cond][i] ->Fit(f_brho_mw3[cond][i], "R","");
    ////
    c->cd(9+i);
    conditions = beamcondition + "1";
    h_beta_mw3[cond][i] = new TH2D(Form("h_beta_mw3%i%i",cond,i),Form("Beta MW3 correlation (Run %i);%s /mm; #Delta #beta",runnum,mw3_dummy[i].Data()), 500, -0.01,0.01, 500, -500, 500);
    if(i==0){
      dummystring = Form("Mw3_X:%s>>h_beta_mw3%i%i",beta_tofw_mod.Data(),cond,i);
    }else if(i==1){
      dummystring = Form("Mw3_Y:%s>>h_beta_mw3%i%i",beta_tofw_mod.Data(),cond,i);
    }else{
      //dummystring = Form("%s:%s>>h_beta_mw3%i%i", Mw3_X_mod.Data(),beta_tofw_mod.Data(), cond, i);
    }
    ch -> Draw(dummystring,conditions,"col");
    ////
    c->cd(13+i);
    prof_beta_mw3[cond][i] = h_beta_mw3[cond][i]->ProfileX();
    prof_beta_mw3[cond][i] ->Draw();
    ////
  }
  

  //
  c->cd(4);
  conditions = "";
  //conditions = "abs(MusicZ-20.)<0.4 && abs(AoQ_S2_Cave-2.45)<0.02"; // 49Ca
  h_brhobrho = new TH2D("hbrhobrho", "Brho correlation in FRS and Cave; Brho FRS /Tm; Brho Cave /Tm", 500, 8.7, 9.3, 500, 8.7, 9.3);
  fragbrhostring = "FragBrho";//Form("(%s-(%f))/(%f)", Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(0), f_brho_mw3[cond][2]->GetParameter(1));
  ch->Draw(Form("(%s):%s>>hbrhobrho", fragbrhostring.Data(), brho.Data()), conditions, "col");
  //
  c->cd(8);
  h_aoqaoq = new TH2D("haoqaoq", "Brho correlation in FRS and Cave; AoQ FRS ; AoQ Cave", 500, 2.2, 2.7, 500, 2.2, 2.7);
  fragaoqstring = "FragAoQ"; //Form("(%s)*sqrt(1-%s*%s)/((%s)*(%f))", fragbrhostring.Data(), beta_tofw_mod.Data(),beta_tofw_mod.Data(),beta_tofw_mod.Data(), mc_e);
  dummystring = Form("%s:AoQ_S2_Cave>>haoqaoq",fragaoqstring.Data());
  cout << dummystring <<endl;
  ch->Draw(dummystring, conditions, "col");
  //
  c->cd(12);
  h_pid = new TH2D("hpid", "PID of fragment; AoQ Cave; TwimZ", 500, 2.2, 2.7, 500, 10, 30);
  ch->Draw(Form("TwimZ:%s>>hpid",fragaoqstring.Data()), conditions, "col");
  //
  c -> Print(recobrho_outpdf);
  //delete c;
  for(int i=0; i<3; i++){
    delete h_brho_mw3[cond][i];
    delete f_brho_mw3[cond][i];
    delete h_beta_mw3[cond][i];
  }
  delete h_brhobrho;
  delete h_aoqaoq;
  delete h_pid;
  //
  return 0;//cond++;
}
void initialise(){
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleFont(42,"T");

  /*
  gStyle->SetLabelSize(25,"XYZ");
  gStyle->SetTitleSize(25,"XYZ");
  */
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");

  gStyle->SetTitleOffset(.9,"X");
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetLabelOffset(0.01,"X");
  gStyle->SetLabelOffset(0.01,"Y");
}
