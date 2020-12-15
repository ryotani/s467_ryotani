#define NUMPADDLE 27

#include "fragmentana.h"
// *
TString infile = "./tofw/output/mktree_tofw_frs_empty.root";
TString outpdf = "./tofw/output/fragment_reco.pdf";
TString brho_outpdf = "./tofw/output/fragment_reco_brho.pdf";
TString recobrho_outpdf = "./tofw/output/fragment_reco_brho_draw.pdf";
/* /
TString infile = "./tofw/output/mktree_tofw_frs_ch2-24mm.root";
TString outpdf = "./tofw/output/fragment_reco_ch2.pdf";
TString brho_outpdf = "./tofw/output/fragment_reco_brho_ch2.pdf";
TString recobrho_outpdf = "./tofw/output/fragment_reco_brho_ch2.pdf";
//TString recobrho_outpdf = "./tofw/output/fragment_reco_brho_ch2_49Ca.pdf";
*/

void initialise();
//int tofw_calib();
int loadtofpara();
//int transfer_mat();
int draw_pidgate(TString conditions);
//int draw_transfer_corr(TString conditions);
//int draw_toflength_corr(TString conditions);
int brho_corr();

int fragment_reco(){
  if (NUMPADDLE>28) return 1;
  initialise();
  //tofw_calib();
  //loadtofpara();
  //transfer_mat();

  brho_corr();
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(2,2);
  c->cd(1);
  h_brho_mw3[cond][2]->Draw("colz");
  prof_brho_mw3[cond][2] ->Draw("same");
  c->cd(3);
  h_aoqaoq->Draw("colz");
  c->cd(2);
  h_pid->Draw("colz");
  c->cd(4);
  h_fragaoq_proj->Draw();
  c -> Print("./tofw/output/fragment_reco_r3bmeeting20.pdf");
  delete c;
  //
  p->Close();
  delete ch;
  return 0;
}

int brho_corr(){
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(4,4);
  Mw3_X_mod = "(Mw3_X - (108.147723) - (Mw2_X-Mw1_X) *(7.365536))+ (- (6.988648) - Mw2_X *(0.800182))+ (- (-13.404625) - (Mw2_Y-Mw1_Y) *(-0.331842))+ (- (-0.308584) - Mw1_Y *(0.052478))";
  cut_mw = Form("(Mw1_X>%f)&&(Mw1_X<%f)&&", -30.,25.);
  for(int i=0;i<4;i++) cut_mw += Form("(%s>%f)&&(%s<%f)&&",axis_mw3[i].Data(),range_cut_mw3_low[i],axis_mw3[i].Data(),range_cut_mw3_high[i]);
  
  TString conditions ="";
  TString dummystring ="";
  TString mw3_dummy[3] = {"Mw3_X", "Mw3_Y", "Corrected Mw3_X"};
  brho = "Brho_S2_Cave";
  beta_tofw_mod = "FragBeta";
  for(int i=2; i<3; i++){
    c->cd(1+i);
    conditions = cut_mw;
    conditions += Zgate;
    conditions += "&& (abs(FragBeta - Beta_S2_Cave)<0.004)";
    if(i<2){
      h_brho_mw3[cond][i] = new TH2D(Form("h_brho_mw3%i%i",cond,i),Form("Brho MW3 correlation;%s /mm; Brho in FRS /Tm",mw3_dummy[i].Data()), 500, 8.7, 9.3, 500, -500, 500);
    }else{
      h_brho_mw3[cond][2] = new TH2D(Form("h_brho_brho%i%i",cond, i),Form("Brho in FRS vs Subtracted X position of MWPC3;Brho in FRS /Tm; Modified X position of MWPC3 /mm"), 500, 8.7, 9.3, 400, -200, 200);
    }
    if(i==0){
      dummystring = Form("Mw3_X:%s>>h_brho_mw3%i%i",brho.Data(),cond,i);
    }else if(i==1){
      dummystring = Form("Mw3_Y:%s>>h_brho_mw3%i%i",brho.Data(),cond,i);
    }else{
      //dummystring = Form(":>>h_brho_mw3%i%i", Mw3_X_mod.Data(),brho.Data(), cond, i);
      //dummystring = Form("Brho_S2_Cave:FragBrho>>h_brho_brho%i%i",cond, i);
      dummystring = Form("%s:Brho_S2_Cave>>h_brho_brho%i%i",Mw3_X_mod.Data(), cond, i);
      //continue;
    }
    ch -> Draw(dummystring,conditions,"col");
    ////
    c->cd(5+i);
    f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"[0]+[1]*x",9.0,9.12);
    //f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"(sqrt([1]*[1]-4.*[0]*([2]-x))-[1])/(2.*[0])",9.0,9.13);
    //f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"[2]+[1]*x+[0]*x*x",9.0,9.13);
    f_brho_mw3[cond][i]->SetLineWidth(1);
    prof_brho_mw3[cond][i] = h_brho_mw3[cond][i]->ProfileX();
    prof_brho_mw3[cond][i] ->Draw();
    prof_brho_mw3[cond][i] ->Fit(f_brho_mw3[cond][i], "R","");
    ////
    c->cd(9+i);
    //conditions = beamcondition;
    h_beta_mw3[cond][i] = new TH2D(Form("h_beta_mw3%i%i",cond,i),Form("Beta MW3 correlation;%s /mm; #Delta #beta",mw3_dummy[i].Data()), 500, -0.01,0.01, 500, -500, 500);
    if(i==0){
      dummystring = Form("Mw3_X:%s>>h_beta_mw3%i%i",beta_tofw_mod.Data(),cond,i);
    }else if(i==1){
      dummystring = Form("Mw3_Y:%s>>h_beta_mw3%i%i",beta_tofw_mod.Data(),cond,i);
    }else{
      //dummystring = Form("%s:%s>>h_beta_mw3%i%i", Mw3_X_mod.Data(),beta_tofw_mod.Data(), cond, i);
      continue;
    }
    ch -> Draw(dummystring,conditions,"col");
    ////
    c->cd(13+i);
    prof_beta_mw3[cond][i] = h_beta_mw3[cond][i]->ProfileX();
    prof_beta_mw3[cond][i] ->Draw();
    ////
  }

  c->cd(11);
  TH2F* h_brho_corr = new TH2F("h_brho_corr", "h_brho_old_new", 500, 8.7, 9.3,500, 8.7, 9.3);
  fragbrhostring = Form("(%s-(%f))/(%f)", Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(0), f_brho_mw3[cond][2]->GetParameter(1));
  //fragbrhostring = Form("(%f)*pow(%s,2.0)+(%f)*%s+(%f)", f_brho_mw3[cond][2]->GetParameter(0), Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(1), Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(2));
  //fragbrhostring = Form("(sqrt(pow(%f,2.0)-4.*%f*(%f-%s)-(%f))/(%f*2.0))", f_brho_mw3[cond][2]->GetParameter(1), f_brho_mw3[cond][2]->GetParameter(0),f_brho_mw3[cond][2]->GetParameter(2), Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(1),f_brho_mw3[cond][2]->GetParameter(0));
  cout << "fragbrhostring: "<< fragbrhostring <<endl;
  ch -> Draw(Form("%s:FragBrho>>h_brho_corr", fragbrhostring.Data()) ,conditions, "colz");
  
  //
  c->cd(4);
  //conditions = "";
  //conditions = "abs(MusicZ-20.)<0.4 && abs(AoQ_S2_Cave-2.45)<0.02"; // 49Ca
  h_brhobrho = new TH2D("hbrhobrho", "Brho correlation in FRS and Cave; Brho FRS /Tm; Brho Cave /Tm", 500, 8.7, 9.3, 500, 8.7, 9.3);
  //fragbrhostring = "FragBrho";//Form("(%s-(%f))/(%f)", Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(0), f_brho_mw3[cond][2]->GetParameter(1));
  ch->Draw(Form("(%s):%s>>hbrhobrho", fragbrhostring.Data(), brho.Data()), conditions, "col");
  //
  c->cd(8);
  h_aoqaoq = new TH2D("haoqaoq", "Brho correlation in FRS and Cave; AoQ FRS ; AoQ Cave", 500, 2.2, 2.7, 500, 2.2, 2.7);
  //fragaoqstring = "FragAoQ";
  fragaoqstring = Form("(%s)*sqrt(1-%s*%s)/((%s)*(%f))", fragbrhostring.Data(), beta_tofw_mod.Data(),beta_tofw_mod.Data(),beta_tofw_mod.Data(), mc_e);
  dummystring = Form("%s:AoQ_S2_Cave>>haoqaoq",fragaoqstring.Data());
  cout << dummystring <<endl;
  ch->Draw(dummystring, conditions, "col");
  //
  c->cd(12);
  h_pid = new TH2D("hpid", "PID of fragment; AoQ Cave; TwimZ", 500, 2.2, 2.7, 500, 10, 30);
  ch->Draw(Form("TwimZ:%s>>hpid",fragaoqstring.Data()), conditions, "colz");
  //
  c->cd(16);
  //h_fragaoq_proj = new TH1D(Form("h_fragaoq_proj%i",20), Form("AoQ projected with Z=%i",20), 500, 2.2, 2.7);
  h_fragaoq_proj = new TH1D("h_fragaoq_proj", "AoQ projected with Z=18", 500, 2.2, 2.7);
  //conditions += Form("&& abs(TwimZ-%f)<0.4",(Double_t)atom);
  conditions += Form("&& abs(TwimZ-%f)<0.4",(Double_t)18);
  //conditions + "&& abs(TwimZ-20.)<0.4"
  ch->Draw(Form("%s>>h_fragaoq_proj",fragaoqstring.Data()), conditions, "");
  /*
  for(int atom=20; atom<=20; atom++){
    h_fragaoq_proj[atom] = new TH1D(Form("h_fragaoq_proj%i",atom), Form("AoQ projected with Z=%i",atom), 500, 2.2, 2.7);
    //conditions += Form("&& abs(TwimZ-%f)<0.4",(Double_t)atom);
    ch->Draw(Form("%s>>h_fragaoq_proj%i",fragaoqstring.Data(),atom), conditions + Form("&& abs(TwimZ-%f)<0.4",(Double_t)atom), "");//atom==18?"":"same");
    //h_fragaoq_proj[atom]->SetLineColor(atom-17);
    //h_fragaoq_proj[atom]->GetXaxis()->SetRangeUser(2.2,2.7);
    }*/
  c -> Print(recobrho_outpdf);
  delete c;
  return 0;//cond++;
}

int draw_pidgate(TString conditions){
  cout<<"\033[1;31m Index: "<<cond<<", Draw: "<<conditions<<"\033[m"<<endl;
  c->cd(1);
  h_frspid_mw[cond] = new TH2D(Form("h_frspid%i",cond), "PID in FRS (S2-CaveC); AoQ; MusicZ", 500, 2.2,2.7, 500,10,30);
  ch->Draw(Form("MusicZ:AoQ_S2_Cave>>h_frspid%i",cond),conditions,"col");
  //
  c->cd(2);
  h_music_twim_mw[cond] = new TH2D(Form("hmusictwim_mw%i",cond),"R3BMusic and Twim with #beta in FRS;TwimZ;MusicZ",200,10,30,200,10,30);
  ch->Draw(Form("MusicZ:TwimZ>>hmusictwim_mw%i",cond),conditions,"col");
  //
  c->cd(3);
  conditionwithbetacut = conditions + "&&(";

  h_beta_beta_mw[cond][NUMPADDLE] = new TH2D(Form("hbetabeta_mw%i",cond), "#beta correlation: Frs vs CaveC; #beta in FRS; #beta in Cave (ns)", 500, 0.7, 0.8, 500, 0.7, 0.8);
  //
  for(int i = 0 ; i<NUMPADDLE; i++){
    TString condition_temp=Form("(Tofw_Paddle==%i && abs(%s - Beta_S2_Cave)<0.004)",i+1, fragbeta[i].Data());
    conditionwithbetacut += condition_temp;
    condition_temp += "&&";
    condition_temp += conditions;
    if(i==NUMPADDLE-1){
      conditionwithbetacut += ")";
    }else{
      conditionwithbetacut += "||";
    }
  }
  cout << "Beta Tofw: "<<beta_tofw<<endl<<" conditionwithbetacut: "<< conditionwithbetacut << endl;
  ch->Draw(Form("%s:Beta_S2_Cave>>hbetabeta_mw%i", beta_tofw.Data(), cond), conditionwithbetacut, "col");
  //
  c->cd(4);
  //
  //cout << endl << conditionwithbetacut << endl;
  return 0;
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
 

  ch = new TChain("Tree");
  ch -> Add(infile);
  ch -> SetProof();
}
