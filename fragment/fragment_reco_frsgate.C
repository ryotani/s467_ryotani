#define MINPADDLE 14
#define NUMPADDLE 16

#include "fragmentana_twim.h"
int zet = 18, mass = 44;
double min_aoq = 1.8, max_aoq = 2.7;
//
bool IsEmpty=true;
TString targ = "empty";
TString infile = "./fragment/output/mktree_fragment_Dec_empty.root";
double beta_offset = 0.0;

//
/*
bool IsEmpty=false;
TString targ = "ch2";
TString infile = "./fragment/output/mktree_fragment_ch2-24mm.root";
//TString infile = "./rootfiles/rootfiletmp/fragment_Sep2021/s467_filltree_Setting13_0354_28Sep.root";//1hr CH2 run
//double beta_offset = 0.006;
double beta_offset = 0.000;
*/
//
/*
bool IsEmpty=false;
TString targ = "carbon";
TString infile = "./fragment/output/mktree_fragment_carbon.root";
double beta_offset = 0.000;
*/
//
/*
bool IsEmpty=false;
TString targ = "PP";
TString infile = "./fragment/output/mktree_fragment_PP.root";
double beta_offset = 0.000;
*/
//
TString outpdf = "./fragment/output/fragment_reco_frsgate_Dec_" + targ + Form("_z%i_a%i",zet,mass) + "_beta_offset_paddle.pdf";
TString outcsv = "./fragment/output/fragment_reco_frsgate_Dec_" + targ + Form("_z%i_a%i",zet,mass) + "_beta_offset_paddle.csv";
TString outroot = "./fragment/output/fragment_reco_frsgate_Dec_" + targ + Form("_z%i_a%i",zet,mass) + "_beta_offset_paddle.root";
TString recobrho_outpdf = "./fragment/output/fragment_reco_frsgate_brho_Dec_" + targ + Form("_z%i_a%i",zet,mass) + "_beta_offset_paddle.pdf";

ofstream fcsv(outcsv, ofstream::out);
TFile *fout = new TFile(outroot,"RECREATE");

//////////
void initialise();
//int tofw_calib();
int loadtofpara();
//int transfer_mat();
int draw_pidgate(TString conditions);
//int draw_transfer_corr(TString conditions);
//int draw_toflength_corr(TString conditions);
int brho_corr();
int draw_tofw();

int fragment_reco_frsgate(){
  if (NUMPADDLE>28) return 1;
  initialise();
  
  //tofw_calib();
  //loadtofpara();
  //transfer_mat();
  c = new TCanvas("c","c",3000,2500);
  c->Print(outpdf + "[");
  
  brho_corr();
  /*
  c = new TCanvas("c","c",3000,2500);
  c -> Divide(2,2);
  c->cd(1);
  h_brho_mw3[cond][2]->Draw("colz");
  //prof_brho_mw3[cond][2] ->Draw("same");
  c->cd(3);
  h_aoqaoq->Draw("colz");
  c->cd(2);
  h_pid->Draw("colz");
  c->cd(4)->SetLogy();
  h_fragaoq_proj->Draw();
  c -> Print(outpdf);
  delete c;
  / */
  draw_tofw();
  //
  c->Print(outpdf + "]");
  //
  /*
  h_brho_mw3[cond][2]->Write();
  h_aoqaoq->Write();
  h_pid->Write();
  h_fragaoq_proj->Write();
  */
  h_zaoq_mod[0]->Write();
  h_zaoq_mod[1]->Write();
  fout->Close();
  //
  //p->Close();
  delete ch;
  return 0;
}

int draw_tofw(){
  TString dummystring = "";
  c = new TCanvas("c","c",3000,2500);
  c -> Divide(6,5);
  //
  if(NUMPADDLE>28) return 1;
  c->cd(NUMPADDLE+1);
  h_paddleaoq = new TH2D("h_paddleaoq","Paddle ID vs AoQ; AoQ; PaddleID", 100, min_aoq, max_aoq, NUMPADDLE+1, -0.5, NUMPADDLE + 0.5);
  ch->Draw(Form("Tofw_Paddle:%s>>h_paddleaoq",fragaoqstring.Data()), "", "colz");
  c->cd(NUMPADDLE+2);
  h_paddleaoq_gated = new TH2D("h_paddleaoq_gated","Paddle ID vs AoQ (gated in FRS); AoQ; PaddleID", 100, min_aoq, max_aoq, NUMPADDLE+1, -0.5, (double)NUMPADDLE + 0.5);
  ch->Draw(Form("Tofw_Paddle:%s>>h_paddleaoq_gated",fragaoqstring.Data()), frspidgate, "colz");
  //
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd(i+1);
    h_zaoq_paddle[i] = new TH2D(Form("h_zaoq_paddle%i",i), Form("Z vs AoQ difference (CaveC - AoQ in gate) in Paddle%i (FRS gated); #Delta AoQ ; TwimZ",i+1), 500, -0.5, 0.5, 500, 10., 30.);
    dummystring = Form("Tofw_Paddle==%i&&",i+1);
    dummystring += frspidgate;
    ch -> Draw(Form("FragZ:%s-%f>>h_zaoq_paddle%i",fragaoqstring.Data(),(double)mass/(double)zet,i),dummystring, "colz");
  }
  c->Print(outpdf);
  //
  fcsv<<"Diff in AoQ";
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd(i+1);
    h_proj_deltaaoq[i] = h_zaoq_paddle[i]->ProjectionX(Form("_px%i",i),187,212);//Z=18
    TSpectrum *tspec = new TSpectrum(10);
    Int_t npeaks = tspec->Search(h_proj_deltaaoq[i]);
    Double_t *peakX = tspec->GetPositionX();
    Double_t x_tmp=-100;
    //for(int j=0; j<npeaks; j++){
    for(int j=0; j<min(4,npeaks); j++){
      cout<<j<<" "<<peakX[j]<<endl;
      if(abs(x_tmp) > abs(peakX[j])) x_tmp = peakX[j];
    }
    cout<<x_tmp<<" is selected"<<endl<<endl;
    /*
    f_aoq_paddle[i] = new TF1(Form("f_aoq_paddle%i",i), "gaus");
    f_aoq_paddle[i] -> SetParameter(1,x_tmp);
    f_aoq_paddle[i] -> SetLineColor(kBlue);
    h_proj_deltaaoq[i]->Fit(Form("f_aoq_paddle%i",i),"L","",x_tmp-0.05,x_tmp+0.05);
    aoq0[i] = f_aoq_paddle[i] -> GetParameter(1);*/
    aoq0[i] = x_tmp;
    fcsv<<", "<<aoq0[i];
    h_proj_deltaaoq[i]->Draw();
  }
  fcsv<<endl;
  c->cd(NUMPADDLE+1);
  h_beta_paddle[0] = new TH2D("h_beta_paddle0","Beta (before correction) in each paddle", 500, 0.5, 0.8, NUMPADDLE+1, -0.5, NUMPADDLE + 0.5);
  ch->Draw(Form("Tofw_Paddle:%s>>h_beta_paddle0",beta_tofw_mod.Data()), frspidgate, "colz");
  //
  c->Print(outpdf);
  //
  fcsv<<"Beta offset";
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd(i+1);
    h_beta_proj[i] = h_beta_paddle[0]->ProjectionX(Form("_px%i",i),i+2,i+2);
    /*
    f_beta_proj[i] = new TF1(Form("f_beta_proj%i",i),"gaus");
    h_beta_proj[i]->Fit(Form("f_beta_proj%i",i),"L","",0.7,0.8);
    ave_beta_org[i] = f_beta_proj[i]->GetParameter(1);
    */
    TSpectrum *tspec = new TSpectrum(5);
    Int_t npeaks = tspec->Search(h_beta_proj[i]);
    Double_t *peakX = tspec->GetPositionX();
    ave_beta_org[i] = peakX[0];
    double beta_offset_tmp = ((double)zet/(double)mass)/(1./ave_beta_org[i]+ave_beta_org[i]/sqrt(1.-pow(ave_beta_org[i],2.))) * aoq0[i];
    beta_tofw_mod +=  Form("+(Tofw_Paddle==%i)*(%f)",i+1,beta_offset_tmp);
    fcsv<<", "<<beta_offset_tmp;
  }
  fcsv<<endl;
  c->cd(NUMPADDLE+1);
  h_beta_paddle[1] = new TH2D("h_beta_paddle1","Beta (after correction) in each paddle", 500, 0.5, 0.8, NUMPADDLE+1, -0.5, NUMPADDLE + 0.5);
  ch->Draw(Form("Tofw_Paddle:%s>>h_beta_paddle1",beta_tofw_mod.Data()), frspidgate, "colz");
  //
  c->cd(NUMPADDLE+2);
  fragaoqstring = Form("(%s)*sqrt(1-(%s)*(%s))/((%s)*(%f))", fragbrhostring.Data(), beta_tofw_mod.Data(),beta_tofw_mod.Data(),beta_tofw_mod.Data(), mc_e);
  h_paddleaoq_gated->SetTitle("Paddle ID vs AoQ after correction (gated in FRS); AoQ; PaddleID");
  ch->Draw(Form("Tofw_Paddle:%s>>h_paddleaoq_gated",fragaoqstring.Data()), frspidgate, "colz");
  c->Print(outpdf);
  //
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd(i+1);
    h_aoqaoq_paddle[i] = new TH2D(Form("h_aoqaoq_paddle%i",i), Form("Reconstructed AoQ vs AoQ difference (CaveC - FRS), gated incoming in Paddle%i; #Delta AoQ; AoQ in Cave",i+1), 500, -1, 1, 500, min_aoq, max_aoq);
    dummystring = Form("Tofw_Paddle==%i",i+1);
    ch -> Draw(Form("%s:%s-FRSAoQ>>h_aoqaoq_paddle%i",fragaoqstring.Data(),fragaoqstring.Data(),i),dummystring, "colz");
  }
  c->Print(outpdf);
  //
  /* / y pos MW3 vs TofW tdiff
  c->cd(NUMPADDLE+1);
  h_mw3y_paddle = new TH2D("h_mw3y_paddle", "MW3_Y vs tdiff in paddles", NUMPADDLE+1, -0.5, NUMPADDLE +0.5, 500, -300, 300);
  ch->Draw("Mw3_Y:....");
  / */  //
  // x pos MW3 vs TofW tdiff
  c->cd(NUMPADDLE+1)->SetLogz();
  h_mw3x_paddle = new TH2D("h_mw3x_paddle", "MW3_X vs paddle id", NUMPADDLE+1, -0.5, NUMPADDLE +0.5, 500, -300, 300);
  ch->Draw("Mw3_X:Tofw_Paddle>>h_mw3x_paddle","","colz");
  //
  c->cd(NUMPADDLE+2)->SetLogz();
  h_mw3x_paddle_gated = new TH2D("h_mw3x_paddle_gated", "MW3_X vs paddle_gated id", NUMPADDLE+1, -0.5, NUMPADDLE +0.5, 500, -300, 300);
  ch->Draw("Mw3_X:Tofw_Paddle>>h_mw3x_paddle_gated",frspidgate,"colz");
  //
  // y pos MW3 vs AoQ (gated)
  for(int i = 0; i<NUMPADDLE; i++){
    c -> cd(i+1);
    h_mw3y_aoq[i] = new TH2D(Form("h_mw3y_aoq%i",i), Form("MW3_Y vs AoQ gated incoming and Z in Padle%i; AoQ in Cave; Y in Mwpc3",i+1), 500, min_aoq, max_aoq, 500, -200, 200);
    TString dummystring = Form("(Tofw_Paddle == %i)&& abs(FragZ-18)<0.5&&", i+1);
    dummystring += frspidgate;
    ch->Draw(Form("Mw3_Y:%s>>h_mw3y_aoq%i",fragaoqstring.Data(),i),dummystring,"colz");
  }
  //
  c->Print(outpdf);
  ////
  for(int i = MINPADDLE; i<NUMPADDLE; i++){
    c -> cd(i+1);
    h_zaoq_paddle_mod[i] = new TH2D(Form("h_zaoq_paddle_mod%i",i), Form("Z vs AoQ in Paddle%i (FRS gated); #Delta AoQ ; TwimZ",i+1), 500, min_aoq, max_aoq, 500, 10., 30.);
    dummystring = Form("Tofw_Paddle==%i&&",i+1);
    dummystring += frspidgate;
    //ch -> Draw(Form("FragZ:%s-FRSAoQ>>h_zaoq_paddle_mod%i",fragaoqstring.Data(),i),dummystring, "colz");
    ch -> Draw(Form("FragZ:%s>>h_zaoq_paddle_mod%i",fragaoqstring.Data(),i),dummystring, "colz");
  }
  //
  c->cd(NUMPADDLE+1);
  h_zaoq_mod[0] = new TH2D("h_zaoq_mod0","PID after beta correction",500,min_aoq,max_aoq,500,10.,30.);
  ch->Draw(Form("FragZ:%s>>h_zaoq_mod0",fragaoqstring.Data()),"","colz");
  //
  c->cd(NUMPADDLE+2);
  h_zaoq_mod[1] = new TH2D("h_zaoq_mod1","PID after beta correction (FRS gated)",500,min_aoq,max_aoq,500,10.,30.);
  ch->Draw(Form("FragZ:%s>>h_zaoq_mod1",fragaoqstring.Data()),frspidgate,"colz");
  //
  c->Print(outpdf);
  return 1;
}

int brho_corr(){
  //  c = new TCanvas("c","c",3000,2500);
  c -> Divide(4,4);
  
  TString conditions ="";
  TString dummystring ="";
  TString mw3_dummy[3] = {"Mw3_X", "Mw3_Y", "Corrected Mw3_X"};
  for(int cond =0; cond<2; cond++){
    if(cond==1)
      conditions += frspidgate;
    c->cd(1+cond*8);
    h_frspid_mw[cond] = new TH2D(Form("h_frspid%i",cond), "PID in FRS (S2-CaveC); AoQ; MusicZ", 500, min_aoq,max_aoq, 500,10,30);
    ch->Draw(Form("MusicZ:FRSAoQ>>h_frspid%i",cond),conditions,"col");
    for(int imw=0; imw<4; imw++){
      TH2F *h = new TH2F(Form("h%i",4+imw+cond*8),Form("MW%i; x[mm]; y[mm]",imw), 500,-200,200,500,-200,200);
      c->cd(5+imw+cond*8);
      ch->Draw(Form("Mw%i_Y:Mw%i_X>>h%i",imw,imw,4+imw+cond*8),conditions,"colz");
    }
    c->cd(2+cond*8);
    TH2F *musz = new TH2F(Form("hmusz%i",cond), "Z correlation; R3BMusic; Twim", 500,10,30,500,10,30);
    ch->Draw(Form("TwimZ:MusicZ>>hmusz%i",cond),conditions,"colz");
    c->cd(3+cond*8);
    TH2F *musangle = new TH2F(Form("musangle%i",cond), "Music Angles; R3BMusic; Twim", 500,-100,100,500,-100,100);
    ch->Draw(Form("TwimTheta*1000.:MusicTheta*1000.>>musangle%i",cond),conditions,"colz"); 
    /*
    TH2F *hmusic = new TH2F(Form("hmusic%i",cond),"R3BMusic; angle [mrad]; z", 500,-100,100,500,10,30);
    ch->Draw(Form("MusicZ:MusicTheta*1000.>>hmusic%i",cond),conditions,"colz");
    c->cd(3+cond*8);
    TH2F* htwim = new TH2F(Form("htwim%i",cond),"Twim; angle [mrad]; z", 500,-100,100,500,10,30);
    ch->Draw(Form("TwimZ:TwimTheta*1000.>>htwim%i",cond),conditions,"colz");
    */
    c->cd(4+cond*8);
    TH2F* hbeta = new TH2F(Form("hbeta%i",cond),"Beta correlations; FRS; fragment", 500,0.5,0.9,500,0.5,0.9);
    ch->Draw(Form("FragBeta:FRSBeta>>hbeta%i",cond),conditions,"colz");
  }
  c->Print(outpdf);
  c->Clear();
  c->Divide(4,4);
  //
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(1 + i%16);
    TH2F* htofw = new TH2F(Form("htofw%i",i),Form("Beta correlations in paddle%i; beta in FRS; ToF",i+1), 500,0.5,0.9,500,0,30);
    ch->Draw(Form("FragTof:FRSBeta>>htofw%i",i),Form("Tofw_Paddle==%i",i+1),"colz");
    if(i%16==15 || i==NUMPADDLE-1){
      c->Print(outpdf);
      c->Clear();
      c->Divide(4,4);
    }
  }
  //
  c->Clear();
  //delete c;
  //c = new TCanvas("c","c",3000,2500);
  c -> Divide(4,4);
  
  conditions += "&&";
  conditions += cut_mw;
  //conditions += Zgate;
  if(IsEmpty){ conditions += " (abs(FragBeta - FRSBeta)<0.004)";}
  else{ conditions +="1";}
  //for(int i=2; i<3; i++){
  for(int i=0; i<3; i++){
    c->cd(1+i);
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
      dummystring = Form("%s:FRSBrho>>h_brho_brho%i%i",Mw3_X_mod.Data(), cond, i);
      //continue;
    }
    ch -> Draw(dummystring,conditions,"col");
    ////
    c->cd(5+i);
    f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"[0]+[1]*x",9.0,9.12);
    //f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"(sqrt([1]*[1]-4.*[0]*([2]-x))-[1])/(2.*[0])",9.0,9.13);
    //f_brho_mw3[cond][i] = new TF1(Form("func_brho_mw3%i%i",cond,i),"[2]+[1]*x+[0]*x*x",9.0,9.13);
    f_brho_mw3[cond][i]->SetLineWidth(1);
    if(IsEmpty){
      prof_brho_mw3[cond][i] = h_brho_mw3[cond][i]->ProfileX();
      prof_brho_mw3[cond][i] ->Draw();
      prof_brho_mw3[cond][i] ->Fit(f_brho_mw3[cond][i], "R","");
    }
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
  TH2F* h_brho_corr = new TH2F("h_brho_corr", "h_brho_old_new", 500, 6, 9.3,500, 6, 9.3);
  if(IsEmpty)
    fragbrhostring = Form("(%s-(%f))/(%f)", Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(0), f_brho_mw3[cond][2]->GetParameter(1));
  //fragbrhostring = Form("(%f)*pow(%s,2.0)+(%f)*%s+(%f)", f_brho_mw3[cond][2]->GetParameter(0), Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(1), Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(2));
  //fragbrhostring = Form("(sqrt(pow(%f,2.0)-4.*%f*(%f-%s)-(%f))/(%f*2.0))", f_brho_mw3[cond][2]->GetParameter(1), f_brho_mw3[cond][2]->GetParameter(0),f_brho_mw3[cond][2]->GetParameter(2), Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(1),f_brho_mw3[cond][2]->GetParameter(0));
  cout << "fragbrhostring: "<< fragbrhostring <<endl;
  fcsv << "fragbrhostring, "<< fragbrhostring <<endl;
  ch -> Draw(Form("%s:FragBrho>>h_brho_corr", fragbrhostring.Data()) ,conditions, "colz");
  
  //
  c->cd(4);
  //conditions = "";
  //conditions = "abs(MusicZ-20.)<0.4 && abs(FRSAoQ-2.45)<0.02"; // 49Ca
  h_brhobrho = new TH2D("hbrhobrho", "Brho correlation in FRS and Cave; Brho FRS /Tm; Brho Cave /Tm", 500, 8.7, 9.3, 500, 8.7, 9.3);
  ch->Draw(Form("(%s):%s>>hbrhobrho", fragbrhostring.Data(), brho.Data()), conditions, "col");
  //
  c->cd(8);
  h_aoqaoq = new TH2D("haoqaoq", "Brho correlation in FRS and Cave; AoQ FRS ; AoQ Cave", 500, min_aoq, max_aoq, 500, min_aoq, max_aoq);
  dummystring = Form("%s:FRSAoQ>>haoqaoq",fragaoqstring.Data());
  cout << dummystring <<endl;
  ch->Draw(dummystring, conditions, "col");
  //
  c->cd(12);
  h_pid = new TH2D("hpid", "PID of fragment; AoQ Cave; TwimZ", 500, min_aoq, max_aoq, 500, 10, 30);
  ch->Draw(Form("FragZ:%s>>hpid",fragaoqstring.Data()), conditions, "colz");
  //
  c->cd(16);
  //h_fragaoq_proj = new TH1D(Form("h_fragaoq_proj%i",20), Form("AoQ projected with Z=%i",20), 500, min_aoq, max_aoq);
  h_fragaoq_proj = new TH1D("h_fragaoq_proj", "AoQ projected with Z=18", 500, min_aoq, max_aoq);
  //conditions += Form("&& abs(TwimZ-%f)<0.4",(Double_t)atom);
  // conditions += Form("&& abs(TwimZ-%f)<0.4",(Double_t)18);
  conditions += Form("&& abs(FragZ-%f)<0.4",(Double_t)18);
  //conditions + "&& abs(TwimZ-20.)<0.4"
  ch->Draw(Form("%s>>h_fragaoq_proj",fragaoqstring.Data()), conditions, "");
  /*
  for(int atom=20; atom<=20; atom++){
    h_fragaoq_proj[atom] = new TH1D(Form("h_fragaoq_proj%i",atom), Form("AoQ projected with Z=%i",atom), 500, min_aoq, max_aoq);
    //conditions += Form("&& abs(TwimZ-%f)<0.4",(Double_t)atom);
    ch->Draw(Form("%s>>h_fragaoq_proj%i",fragaoqstring.Data(),atom), conditions + Form("&& abs(TwimZ-%f)<0.4",(Double_t)atom), "");//atom==18?"":"same");
    //h_fragaoq_proj[atom]->SetLineColor(atom-17);
    //h_fragaoq_proj[atom]->GetXaxis()->SetRangeUser(min_aoq,max_aoq);
    }*/
  c -> Print(recobrho_outpdf);
  delete c;
  return 0;//cond++;
}

int draw_pidgate(TString conditions){
  cout<<"\033[1;31m Index: "<<cond<<", Draw: "<<conditions<<"\033[m"<<endl;
  c->cd(1);
  h_frspid_mw[cond] = new TH2D(Form("h_frspid%i",cond), "PID in FRS (S2-CaveC); AoQ; MusicZ", 500, min_aoq,max_aoq, 500,10,30);
  ch->Draw(Form("MusicZ:FRSAoQ>>h_frspid%i",cond),conditions,"col");
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
    TString condition_temp=Form("(Tofw_Paddle==%i && abs(%s - FRSBeta)<0.004)",i+1, fragbeta[i].Data());
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
  ch->Draw(Form("%s:FRSBeta>>hbetabeta_mw%i", beta_tofw.Data(), cond), conditionwithbetacut, "col");
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
  //ch -> SetProof();
  //
  frspidgate = Form("abs(MusicZ-%i)<0.4 && abs(FRSAoQ -%f)<%f",zet, (double)mass/(double)zet, 0.22/(double)zet);
  //
  //Mw3_X_mod = "(Mw3_X - (120.606710) - TwimTheta *(5753.779516))+ (- (4.469672) - Mw2_X *(0.433325))+ (- (4.917787) - (Mw2_Y-Mw1_Y) *(-0.131103))+ (- (18.555400) - Mw1_Y *(0.542263))";
  //Sep2021 fit
  //Mw3_X_mod = "(Mw3_X - (120.342308) - TwimTheta *(5685.146295))+ (- (4.625269) - (Mw1_X+Mw2_X)/2. *(0.571696))+ (- (1.241102) - (Mw3_Y-Mw1_Y) *(-0.056444))+ (- (10.569345) - (Mw1_Y+Mw2_Y)/2. *(0.764638))";
  // Test 2021
  //Mw3_X_mod = "(Mw3_X - (120.342308) - TwimTheta *(5685.146295))+ (- (4.625269) - (Mw1_X+Mw2_X)/2. *(0.571696))+ (- (1.241102))";// - (Mw3_Y-Mw1_Y) *(-0.056444))+ (- (10.569345) - (Mw1_Y+Mw2_Y)/2. *(0.764638))";
  //
  // 1 Dec new brho in FRS
  //Mw3_X_mod = "(Mw3_X - (122.987775) - TwimTheta *(6092.815790))+ (- (2.039319) - (Mw1_X+Mw2_X)/2. *(0.306839))";//+ (- (1.011031) - (Mw3_Y-Mw1_Y) *(-0.051042))+ (- (10.434147) - (Mw1_Y+Mw2_Y)/2. *(0.744019))
  // 3 Dec
  Mw3_X_mod = "(Mw3_X - (122.993625) - TwimTheta *(6093.728091))+ (- (2.049605) - (Mw1_X+Mw2_X)/2. *(0.307520))";//+ (- (1.004413) - (Mw3_Y-Mw1_Y) *(-0.050580))+ (- (10.419900) - (Mw1_Y+Mw2_Y)/2. *(0.743001))
  fragbrhostring = "((" + Mw3_X_mod + "+1282.411556)/141.549357)"; 
    //Form"((%s-(%f))/%f)",Mw3_X_mod.Data(),-1282.411556,141.549357);
    //  "((Mw3_X - (120.342308) - TwimTheta *(5685.146295))+ (- (4.625269) - (Mw1_X+Mw2_X)/2. *(0.571696))+ (- (1.241102) - (Mw3_Y-Mw1_Y) *(-0.056444))+ (- (10.569345) - (Mw1_Y+Mw2_Y)/2. *(0.764638))-(-1282.411556))/(141.549357)";
  ///
  brho = "FRSBrho";
  beta_tofw_mod = Form("(FragBeta - %f)",beta_offset);
  //fragaoqstring = "FragAoQ";
  fragaoqstring = Form("(%s)*sqrt(1-%s*%s)/((%s)*(%f))", fragbrhostring.Data(), beta_tofw_mod.Data(),beta_tofw_mod.Data(),beta_tofw_mod.Data(), mc_e);
  //
  cut_mw = Form("(Mw1_X>%f)&&(Mw1_X<%f)&&", -30.,25.);
  for(int i=0;i<(IsEmpty?4:3);i++) cut_mw += Form("(%s>%f)&&(%s<%f)&&",axis_mw3[i].Data(),range_cut_mw3_low[i],axis_mw3[i].Data(),range_cut_mw3_high[i]);

}
