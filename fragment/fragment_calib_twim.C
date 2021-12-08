#define NUMPADDLE 28

#include "fragmentana_twim.h"

TString infile = "./fragment/output/mktree_fragment_Dec_empty.root";
//TString infile = "./rootfiles/rootfiletmp/MUSIC/s467_filltree_Setting13_0340_22Sep.root";
TString outpdf = "./fragment/output/fragment_calib_twim_Dec.pdf";
TString brho_outpdf = "./fragment/output/fragment_calib_twim_brho_Dec.pdf";
TString recobrho_outpdf = "./fragment/output/fragment_reco_twim_brho_Dec.pdf";
/*
TString infile = "./fragment/output/mktree_tofw_frs_ch2-24mm.root";
TString outpdf = "./fragment/output/fragment_calib_twim_ch2_test.pdf";
TString brho_outpdf = "./fragment/output/fragment_calib_twim_brho_ch2_test.pdf";
TString recobrho_outpdf = "./fragment/output/fragment_reco_brho_ch2_test.pdf";
*/

void initialise();
int tofw_calib();
int loadtofpara();
int transfer_mat();
int draw_pidgate(TString conditions);
int draw_transfer_corr(TString conditions);
int draw_toflength_corr(TString conditions);
int brho_corr();

int fragment_calib_twim(){
  if (NUMPADDLE>28) return 1;
  initialise();
  //tofw_calib();
  //loadtofpara();
  beta_tofw = "FragBeta";
  cut_mw = Form("(Mw1_X>%f)&&(Mw1_X<%f)&&", -30.,25.);
  for(int i=0;i<4;i++) cut_mw += Form("(%s>%f)&&(%s<%f)&&",axis_mw3[i].Data(),range_cut_mw3_low[i],axis_mw3[i].Data(),range_cut_mw3_high[i]);
  //cut_mw += Form("(Mw2_X>%f)&&(Mw2_X<%f)&&", -40, 40);
  transfer_mat();
  //brho_corr(); // To be done fragment_reco...
  //
  //p->Close();
  delete ch;
  return 0;
}

int brho_corr(){
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(4,4);
  TString conditions ="";
  TString dummystring ="";
  TString mw3_dummy[3] = {"Mw3_X", "Mw3_Y", "Corrected Mw3_X"};
  for(int i=0; i<3; i++){
    c->cd(1+i);
    conditions = Zgate;
    h_brho_mw3[cond][i] = new TH2D(Form("h_brho_mw3%i%i",cond,i),Form("Brho MW3 correlation;%s /mm; Brho in FRS /Tm",mw3_dummy[i].Data()), 500, 8.7, 9.3, 500, -500, 500);
    if(i==0){
      dummystring = Form("Mw3_X:%s>>h_brho_mw3%i%i",brho.Data(),cond,i);
    }else if(i==1){
      dummystring = Form("Mw3_Y:%s>>h_brho_mw3%i%i",brho.Data(),cond,i);
    }else{
      dummystring = Form("%s:%s>>h_brho_mw3%i%i", Mw3_X_mod.Data(),brho.Data(), cond, i);
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
    conditions = beamcondition;
    h_beta_mw3[cond][i] = new TH2D(Form("h_beta_mw3%i%i",cond,i),Form("Beta MW3 correlation;%s /mm; #Delta #beta",mw3_dummy[i].Data()), 500, -0.01,0.01, 500, -500, 500);
    if(i==0){
      dummystring = Form("Mw3_X:%s>>h_beta_mw3%i%i",beta_tofw_mod.Data(),cond,i);
    }else if(i==1){
      dummystring = Form("Mw3_Y:%s>>h_beta_mw3%i%i",beta_tofw_mod.Data(),cond,i);
    }else{
      dummystring = Form("%s:%s>>h_beta_mw3%i%i", Mw3_X_mod.Data(),beta_tofw_mod.Data(), cond, i);
    }
    ch -> Draw(dummystring,conditions,"col");
    ////
    c->cd(13+i);
    prof_beta_mw3[cond][i] = h_beta_mw3[cond][i]->ProfileX();
    prof_beta_mw3[cond][i] ->Draw();
    ////
    /*
    c ->cd(1);
  // FRS PID with Z=20, R3BMusic && Twim, Beta corr
  c->cd(2);
  // FRS Brho vs MW3X_mod
  c->cd(3);
  // TProfile and fit
  c->cd(4);*/
  // PID in CaveC
  //  brhocave = Form("((%s-(%f))/(%f))",Mw3_X_mod.Data(), f_mw3brho->GetParameter(0), f_mw3brho->GetParameter(1));
  }
  //
  c->cd(4);
  conditions = "abs(MusicZ-20.)<0.4 && abs(FRSAoQ-2.45)<0.02"; // 49Ca
  h_brhobrho = new TH2D("hbrhobrho", "Brho correlation in FRS and Cave; Brho FRS /Tm; Brho Cave /Tm", 500, 8.7, 9.3, 500, 8.7, 9.3);
  fragbrhostring = Form("(%s-(%f))/(%f)", Mw3_X_mod.Data(), f_brho_mw3[cond][2]->GetParameter(0), f_brho_mw3[cond][2]->GetParameter(1));
  ch->Draw(Form("(%s):%s>>hbrhobrho", fragbrhostring.Data(), brho.Data()), conditions, "col");
  //
  c->cd(8);
  h_aoqaoq = new TH2D("haoqaoq", "Brho correlation in FRS and Cave; AoQ FRS ; AoQ Cave", 500, 2.2, 2.7, 500, 2.2, 2.7);
  fragaoqstring = Form("(%s)*sqrt(1-%s*%s)/((%s)*(%f))", fragbrhostring.Data(), beta_tofw_mod.Data(),beta_tofw_mod.Data(),beta_tofw_mod.Data(), mc_e);
  dummystring = Form("%s:FRSAoQ>>haoqaoq",fragaoqstring.Data());
  cout << dummystring <<endl;
  ch->Draw(dummystring, conditions, "col");
  //
  c->cd(12);
  h_pid = new TH2D("hpid", "PID of fragment; AoQ Cave; TwimZ", 500, 2.2, 2.7, 500, 10, 30);
  ch->Draw(Form("TwimZ:%s>>hpid",fragaoqstring.Data()), conditions, "col");
  //
  c -> Print(recobrho_outpdf);
  delete c;
  return cond++;
}

int transfer_mat(){
  ofstream fout("./fragment/output/fragment_fit_brho.txt", ofstream::out);
  //
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(4,4);
  c -> Print(brho_outpdf + "[");
  ////
  TString conditions = "1";
  draw_transfer_corr(conditions);
  //
  beamcondition += cut_mw + Zgate;
  conditions = "(";
  conditions += beamcondition;
  conditions += ")";
  draw_transfer_corr(conditions);
  //
  mwcondition = Form("&& abs(%s-(%f))<2.", axis_mw_h[0].Data(), max_mw1[0][0]);
  mwcondition += Form("&& abs(%s-(%f))<2.", axis_mw_h[2].Data(), max_mw1[2][0]);
  mwcondition += Form("&& abs(%s-(%f))<2.", axis_mw_v[3].Data(), max_mw1[3][1]);
  //mwcondition += "&& Mw3_X < 300.";
  
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  draw_transfer_corr(conditions);
  //fout << f_mw3[cond-1][2]->GetParameter(0) <<", " <<f_mw3[cond-1][2]->GetParameter(1) <<", "<<axis_mw3[2]<<endl;
  //
  // Reset the mwcondition with y and b only
  mwcondition = Form("&& abs(%s-(%f))<2", axis_mw_h[2].Data(), max_mw1[2][0]);
  mwcondition += Form("&& abs(%s-(%f))<2", axis_mw_v[3].Data(), max_mw1[3][1]);
  //mwcondition += "&& Mw3_X < 300.";
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  Mw3_X_mod = Form("(Mw3_X - (%f) - %s *(%f))", f_mw3[cond-1][2]->GetParameter(0), axis_mw3[2].Data(), f_mw3[cond-1][2]->GetParameter(1));
  cout << "Newly Corrected Mw3_X = " << Mw3_X_mod <<endl;
  draw_transfer_corr(conditions);
  //
  // Reset the mwcondition with y  only
  mwcondition = Form("&& abs(%s-(%f))<2", axis_mw_h[2].Data(), max_mw1[2][0]);
  //mwcondition += Form("&& abs(%s-(%f))<2", axis_mw_v[3].Data(), max_mw1[3][1]);
  //mwcondition += "&& Mw3_X < 300.";
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  Mw3_X_mod += Form("+ (- (%f) - %s *(%f))", f_mw3[cond-1][0]->GetParameter(0), axis_mw3[0].Data(), f_mw3[cond-1][0]->GetParameter(1));
  cout << "Newly Corrected Mw3_X = " << Mw3_X_mod <<endl;
  draw_transfer_corr(conditions);
  //
  // Reset the mwcondition
  //mwcondition = "&& Mw3_X < 300.";
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  Mw3_X_mod += Form("+ (- (%f) - %s *(%f))", f_mw3[cond-1][3]->GetParameter(0), axis_mw3[3].Data(), f_mw3[cond-1][3]->GetParameter(1));
  cout << "Newly Corrected Mw3_X = " << Mw3_X_mod <<endl;
  draw_transfer_corr(conditions);
  //
  //mwcondition = "&& Mw3_X < 300.";
  conditions = "(";
  //conditions += beamcondition; // remove beta condition
  conditions += cut_mw;
  conditions += "abs(MusicZ-20.)<0.4";
  conditions += mwcondition;
  conditions += ")";
  Mw3_X_mod += Form("+ (- (%f) - %s *(%f))", f_mw3[cond-1][1]->GetParameter(0), axis_mw3[1].Data(), f_mw3[cond-1][1]->GetParameter(1));
  cout << "Newly Corrected Mw3_X = " << Mw3_X_mod <<endl;
  draw_transfer_corr(conditions);
  //
  // Final figures
  conditions = "(";
  //conditions += beamcondition; // remove beta condition
  conditions += cut_mw;
  conditions += Zgate;
  //conditions += "abs(MusicZ-20.)<0.4";
  //conditions += mwcondition;
  conditions += ")";
  draw_transfer_corr(conditions);
  ///
  fout << "Mw3_X_mod: "<<endl <<Mw3_X_mod <<endl<<endl;
  ////
  /////////////////////////////////////////////////////////////
  /*
  // Correction of ToF length
  conditions = "1";
  draw_toflength_corr(conditions);
  // Gate on 49Ca
  //beamcondition += "&& abs(FRSAoQ-2.45)<0.02";
  mwcondition = Form("&& abs(%s-(%f))<2.", axis_mw_h[0].Data(), max_mw1[0][0]);
  mwcondition += Form("&& abs(%s-(%f))<2.", axis_mw_h[2].Data(), max_mw1[2][0]);
  mwcondition += Form("&& abs(%s-(%f))<2.", axis_mw_v[3].Data(), max_mw1[3][1]);
  //mwcondition += "&& Mw3_X < 300.";
  //range_fit_mw3_high[0] = 0.;
  beta_tofw_mod = beta_tofw;
  //
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  draw_toflength_corr(conditions);
  ////
  mwcondition = Form("&& abs(%s-(%f))<2.", axis_mw_h[2].Data(), max_mw1[2][0]); //y
  mwcondition += Form("&& abs(%s-(%f))<2.", axis_mw_v[3].Data(), max_mw1[3][1]); //b
  //mwcondition += "&& Mw3_X < 300.";
  beta_tofw_mod = Form("(%s-(%f+(%f)*%s))", beta_tofw_mod.Data(), f_mwbeta[cond-1][2]->GetParameter(0), f_mwbeta[cond-1][2]->GetParameter(1), axis_mw3[2].Data()); // Subtract a dependence
  //
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  draw_toflength_corr(conditions);
  ////
  mwcondition = Form("&& abs(%s-(%f))<2.", axis_mw_h[2].Data(), max_mw1[2][0]); //y
  //mwcondition += "&& Mw3_X < 300.";
  beta_tofw_mod = Form("(%s-(%f+(%f)*%s))", beta_tofw_mod.Data(), f_mwbeta[cond-1][0]->GetParameter(0), f_mwbeta[cond-1][0]->GetParameter(1), axis_mw3[0].Data()); // Subtract x dependence
  //
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  draw_toflength_corr(conditions);
  ////
  //mwcondition = "&& Mw3_X < 300.";
  beta_tofw_mod = Form("(%s-(%f+(%f)*%s))", beta_tofw_mod.Data(), f_mwbeta[cond-1][3]->GetParameter(0), f_mwbeta[cond-1][3]->GetParameter(1), axis_mw3[3].Data()); // Subtract b dependence
  //
  conditions = "(";
  conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  draw_toflength_corr(conditions);
  ////
  //mwcondition = "&& Mw3_X < 300.";
  beta_tofw_mod = Form("(%s-(%f+(%f)*%s))", beta_tofw_mod.Data(), f_mwbeta[cond-1][1]->GetParameter(0), f_mwbeta[cond-1][1]->GetParameter(1), axis_mw3[1].Data()); // Subtract y dependence
  //
  conditions = "(1";
  //conditions += beamcondition;
  conditions += mwcondition;
  conditions += ")";
  draw_toflength_corr(conditions);
  ////
  fout << "beta_tofw_mod: "<<endl <<beta_tofw_mod <<endl<<endl;
  / */
  c -> Print(brho_outpdf + "]");
  delete c;
  fout.close();
  return 0;
}

int draw_toflength_corr(TString conditions){
  if(cond >= NUMCOND){
    cerr<<"Too many conditions"<<endl;
    return -1;
  }
  draw_pidgate(conditions); // cd(1-4) will be used.
  ////
  for(int i = 0; i < 4 ; i++){
      c->cd(5+i);
      h_mw12[cond][i] = new TH2D(Form("h_mw12_%i%i",cond,i), Form("Phase space before GLAD;%s / mm;%s / mm",axis_mw_h[i].Data(),axis_mw_v[i].Data()), 500, -1.*range_mw, range_mw, 500, -1.*range_mw, range_mw);
      ch->Draw(Form("%s:%s>>h_mw12_%i%i",axis_mw_v[i].Data(),axis_mw_h[i].Data(),cond,i),conditionwithbetacut,"colz");
      if(cond==0){
	cerr<<"Get maximum of the MWPC distribution first"<<endl;
	return -1;
      }
      //
      c->cd(9+i);
      h_mwbeta[cond][i] = new TH2D(Form("h_mwbeta_%i%i",cond,i), Form("Flight length correction;%s / mm;#beta_{cave}-#beta_{frs}", axis_mw3[i].Data()), 500, i==2?-0.04:-1.*range_mw, i==2?0.04:range_mw, 500, -0.01, 0.01);
      ch->Draw(Form("%s-Beta_S2_Cave:%s>>h_mwbeta_%i%i", beta_tofw_mod.Data(), axis_mw3[i].Data(), cond, i), conditionwithbetacut,"colz");
      //
      c->cd(9+i);
      //prof_mwbeta[cond][i] = h_mwbeta[cond][i]->ProfileX();
      //prof_mwbeta[cond][i] ->Draw("same");
      h_median_mwbeta[cond][i] = h_mwbeta[cond][i]->QuantilesX();
      h_median_mwbeta[cond][i] ->GetYaxis() ->SetRangeUser(-0.01,0.01);
      h_median_mwbeta[cond][i] ->Draw("same");

      f_mwbeta[cond][i] = new TF1(Form("func_mwbeta_%i_%i",cond,i),"[0]+[1]*x", range_fit_mw3_low[i], range_fit_mw3_high[i]);
      f_mwbeta[cond][i]->SetLineWidth(1);
      //prof_mwbeta[cond][i] -> Fit(f_mwbeta[cond][i], "R","");
      h_median_mwbeta[cond][i] -> Fit(f_mwbeta[cond][i], "R","");
      //
  }
  //  
  c -> Print(brho_outpdf);
  return cond++;
}


int draw_transfer_corr(TString conditions){
  if(cond >= NUMCOND){
    cerr<<"Too many conditions"<<endl;
    return -1;
  }
  draw_pidgate(conditions);
  //
  for(int i = 0; i < 4 ; i++){
      c->cd(5+i);
      h_mw12[cond][i] = new TH2D(Form("h_mw12_%i%i",cond,i), Form("Phase space before GLAD;%s / mm;%s / mm",axis_mw_h[i].Data(),axis_mw_v[i].Data()), 500, -1.*range_mw, range_mw, 500, -1.*range_mw, range_mw);
      ch->Draw(Form("%s:%s>>h_mw12_%i%i",axis_mw_v[i].Data(),axis_mw_h[i].Data(),cond,i),conditionwithbetacut,"colz");
      if(cond==0){
	Int_t xbin, ybin, zbin;
	Int_t bin = h_mw12[cond][i]->GetMaximumBin(xbin, ybin, zbin);
	//cout<<"Xaxis: "<<h_mw12[cond][i]->GetXaxis()->GetBinCenter(xbin)<<", Yaxis: "<< h_mw12[cond][i]->GetYaxis()->GetBinCenter(ybin)<<endl;
	max_mw1[i][0]=h_mw12[cond][i]->GetXaxis()->GetBinCenter(xbin);
	max_mw1[i][1]=h_mw12[cond][i]->GetYaxis()->GetBinCenter(ybin);
	line_mw[i][0][0]= new TLine(max_mw1[i][0]-2.,-1.*range_mw,max_mw1[i][0]-2.,range_mw);
	line_mw[i][0][1]= new TLine(max_mw1[i][0]+2.,-1.*range_mw,max_mw1[i][0]+2.,range_mw);
	line_mw[i][1][0]= new TLine(-1.*range_mw,max_mw1[i][1]-2.,range_mw,max_mw1[i][1]-2.);
	line_mw[i][1][1]= new TLine(-1.*range_mw,max_mw1[i][1]+2.,range_mw,max_mw1[i][1]+2.);
	for(int j=0; j<4; j++){
	  //if((i==1&&j/2==1)||(i==3&&j/2==0)) continue;
	  line_mw[i][j/2][j%2]->Draw("same");
	}
      }
      //
      c->cd(9+i);
      h_mw3[cond][i] = new TH2D(Form("h_mw3_%i%i",cond,i), Form("Ion transfer through GLAD;%s / mm;Mw3_X_mod / mm", axis_mw3[i].Data()), 500, i==2?-0.04:-1.*range_mw, i==2?0.04:range_mw, 500, -500, 500);
      ch->Draw(Form("%s:%s>>h_mw3_%i%i", Mw3_X_mod.Data(), axis_mw3[i].Data(), cond, i), conditionwithbetacut,"colz");
      //
      //c->cd(13+i);
      c->cd(9+i);
      prof_mw3[i] = h_mw3[cond][i]->ProfileX();
      prof_mw3[i] ->Draw("same");
      //if(mwcondition.Sizeof()>1){	      }
      f_mw3[cond][i] = new TF1(Form("func_mw3_%i_%i",cond,i),"[0]+[1]*x", range_fit_mw3_low[i], range_fit_mw3_high[i]);
      f_mw3[cond][i]->SetLineWidth(1);
      prof_mw3[i] -> Fit(f_mw3[cond][i], "R","");
  }
  //
  c -> Print(brho_outpdf);
  return cond++;
}

int draw_pidgate(TString conditions){
  cout<<"\033[1;31m Index: "<<cond<<", Draw: "<<conditions<<"\033[m"<<endl;
  c->cd(1);
  h_frspid_mw[cond] = new TH2D(Form("h_frspid%i",cond), "PID in FRS (S2-CaveC); AoQ; MusicZ", 500, 2.2,2.7, 500,10,30);
  ch->Draw(Form("MusicZ:FRSAoQ>>h_frspid%i",cond),conditions,"col");
  //
  c->cd(2);
  h_music_twim_mw[cond] = new TH2D(Form("hmusictwim_mw%i",cond),"R3BMusic and Twim with #beta in FRS;TwimZ;MusicZ",200,10,30,200,10,30);
  ch->Draw(Form("MusicZ:TwimZ>>hmusictwim_mw%i",cond),conditions,"col");
  //
  c->cd(3);

  h_beta_beta_mw[cond][NUMPADDLE] = new TH2D(Form("hbetabeta_mw%i",cond), "#beta correlation: Frs vs CaveC; #beta in FRS; #beta in Cave (ns)", 500, 0.7, 0.8, 500, 0.7, 0.8);
  //
  conditionwithbetacut = conditions + "&& (abs(FragBeta - Beta_S2_Cave)<0.004)";
  /*
  conditionwithbetacut = conditions + "&&(";
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
  }*/
  cout << "Beta Tofw: "<<beta_tofw<<endl<<" conditionwithbetacut: "<< conditionwithbetacut << endl;
  ch->Draw(Form("%s:Beta_S2_Cave>>hbetabeta_mw%i", beta_tofw.Data(), cond), conditionwithbetacut, "col");
  //
  c->cd(4);
  //
  //cout << endl << conditionwithbetacut << endl;
  return 0;
}

int tofw_calib(){
  ofstream fout("./fragment/output/fragment_fit_tof_Nov.csv", ofstream::out);
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(6,5);
  c->Print(outpdf + "[");
  c->cd(1);
  h_frspid = new TH2D("h_frspid", "PID in FRS (S2-CaveC); AoQ; MusicZ", 500, 2.2,2.7, 500,10,30);
  ch->Draw("MusicZ:FRSAoQ>>h_frspid","","col");
  c->cd(2);
  h_tofw_paddle = new TH1D("h_tofw_paddle", "Counts in TofW paddles; Paddle ID; Counts", NUMPADDLE+1, -0.5, NUMPADDLE + 0.5);
  ch->Draw("Tofw_Paddle>>h_tofw_paddle","","");
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    h_beta_tof[i] = new TH2D(Form("hbetatof%i",i+1), Form("ToF vs #beta in Frs: Paddle %i; #beta in FRS; ToF between SofStart and TofW (ns)",i+1), 500, 0.7, 0.8, 500, 40,50);
    ch->Draw(Form("FragTof:Beta_S2_Cave>>hbetatof%i",i+1),Form("Tofw_Paddle==%i",i+1),"col");
  }
  c->Print(outpdf);
  //
  for(int i = 0 ; i<2; i++){
    c->cd(i+1);
    h_music_twim[i] = new TH2D(Form("hmusictwim%i",i),Form("R3BMusic and Twim with #beta in FRS%s",i==0?"":" (gated)"),200,10,30,200,10,30);
    ch->Draw(Form("MusicZ:TwimZ>>hmusictwim%i",i),i==1?Zgate:"","col");
  }
  //
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    h_beta_tof_cut[i] = new TH2D(Form("hbetatofcut%i",i+1), Form("ToF vs #beta in Frs: Paddle %i (Z gated); #beta in FRS; ToF between SofStart and TofW (ns)",i+1), 100, 0.7, 0.8, 100, 40,50);
    ch->Draw(Form("FragTof:Beta_S2_Cave>>hbetatofcut%i",i+1),Form("Tofw_Paddle==%i&&%s",i+1,Zgate.Data()),"col");
  }
  c->Print(outpdf);
  //
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    f_tof[i] = new TF1(Form("func%i",i),"[0]+[1]/x",0.75,0.78);
    f_tof[i] ->SetLineWidth(1);
    h_median_tof[i] = h_beta_tof_cut[i]->QuantilesX();
    h_median_tof[i] ->GetYaxis() ->SetRangeUser(40,50);
    h_median_tof[i] ->Draw();
    if(i > 21){
      f_tof[i] ->SetLineColor(4);
      f_tof[i] ->FixParameter(0, f_tof[21]->GetParameter(0));
      f_tof[i] ->FixParameter(1, f_tof[21]->GetParameter(1));
    }
    h_median_tof[i] ->Fit(f_tof[i], "RQ", "");
    //prof_tof[i] = h_beta_tof_cut[i]->ProfileX();
    //prof_tof[i] ->Draw("same");
    fragbeta[i] += Form("(1./(FragTof-(%f))*(%f))", f_tof[i]->GetParameter(0), f_tof[i]->GetParameter(1));
    fout << i+1 << ", " << f_tof[i]->GetParameter(0) <<", " <<f_tof[i]->GetParameter(1) <<endl;
  }
  c->Print(outpdf);
  //
  for(int i = 0 ; i<2; i++){
    c->cd(i+1)->Clear();
  }
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    h_beta_beta[i] = new TH2D(Form("hbetabeta%i",i+1), Form("#beta correlation: Frs vs Paddle %i; #beta in FRS; #beta in Cave (ns)",i+1), 500, 0.7, 0.8, 500, 0.7, 0.8);
    ch->Draw(Form("%s:Beta_S2_Cave>>hbetabeta%i", fragbeta[i].Data(), i+1),Form("Tofw_Paddle==%i",i+1),"col");
  }
  c->Print(outpdf);
  //
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    h_beta_beta_cut[i] = new TH2D(Form("hbetabetacut%i",i+1), Form("#beta correlation: Frs vs Paddle %i (beta-gated); #beta in FRS; #beta in Cave (ns)",i+1), 500, 0.7, 0.8, 500, 0.7, 0.8);
    ch->Draw(Form("%s:Beta_S2_Cave>>hbetabetacut%i", fragbeta[i].Data(), i+1),Form("Tofw_Paddle==%i && abs(%s - Beta_S2_Cave)<0.003",i+1, fragbeta[i].Data()),"col");
  }
  c->Print(outpdf);
  //
  h_mw_brho[NUMPADDLE] = new TH2D("hmwbrhototal", "Rigidity vs MWPC3-X; Brho-frs (#beta#gamma AoQ); MWPC3-X (mm)", 200, min_brho, max_brho, 200, -500, 500);
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    h_mw_brho[i] = new TH2D(Form("hmwbrho%i", i+1), Form("Rigidity vs MWPC3-X, Paddle %i; #beta#gamma AoQ; MWPC3-X (mm)", i+1), 200, min_brho, max_brho, 200, -500, 500);
    TString condition =  Form("Tofw_Paddle==%i && abs(%s - Beta_S2_Cave)<0.003 && %s && Beta_S2_Cave<1 && Beta_S2_Cave>0.5", i+1, fragbeta[i].Data(), Zgate.Data());
    cout << condition <<endl;
    ch->Draw(Form("Mw3_X:%s>>hmwbrho%i", brho.Data(), i+1), condition,"col");
    //ch->Draw(Form("Mw3_X:Beta_S2_Cave*TheGamma*TheAoQ>>hmwbrho%i", i+1, brho.Data()), condition,"col");
    if(i<2) continue; // As there's many events
    h_mw_brho[NUMPADDLE] -> Add(h_mw_brho[i]);
  }
  c->cd(1);
  h_mw_brho[NUMPADDLE]->Draw("col");
  c->cd(2);
  prof_mw_brho = h_mw_brho[NUMPADDLE]->ProfileX();
  prof_mw_brho->Draw();
  f_mw3brho = new TF1("f_mw3brho","[0]+[1]*x",9.04,9.13); // MW3 = [0] + [1] * Brho --> Brho = (MW3-[0])/[1]
  f_mw3brho->SetLineWidth(1);
  prof_mw_brho->Fit(f_mw3brho, "R", "");
  //cout<<"Fit brho: "<<f_mw3brho->GetParameter(0)<<" + x*"<<f_mw3brho->GetParameter(1)<<endl;
  fout << 100 << ", " << f_mw3brho->GetParameter(0) <<", " <<f_mw3brho->GetParameter(1) <<endl;
  brhocave += Form("((Mw3_X-(%f))/(%f))",f_mw3brho->GetParameter(0), f_mw3brho->GetParameter(1));
  c->Print(outpdf);
  //
  //p->ClearCache();
  cout<<"CLEARED PROOF CACHE"<<endl;
  //
  h_fragaoq[NUMPADDLE] = new TH2D("h_fragaoq", "AoQ of FRS and Fragment; AoQ in FRS; AoQ in CaveC", 500, 2.2,2.7, 500,2.2,2.7);
  for(int i=0; i<NUMPADDLE; i++){
    c->cd(i+3);
    fragaoq[i] += Form("((%s)/(%s*%f)*sqrt(1-%s*%s))", brhocave.Data(), fragbeta[i].Data(), mc_e, fragbeta[i].Data(), fragbeta[i].Data());
    //Form("((%s)/%s*sqrt(1-%s*%s))", brho.Data(), fragbeta[i].Data(), fragbeta[i].Data(), fragbeta[i].Data());
    cout<<"Paddle"<<i+1<<": "<<fragaoq[i]<<endl;
    //
    h_fragaoq[i] = new TH2D(Form("h_fragaoq%i",i), Form("AoQ of FRS and Fragment (Paddle %i); AoQ in FRS; AoQ in CaveC", i+1), 500, 2.2, 2.7, 500,2.2,2.7);
    ch->Draw(Form("%s:FRSAoQ>>h_fragaoq%i", fragaoq[i].Data(),i),Form("Tofw_Paddle==%i",i+1),"col");
    h_fragaoq[NUMPADDLE]->Add(h_fragaoq[i]);
  }
  c->cd(1);
  h_frspid->Draw("col");
  c->cd(2);
  h_fragaoq[NUMPADDLE]->Draw("col");
  c->Print(outpdf);
  //
  h_fragpid[NUMPADDLE] = new TH2D("h_fragpid", "Fragment PID; AoQ; Twim Z", 500, 2.2, 2.7, 500, 10,30);
  for(int i=0; i<NUMPADDLE; i++){
    c->cd(i+3);
    h_fragpid[i] = new TH2D(Form("h_fragpid%i",i), Form("Fragment PID (Paddle %i); AoQ; Twim Z", i+1), 500, 2.2,2.7, 500,10,30);
    ch->Draw(Form("TwimZ:%s>>h_fragpid%i",fragaoq[i].Data(),i),Form("Tofw_Paddle==%i",i+1),"col");
    h_fragpid[NUMPADDLE]->Add(h_fragpid[i]);
  }
  c->cd(2);
  h_fragpid[NUMPADDLE]->Draw("col");
  c->Print(outpdf);
  //
  c->Print(outpdf + "]");
  delete c;
  fout.close();
  return 0;
}

int loadtofpara(){
  ifstream fin("./fragment/output/fragment_fit_tof_Nov.csv", ofstream::in);  if(!fin.is_open()||!fin.good()){
    cerr<<"No CSV file found: "<<endl;
    return 1;
  }
  string line;
  int dummyindex = 0;
  char dummychar;
  for(int i = 0; i < NUMPADDLE; i++){ 
    fin >> dummyindex >> dummychar >> par[i][0] >> dummychar >> par[i][1];
    if(dummyindex != i+1) cerr<<"Index of parameters does not match"<<endl;
    //cout<< dummyindex << dummychar << par[i][0] << dummychar << par[i][1] <<endl;
    fragbeta[i] += Form("(1./(FragTof-(%f))*(%f))", par[i][0], par[i][1]);
    fragaoq[i] += Form("((%s)/%s*sqrt(1-%s*%s))", brho.Data(), fragbeta[i].Data(), fragbeta[i].Data(), fragbeta[i].Data());
  }
  fin.close();

  beta_tofw = "(";
  for(int i = 0 ; i<NUMPADDLE; i++){
    beta_tofw += Form("(1.0*(Tofw_Paddle==%i))* %s", i+1, fragbeta[i].Data());
    if(i==NUMPADDLE-1){
      beta_tofw += Form("+(1.0*(Tofw_Paddle>%i))* (1./(FragTof-(2.773780))*(31.234900))", NUMPADDLE); 
      beta_tofw += ")";
    }else{
      beta_tofw += "+";
    }
  }

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
}
