#define NUMPADDLE 28
TProof *p =  TProof::Open("");

TCanvas *c;
TString infile = "./tofw/output/mktree_tofw_frs.root";
TString outpdf = "./tofw/output/fragment_calib.pdf";
TString Zgate = "";
TString fragbeta[NUMPADDLE] = {""};
TString fragaoq[NUMPADDLE] = {""};
TString brho = "Beta_S2_Cave / sqrt(1-Beta_S2_Cave*Beta_S2_Cave) * AoQ_S2_Cave";
TChain *ch;
TH2D *h_beta_tof[NUMPADDLE], *h_beta_tof_cut[NUMPADDLE], *h_music_twim[2], *h_beta_beta[NUMPADDLE], *h_beta_beta_cut[NUMPADDLE];
TH2D *h_mw_brho[NUMPADDLE+1], *h_fragaoq[NUMPADDLE+1], *h_fragpid[NUMPADDLE+1], *h_frspid;
TProfile *prof_tof[NUMPADDLE], *prof_mw_brho;
TH1D *h_tofw_paddle, *h_median_tof[NUMPADDLE];
TF1 *f_tof[NUMPADDLE];
ofstream fout("./tofw/output/fragment_fit_tof.csv", ofstream::out);

void initialise();

int fragment_calib(){
  initialise();
  if (NUMPADDLE>28) return 1;
  c->cd(1);
  h_frspid = new TH2D("h_frspid", "PID in FRS (S2-CaveC); AoQ; MusicZ", 500, 2.2,2.7, 500,10,30);
  ch->Draw("MusicZ:AoQ_S2_Cave>>h_frspid","","col");
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
    if(i==1)
      Zgate += "abs(MusicZ-TwimZ)<0.5";
    ch->Draw(Form("MusicZ:TwimZ>>hmusictwim%i",i),Zgate,"col");
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
    fragbeta[i] += Form("(1./(FragTof-%f)*%f)", f_tof[i]->GetParameter(0), f_tof[i]->GetParameter(1));
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
  h_mw_brho[NUMPADDLE] = new TH2D("hmwbrhototal", "Rigidity vs MWPC3-X; #beta#gamma AoQ; MWPC3-X (mm)", 200, 2.8, 3., 200, -500, 500);
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1+2);
    h_mw_brho[i] = new TH2D(Form("hmwbrho%i", i+1), Form("Rigidity vs MWPC3-X, Paddle %i; #beta#gamma AoQ; MWPC3-X (mm)", i+1), 200, 2.8, 3., 200, -500, 500);
    TString condition =  Form("Tofw_Paddle==%i && abs(%s - Beta_S2_Cave)<0.003 && %s && Beta_S2_Cave<1 && Beta_S2_Cave>0.5", i+1, fragbeta[i].Data(), Zgate.Data());
    cout << condition <<endl;
    ch->Draw(Form("Mw3_X:%s>>hmwbrho%i", brho.Data(), i+1), condition,"col");
    //ch->Draw(Form("Mw3_X:Beta_S2_Cave*TheGamma*TheAoQ>>hmwbrho%i", i+1, brho.Data()), condition,"col");
    h_mw_brho[NUMPADDLE] -> Add(h_mw_brho[i]);
  }
  c->cd(1);
  h_mw_brho[NUMPADDLE]->Draw("col");
  c->cd(2);
  prof_mw_brho = h_mw_brho[NUMPADDLE]->ProfileX();
  prof_mw_brho->Draw();
      //prof_tof[i] = h_beta_tof_cut[i]->ProfileX();
    //prof_tof[i] ->Draw("same");
  c->Print(outpdf);
  //
  h_fragaoq[NUMPADDLE] = new TH2D("h_fragaoq", "AoQ of FRS and Fragment; AoQ in FRS; AoQ in CaveC", 500, 2.2,2.7, 500,2.2,2.7);
  for(int i=0; i<NUMPADDLE; i++){
    c->cd(i+3);
    fragaoq[i] += Form("((%s)/%s*sqrt(1-%s*%s))", brho.Data(), fragbeta[i].Data(), fragbeta[i].Data(), fragbeta[i].Data());
    cout<<fragaoq[i]<<endl;
    //
    h_fragaoq[i] = new TH2D(Form("h_fragaoq%i",i), Form("AoQ of FRS and Fragment (Paddle %i); AoQ in FRS; AoQ in CaveC", i+1), 500, 2.2, 2.7, 500,2.2,2.7);
    ch->Draw(Form("%s:AoQ_S2_Cave>>h_fragaoq%i", fragaoq[i].Data(),i),Form("Tofw_Paddle==%i",i+1),"col");
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
    //ch->Draw(Form("TwimZ:%s/%s*sqrt(1-%s*%s)>>h_fragpid%i",brho.Data(), fragbeta[i].Data(), fragbeta[i].Data(), fragbeta[i].Data(),i),"","col");
    h_fragpid[NUMPADDLE]->Add(h_fragpid[i]);
  }
  c->cd(2);
  h_fragpid[NUMPADDLE]->Draw("col");
  c->Print(outpdf);
  //
  c->Print(outpdf + "]");
  delete ch;
  p->Close();
  fout.close();
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
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(6,5);
  c->Print(outpdf + "[");

}
