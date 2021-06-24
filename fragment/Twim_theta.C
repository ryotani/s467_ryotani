auto runnum = 340, settings=13;
//auto runnum = 238, settings=6;
auto *f = new TFile(Form("./rootfiles/rootfiletmp/TofW/s467_FRSTree_Setting%i_%04d_FragmentTree_june2021.root",settings,runnum));
//auto *f = new TFile("./rootfiles/rootfiletmp/TofW/s467_FRSTree_Setting13_0353_FragmentTree_june2021.root"); // Empty target run
auto Tree = (TTree*)( f->Get("Tree"))->Clone();
auto c = new TCanvas();
TString outfile = Form("./fragment/Twim_theta%04d.pdf",runnum);
TH1F* hmw[4][2];
TH2F* hmwcor[4][4][2], *hmwcorres[4][4][2], *hmwmusic[4][2][1], *hmwmusicres[4][2][1];
//TProfile *pmwcor[4][4][2];
TF1* fmw[4][2], *fmwcor[4][4][2], *fmwmusic[4][2][1];
TString axis[3] = {'X','Y','Z'}, smw[4][2];
double mean[4][2]={NAN}, deviation[4][2]={NAN};
double zpos[4] = {-1.*(49.5+500+336.5+141.25+360+444.5), 1408.-444.5 , 1408.-444.5 +81.+550.+26., 1408.-444.5 -278. + 5925+1050-900.};
double res[4][4][2]={NAN};
// https://elog.gsi.de/land/s444_s467/666
// https://elog.gsi.de/land/s444_s467/1130
// https://elog.gsi.de/land/s444_s467/35

double rangediffmw(int i, int i_axis, double sigma){
  //return (mean[i][i_axis] - mean[0][i_axis]) + sigma * sqrt(pow(deviation[i][i_axis],2) + pow(deviation[0][i_axis],2));
  return ((mean[i][i_axis] - mean[0][i_axis]) + sigma * sqrt(pow(deviation[i][i_axis],2) + pow(deviation[0][i_axis],2)))/(zpos[i]-zpos[0])*1000.;
}

void drawfit(int i_v, int i_h, int i_axis){
  TString str_v =  "(" + smw[i_v][i_axis] + "-" + smw[0][i_axis] + ")/("+Form("%f",(zpos[i_v]-zpos[0])/1000.)+")";
  TString str_h =  "(" + smw[i_h][i_axis] + "-" + smw[0][i_axis] + ")/("+Form("%f",(zpos[i_h]-zpos[0])/1000.)+")";
  TString str_range = Form("(500,%f,%f,500,%f,%f)", rangediffmw(i_h,i_axis,-5.), rangediffmw(i_h,i_axis,5.) , rangediffmw(i_v,i_axis,-5.), rangediffmw(i_v,i_axis,5.));
  //cout<< str_h << " "<<str_range<<endl;
  Tree->Draw(Form("%s:%s>>hmwcor%i%i%i%s", str_v.Data(), str_h.Data(), i_v,i_h,i_axis, str_range.Data()),"","colz");
  hmwcor[i_v][i_h][i_axis] = (TH2F*) gDirectory->Get(Form("hmwcor%i%i%i",i_v,i_h,i_axis));
  auto profx = hmwcor[i_v][i_h][i_axis]->ProfileX();
  fmwcor[i_v][i_h][i_axis] = new TF1(Form("fwmcor%i%i%i",i_v,i_h,i_axis),"pol1",rangediffmw(i_h,i_axis,-2.),rangediffmw(i_h,i_axis,2.));
  profx->Fit(fmwcor[i_v][i_h][i_axis], "R","");
  hmwcor[i_v][i_h][i_axis] ->Draw("colz");
  profx->Draw("same");
  fmwcor[i_v][i_h][i_axis] ->Draw("same");
  c->Print(outfile);
  delete profx;
  //
  str_range = Form("(500,%f,%f,500,%f,%f)", rangediffmw(i_h,i_axis,-5.), rangediffmw(i_h,i_axis,5.) , -10., 10.);
  Tree->Draw(Form("%s-(%s*(%f) + (%f)):%s>>hmwcorres%i%i%i%s", str_v.Data(), str_h.Data(), fmwcor[i_v][i_h][i_axis]->GetParameter(1), fmwcor[i_v][i_h][i_axis]->GetParameter(0),
		  str_h.Data(), i_v,i_h,i_axis, str_range.Data())
	     ,"","colz");
  hmwcorres[i_v][i_h][i_axis] = (TH2F*) gDirectory->Get(Form("hmwcorres%i%i%i",i_v,i_h,i_axis));
  c->Print(outfile);
}

void drawfitmusic(int i_v, int i_h){
  int i_axis=0;
  //TString str_v =  "(" + smw[i_v][i_axis] + "-" + smw[0][i_axis] + ")/("+Form("%f",(zpos[i_v]-zpos[0])/1000.)+")";
  TString str_v =  (i_v==0?"MusicTheta":"TwimTheta");
  TString str_h =  "(" + smw[i_h][i_axis] + "-" + smw[0][i_axis] + ")/("+Form("%f",(zpos[i_h]-zpos[0])/1000.)+")";
  TString str_range = Form("(500,%f,%f,500,%f,%f)", rangediffmw(i_h,i_axis,-5.), rangediffmw(i_h,i_axis,5.) , -0.01, 0.01);
  Tree->Draw(Form("%s:%s>>hmwmusic%i%i%i%s", str_v.Data(), str_h.Data(), i_v,i_h,i_axis, str_range.Data()),"","colz");
  hmwmusic[i_v][i_h][i_axis] = (TH2F*) gDirectory->Get(Form("hmwmusic%i%i%i",i_v,i_h,i_axis));
  auto profx = hmwmusic[i_v][i_h][i_axis]->ProfileX();
  fmwmusic[i_v][i_h][i_axis] = new TF1(Form("fwmmusic%i%i%i",i_v,i_h,i_axis),"pol1",rangediffmw(i_h,i_axis,-2.),rangediffmw(i_h,i_axis,2.));
  profx->Fit(fmwmusic[i_v][i_h][i_axis], "R","");
  hmwmusic[i_v][i_h][i_axis] ->Draw("colz");
  profx->Draw("same");
  fmwmusic[i_v][i_h][i_axis] ->Draw("same");
  c->Print(outfile);
  delete profx;
  //
  Tree->Draw(Form("%s-(%s*(%f) + (%f)):%s>>hmwmusicres%i%i%i%s", str_v.Data(), str_h.Data(), fmwmusic[i_v][i_h][i_axis]->GetParameter(1), fmwmusic[i_v][i_h][i_axis]->GetParameter(0),
		  str_h.Data(), i_v,i_h,i_axis, str_range.Data())
	     ,"","colz");
  hmwmusicres[i_v][i_h][i_axis] = (TH2F*) gDirectory->Get(Form("hmwmusicres%i%i%i",i_v,i_h,i_axis));
  c->Print(outfile);
}



void Twim_theta(void) {
  //c->Divide(2,2);
  c->Print(outfile+"[");
  //Tree->Scan("Mw1_X-Mw2_X: (Mw1_X-Mw2_X)/TwimTheta","","",10);
  //Tree->Scan("Mw1_X-Mw2_X: TwimTheta","","",10);
  Tree->Scan("Mw0_X:Mw1_X:Mw2_X: TwimTheta","","",10);
  //c->cd(1);
  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      smw[i][j]=Form("Mw%i_%s",i,axis[j].Data());
      Tree->Draw(Form("%s>>hmw%i_%s(600,-200,200)",smw[i][j].Data(),i,axis[j].Data()));
      hmw[i][j] = (TH1F*) gDirectory->Get(Form("hmw%i_%s",i,axis[j].Data()));
      hmw[i][j] -> GetXaxis()->SetTitle("Position / mm");
      fmw[i][j] = new TF1(Form("fwm%i%i",i,j),"gaus");
      hmw[i][j] -> Fit(fmw[i][j]);
      hmw[i][j] -> Fit(Form("fwm%i%i",i,j),"","",-150,150);
      mean[i][j] = fmw[i][j]->GetParameter(1);
      deviation[i][j] = fmw[i][j]->GetParameter(2);
      //smw[i][j] += Form("- (%.2f)",mean[i][j]);
      hmw[i][j] -> GetXaxis() -> SetRangeUser(mean[i][j]-5*deviation[i][j], mean[i][j]+5*deviation[i][j]);
      hmw[i][j] ->Draw();
      fmw[i][j] ->Draw("same");
      //c->Print(outfile);
    }
  }
  //
  for(int j=0; j<2; j++){
    for(int i=1; i<4; i++){
      if (i!=1) drawfit(i,1,j);
      if (i!=2) drawfit(i,2,j);
    }
  }
  Tree->Draw("MusicTheta:TwimTheta>>(500,-.01,0.01,500,-0.01,0.01)","","colz");
  //c->Print(outfile);
  for(int j=0; j<2; j++){
    for(int i=1; i<4; i++){
      drawfitmusic(j,i);
    }
  }
  //
  for(int j=0; j<2; j++){
    for(int i=1; i<4; i++){
      for(int k=1; k<3; k++){
      //res[i][j] = ((hmwcorres[i][1][j]->ProjectionY())->Fit("gaus"))->GetParameter(1);
	if(i==k)
	  continue;
	//res[i][k][j] = hmwcorres[i][k][j]->GetStdDev(2) / ((k!=1)?fmwcor[i][1][j]->GetParameter(1):(i==1?1.:fmwcor[i][k][j]->GetParameter(1)/fmwcor[i][1][j]->GetParameter(1)));
	//res[i][k][j] = hmwcorres[i][k][j]->GetStdDev(2) / ((i==1)?1.:fmwcor[i][k][j]->GetParameter(1)/((k!=1)?fmwcor[1][k][j]->GetParameter(1):1.));
	res[i][k][j] = hmwcorres[i][k][j]->GetStdDev(2) / ((i==1)?1.:fmwcor[i][1][j]->GetParameter(1));
	//double err = (i==1)?0.:(hmwcorres[i][k][j]->GetStdDev(2) / pow(fmwcor[i][1][j]->GetParameter(1),2) *fmwcor[i][1][j]->GetParError(1));
	cout<<"Resolution in Mw"<<k<<"-Mw0 and Mw"<<i<<"-Mw0 in "<<axis[j]<<"std: "<< hmwcorres[i][k][j]->GetStdDev(2)<<" axis:" << res[i][k][j] << endl;
      }
    }
  }
  //
  TString str_angle[2];
  for(int i_axis=0; i_axis<2; i_axis++){
    int num_measurements = (i_axis==0?5:3);
    str_angle[i_axis]= "((" + smw[1][i_axis] + "-" + smw[0][i_axis] + Form(")/(%f)",(zpos[1]-zpos[0])/1000.);
    //for(int i=2; i<4; i++){// incl Mw3
    num_measurements += -1;
    for(int i=2; i<3; i++){// excl Mw3
      str_angle[i_axis]+= "+(" + smw[i][i_axis] + "-" + smw[0][i_axis] + Form(")/(%f)",fmwcor[i][1][i_axis]->GetParameter(1)*(zpos[i]-zpos[0])/1000.);
    }
    if(i_axis==0){
      str_angle[i_axis]+= Form("+(MusicTheta)/(%f)",fmwmusic[0][1][i_axis]->GetParameter(1));
      str_angle[i_axis]+= Form("+(TwimTheta)/(%f)",fmwmusic[1][1][i_axis]->GetParameter(1));
    }
    str_angle[i_axis]+=Form(")/%i.",num_measurements);
    //cout<<str_angle[i_axis]<<endl;
    ////////////
    for (int i_h = 1; i_h<6; i_h++){
      TString str_h;
      if(i_h<3){
	str_h=  "(" + smw[i_h][i_axis] + "-" + smw[0][i_axis] + ")/("+Form("%f",(i_h==1?1.:fmwcor[i_h][1][i_axis]->GetParameter(1))*(zpos[i_h]-zpos[0])/1000.)+")";
      }else if(i_h==3){
	str_h=  "(" + smw[2][i_axis] + "-" + smw[1][i_axis] + ")/("+Form("%f",(2==1?1.:fmwcor[2][1][i_axis]->GetParameter(1))*(zpos[2]-zpos[1])/1000.)+")";
      }else if(i_h==4){
	str_h=  Form("MusicTheta/(%f)",fmwmusic[0][1][i_axis]->GetParameter(1));
      }else{
	str_h=  Form("TwimTheta/(%f)",fmwmusic[1][1][i_axis]->GetParameter(1));
      }
      Tree->Draw(Form("%s:%s>>(500,-20,20,500,-20,20)",str_h.Data(),str_angle[i_axis].Data()),"","colz");
      c->Print(outfile);
    }
  }

  double angle[2]={1.};
  for(int j=0; j<2; j++){
    //  angle[j] = ;
  }

  
  // Old histograms
  //
  Tree->Draw("Mw2_X-Mw1_X:Mw2_X-Mw0_X>>(500,-20,20,500,-20,20)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw2_X-Mw1_X:Mw2_X+Mw0_X>>(500,-40,40,500,-20,20)","","colz");
  c->Print(outfile);
  //
  //  c->cd(2);
  Tree->Draw("Mw2_X-Mw1_X:TwimTheta>>(500,-.01,.01,500,-20,20)","","colz");
  c->Print(outfile);
  //
  //  c->cd(3);
  Tree->Draw("Mw2_X-Mw0_X:TwimTheta>>(500,-.01,.01,500,-20,20)","","colz");
  c->Print(outfile);
  //
  //  c->cd(3);
  Tree->Draw("Mw2_X+Mw0_X:TwimTheta>>(500,-.01,.01,500,-40,40)","","colz");
  c->Print(outfile);
  //
  //  c->cd(4);
  Tree->Draw("Mw1_X-Mw0_X:TwimTheta>>(500,-.01,.01,500,-20,20)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw1_X+Mw0_X:TwimTheta>>(500,-.01,.01,500,-40,40)","","colz");
  c->Print(outfile);
  //
  //////////
  Tree->Draw("Mw2_Y-Mw1_Y:Mw2_Y-Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw2_Y-Mw1_Y:Mw2_Y+Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw2_Y+Mw1_Y:Mw2_Y+Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw2_Y+Mw1_Y:Mw2_Y-Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw3_Y-Mw2_Y:Mw2_Y-Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw3_Y-Mw2_Y:Mw2_Y+Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw3_Y+Mw2_Y:Mw2_Y-Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw3_Y+Mw2_Y:Mw2_Y+Mw0_Y>>(500,-150,150,500,-150,150)","","colz");
  c->Print(outfile);
  //
  ////////
  /*
  Tree->Draw("Mw2_Y-Mw0_Y:Mw1_Y-Mw0_Y>>hy2010(500,-150,150,500,-150,150)","","colz");
  auto * hy2010 = (TH2D*)(gDirectory->Get("hy2010"));
  auto * py2010 = hy2010->ProfileX();
  auto * fy2010 = new TF1("fy2010","pol1");
  py2010->Fit(fy2010,"","",-30,-10);
  hy2010->Draw("colz");
  py2010->Draw("same");
  fy2010->Draw("same");
  c->Print(outfile);
  //
  Tree->Draw(Form("(Mw2_X-Mw0_X):(Mw1_X-Mw0_X)*%f>>(500,-150,150,500,-150,150)",fy2010->GetParameter(1)),"","colz");
  c->Print(outfile);
  //
  / *
  //  c->cd(2);
  Tree->Draw("Mw2_Y-Mw1_Y:TwimTheta>>(500,-.01,.01,500,-150,150)","","colz");
  c->Print(outfile);
  //
  //  c->cd(3);
  Tree->Draw("Mw2_Y-Mw0_Y:TwimTheta>>(500,-.01,.01,500,-150,150)","","colz");
  c->Print(outfile);
  //
  //  c->cd(3);
  Tree->Draw("Mw2_Y+Mw0_Y:TwimTheta>>(500,-.01,.01,500,-150,150)","","colz");
  c->Print(outfile);
  //
  //  c->cd(4);
  Tree->Draw("Mw1_Y-Mw0_Y:TwimTheta>>(500,-.01,.01,500,-150,150)","","colz");
  c->Print(outfile);
  //
  Tree->Draw("Mw1_Y+Mw0_Y:TwimTheta>>(500,-.01,.01,500,-150,150)","","colz");
  c->Print(outfile);
  //
  */

  c->Print(outfile+"]");
}
