TString strfilein[2] = {"rootfiles_local/incl_out_test_indivfit_12-2Apr_38Ca_all.root", "rootfiles_local/incl_out_test_indivfit_12-2Apr_50Ca_all.root"};
TFile *filein[2];
TGraphErrors *g[2][3][2]; // file name, target, reaction
TCanvas *c, *c2;
void incl_49Ca1n();
double ch2_1n =0.;
double Ech2_1n=0.;
double carbon_1n=0.;
double Ecarbon_1n=0.;
int color[3]={kBlack, kBlue, kRed};

void draw_incl(){
  incl_49Ca1n();
  c=new TCanvas("c","c",1000,1000);
  c->Divide(2,3);
  for(int i=0; i<2; i++){//50Ca and 38Ca
    filein[i] = TFile::Open(strfilein[i]);
    filein[i]->cd();
    for(int targ=0; targ<3; targ++){
      for(int reaction=0; reaction<2; reaction++){
	g[i][targ][reaction] = (TGraphErrors*)filein[i]->FindObjectAny(Form("g_incl_20_%i_%i",targ,reaction));
	c->cd(targ*2+reaction+1);
	/*
	if(i==0){
	  g[i][targ][reaction]->Draw("ape");
	  //g[i][targ][reaction]->GetXaxis()->SetRange(38,51);
	}else{
	  g[i][targ][reaction]->Draw("ape");
	}
	*/
	if(i==1){
	  for(int N=0; N<g[i][targ][reaction]->GetN(); N++){
	    int N0=g[0][targ][reaction]->GetN();
	    double X, Y;
	    g[1][targ][reaction]->GetPoint(N, X, Y);
	    double eX = g[1][targ][reaction]->GetErrorX(N);
	    double eY = g[1][targ][reaction]->GetErrorY(N);
	    if(abs(X-49.)<0.2 && reaction ==0){
	      if(targ==0){
		carbon_1n = Y;
		Ecarbon_1n = eY;
	      }else if(targ==1){
		Y = ch2_1n;
		eY = Ech2_1n;
	      }else if(targ==2){
		Y = (ch2_1n-carbon_1n)/2.;
		eY = Y * TMath::Sqrt(TMath::Power(Ecarbon_1n/carbon_1n,2.)+TMath::Power(Ech2_1n/ch2_1n,2.));
	      }
	    }
	    g[0][targ][reaction]->SetPoint(N0, X, Y);
	    g[0][targ][reaction]->SetPointError(N0, eX, eY);	    
	  }
	  for(int N=0; N<g[0][targ][reaction]->GetN(); N++){
	    double X, Y;
	    g[0][targ][reaction]->GetPoint(N, X, Y);
	    g[0][targ][reaction]->SetPoint(N, X + (targ==0?0.05:-0.05), Y);
	  }
	}
      }
    }
  }
  for(int targ=0; targ<3; targ++){
    for(int reaction=0; reaction<2; reaction++){
      c->cd(targ*2+reaction+1);
      g[0][targ][reaction]->SetMarkerStyle(20+targ);//kCircle);
      g[0][targ][reaction]->SetMarkerSize(1);
      g[0][targ][reaction]->SetMarkerColor(color[targ]);
      g[0][targ][reaction]->SetLineColor(color[targ]);
      g[0][targ][reaction]->GetXaxis()->SetRangeUser(37.8,51.2);
      //g[0][targ][reaction]->GetXaxis()->SetRange(37.5,51.5);
      //g[0][targ][1]->GetXaxis()->SetRangeUser(37.5,51.5);
      //g[0][targ][reaction]->GetYaxis()->SetRangeUser(0,g[0][targ][reaction]->GetYaxis()->GetXmax());
      g[0][targ][reaction]->GetYaxis()->SetRangeUser(0,targ==1?200:100);
      g[0][targ][reaction]->GetXaxis()->SetTitle("Mass / A");
      g[0][targ][reaction]->GetYaxis()->SetTitle("Inclusive cross sections / mb");
      g[0][targ][reaction]->GetXaxis()->SetLabelSize(30);
      g[0][targ][reaction]->GetXaxis()->SetTitleSize(30);
      g[0][targ][reaction]->GetYaxis()->SetLabelSize(30);
      g[0][targ][reaction]->GetYaxis()->SetTitleSize(30);
      g[0][targ][reaction]->Draw("ape");
    }
  }
  c2 = new TCanvas("c2","c2",1000,1000);
  //c2->Divide(2,1);
  for(int reaction=0; reaction<2; reaction++){
    //c2->cd(2-reaction);
    /*    for(int targ=0; targ<2; targ++){
      for(int N=0; N<g[0][2*targ][reaction]->GetN(); N++){
	double X, Y;
	g[0][2*targ][reaction]->GetPoint(N, X, Y);
	g[0][2*targ][reaction]->SetPoint(N, X + TMath::Power(-1.,targ+1) * 0.1, Y);
      }
      }*/
    g[0][2][reaction]->Draw("ape");
    g[0][0][reaction]->Draw("pe");
    //
    auto legend = new TLegend(0.55,0.75,0.89,0.89);
    //legend->SetHeader("","C"); // option "C" allows to center the header
    legend->AddEntry(g[0][0][reaction],"Carbon contribution","pe");
    legend->AddEntry(g[0][2][reaction],"Proton contribution","pe");
    legend->SetTextSize(0.03);
    legend->SetTextAlign(12);
    legend->SetLineWidth(0);
    legend->Draw();
    c2->SaveAs(Form("incl_%s.png",reaction==0?"-1n":"-1p"));
  }
}
/*
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           4.53230e+02   6.96931e+00   2.64740e-02   2.05810e-06
   2  p1           2.39546e+00   1.69119e-04   1.14225e-06  -7.22916e-02
   3  p2           1.22674e-02   1.67631e-04   5.80319e-07   4.69452e-02
   4  p3           2.16636e+04   5.64625e+02   3.89441e-01  -1.64059e-08
   5  p4           2.45080e+00   3.86447e-04   1.16863e-06   8.84813e-02
   6  p5           1.02346e-02   1.27729e-04   7.88303e-08   4.54834e-01
   Transmission for Ch2: 0.778336	 +/-	0.00124052

For empty
60.7088, 0.0151803 for 49Ca->48Ca

*/
void incl_49Ca1n(){
  double in_total[2] = {199837,418775};
  double intensity[2] = {60.7088, 4.53230e+02};
  double width[2] = {1.54420e-02, 1.22674e-02};
  double Eintensity[2] = {2.60879, 6.96931};
  double Ewidth[2] = {7.65359e-04, 1.67631e-04};
  double Nreactions[2]={0.};
  double ENreactions[2]={0.};
  for(int i=0; i<2; i++){
    Nreactions[i] = intensity[i]* width[i] * TMath::Sqrt(2.*TMath::Pi()) * 500. / (2.75-1.9);
    ENreactions[i]= Nreactions[i] * TMath::Sqrt(TMath::Power((Eintensity[i]/intensity[i]),2)+TMath::Power(Ewidth[i]/width[i],2));
  }
  double transmission[2] = {0.863985,0.778336};
  double Etransmission[2] = {0.00201588, 0.00124052};
  //double density = 1.;
  //double thickness = TMath::Na() * 2.4 / density / 14.; // 14 = CH2 mass, 2.4 for cm
  double to_mb = 1.e27;
  const Double_t d_targ[3] = {0, 1.23E+23, 9.85E+22}; // Num / cm^2
  for (int i=0; i<2; i++){
    double tmp= Nreactions[i] / in_total[i] /transmission[i] / d_targ[2] * to_mb;
    ch2_1n += TMath::Power(-1.,i+1) *tmp;
    Ech2_1n = TMath::Sqrt(Ech2_1n*Ech2_1n
			  + tmp * tmp * (TMath::Power(ENreactions[i]/Nreactions[i],2.) + TMath::Power(Etransmission[i]/transmission[i],2)));
  }
  cout<<"49Ca intensity: "<<Nreactions[1]<<"+/-"<<ENreactions[1]<<", transmission:"<<transmission[1]<<"+/-"<<Etransmission[1]<<", incl:"<<ch2_1n<<"+/-"<<Ech2_1n<<endl;

  /*
+const int MINZ=18, MAXZ=21, MAXA=52, MININCLZ=19, MAXINCLZ=21, BINAOQ=500, BINZ=500;
 const double MINAOQ=1.9, MAXAOQ=2.75;
+TString frsname = "50Ca";
+double zoff =0;
+
+const int MINZ=18, MAXZ=21, MAXA=52, MININCLZ=19, MAXINCLZ=21, BINAOQ=500, BINZ=500;
+const double MINAOQ=1.65, MAXAOQ=2.5;
+TString frsname = "38Ca";
+double zoff=0.3; // tentative.
+
*/
}
