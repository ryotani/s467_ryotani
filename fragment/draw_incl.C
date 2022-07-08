TString string_Date ="21Jun";
//TString strfilein[2] = {"fragment/output/incl_out_indivfit_30May_38Ca_all.root", "fragment/output/incl_out_indivfit_30May_50Ca_all.root"};
//TString strfilein[2] = {"fragment/output/incl_out_indivfit_14Jun_38Ca_all.root", "fragment/output/incl_out_indivfit_14Jun_50Ca_all.root"};
TString strfilein[2] = {"fragment/output/incl_out_indivfit_21Jun_38Ca_all.root", "fragment/output/incl_out_indivfit_21Jun_50Ca_all.root"};
TFile *filein[2];
TGraphErrors *g[3][3][2]; // file name, target, reaction
TCanvas *c, *c2, *cfit;;
void incl_49Ca1n(), draw_aoq_fit();;
double ch2_1n =0.;
double Ech2_1n=0.;
double carbon_1n=0.;
double Ecarbon_1n=0.;
double sigma[100][3][3]={0.}; // [target][reaction]
double sigma_err[100][3][3]={0.}; // 
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
	if(i==0){
	  g[2][targ][reaction] = new TGraphErrors();
	  for(int N=0; N<g[0][targ][reaction]->GetN(); N++){
	    int N0 = g[2][targ][reaction]->GetN();
	    double X, Y;
	    g[0][targ][reaction]->GetPoint(N, X, Y);
	    double eX = g[0][targ][reaction]->GetErrorX(N);
	    double eY = g[0][targ][reaction]->GetErrorY(N);
	    if(X < 38.5 ||(42.5 < X && X < 46.5) || (41.5 < X && X<42.5)) continue;
	    g[2][targ][reaction]->SetPoint(N0, X, Y);
	    g[2][targ][reaction]->SetPointError(N0, eX, eY);
	    sigma[(int)X][targ][reaction] = Y;
	    sigma_err[(int)X][targ][reaction] = eY;
	  }
	} else if(i==1){
	  for(int N=0; N < g[1][targ][reaction]->GetN(); N++){
	    int N0=g[2][targ][reaction]->GetN();
	    double X, Y;
	    g[1][targ][reaction]->GetPoint(N, X, Y);
	    double eX = g[1][targ][reaction]->GetErrorX(N);
	    double eY = g[1][targ][reaction]->GetErrorY(N);
	    /*
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
	    }*/
	    if((40.5 < X && X < 46.5) || (50 < X && reaction == 0)) continue;
	    //if(X < 38.5 ||(40.5 < X && X < 46.5)) continue;
	    g[2][targ][reaction]->SetPoint(N0, X, Y);
	    g[2][targ][reaction]->SetPointError(N0, eX, eY);
	    sigma[(int)X][targ][reaction] = Y;
	    sigma_err[(int)X][targ][reaction] = eY;
	  }
	  /*
	  for(int N=0; N<g[0][targ][reaction]->GetN(); N++){
	    double X, Y;
	    g[0][targ][reaction]->GetPoint(N, X, Y);
	    g[0][targ][reaction]->SetPoint(N, X + (targ==0?0.05:-0.05), Y);
	    }*/
	}
      }
    }
  } 
 for(int targ=0; targ<3; targ++){
    for(int reaction=0; reaction<2; reaction++){
      for(int N=0; N<g[2][targ][reaction]->GetN(); N++){
	double X, Y;
	g[2][targ][reaction]->GetPoint(N, X, Y);
	double eX = g[2][targ][reaction]->GetErrorX(N);
	double eY = g[2][targ][reaction]->GetErrorY(N);
	/*
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
	  g[1][targ][reaction]->SetPointError(N, eX, eY);
	}
	*/
	
	g[2][targ][reaction]->SetPoint(N, X + (targ==0?0.05:-0.05), Y);
	}
      c->cd(targ*2+reaction+1);
      g[2][targ][reaction]->SetMarkerStyle(20+targ);//kCircle);
      g[2][targ][reaction]->SetMarkerSize(1);
      g[2][targ][reaction]->SetMarkerColor(color[targ]);
      g[2][targ][reaction]->SetLineColor(color[targ]);
      g[2][targ][reaction]->SetTitle(Form("Inclusive cross sections for %s channel",reaction==0?"-1n":"-1p"));
      g[2][targ][reaction]->GetXaxis()->SetRangeUser(37.8,51.2);
      //g[2][targ][reaction]->GetXaxis()->SetRange(37.5,51.5);
      //g[2][targ][2]->GetXaxis()->SetRangeUser(37.5,51.5);
      //g[2][targ][reaction]->GetYaxis()->SetRangeUser(0,g[2][targ][reaction]->GetYaxis()->GetXmax());
      g[2][targ][reaction]->GetYaxis()->SetRangeUser(0,targ==1?200:100);
      g[2][targ][reaction]->GetXaxis()->SetTitle("Mass of parent nuclei / A");
      g[2][targ][reaction]->GetYaxis()->SetTitle("Inclusive cross sections / mb");
      g[2][targ][reaction]->GetXaxis()->SetLabelSize(30);
      g[2][targ][reaction]->GetXaxis()->SetTitleSize(30);
      g[2][targ][reaction]->GetYaxis()->SetLabelSize(30);
      g[2][targ][reaction]->GetYaxis()->SetTitleSize(30);
      g[2][targ][reaction]->Draw("ape");
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
    g[2][2][reaction]->Draw("ape");
    g[2][0][reaction]->Draw("pe");
    //
    auto legend = new TLegend(0.55,0.75,0.89,0.89);
    //legend->SetHeader("","C"); // option "C" allows to center the header
    //legend->AddEntry(g[2][0][reaction],"Carbon contribution","pe");
    if(reaction==0){
      legend->AddEntry(g[2][0][reaction],"#sigma_{(C,C*n)}","pe");
      legend->AddEntry(g[2][2][reaction],"#sigma_{(p,pn)}","pe");
    }else{
      legend->AddEntry(g[2][0][reaction],"#sigma_{(C,C*p)}","pe");
      legend->AddEntry(g[2][2][reaction],"#sigma_{(p,2p)}","pe");
    }
    //legend->AddEntry(g[2][2][reaction],"Proton contribution","pe");
    //legend->AddEntry(g[2][0][reaction],"#sigma_{C,}","pe");
    //legend->AddEntry(g[2][2][reaction],"Proton contribution","pe");
    legend->SetTextSize(0.03);
    legend->SetTextAlign(12);
    legend->SetLineWidth(0);
    legend->Draw();
    c2->SaveAs(Form("./fragment/output/incl_%s_50Ca_%s.png",reaction==0?"-1n":"-1p",string_Date.Data()));
    c2->SaveAs(Form("./fragment/output/incl_%s_50Ca_%s.pdf",reaction==0?"-1n":"-1p",string_Date.Data()));
  }
  draw_aoq_fit();
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


void draw_aoq_fit(){
  cfit = new TCanvas("aoqfit","aoqfit",1000,1000);
  cfit->Divide(3,2);
  TH2F* htmp[6];
  //for(int i =0; i < 6 ;i++) cfit->cd(1+3*(i%2)+i/2)->SetLogy();
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf[");
  //
  int mass_list[] = {39, 40, 41, 47, 48, 49, 50, 51};
  //
  for(int l = 0; l < (sizeof mass_list)/(sizeof mass_list[0]); l++){
    TFile *tmpfile;
    if(l<3){
      tmpfile = filein[0];
    }else{
      tmpfile = filein[1];
    }
    htmp[0] = (TH2F*) tmpfile->FindObjectAny(Form("h_aoq_gated_all_Z20_A%i_Z20_empty", mass_list[l]))->Clone();
    htmp[3] = (TH2F*) tmpfile->FindObjectAny(Form("h_aoq_gated_all_Z20_A%i_Z19_empty", mass_list[l]))->Clone();
    htmp[1] = (TH2F*) tmpfile->FindObjectAny(Form("h_aoq_gated_all_Z20_A%i_Z20_carbon", mass_list[l]))->Clone();
    htmp[4] = (TH2F*) tmpfile->FindObjectAny(Form("h_aoq_gated_all_Z20_A%i_Z19_carbon", mass_list[l]))->Clone();
    htmp[2] = (TH2F*) tmpfile->FindObjectAny(Form("h_aoq_gated_all_Z20_A%i_Z20_ch2", mass_list[l]))->Clone();
    htmp[5] = (TH2F*) tmpfile->FindObjectAny(Form("h_aoq_gated_all_Z20_A%i_Z19_ch2", mass_list[l]))->Clone();
    //    for(int i =0; i < 6 ; i++){
    //cfit->cd(1+3*(i%2)+i/2);
    //htmp[i]->Draw();
    //}
    //TLatex latex;
    //latex.SetTextSize(0.025);
    for(int i =0; i < 6 ; i++){
      cfit->cd(i+1)->SetLogy();
      htmp[i]->Draw();
      if(i%3==0) continue;
      if(i%3==1){
	auto *t = new TLatex((double)(mass_list[l]+1)/20.,htmp[i]->GetMaximum(),Form("#sigma_{c} = %.1f #pm %.1f", sigma[mass_list[l]][0][i/3], sigma_err[mass_list[l]][0][i/3]) );
	t->Draw();
	//latex.DrawLatex(2.2,.6,Form("#sigma_{c} = %.1f #pm %.1f", sigma[mass_list[l]][0][i/3], sigma_err[mass_list[l]][0][i/3]) );
      }
      if(i%3==2){
	auto *t = new TLatex( (double)(mass_list[l]+1)/20.,1.2*htmp[i]->GetMaximum(),Form("#sigma_{ch2} = %.1f #pm %.1f", sigma[mass_list[l]][1][i/3], sigma_err[mass_list[l]][1][i/3]) );
	auto *t2 = new TLatex((double)(mass_list[l]+1)/20.,0.9*htmp[i]->GetMaximum(),Form("#sigma_{p} = %.1f #pm %.1f", sigma[mass_list[l]][2][i/3], sigma_err[mass_list[l]][2][i/3]) );
	t->Draw();
	t2->Draw();
      }
    }
    cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
    //
    for(int i =0; i < 6 ; i++){
      cfit->cd(i+1)->SetLogy(0);
      //htmp[i]->Draw();
    }
    cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  }
   
  /*
  htmp[0] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A40_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A40_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A40_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A40_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A40_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A40_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  //
  htmp[0] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A41_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A41_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A41_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A41_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A41_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[0]->FindObjectAny("h_aoq_gated_all_Z20_A41_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  //
  htmp[0] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A47_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A47_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A47_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A47_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A47_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A47_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  //
  htmp[0] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A48_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A48_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A48_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A48_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A48_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A48_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  //
  htmp[0] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A49_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A49_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A49_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A49_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A49_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A49_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  //
  htmp[0] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A50_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A50_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A50_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A50_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A50_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A50_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  //
  htmp[0] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A51_Z20_empty")->Clone();
  htmp[1] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A51_Z19_empty")->Clone();
  htmp[2] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A51_Z20_carbon")->Clone();
  htmp[3] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A51_Z19_carbon")->Clone();
  htmp[4] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A51_Z20_ch2")->Clone();
  htmp[5] = (TH2F*) filein[1]->FindObjectAny("h_aoq_gated_all_Z20_A51_Z19_ch2")->Clone();
  for(int i =0; i < 6 ; i++){
    cfit->cd(1+3*(i%2)+i/2);
    htmp[i]->Draw();
  }
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf");
  / */
  cfit->SaveAs("fragment/output/frag_aoq_fit.pdf]");
}
