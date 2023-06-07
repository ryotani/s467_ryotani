void plot_rescattering()
{
  TCanvas* canvas = new TCanvas("canvas", "Function Plot", 800, 600);
  double mb = 1e-27;//cm^-2
  double sigma = TMath::Pi() * pow(1.2e-13 * pow(50.,1./3.),2.);
  double ntarg = 1.e23; // for carbon and ch2
  double nsigma = ntarg * sigma;// / mb;
  cout<<"sigma="<<sigma<<", nsigma="<<nsigma<<endl;

  TF1 *f[3];
  int colour[3]={kBlack,kRed,kBlue};
  double sigman[3]={3*nsigma,nsigma,1};
  for(int i =0; i<2; i++)
  //for (double p = 1; p<=1; p++)
      {
	//func->Draw();
	//TF1* f_tmp = (TF1*)func->Clone();
	TF1* f_tmp = new TF1("Reaction loss deviation;Total reaction cross section asymmetry(%);Difference in cross section (\%) ", [](double* x, double* p) {
	    double pVal = p[0];
	    //η=2(σ’R - σR)/(σ’R + σR),
	    //σ-xn=N’/(NNtexp(-(σR)Nt)) * (-η/)/(exp(-η(σR+σ’R)/2.)-1)* (σR+σ’R)Nt/2.
	    return ((-x[0]/100.) / (exp(-pVal/2. * (x[0])/100.) - 1) * pVal/2.)*100. - 100.;
	  }, -50, 50, 1);
	f_tmp->SetParameter(0, sigman[i]);
	f_tmp->SetLineColor(colour[i]);
	f_tmp->GetYaxis()->SetRangeUser(-10,10);
	f[i] = f_tmp;
	f[i]->Draw(i==0?"":"same");
    }

  TLegend *ld = new TLegend(0.2,0.65,0.4,0.8);
  ld->AddEntry(f[1], "Carbon target", "l");
  ld->AddEntry(f[0], "CH_{2} target", "l");
  ld->Draw();
    //canvas->SetGrid();
    //canvas->Draw();
  canvas->SaveAs("rescattering.png");
}
