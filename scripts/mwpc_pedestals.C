#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>

auto*c = new TCanvas();
TH1D *hist_tmp;
TH2D *hist2_tmp;
constexpr int E_RANGE = 4000;
TString outdir="/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/scripts/";
TString outpdf=outdir + "mwpc_pedestals_ID.pdf";
TString outroot = outdir + "mwpc_pedestals_ID.root";
TString in_paramfile = outdir + "../parameters/common.par";
TString out_paramfile = outdir + "common_new_ID.par";
vector<TH1D*> h1_mw0, h1_mw1, h1_mw2, h1_mw3;
// Channels to be excluded (id starts from 0)
vector<int> exclude0 = {};
vector<int> exclude1 = {};//23, 42, 23+64, 42+64};
vector<int> exclude2 = {};
vector<int> exclude3 = {};//16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36,38,40,42,44,46,47,48};
//,49,50,51,52,53,54,55,56,57,58,59,60,61,62};
//
vector<TString> mwname = {"mwpc0", "mwpc1", "mwpc2", "mwpc3"};
vector<int> xpads = {64, 128, 128, 288};
vector<int> ypads = {64, 40, 40, 120};
vector<Double_t> pedestal0, pedestal1, pedestal2, pedestal3;
vector<TString> par_template=
  {"##############################################################################",
   "# Class:   MWPCCalPar",
   "# Context: MWPCCalParContext",
   "##############################################################################",
   "[MWPCCalPar]",
   "//----------------------------------------------------------------------------",
   "MWPCCalPar:  Float_t				\\",
   "MWPCPadXNumberPar:  Int_t  XPADS ",
   "MWPCPadYNumberPar:  Int_t  YPADS ",
   "MWPCParamsFitPar:  Int_t  2 ",
   "##############################################################################"
  };

  
int mwpc_pedestals(int set_id = 13) {
  auto* fout = TFile::Open(outroot,"RECREATE");
  // Set the number of threads for multi-threading
  ROOT::EnableImplicitMT(10);
    
  // Input file pattern using wildcard
  std::string inputFilePattern = Form("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/mwpc_cal_Jul2023/s467_filltree_Setting%i*_20Jul.root", set_id);
  outpdf.ReplaceAll("ID",TString::Itoa(set_id,10));
  outroot.ReplaceAll("ID",TString::Itoa(set_id,10));
  out_paramfile.ReplaceAll("ID",TString::Itoa(set_id,10));

  // Create an RDataFrame for reading the input trees
  ROOT::RDataFrame rdf("evt", inputFilePattern);

  auto hist0 = rdf.Define("i", "Mwpc0MappedData.fPad - 1 + 64* ((Mwpc0MappedData.fPlane+1)/2 - 1)").Histo2D({"h2_mw0", "Mwpc0", E_RANGE, -0.5, E_RANGE-0.5, 64*2, 0, 64*2},"Mwpc0MappedData.fCharge", "i");
  auto hist1 = rdf.Define("i", "Mwpc1MappedData.fPad - 1 + 64* (Mwpc1MappedData.fPlane-1)").Histo2D({"h2_mw1", "Mwpc1", E_RANGE, -0.5, E_RANGE-0.5, 64*2 + 40, 0, 64*2 + 40},"Mwpc1MappedData.fCharge", "i");
  auto hist2 = rdf.Define("i", "Mwpc2MappedData.fPad - 1 + 64* (Mwpc2MappedData.fPlane-1)").Histo2D({"h2_mw2", "Mwpc2", E_RANGE, -0.5, E_RANGE-0.5, 64*2 + 40, 0, 64*2 + 40},"Mwpc2MappedData.fCharge", "i");
  auto hist3 = rdf.Define("i", "Mwpc3MappedData.fPad - 1 + 288* ((Mwpc3MappedData.fPlane+1)/2 - 1)").Histo2D({"h2_mw3", "Mwpc3", E_RANGE, 0, E_RANGE-0.5, 288 + 120, 0, 288 + 120},"Mwpc3MappedData.fCharge", "i");
  //
  auto hcal0 = rdf.Define("i", "Mwpc0CalData.fPad - 1 + 64* ((Mwpc0CalData.fPlane+1)/2 - 1)").Histo2D({"hcal_mw0", "Mwpc0", E_RANGE, -0.5, E_RANGE-0.5, 64*2, 0, 64*2},"Mwpc0CalData.fCharge", "i");
  auto hcal1 = rdf.Define("i", "Mwpc1CalData.fPad - 1 + 64* (Mwpc1CalData.fPlane-1)").Histo2D({"hcal_mw1", "Mwpc1", E_RANGE, -0.5, E_RANGE-0.5, 64*2 + 40, 0, 64*2 + 40},"Mwpc1CalData.fCharge", "i");
  auto hcal2 = rdf.Define("i", "Mwpc2CalData.fPad - 1 + 64* (Mwpc2CalData.fPlane-1)").Histo2D({"hcal_mw2", "Mwpc2", E_RANGE, -0.5, E_RANGE-0.5, 64*2 + 40, 0, 64*2 + 40},"Mwpc2CalData.fCharge", "i");
  auto hcal3 = rdf.Define("i", "Mwpc3CalData.fPad - 1 + 288* ((Mwpc3CalData.fPlane+1)/2 - 1)").Histo2D({"hcal_mw3", "Mwpc3", E_RANGE, 0, E_RANGE-0.5, 288 + 120, 0, 288 + 120},"Mwpc3CalData.fCharge", "i");
  //
  auto hhit0 = rdf.Histo2D({"hhit0","hhit0", 240, -120, 120, 240, -120, 120}, "Mwpc0HitData.fX", "Mwpc0HitData.fY");
  auto hhit1 = rdf.Histo2D({"hhit1","hhit1", 240, -120, 120, 240, -120, 120}, "Mwpc1HitData.fX", "Mwpc1HitData.fY");
  auto hhit2 = rdf.Histo2D({"hhit2","hhit2", 240, -120, 120, 240, -120, 120}, "Mwpc2HitData.fX", "Mwpc2HitData.fY");
  auto hhit3 = rdf.Histo2D({"hhit3","hhit3", 420, -420, 420, 240, -120, 120}, "Mwpc3HitData.fX", "Mwpc3HitData.fY");
  //
  /*
  auto hrolu0x = rdf.Histo2D({"hrolu0x","hrolu0x", 240, -120, 120, 240, -120, 120}, "RoluPosData.fX", "Mwpc0HitData.fX");
  auto hrolu1x = rdf.Histo2D({"hrolu1x","hrolu1x", 240, -120, 120, 240, -120, 120}, "RoluPosData.fX", "Mwpc1HitData.fX");
  auto hrolu2x = rdf.Histo2D({"hrolu2x","hrolu2x", 240, -120, 120, 240, -120, 120}, "RoluPosData.fX", "Mwpc2HitData.fX");
  auto hrolu3x = rdf.Histo2D({"hrolu3x","hrolu3x", 420, -420, 420, 240, -120, 120}, "RoluPosData.fX", "Mwpc3HitData.fX");
  auto hrolu0y = rdf.Histo2D({"hrolu0y","hrolu0y", 240, -120, 120, 240, -120, 120}, "RoluPosData.fX", "Mwpc0HitData.fY");
  auto hrolu1y = rdf.Histo2D({"hrolu1y","hrolu1y", 240, -120, 120, 240, -120, 120}, "RoluPosData.fX", "Mwpc1HitData.fY");
  auto hrolu2y = rdf.Histo2D({"hrolu2y","hrolu2y", 240, -120, 120, 240, -120, 120}, "RoluPosData.fX", "Mwpc2HitData.fY");
  auto hrolu3y = rdf.Histo2D({"hrolu3y","hrolu3y", 420, -420, 420, 240, -120, 120}, "RoluPosData.fX", "Mwpc3HitData.fY");
  */
  /*
  // Plot max pad charges // Hmm difficult...
  // Define a lambda function to get the maximum fCharge for each event
  auto getMaxCharge = [](TClonesArray* fChargeArray) {
			double maxCharge = -1.0;
			for (Int_t i = 0; i < fChargeArray->GetEntries(); ++i) {
			  Double_t charge = dynamic_cast<R3BMwpcCalData*>(fChargeArray->At(i))->GetQ();
			  if (charge > maxCharge) {
			    maxCharge = charge;
			  }
			}
			return maxCharge;
		      };

  // Define a new column with the maximum fCharge for each event
  auto rdfMaxCharge = rdf.Define("maxCharge", getMaxCharge, {"Mwpc1CalData"});

  */  
  /*
  auto rdfModified = rdf.Define("maxCharge", [](const std::vector<double>& v_charge) { //[](const std::vector<R3BMwpcCalData>& calData) {
					       Double_t maxCharge = 0;
					       Int_t maxID = -1;
					       for (const auto& charge : v_charge) {
						 if (charge > maxCharge) {
						   maxCharge = charge;
						   //maxID = data.GetPad();
						 }
					       }
					       return maxCharge; //std::make_pair(maxCharge, maxID);
					     }, {"Mwpc0CalData.fCharge"});

  // Create a 2D histogram
  //TH2D hist("histogram", "Maximum Charge vs ID", 100, 0, 100, 100, 0, 1000);

  // Fill the histogram using RDataFrame with the modified column
  auto hist = rdfModified.Histo2D({"hcal_mw0", "hcal_mw0",  E_RANGE, -0.5, E_RANGE-0.5, 64*2, 0, 64*2}, "maxCharge.first", "maxCharge.second");
  hist->Draw();
  */
  //
  c->Divide(2,2);
  c->cd(1); hist0->Draw("colz");
  c->cd(2); hist1->Draw("colz");
  c->cd(3); hist2->Draw("colz");
  c->cd(4); hist3->Draw("colz");
  hist0->Write();
  hist1->Write();
  hist2->Write();
  hist3->Write();
  //c->Print(outpdf + "[");
  c->Print(outpdf + "(");
  //
  c->cd(1); hcal0->Draw("colz");
  c->cd(2); hcal1->Draw("colz");
  c->cd(3); hcal2->Draw("colz");
  c->cd(4); hcal3->Draw("colz");
  hcal0->Write();
  hcal1->Write();
  hcal2->Write();
  hcal3->Write();
  c->Print(outpdf);
  //
  c->cd(1); hhit0->Draw("colz");
  c->cd(2); hhit1->Draw("colz");
  c->cd(3); hhit2->Draw("colz");
  c->cd(4); hhit3->Draw("colz");
  hhit0->Write();
  hhit1->Write();
  hhit2->Write();
  hhit3->Write();
  c->Print(outpdf);
  /*
  c->cd(1); hrolu0x->Draw("colz");
  c->cd(2); hrolu1x->Draw("colz");
  c->cd(3); hrolu2x->Draw("colz");
  c->cd(4); hrolu3x->Draw("colz");
  hrolu0x->Write();
  hrolu1x->Write();
  hrolu2x->Write();
  hrolu3x->Write();
  c->Print(outpdf);
  //
  c->cd(1); hrolu0y->Draw("colz");
  c->cd(2); hrolu1y->Draw("colz");
  c->cd(3); hrolu2y->Draw("colz");
  c->cd(4); hrolu3y->Draw("colz");
  hrolu0y->Write();
  hrolu1y->Write();
  hrolu2y->Write();
  hrolu3y->Write();
  c->Print(outpdf);
  / */
  auto get_slice = [&](auto hist, vector<TH1D*> &v) {
		     for(int p=0; p<hist->GetNbinsY(); p++){
		       //cout<<p<<" "<<endl;
		       auto h = hist->ProjectionX(Form("%s_%i",hist->GetName(),p), p+1, p+1);
		       //h->GetXaxis()->SetRangeUser(0,60);
		       h->SetTitle(Form("%s pad%i",hist->GetTitle(),p));
		       v.push_back(h);
		     }
		   };
  get_slice(hist0, h1_mw0);
  get_slice(hist1, h1_mw1);
  get_slice(hist2, h1_mw2);
  get_slice(hist3, h1_mw3);
  //
  auto fit_slice = [&](vector<TH1D*> &v, vector<Double_t> &v_ped) {
		     auto p=0;
		     for(auto h: v){
		       c->cd(p%4+1);
		       auto max=-1;
		       auto pedestal= 0.;
		       auto fun_name = Form("f_%s",h->GetName());
		       unique_ptr<TF1> f  = make_unique<TF1>(fun_name, "gaus", 0, 50);
		       f->SetParLimits(2,0.,3.);
		       for(int b=0; b<h->GetNbinsX(); b++){
			 auto val = h->GetBinContent(b);
			 if(val < h->GetMaximum()/10.) continue;
			 if(val > max){
			   max = val;
			 }else{
			   pedestal = h->GetBinCenter(b-1); // The last one was the largest
			   //cout<<"Pedestal"<<pedestal<<" bin"<<b<<" ";
			   h->GetXaxis()->SetRangeUser(0,200);//2.5*b);
			   h->Fit(f.get(),"Q WW","",pedestal-8.,pedestal+1.5);
			   pedestal = f->GetParameter(1);
			   break;
			 }
		       }
		       v_ped.push_back(pedestal);
		       //cout<<h->GetTitle()<<endl;
		       h->Write();
		       h->Draw();
		       if(p%4==3 || h == v.back())
			 c->Print(outpdf);
		       //
		       p++;
		     }
		   };
  fit_slice(h1_mw0, pedestal0);
  fit_slice(h1_mw1, pedestal1);
  fit_slice(h1_mw2, pedestal2);
  fit_slice(h1_mw3, pedestal3);
  //
  //
  ofstream pout;
  pout.open(out_paramfile, ios::out);
  auto write_params = [&](int id, vector <Double_t> &v_ped, vector <Int_t> &v_ex) {
			auto replace_string = [&](TString &strline){
						strline.ReplaceAll("MWPC", mwname.at(id));
						strline.ReplaceAll("XPADS", TString::Itoa(xpads.at(id),10));
						strline.ReplaceAll("YPADS", TString::Itoa(ypads.at(id),10));
					      };
			for(int i = 0; i<7; i++){
			  TString dummy = par_template.at(i);
			  replace_string(dummy);
			  pout<<dummy<<endl;
			  cout<<dummy<<endl;
			}
			int p = -1;
			cout<<"  ";
			for(auto ped: v_ped){
			  p++;
			  if( std::find(v_ex.begin(), v_ex.end(), p) != v_ex.end())
			    {	
			      pout<<"9999 0 ";
			      cout<<"9999 0 ";
			    }
			  else
			    {
			      pout<<v_ped.at(p) << " 0 ";
			      cout<<v_ped.at(p) << " 0 ";
			    }
			  if(p%5==4){
			    pout<<" \\"<<endl<<"  ";
			    cout<<" \\"<<endl<<"  ";			    
			  }
			}
			pout<<endl;
			cout<<endl;
			for(int i = 7; i < 11; i++)
			  {
			    TString dummy = par_template.at(i);
			    replace_string(dummy);
			    pout<<dummy<<endl;
			    cout<<dummy<<endl;
			  }
		      };
  write_params(0, pedestal0, exclude0);
  write_params(1, pedestal1, exclude1);
  write_params(2, pedestal2, exclude2);
  write_params(3, pedestal3, exclude3);
  pout.close();
  cout<<"Parameter: " <<out_paramfile<<endl;
  c->Print(outpdf + "]");
  return 0;
}



//evt->Draw("Mwpc1MappedData[].fPad:Mwpc1MappedData[].fCharge>>(4000,0,4000,64,0,64)","Mwpc1MappedData[].fPlane==1 ","colz")


