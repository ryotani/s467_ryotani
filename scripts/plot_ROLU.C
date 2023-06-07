TString in_dir = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/mktree_2023/";
TString in_filename = "mktree_califa_May23_noROLUcut_FRS_TARGET.root";
TString out_filename = "plot_ROLU.png";
TString FRS[2]={"38Ca", "50Ca"};
TString TARGET[3]={"empty","carbon","ch2"};
TFile *in_file[6];
TCanvas *c;

void plot_ROLU(){
  c = new TCanvas("c","",1200,1000);
  c->Divide(3,2);
  for(int i=0; i<6; i++){
    TString dummy(in_dir + in_filename);
    dummy.ReplaceAll("FRS",FRS[i/3]);
    dummy.ReplaceAll("TARGET",TARGET[i%3]);
    cout<<"Opening: " << dummy<<endl;
    in_file[i] = TFile::Open(dummy,"READ");
    c->cd(i+1);
    //
    auto* Tree = dynamic_cast<TTree*>(in_file[i]->Get("Tree"));
    TString title = "Z=20 fragments for "+ FRS[i/3] + " setting with " + TARGET[i%3] + " target;"
      + "Time of Flight of fragments / ns;ROLU X position / mm";
    TH2F* hist = new TH2F(Form("h%i",i),title, 500,35,65,500,-10,20);
    hist->SetTitleOffset(1.,"Y");
      //
    gStyle->SetTitleFont(43,"T");
    gStyle->SetTitleSize(15,"T");
    gStyle->SetTitleAlign(13);
    Tree->Draw(Form("ROLU_X:Frag_Tof>>h%i",i), "pow(FRS_Z-20.,2.)<pow(3.*0.111,2.) && abs(Frag_Z-20.)<0.4","colz");
  }
  
  cout<<"End"<<endl;
}
