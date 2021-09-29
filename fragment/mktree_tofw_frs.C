#define NUMPADDLE 28
#define NUMHIST 5
int posmin=1325, posmax=1424; TString targetname = "empty";
//int posmin=539, posmax=539; TString targetname = "ch2-24mm";
TProof *p =  TProof::Open("");

TCanvas *c;
TH1F *htime[NUMPADDLE], *hpos[NUMPADDLE];
TH2F *htimeoffset;
TF1 *func[NUMPADDLE], *funcpos[NUMPADDLE];

TChain *ch;//*chout;
//TTree *tree[500];
TEventList *el;
TString filename;
TFile* fout;
TH2F* hfrs[NUMHIST], *htofw[NUMHIST];

void mktree_tofw_frs(int runnum){
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle();
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat();
  //TString pdfout;
  ch = new TChain("Tree");
  Int_t firstrun = runnum>500?274:runnum;
  Int_t lastrun = (runnum>500?358+1:runnum);

  std::ifstream RunList("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/RunSummary.csv", std::ios::in);
  if(!RunList.is_open()) std::cerr <<"No run summary found\n";
  int runnumcsv[500], targetpos[500], musicgain[500], junk[500];
  int FRSsetting[500]; // calib:0, ToFCalib:6-8, 40Ca:9, 39Ca:10, 38Ca:11,12, 50Ca:13, ToFWcalib:14
  string dummyline;
  char dumchar;
  double brhocsv[500];
  std::getline (RunList, dummyline);
  //std::cout << dummyline << std::endl;
  Int_t i=0;
      
  while(true){
    i++;
    if(i > 400 || !RunList.good()){
      //std::cerr << "No info for run found" <<std::endl;
      break;
    }

    RunList>>runnumcsv[i]>>dumchar>>FRSsetting[i]>>dumchar>>brhocsv[i]>>dumchar>>targetpos[i]>>dumchar>>musicgain[i]>>dumchar>>junk[i];

    if(FRSsetting[i]!=13) continue;
    if(junk[i]!=0) continue;
    if(targetpos[i]<posmin || targetpos[i]>posmax) continue;
    cout<<runnumcsv[i]<<" "<<dumchar<<" "<<FRSsetting[i]<<" "<<dumchar<<" "<<brhocsv[i]<<" "<<dumchar<<" "<<targetpos[i]<<" "<<dumchar<<" "<<musicgain[i]<<" "<<dumchar<<" "<<junk[i]<<endl;

    filename = Form("/u/taniuchi/s467/rootfiles/rootfiletmp/fragment_Sep2021/s467_filltree_Setting13_%04d_28Sep.root", runnumcsv[i]);
    //ch -> Add(filename);
    ch -> AddFile(filename);
  }
  ch->Merge(Form("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/fragment/output/mktree_fragment_%s.root",targetname.Data()));
  /*
  for(int i = firstrun; i < lastrun+1; i++){
    filename = Form("/u/taniuchi/s467/rootfiles/rootfiletmp/TofW/s467_FRSTree_Setting13_%04d_ToFWhitpar.root", i);
    //ch -> Add(filename);
    ch -> AddFile(filename);
    //cout << "Filename: " << filename <<endl;
  }
  */
  //fout = new TFile("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/output/mktree_tofw_frs.root","recreate");
  /*
  ch -> SetProof();
  ch->SetBranchStatus("*",0);
  ch->SetBranchStatus("SofFrsData.fZ",1);
  ch->SetBranchStatus("SofFrsData.fAq",1);
  ch->SetBranchStatus("SofFrsData.fBeta",1);
  ch->SetBranchStatus("SofFrsData.fBrho",1);
  ch->SetBranchStatus("TofWHitData.*",1);
  //
  Double_t fZ[3] = {0.};
  ch->SetBranchAddress("SofFrsData.fZ", fZ);

  / *
  for(int i=0; i<NUMHIST; i++){
    hfrs[i] = new TH2F(Form("hfrs%i",i), Form("hfrs%i",i) , 500, 2.2, 2.8, 500, 10, 35);
    htofw[i] = new TH2F(Form("htofw%i",i), Form("htofw%i",i) , 30, -0.5, 29.5, 1000, -5, 5);
  }
  ch->Draw("fZ:fAq>>hfrs0","SofFrsData.fZ>0", "colz");
  ch->Draw("TofWHitData.fTime:TofWHitData.fPaddleId>>htofw0", "SofFrsData.fZ>0", "colz");
  ch->Draw("fZ:fAq>>hfrs1","SofFrsData.fZ>0 && abs(SofFrsData.fZ-20)<0.5", "colz");
  ch->Draw("TofWHitData.fTime:TofWHitData.fPaddleId>>htofw1", "SofFrsData.fZ>0 && abs(SofFrsData.fZ-20)<0.5", "colz");
  ch->Draw("fZ:fAq>>hfrs2","SofFrsData.fZ>0 && abs(SofFrsData.fZ-19)<0.5", "colz");
  ch->Draw("TofWHitData.fTime:TofWHitData.fPaddleId>>htofw2", "SofFrsData.fZ>0 && abs(SofFrsData.fZ-19)<0.5", "colz");
  ch->Draw("fZ:fAq>>hfrs3","SofFrsData.fZ>0 && abs(SofFrsData.fZ-21)<0.5", "colz");
  ch->Draw("TofWHitData.fTime:TofWHitData.fPaddleId>>htofw3", "SofFrsData.fZ>0 && abs(SofFrsData.fZ-21)<0.5", "colz");
  */
  /* /
  el = new TEventList("elist","elist");  
  ch -> Draw(">>elist","SofFrsData.fZ>0");
  cout << "Num Tree = " << ch->GetNtrees();//ch->GetTreeNumber();
  cout << ", Num Entry = " << ch->GetEntries()<<", Num Valid Entries = " <<el->GetN()<<", Valid ratio (%) = "<<(double)(el->GetN())/(double)(ch->GetEntries()) *100.<< endl;
  ch ->SetEventList(el);
  */
    
  /*
  for(int i = firstrun; i < lastrun+1; i++){
    tree[i] = new TTree(Form("tree%i",i), Form("tree%i",i));
    tree[i] = (ch->GetTree())->CopyTree("SofFrsData.fZ>0");
    //tree[i] = (ch->GetTree())->CloneTree(-1,"fast");
    cout<<tree[i] -> GetEntries() <<endl;
  }
  //ch->Merge("test.root");
  / * /
  for(Long64_t i = 0; i < 100; i++){
    //Long64_t ientry = LoadTree(el->GetEntry(i));
    //Long64_t ientry = ch->LoadTree(el->GetEntry(i));
    //ch->GetEntry(el->GetEntry(i));
    ch->GetEntry(i);
    //cout <<el->GetEntry(i) << " "<< fZ <<endl;
    cout <<fZ[0] <<endl;
  }
  * /
  //tree->Write();
  for(int i =0; i<NUMHIST; i++){
    hfrs[i] ->Write();
    htofw[i] ->Write();
  }
  */
  
  //ch->CloneTree()->Write();
  //fout->Close();


  
  //  c->Print(pdfout);

}


void mktree_tofw_frs(){
  mktree_tofw_frs(1000);
}
