TString indir = "/u/taniuchi/s467/rootfiles/Feb2023_FSnov22/";
TString outdir = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/mktree_2023/";
Bool_t verbose =false;

TChain *ch;// get evt as TChain
TTree *tree;
TString filename;
TFile* fout;
//Define output tree
Int_t Tpat;
TClonesArray* FrsDataCA, *MusicHitCA, *TwimHitCA, *SofTrackingCA, *Mwpc0CA, *Mwpc1CA, *Mwpc2CA, *Mwpc3CA, *RoluPosCA, *TofWHitCA;
Double_t FRSAoQ, FRSBeta, FRSBrho, FRSGamma, MusicZ;
Double_t MusicE, TwimE, TwimZ, FragZ, MusicTheta, TwimTheta, Mw0_X, Mw1_X, Mw2_X, Mw3_X, Mw0_Y, Mw1_Y, Mw2_Y, Mw3_Y, ROLU_X,ROLU_Y, Tofw_Y, FragTof, FragTof_corr, FragAoQ, FragAoQ_corr, FragBrho, FragBeta, FragGamma;
Int_t Tofw_Paddle;

Long64_t total_entry=0, cnt_tpat=0, cnt_frs=0, cnt_beta=0, cnt_rolu=0, cnt_rolu0=0;

void mktree_evt(int FRSset1 = 13, int FRSset2 = 13, int i_target=0, TString suffix = "12Feb");
Long64_t loadtree(TString infilename){
  
  TFile* fin = TFile::Open(infilename,"READ");
  TTree* ch = (TTree*)fin->Get("evt");
  //ch = new TChain("evt");
  Long64_t nentry = ch -> GetEntries();
  cout<<infilename << ", Entry: "<<nentry<<endl;

  ch->SetBranchAddress("EventHeader.fTpat",&Tpat);
  // Frs has array
  ch->SetBranchAddress("FrsData",&FrsDataCA);
  
  //  ch->SetBranchAddress("FrsData.fAq",FRSAoQ);
  //  ch->SetBranchAddress("FrsData.fBrho",FRSBrho);
  //  ch->SetBranchAddress("FrsData.fBeta",FRSBeta);
  //  ch->SetBranchAddress("FrsData.fZ", MusicZ);
  //  ch->SetBranchAddress("FrsData(0).Gamma",&FRSGamma);
  /*
  ch->SetBranchAddress("MusicHitData.fE", &MusicE);
  ch->SetBranchAddress("MusicHitData.fTheta", &MusicTheta);
  ch->SetBranchAddress("TwimHitData.fE", &TwimE);
  ch->SetBranchAddress("TwimHitData.fTheta", &TwimTheta);
  ch->SetBranchAddress("SofTrackingData.fZ", &FragZ);
  ch->SetBranchAddress("SofTrackingData.fBrho",&FragBrho);
  ch->SetBranchAddress("SofTrackingData.fAq", &FragAoQ);
  ch->SetBranchAddress("Mwpc0HitData.fX", &Mw0_X);
  ch->SetBranchAddress("Mwpc1HitData.fX", &Mw1_X);
  ch->SetBranchAddress("Mwpc2HitData.fX", &Mw2_X);
  ch->SetBranchAddress("Mwpc3HitData.fX", &Mw3_X);
  ch->SetBranchAddress("Mwpc0HitData.fY", &Mw0_Y);
  ch->SetBranchAddress("Mwpc1HitData.fY", &Mw1_Y);
  ch->SetBranchAddress("Mwpc2HitData.fY", &Mw2_Y);
  ch->SetBranchAddress("Mwpc3HitData.fY", &Mw3_Y);
  ch->SetBranchAddress("RoluPosData.fX", &ROLU_X);
  ch->SetBranchAddress("RoluPosData.fY", &ROLU_Y);
  ch->SetBranchAddress("TofWHitData.fY", &Tofw_Y);
  ch->SetBranchAddress("TofWHitData.fTof", &FragTof);
  ch->SetBranchAddress("TofWHitData.fPaddleId", &Tofw_Paddle);
  */
  ch->SetBranchAddress("MusicHitData",&MusicHitCA);
  ch->SetBranchAddress("TwimHitData",&TwimHitCA);
  ch->SetBranchAddress("SofTrackingData",&SofTrackingCA);
  ch->SetBranchAddress("Mwpc0HitData",&Mwpc0CA);
  ch->SetBranchAddress("Mwpc1HitData",&Mwpc1CA);
  ch->SetBranchAddress("Mwpc2HitData",&Mwpc2CA);
  ch->SetBranchAddress("Mwpc3HitData",&Mwpc3CA);
  ch->SetBranchAddress("RoluPosData",&RoluPosCA);
  ch->SetBranchAddress("TofWHitData",&TofWHitCA);
  
  cout << "Start Merging" <<endl;
  for(Long64_t n=0; n<nentry; n++){
    if((n+1)%10000==0 || verbose) {
      cout<<n+1<<" entry/"<<nentry<<endl;
      //ch->Show(n);
    }
    ch->GetEntry(n);
    if(verbose) cout<<n<<" tpat:"<<Tpat <<" frs:"<<FrsDataCA->GetEntriesFast()<<endl;
    if(Tpat!=1){
      if(verbose) cout<<"Tpat"<<endl;
      cnt_tpat++;
      continue; //only select non reacted
    }
    if(FrsDataCA->GetEntriesFast()!=4){
      if(verbose) cout<<"FRS"<<endl;
      cnt_frs++;
      continue;
    }
    auto frsdata = (R3BFrsData*)FrsDataCA->At(0);
    if(isnan(frsdata->GetZ())){
      if(verbose) cout<<"Beta"<<endl;
      cnt_beta++;
      continue;
    }
    if(RoluPosCA->GetEntriesFast()!=1){
      cnt_rolu++;
      continue;
    }
    auto roludata = (R3BMwpcHitData*)RoluPosCA->At(0);
    ROLU_X = roludata ->GetX();
    ROLU_Y = roludata ->GetY();
    if(ROLU_X < -3.){
      if(verbose) cout<<"ROLU"<<endl;
      cnt_rolu++;
      continue;
    }
    if(ROLU_X<0.){
      cnt_rolu0++;
    }
    //
    FRSAoQ = frsdata -> GetAq();
    FRSBrho= frsdata -> GetBrho();
    FRSBeta= frsdata -> GetBeta();
    MusicZ = frsdata -> GetZ();
    //
    if(MusicHitCA->GetEntriesFast()==1){
      auto musicdata = (R3BMusicHitData*)MusicHitCA->At(0);
      MusicE = musicdata->GetEave();
      MusicTheta = musicdata->GetTheta();
    }else{
      MusicE = NAN;
      MusicTheta = NAN;
    }
    if(TwimHitCA->GetEntriesFast()==1){
      auto twimdata = (R3BTwimHitData*)TwimHitCA->At(0);
      TwimE = twimdata->GetEave();
      TwimZ = twimdata->GetZcharge();
      TwimTheta = twimdata->GetTheta();
    }else{
      TwimE = NAN;
      TwimZ = NAN;
      TwimTheta = NAN;
    }
    //
    if(SofTrackingCA->GetEntriesFast()==1){
      auto trackingdata = (R3BSofTrackingData*)SofTrackingCA->At(0);
      FragZ = trackingdata->GetZ();
      FragBrho = trackingdata->GetBrho();
      FragAoQ = trackingdata->GetAq();
    }else{
      FragZ = NAN;
      FragBrho = NAN;
      FragAoQ = NAN;
    }
    //
    if(Mwpc0CA->GetEntriesFast()==1){
      auto mwpcdata = (R3BMwpcHitData*)Mwpc0CA->At(0);
      Mw0_X = mwpcdata->GetX();
      Mw0_Y = mwpcdata->GetY();
    }else{
      Mw0_X = NAN;
      Mw0_Y = NAN;
    }
    if(Mwpc1CA->GetEntriesFast()==1){
      auto mwpcdata = (R3BMwpcHitData*)Mwpc1CA->At(0);
      Mw1_X = mwpcdata->GetX();
      Mw1_Y = mwpcdata->GetY();
    }else{
      Mw1_X = NAN;
      Mw1_Y = NAN;
    }
    if(Mwpc2CA->GetEntriesFast()==1){
      auto mwpcdata = (R3BMwpcHitData*)Mwpc2CA->At(0);
      Mw2_X = mwpcdata->GetX();
      Mw2_Y = mwpcdata->GetY();
    }else{
      Mw2_X = NAN;
      Mw2_Y = NAN;
    }
    if(Mwpc3CA->GetEntriesFast()==1){
      auto mwpcdata = (R3BMwpcHitData*)Mwpc3CA->At(0);
      Mw3_X = mwpcdata->GetX();
      Mw3_Y = mwpcdata->GetY();
    }else{
      Mw3_X = NAN;
      Mw3_Y = NAN;
    }
    //
    if(TofWHitCA->GetEntriesFast()==1){
      auto tofwdata = (R3BSofTofWHitData*)TofWHitCA->At(0);
      Tofw_Y = tofwdata->GetY();
      FragTof = tofwdata->GetTof();
      Tofw_Paddle = tofwdata->GetPaddle();
    }else{
      Tofw_Y = NAN;
      FragTof = NAN;
      Tofw_Paddle = -1;
    }
    tree->Fill();
  }
  return nentry;
}

void mktree_evt(int FRSset1, int FRSset2, int i_target, TString suffix){
  TString FRSname, targetname;
  int posmin, posmax;
  int FRSsettingRange[2] = {FRSset1, FRSset2};
  if(FRSset1==13){
    FRSname = "50Ca";
  }else if(FRSset1==122){
    FRSname = "38Ca_122"; // Only the latter half
  }else if(FRSset1==11||FRSset1==12){
    FRSname = "38Ca";
  }else{
    cerr<<"Check if the FRSsetting number is correct." <<endl;
    return;
  }
  //
  if(i_target==0){
    targetname = "empty";
    posmin=1325;
    posmax=1424;
  }else if(i_target==1){
    targetname = "ch2-24mm";
    posmin=539;
    posmax=539;
  }else if(i_target==2){
    targetname = "carbon";
    posmin=362;
    posmax=362;
  }else if(i_target==2){
    targetname = "PP";
    posmin=893;
    posmax=893;
  }else{
    cerr<<"No target defined"<<endl;
    return;
  }
  //
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle();
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat();
  //TString pdfout;
  
  //
  FrsDataCA = new TClonesArray("R3BFrsData", 5);
  MusicHitCA = new TClonesArray("R3BMusicHitData",1);
  TwimHitCA = new TClonesArray("R3BTwimHitData",1);
  SofTrackingCA = new TClonesArray("R3BSofTrackingData",1);
  Mwpc0CA = new TClonesArray("R3BMwpcHitData",1);
  Mwpc1CA = new TClonesArray("R3BMwpcHitData",1);
  Mwpc2CA = new TClonesArray("R3BMwpcHitData",1);
  Mwpc3CA = new TClonesArray("R3BMwpcHitData",1);
  RoluPosCA = new TClonesArray("R3BMwpcHitData",1);
  TofWHitCA = new TClonesArray("R3BSofTofWHitData",1);

  TString outfilename=Form("%smktree_fragment_%s_%s.root",outdir.Data(),FRSname.Data(),targetname.Data());
  cout<<"Output file: "<<outfilename<<endl;
  TFile *fout = TFile::Open(outfilename,"RECREATE");
  //ch = new TChain("evt");
  tree = new TTree("Tree", "Tree");

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
  //
  tree->Branch("tpat",&Tpat);
  tree->Branch("FRSAoQ",&FRSAoQ);
  tree->Branch("FRSBrho",&FRSBrho);
  tree->Branch("FRSBeta",&FRSBeta);
  tree->Branch("FRS_Z", &MusicZ);
  //tree->Branch("FRSGamma",&FRSGamma);
  //tree->Branch("Frag_AoQ", &FragAoQ_corr);
  tree->Branch("Frag_AoQ", &FragAoQ);
  //tree->Branch("Frag_AoQ_original", &FragAoQ);
  tree->Branch("Frag_Brho",&FragBrho);
  tree->Branch("Frag_Z", &FragZ);
  //  tree->Branch("Frag_Tof", &FragTof_corr);
  //tree->Branch("Frag_Tof_original", &FragTof);
  tree->Branch("Frag_Tof", &FragTof);
  tree->Branch("Frag_Beta", &FragBeta);
  //tree->Branch("Frag_Gamma", &FragGamma);
  tree->Branch("MusicE", &MusicE);
  tree->Branch("MusicTheta", &MusicTheta);
  tree->Branch("TwimE", &TwimE);
  tree->Branch("TwimZ", &TwimZ);
  tree->Branch("TwimTheta", &TwimTheta);
  tree->Branch("Mw0_X", &Mw0_X);
  tree->Branch("Mw1_X", &Mw1_X);
  tree->Branch("Mw2_X", &Mw2_X);
  tree->Branch("Mw3_X", &Mw3_X);
  tree->Branch("Mw0_Y", &Mw0_Y);
  tree->Branch("Mw1_Y", &Mw1_Y);
  tree->Branch("Mw2_Y", &Mw2_Y);
  tree->Branch("Mw3_Y", &Mw3_Y);
  tree->Branch("ROLU_X", &ROLU_X);
  tree->Branch("Tofw_Y", &Tofw_Y);
  tree->Branch("Tofw_Paddle", &Tofw_Paddle);
  
  while(true){
    i++;
    if(i > 400 || !RunList.good()){
      //std::cerr << "No info for run found" <<std::endl;
      break;
    }

    RunList>>runnumcsv[i]>>dumchar>>FRSsetting[i]>>dumchar>>brhocsv[i]>>dumchar>>targetpos[i]>>dumchar>>musicgain[i]>>dumchar>>junk[i];

    if(FRSsetting[i]<FRSsettingRange[0]||FRSsetting[i]>FRSsettingRange[1]) continue;
    if(FRSsettingRange[1]==122 && FRSsetting[i] > 12 && FRSsetting[i] < 122) continue; // after run 262 removed, it should be fine to use both.
    if(junk[i]!=0) continue;
    if(targetpos[i]<posmin || targetpos[i]>posmax) continue;
    cout<<runnumcsv[i]<<" "<<dumchar<<" "<<FRSsetting[i]<<" "<<dumchar<<" "<<brhocsv[i]<<" "<<dumchar<<" "<<targetpos[i]<<" "<<dumchar<<" "<<musicgain[i]<<" "<<dumchar<<" "<<junk[i]<<endl;
    //if(FRSsetting[i]==122) FRSsetting[i]=12;
    filename = Form("%ss467_filltree_Setting%i_%04d_%s.root", indir.Data(), FRSsetting[i], runnumcsv[i], suffix.Data());
    cout<<"Input file: "<<filename<<endl;
    ///
    /// Load tree
    //TFile *fin = TFile::Open(filename,"READ");
    //ch = new TChain("evt");
    //ch -> AddFile(filename);
    total_entry += loadtree(filename);
    //cout<< ch->GetEntries() <<endl;
  }
  //
  cout<<"Finished:"<<total_entry<<endl;
  cout<<"Tpat invalid events (%)"<<100.*(double)cnt_tpat/(double)total_entry<<endl;
  cout<<"FRS invalid events (%)"<<100.*(double)cnt_frs/(double)total_entry<<endl;
  cout<<"Beta invalid events (%)"<<100.*(double)cnt_beta/(double)total_entry<<endl;
  cout<<"ROLU invalid events (%)"<<100.*(double)cnt_rolu/(double)total_entry<<endl;
  cout<<"ROLU < 0. events (%)"<<100.*(double)cnt_rolu0/(double)total_entry<<endl;
  fout->cd();
  cout<<"Output file: "<<outfilename<<endl;
  tree->Write();
  fout->Close();
  //ch->Merge(outfilename);
}
