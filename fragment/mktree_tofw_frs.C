TString indir = "/u/taniuchi/s467/rootfiles/rootfile_land/202002_s467/";
TString outdir = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/rootfiles/rootfile_land/mktree/";

TChain *ch;
TString filename;
TFile* fout;

void mktree_tofw_frs(int FRSset1 = 13, int FRSset2 = 13, int i_target=0, TString suffix = "8Feb");

void mktree_tofw_frs(int FRSset1, int FRSset2, int i_target, TString suffix){
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
  ch = new TChain("Tree");

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

    if(FRSsetting[i]<FRSsettingRange[0]||FRSsetting[i]>FRSsettingRange[1]) continue;
    if(junk[i]!=0) continue;
    if(targetpos[i]<posmin || targetpos[i]>posmax) continue;
    cout<<runnumcsv[i]<<" "<<dumchar<<" "<<FRSsetting[i]<<" "<<dumchar<<" "<<brhocsv[i]<<" "<<dumchar<<" "<<targetpos[i]<<" "<<dumchar<<" "<<musicgain[i]<<" "<<dumchar<<" "<<junk[i]<<endl;
    if(FRSsetting[i]==122) FRSsetting[i]=12;
    filename = Form("%ss467_filltree_Setting%i_%04d_%s.root", indir.Data(), FRSsetting[i], runnumcsv[i], suffix.Data());
    cout<<"Input file: "<<filename<<endl;
    ch -> AddFile(filename);
  }
  TString outfilename=Form("%smktree_fragment_%s_%s.root",outdir.Data(),FRSname.Data(),targetname.Data());
  cout<<"Output file: "<<outfilename<<endl;
  ch->Merge(outfilename);
}
