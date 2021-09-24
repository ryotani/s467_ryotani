#define NUMPADDLE 28
TProof *p =  TProof::Open("");

TCanvas *c;
TH1F *htime[NUMPADDLE], *hpos[NUMPADDLE];
TH2F *htimeoffset;
TF1 *func[NUMPADDLE], *funcpos[NUMPADDLE];

//TTree *tree;
//TFile *f;
TChain *ch;

void draw_tofw_ftime(int runnum){
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle();
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat();
  TString filename = runnum<1000?Form("/u/taniuchi/s467/rootfiles/rootfiletmp/s467_FRSTree_Setting*%04d_ToFWVFTXpar.root",  runnum): "/u/taniuchi/s467/rootfiles/rootfiletmp/s467_FRSTree_Setting*_ToFWVFTXpar.root";
  TString pdfout =   runnum<1000?Form("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/output/tofw_ftime_%04d.pdf", runnum): "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/output/tofw_ftime.pdf";
  
  //TFile *f= TFile::Open(filename, "READ");
  //tree = (TTree*)(f->Get("evt"));
  ch = new TChain("evt");
  ch -> Add(filename);
  ch -> SetProof();
  cout << "Filename: " << filename <<endl;
  cout << "Num Entry = " << ch->GetEntries()<<endl;
  c = new TCanvas("c","c",1200,1000);
  c -> Divide(6,5);
  for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1);
    htime[i] = new TH1F(Form("h%i",i), Form("Raw Tof in paddle ID%i", i+1), 20000, -200, 200);
    //ch->Draw(Form("TofWHitData.fTime>>h%i",i),Form("TofWHitData.fPaddleId==%i",i+1));
    ch->Draw(Form("SofTofWSingleTcalData.fRawTofNs>>h%i",i),Form("SofTofWSingleTcalData.fDetector==%i",i+1));
    func[i] = new TF1(Form("func%i",i), "gaus", -20, 20);
    //func[i] = new TF1(Form("func%i",i), "[0]*exp(-0.5*((x-[1])/[2])**2)", -50, 50);
    func[i]->SetParameter(0, htime[i]->GetMaximum());
    func[i]->SetParameter(1, htime[i]->GetBinCenter(htime[i]->GetMaximumBin()));
    func[i]->SetParameter(2, 0.01);
    htime[i]->Fit(func[i], "LL", "", htime[i]->GetBinCenter(htime[i]->GetMaximumBin()-5), htime[i]->GetBinCenter(htime[i]->GetMaximumBin()+5));
    //htime[i]->Fit(func[i], "");
    htime[i]->GetXaxis()->SetRangeUser(func[i]->GetParameter(1)-10.*func[i]->GetParameter(2),func[i]->GetParameter(1)+10.*func[i]->GetParameter(2));
    htime[i]->Fit(func[i], "LL", "", func[i]->GetParameter(1)-5.*func[i]->GetParameter(2), func[i]->GetParameter(1)+5.*func[i]->GetParameter(2));
  }    
  c->Print(pdfout+"[");
  c->Print(pdfout);

  /////// Same for positions
    for(int i = 0 ; i<NUMPADDLE; i++){
    c->cd(i+1);
    hpos[i] = new TH1F(Form("hpos%i",i), Form("Raw position in paddle ID%i", i+1), 20000, -200, 200);
    //ch->Draw(Form("TofWHitData.fTime>>h%i",i),Form("TofWHitData.fPaddleId==%i",i+1));
    ch->Draw(Form("SofTofWSingleTcalData.fRawPosNs>>hpos%i",i),Form("SofTofWSingleTcalData.fDetector==%i",i+1));
    funcpos[i] = new TF1(Form("func%i",i), "gaus", -20, 20);
    //funcpos[i] = new TF1(Form("func%i",i), "[0]*exp(-0.5*((x-[1])/[2])**2)", -50, 50);
    funcpos[i]->SetParameter(0, hpos[i]->GetMaximum());
    funcpos[i]->SetParameter(1, hpos[i]->GetBinCenter(hpos[i]->GetMaximumBin()));
    funcpos[i]->SetParameter(2, 0.01);
    hpos[i]->Fit(funcpos[i], "LL", "", hpos[i]->GetBinCenter(hpos[i]->GetMaximumBin()-5), hpos[i]->GetBinCenter(hpos[i]->GetMaximumBin()+5));
    //hpos[i]->Fit(funcpos[i], "");
    hpos[i]->GetXaxis()->SetRangeUser(funcpos[i]->GetParameter(1)-10.*funcpos[i]->GetParameter(2),funcpos[i]->GetParameter(1)+10.*funcpos[i]->GetParameter(2));
    hpos[i]->Fit(funcpos[i], "LL", "", funcpos[i]->GetParameter(1)-5.*funcpos[i]->GetParameter(2), funcpos[i]->GetParameter(1)+5.*funcpos[i]->GetParameter(2));
  }

  //
  for(int i = 0 ; i<NUMPADDLE; i++){ 
    cout << "Paddle" << i+1 << ", Entries:" << htime[i]->GetEntries() <<", Mean Tof:" << func[i]->GetParameter(1) ;
    cout << ", Mean Pos:" << funcpos[i]->GetParameter(1) << endl;
  }
  c->Print(pdfout);
  c->Print(pdfout+"]");
  
  ////////
  TString outputFileNamePar = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/tofw_hit.par";
  ofstream fout;
  fout.open(outputFileNamePar);
  fout<<      "##############################################################################"
      <<endl<<"# Class:   tofwHitPar"
      <<endl<<"# Context: TofWHitParContext"
      <<endl<<"##############################################################################"
      <<endl<<"[tofwHitPar]"
      <<endl<<"//----------------------------------------------------------------------------"
      <<endl;
  fout<<"tofwHitSciPar:  Int_t  "<<(Int_t)NUMPADDLE << endl;
  fout<<"tofwInUsePar:  Int_t  \\"<<endl;
  for(int i = 0; i<NUMPADDLE; i++)  fout<<" "<<((i!=27)?1:0);
  fout<<endl;
  fout<<"tofwPosPar:  Float_t  \\"<<endl;
  for(int i = 0; i<NUMPADDLE; i++)  fout<<" "<<funcpos[i]->GetParameter(1)<<(((i%10)==9&&i+1!=NUMPADDLE)?" \\\n":"");
  fout<<endl;
  fout<<"tofwTofPar:  Float_t  \\"<<endl;
  for(int i = 0; i<NUMPADDLE; i++)  fout<<" "<<func[i]->GetParameter(1)<<(((i%10)==9&&i+1!=NUMPADDLE)?" \\\n":"");
  fout<<endl;
  fout<<"##############################################################################"<<endl;
  /*
  //FairRunOnline* run = new FairRunOnline();
  //run->SetRunId(1);
  //run->Init();
  
  FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
  /// * Runtime data base ------------------------------------ * /
  //FairRuntimeDb* rtdb = run->GetRuntimeDb();
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
  
  //
    rtdb->setContainersStatic(kFALSE);//initContainers();

  if (!rtdb)
    {
      cerr << "No RTDB" <<endl;
      return;
    }
  
  auto fHit_Par = (R3BSofTofWHitPar*)rtdb->getContainer("tofwHitPar");
  if (!fHit_Par)
    {
      cerr << "R3BSofTofWTcal2HitPar:: Couldn't get handle on tofwHitPar container" <<endl;
      return ;
    }
// Ascii file with the Calibration Parameters
  FairParAsciiFileIo* parOut = new FairParAsciiFileIo();
  TString outputFileNamePar = "/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/tofw_test.par";
  parOut->open(outputFileNamePar,"out");
  rtdb->setOutput(parOut);
  //rtdb->addContFactory();
  
  //
  fHit_Par->SetNumSci(NUMPADDLE);
  fHit_Par->SetInUse(1,1);
  
  fHit_Par->setChanged();
  fHit_Par->printParams();
  //FairParamList* list;
  //fHit_Par->putParams(list);
  //

  rtdb->getContFactory("R3BSofTofWContFact");
  rtdb->print();  
  rtdb->saveOutput();
  rtdb->writeContainers();
  rtdb->writeVersions();
  rtdb->saveOutput();
  rtdb->closeOutput();
  */
}


void draw_tofw_ftime(){
  draw_tofw_ftime(366);
}
