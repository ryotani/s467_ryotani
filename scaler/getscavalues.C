void getscavalues() {
  std::ifstream RunList("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/RunSummary.csv", std::ios::in);
  if(!RunList.is_open()) std::cerr <<"No run summary found\n";
  int runnumcsv[500], targetpos[500], musicgain[500], junk[500];
  int FRSsetting[500]; // calib:0, ToFCalib:6-8, 40Ca:9, 39Ca:10, 38Ca:11,12, 50Ca:13, ToFWcalib:14
  string dummyline;
  char dumchar;
  double brhocsv[500];
  Int_t i=-1;
  //
  std::ofstream Outcsv("/u/taniuchi/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/scaler/result.csv", ios::out);//|ios::app);
  Outcsv << "runnumcsv[i], FRSsetting[i], numoffspill, numonspill, numtpat1, numtpat0musicz, numtpat0scicave, numtpat0mastercave, "
	 << "numtpat0trig1, numtpat0trig2, numtpat0trig3, numtpat0trig1mastercave, numtpat0trig2mastercave, numtpat0trig3mastercave, SciL, SciR, Master" << std::endl;
  //
  std::getline (RunList, dummyline);      
  while(true){
    i++;
    RunList>>runnumcsv[i]>>dumchar>>FRSsetting[i]>>dumchar>>brhocsv[i]>>dumchar>>targetpos[i]>>dumchar>>musicgain[i]>>dumchar>>junk[i];
    //std::cout<<runnumcsv[i]<<dumchar<<FRSsetting[i]<<dumchar<<brhocsv[i]<<dumchar<<targetpos[i]<<dumchar<<musicgain[i]<<dumchar<<junk[i]<<std::endl;
    //if(FRSsetting[i] == 13 && junk[i] == 0){
    if(FRSsetting[i] != 13 && junk[i] == 0){
      TString RootFilename = Form("./rootfiles/rootfiletmp/s467_FRSTree_Setting%i_%04d_FRSTree.root", FRSsetting[i], runnumcsv[i]);
      TFile *file = new TFile(RootFilename, "read");
      if(!file->IsOpen()){
	std::cerr<<"File not found: "<<RootFilename <<std::endl;
	continue;
      }
      //
      std::cout<<"Opening root file: "<<RootFilename<<std::endl;
      TTree *FrsTree = (TTree*)file->Get("FrsTree");
      TTree *EvtTree = (TTree*)file->Get("evt");
      if(!FrsTree || !EvtTree)
	continue;
      TH1I *ScaHist = (TH1I*)file->Get("SofScalers1_GeneralView");
      Int_t numoffspill = EvtTree->Draw(">>test","R3BEventHeader.fTpat > 100");
      Int_t numonspill = FrsTree->Draw(">>test","tpat < 100");
      if(numonspill != EvtTree->GetEntries() - numoffspill)//EvtTree->Draw(">>test","R3BEventHeader.fTpat < 100"))
	{
	  std::cerr<<"Mismatch btw FrsTree and evt tree" <<endl;
	  continue;
	}
      Int_t numtpat0 = FrsTree->Draw(">>tpat0","tpat == 0");
      Int_t numtpat1 = numonspill - numtpat0;//= FrsTree->Draw(">>test","tpat&1");
      //auto *evtlist = (TEventList*)gROOT->FindObject("tpat0");
      //FrsTree->SetEventList(evtlist);
      Int_t numtpat0musicz = FrsTree->Draw(">>test2","tpat == 0 && MusicZ>0");
      Int_t numtpat0scicave = FrsTree->Draw(">>test2","tpat == 0 && multMapSci[10]* multMapSci[9]>0");
      Int_t numtpat0mastercave = FrsTree->Draw(">>test2","tpat == 0 && multMapSci[11]>0");
      Int_t numtpat0trig1 = FrsTree->Draw(">>test2","tpat == 0 && trigger==1");
      Int_t numtpat0trig2 = FrsTree->Draw(">>test2","tpat == 0 && trigger==2");
      Int_t numtpat0trig3 = FrsTree->Draw(">>test2","tpat == 0 && trigger>2");
      Int_t numtpat0trig1mastercave = FrsTree->Draw(">>test2","tpat == 0 && multMapSci[11]>0 && trigger==1");
      Int_t numtpat0trig2mastercave = FrsTree->Draw(">>test2","tpat == 0 && multMapSci[11]>0 && trigger==2");
      Int_t numtpat0trig3mastercave = FrsTree->Draw(">>test2","tpat == 0 && multMapSci[11]>0 && trigger>2");
      //
      Int_t SciL = ScaHist->GetBinContent(1);
      Int_t SciR = ScaHist->GetBinContent(2);
      Int_t Master = ScaHist->GetBinContent(3);
      //
      //std::cout<<runnumcsv[i]<<" "<<numoffspill<<" "<<numonspill<<" "<<numtpat1<<" "<<numtpat0musicz<<" "<<numtpat0scicave<<" "<<numtpat0mastercave<<" "<<
      //	numtpat0trig1<<" "<<numtpat0trig2<<" "<<numtpat0trig3<<" "<<numtpat0trig1mastercave<<" "<<numtpat0trig2mastercave<<" "<<numtpat0trig3mastercave<<" "<<SciL<<" "<<SciR<<" "<<Master<<endl;
      Outcsv<<runnumcsv[i]<<", "<<FRSsetting[i]<<", "<<numoffspill<<", "<<numonspill<<", "<<numtpat1<<", "<<numtpat0musicz<<", "<<numtpat0scicave<<", "<<numtpat0mastercave<<", "<<
	numtpat0trig1<<", "<<numtpat0trig2<<", "<<numtpat0trig3<<", "<<numtpat0trig1mastercave<<", "<<numtpat0trig2mastercave<<", "<<numtpat0trig3mastercave<<", "<<SciL<<", "<<SciR<<", "<<Master<<endl;
      //
      delete file;
    } else {
      //std::cout<<"Skipping the number: "<<runnumcsv[i] << std::endl;
    }


    if(i > 400 || !RunList.good()){
      std::cout << "Finish analysis" <<std::endl;
      return;
    }
  }
}

