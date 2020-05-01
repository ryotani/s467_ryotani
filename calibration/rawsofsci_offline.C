struct EXT_STR_h101_t
{
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_SOFSCI_onion_t sci;
};

void rawsofsci_offline(int runnum)
{
  TStopwatch timer;
  timer.Start();

  const Int_t nev = -1; /* number of events to read, -1 - until CTRL+C */

  // *********************************** //
  // PLEASE COMPLETE THE FOLLOWING LINES //
  // *********************************** //
  const Int_t expId = 467; // select experiment: 444 or 467
  // *********************************** //

  // NumSofSci -----------------------------------------
  UShort_t NumSofSci, IdCaveC, IdS2, IdS8;
  if (expId==444){// s444: PRIMARY BEAM EXP, 1 SofSci at CAVE C ONLY
    NumSofSci = 1; 
    IdS2 = 0;
    IdS8 = 0;
  }
  else if(expId==467){
    NumSofSci = 4; // s467: SECONDARY BEAM EXP, 2 at S2, 1 at S8, 1 at CAVE C
    IdS2 = 2;
    IdS8 = 3;
  }
  else{
    NumSofSci = 2; // default: SECONDARY BEAM EXP, 1 at S2, 1 at CAVE C
  }
  IdCaveC = NumSofSci; // cave C is the last SofSci detector

  // --- ----------------------------------- --- //
  // --- Create source using ucesb for input --- //
  // --- ----------------------------------- --- //

  // Create input -----------------------------------------
  TString filename;
  if (expId==444)      filename = "/media/audrey/COURGE/SOFIA/ANALYSE/SOFIA3/data/202002_eng/main*.lmd";
  else if (expId==467)  filename = "/u/taniuchi/s467/s467_lustertmp/main0" + to_string(runnum) + "_*.lmd";
  //filename = "/media/audrey/COURGE/SOFIA/ANALYSE/SOFIA3/data/202002_s467/main0341_0001.lmd";
  else                 filename = "--stream=lxlanddaq01:9000";

  // Output file ----------------------------------------------------
  TString outputFileName = "rootfiles/SofSciRawPosRawTof" + to_string(runnum) + ".root";
 
  // UCESB configuration --------------------------------------------
  TString ntuple_options = "RAW";
  TString ucesb_dir = getenv("UCESB_DIR");
  //TString upexps_dir = "/u/land/fake_cvmfs/upexps";
  TString upexps_dir = ucesb_dir + "/../upexps";
  TString ucesb_path;
  if (expId == 444){
    ucesb_path = upexps_dir + "/202002_s444/202002_s444 --allow-errors --input-buffer=100Mi";
  }
  else if (expId == 467){
    ucesb_path = upexps_dir + "/202002_s467/202002_s467 --allow-errors --input-buffer=100Mi";
  }
  else{
    std::cout << "Experiment was not selected!" << std::endl;
    gApplication->Terminate();
  }
  ucesb_path.ReplaceAll("//", "/");
  
  // Definition of reader --------------------------------------------
  EXT_STR_h101 ucesb_struct;
  R3BUcesbSource* source = new R3BUcesbSource(filename, ntuple_options,ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
  source->AddReader(new R3BUnpackReader((EXT_STR_h101_unpack_t *)&ucesb_struct,offsetof(EXT_STR_h101, unpack)));
  source->AddReader(new R3BSofSciReader((EXT_STR_h101_SOFSCI_t *)&ucesb_struct.sci,offsetof(EXT_STR_h101, sci),NumSofSci));
    
  // --- ----------------- --- //
  // --- Create online run --- //
  // --- ----------------- --- //
  FairRunOnline* run = new FairRunOnline(source);
  run->SetRunId(1);
  run->SetOutputFile(outputFileName);

  // --- ----------------- --- //
  // --- Runtime data base --- //
  // --- ----------------- --- //
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
 
  // --- ------------------------------------------- --- //
  // --- Input parameters :                          --- //
  // ---  ascii file with the calibration parameters --- //
  // --- ------------------------------------------- --- //
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();//Ascii
  parIo1->open("../parameters/CalibParam.par","in");
  rtdb->setFirstInput(parIo1);
  rtdb->print();

  // --- ----------------- --- //
  // --- Add analysis task --- //
  // --- ----------------- --- //

  // === Mapped 2 Tcal for SofSci === //
  R3BSofSciMapped2Tcal* SofSciMap2Tcal = new R3BSofSciMapped2Tcal();
  run->AddTask(SofSciMap2Tcal);

  // === Tcal2RawPos for SofSci === //
  R3BSofSciTcal2RawPosPar* sci_poscalibrator = new R3BSofSciTcal2RawPosPar("R3BSofSciTcal2RawPosPar");
  sci_poscalibrator->SetNumDets(NumSofSci); 
  sci_poscalibrator->SetNumPmts(3); // number of Pmts (2) + reference signal (1)
  sci_poscalibrator->SetNumSignals();
  sci_poscalibrator->SetNumParsPerSignal(2);
  sci_poscalibrator->SetMinStatistics(10000);
  run->AddTask(sci_poscalibrator);

  // === Tcal2RawTof for SofSci === //
  R3BSofSciTcal2RawTofPar* sci_tofcalibrator = new R3BSofSciTcal2RawTofPar("R3BSofSciTcal2RawTofPar");
  sci_tofcalibrator->SetNumDets(NumSofSci); 
  sci_tofcalibrator->SetNumChannels(3); // number of Pmts (2) + reference signal (1)
  sci_tofcalibrator->SetDetIdS2(IdS2);
  sci_tofcalibrator->SetDetIdS8(IdS8);
  sci_tofcalibrator->SetDetIdCaveC(IdCaveC);
  sci_tofcalibrator->SetNumSignals();
  sci_tofcalibrator->SetNumParsPerSignal(2);
  sci_tofcalibrator->SetMinStatistics(1000);
  run->AddTask(sci_tofcalibrator);

  // --- ---------- --- //
  // --- Initialize --- //
  // --- ---------- --- //
  run->Init();
  FairLogger::GetLogger()->SetLogScreenLevel("INFO");
  //FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

  // --- ------------------------------------------- --- //
  // --- output parameters :                         --- //
  // ---  ascii file with the calibration parameters --- //
  // --- ------------------------------------------- --- //
  FairParAsciiFileIo* parOut = new FairParAsciiFileIo();
  TString outputFileNamePar = "parameters/out_sofsci" + to_string(runnum) + ".par";
  parOut->open(outputFileNamePar,"out");
  rtdb->setOutput(parOut);
  rtdb->print();

  // --- --- --- //
  // --- Run --- //
  // --- --- --- //
  run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);
  rtdb->saveOutput();

  // --- ------- --- //
  // --- Cleanup --- //
  // --- ------- --- //
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file is " << outputFileName << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s"
       << endl << endl;

}
