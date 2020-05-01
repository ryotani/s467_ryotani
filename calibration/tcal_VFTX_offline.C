typedef struct EXT_STR_h101_t
{
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_SOFSCI_onion_t sci;
  EXT_STR_h101_SOFTOFW_onion_t tofw;
}EXT_STR_h101;

void tcal_VFTX_offline(int runnum)
{
  TStopwatch timer;
  timer.Start();

  const Int_t nev = -1; /* number of events to read, -1 - until CTRL+C */

  // *********************************** //
  // PLEASE COMPLETE THE FOLLOWING LINES //
  // *********************************** //
  const Int_t expId = 467;               // select experiment: 444 or 467
  // *********************************** //

  // NumSofSci -----------------------------------------
  UShort_t NumSofSci;
  if (expId==444)     NumSofSci = 1; // s444: PRIMARY BEAM EXP, 1 SofSci at CAVE C ONLY
  else if(expId==467) NumSofSci = 4; // s467: SECONDARY BEAM EXP, 2 at S2, 1 at S8, 1 at CAVE C
  else                NumSofSci = 2; // default: SECONDARY BEAM EXP, 1 at S2, 1 at CAVE C

  // Create input -----------------------------------------
  TString filename;
  if (expId==444)      filename = "/media/audrey/COURGE/SOFIA/ANALYSE/SOFIA3/data/202002_eng/main*.lmd";
  else if (expId==467)   filename = "/u/taniuchi/s467/s467_lustertmp/main0" + to_string(runnum) + "_*.lmd";
  //filename = "/media/audrey/COURGE/SOFIA/ANALYSE/SOFIA3/data/202002_s467/main0341_0001.lmd";
  else                 filename = "--stream=lxlanddaq01:9000";
    
  // Output file ------------------------------------------
  TString outputFileName = "rootfiles/calibVftx" + to_string(runnum) + ".root";

  // UCESB configuration ----------------------------------
  TString ntuple_options = "RAW";
  TString ucesb_dir = getenv("UCESB_DIR");
  TString upexps_dir = ucesb_dir + "/../upexps/";
  //TString upexps_dir = "/u/land/fake_cvmfs/upexps";
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
 
  // Setup: Selection of detectors ------------------------
  Bool_t fSci = true;      // Start: Plastic scintillator for ToF
  Bool_t fTofW = true;     // ToF-Wall for time-of-flight of fragments behind GLAD
  
  // Create source using ucesb for input ------------------
  EXT_STR_h101 ucesb_struct;
  R3BUcesbSource* source = new R3BUcesbSource(filename, ntuple_options,ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
  source->SetMaxEvents(nev);
  
  // Definition of reader --------------------------------------------
  R3BUnpackReader* unpackreader;
  unpackreader = new R3BUnpackReader((EXT_STR_h101_unpack*)&ucesb_struct, offsetof(EXT_STR_h101, unpack));
  source->AddReader(unpackreader);

  if (fSci){
    R3BSofSciReader* unpacksci;
    unpacksci = new R3BSofSciReader((EXT_STR_h101_SOFSCI_t*)&ucesb_struct.sci, offsetof(EXT_STR_h101, sci),NumSofSci);
    source->AddReader(unpacksci);
  }

  if (fTofW){
    R3BSofTofWReader* unpacktofw;
    unpacktofw = new R3BSofTofWReader((EXT_STR_h101_SOFTOFW_t*)&ucesb_struct.tofw, offsetof(EXT_STR_h101, tofw));
    source->AddReader(unpacktofw);
  }

    
  /* Create online run ------------------------------------ */
  FairRunOnline* run = new FairRunOnline(source);
  run->SetRunId(1);
  run->SetOutputFile(outputFileName);

  R3BSofTcalContFact needToConstructTcalContFact;
  
  /* Runtime data base ------------------------------------ */
  FairRuntimeDb* rtdb = run->GetRuntimeDb();

  /* Calibrate SofSci ---------------------------------------- */
  if(fSci){
    R3BSofSciMapped2TcalPar* sciTcalibrator = new R3BSofSciMapped2TcalPar("R3BSofSciMapped2TcalPar");
    sciTcalibrator->SetNumSci(NumSofSci); 
    sciTcalibrator->SetNumChannels(3);
    sciTcalibrator->SetNumSignals(NumSofSci,3);
    //sciTcalibrator->SetMinStatistics(1000000);
    sciTcalibrator->SetMinStatistics(100000);
    run->AddTask(sciTcalibrator);
  }

  /* Calibrate time-of-flight wall  ---------------------------------------- */
  if(fTofW){
    R3BSofTofWMapped2TcalPar* tofwTcalibrator = new R3BSofTofWMapped2TcalPar("R3BSofTofWMapped2TcalPar");
    tofwTcalibrator->SetNumPaddles(28);
    tofwTcalibrator->SetNumPmtsPerPaddle(2);
    tofwTcalibrator->SetNumSignals(28,2);
    tofwTcalibrator->SetMinStatistics(50000);
    run->AddTask(tofwTcalibrator);
  }

  /* Initialize ------------------------------------------- */
  run->Init();
  FairLogger::GetLogger()->SetLogScreenLevel("INFO");
  //FairLogger::GetLogger()->SetLogScreenLevel("WARNING");
  /* ------------------------------------------------------ */
  
  // Ascii file with the Calibration Parameters
  FairParAsciiFileIo* parOut = new FairParAsciiFileIo();
  TString outputFileNamePar = "parameters/tcal_VFTX" + to_string(runnum) + ".par";
  parOut->open(outputFileNamePar,"out");
  rtdb->setOutput(parOut);
  
  /* Run -------------------------------------------------- */
  run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);
  rtdb->saveOutput();
  /* ------------------------------------------------------ */

  
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file is " << outputFileName << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s"
       << endl << endl;
}
