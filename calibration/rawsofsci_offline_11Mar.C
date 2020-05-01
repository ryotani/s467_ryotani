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

  // --- ----------------------------------- --- //
  // --- Create source using ucesb for input --- //
  // --- ----------------------------------- --- //
  //TString filename = " --stream=lxlanddaq01:9000";
  //TString filename = "/lustre/land/202002_s467/lustre/r3b/202002_s467/main0300_0001.lmd";
  TString filename = "/u/taniuchi/s467/s467_lustertmp/main0" + to_string(runnum) + "_*.lmd";
  cout<<"Opening "<< filename<<endl;
  
  // Output file ----------------------------------------------------
  TString outputFileName = "rootfiles/SofSciRawPosRawTof" + to_string(runnum) + ".root";
  // UCESB configuration --------------------------------------------
  TString ntuple_options = "RAW";
  TString ucesb_dir = getenv("UCESB_DIR");
  TString upexps_dir = ucesb_dir + "/../upexps/"; //Copied from fake_cvmf and re-compiled
  //TString upexps_dir = "/u/land/fake_cvmfs/upexps";
  TString ucesb_path;
  ucesb_path = upexps_dir + "/202002_s467/202002_s467 --allow-errors --input-buffer=100Mi";
  ucesb_path.ReplaceAll("//","/");
  
  // Definition of reader --------------------------------------------
  EXT_STR_h101 ucesb_struct;
  R3BUcesbSource* source = new R3BUcesbSource(filename, ntuple_options,ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
  source->AddReader(new R3BUnpackReader((EXT_STR_h101_unpack_t *)&ucesb_struct,offsetof(EXT_STR_h101, unpack)));
  source->AddReader(new R3BSofSciReader((EXT_STR_h101_SOFSCI_t *)&ucesb_struct.sci,offsetof(EXT_STR_h101, sci)));
    
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
  sci_poscalibrator->SetNumDetectors(4); 
  sci_poscalibrator->SetNumChannels(3);
  sci_poscalibrator->SetNumSignals();
  sci_poscalibrator->SetNumParsPerSignal(2);
  sci_poscalibrator->SetMinStatistics(10000);
  run->AddTask(sci_poscalibrator);

  // === Tcal2RawTof for SofSci === //
  R3BSofSciTcal2RawTofPar* sci_tofcalibrator = new R3BSofSciTcal2RawTofPar("R3BSofSciTcal2RawTofPar");
  sci_tofcalibrator->SetFirstStaSci(2);
  //sci_tofcalibrator->SetFirstStoSci(4);
  sci_tofcalibrator->SetNumSignals();
  sci_tofcalibrator->SetNumParsPerSignal(2);
  sci_tofcalibrator->SetMinStatistics(10000);
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
  TString outputFileNamePar = "outfile/out_sofsci" + to_string(runnum) + ".par";
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
