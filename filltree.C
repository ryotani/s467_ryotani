/*
 *  Macro to run the online for all the detectors simultaneously
 *
 *  One needs to set up the 2020 experiments: s444 or s467, the unpackers are:
 *
 *  at $UCESB_DIR/../upexps/202002_s444 and $UCESB_DIR/../upexps/202002_s467
 *
 *
 *  Author: Jose Luis <joseluis.rodriguez.sanchez@usc.es>
 *  @since Feb 20th, 2020
 *
 */

#include "profile.h"

typedef struct EXT_STR_h101_t
{
    EXT_STR_h101_unpack_t unpack;
    EXT_STR_h101_TPAT_t unpacktpat;
    EXT_STR_h101_SOFMWPC_onion_t mwpc;
    EXT_STR_h101_MUSIC_onion_t music;
    EXT_STR_h101_SOFSCI_onion_t sci;
    EXT_STR_h101_AMS_t ams;
    EXT_STR_h101_WRMASTER_t wrmaster;
    EXT_STR_h101_WRSOFIA_t wrsofia;
    EXT_STR_h101_WRS2_t wrs2;
    EXT_STR_h101_WRS8_t wrs8;
    EXT_STR_h101_CALIFA_t califa;
    //EXT_STR_h101_WRCALIFA_t wrcalifa;
    EXT_STR_h101_SOFTWIM_onion_t twim;
    EXT_STR_h101_SOFTOFW_onion_t tofw;
    EXT_STR_h101_SOFSCALERS_onion_t scalers;
    EXT_STR_h101_raw_nnp_tamex_t raw_nnp;
    EXT_STR_h101_WRNEULAND_t wrneuland;
    EXT_STR_h101_FRS_t frs;
} EXT_STR_h101;

void filltree(int runnum);

void filltree(){
  filltree(340);
}

void filltree(int runnum)
{
    TStopwatch timer;
    timer.Start();
    
    // NumSofSci, file names and paths -----------------------------
    Int_t sofiaWR, NumSofSci, IdS2, IdS8;
    TString ntuple_options = "RAW";//"RAW,time-stitch=1000" // For no stitched data
    TString filename, outputFilename, sofiacalfilename, vftxcalfilename, tofwhitfilename, musiccalfilename;
    Double_t brho28;
    
    if(runnum==0){
      //filename = "--stream=lxlanddaq01:9000";
      cerr<<"No online analysis available"<<endl;
      return 1;
    }else if (expId==444){ // not modified
      NumSofSci = 1; // s444: PRIMARY BEAM EXP, 1 SofSci at CAVE C ONLY
      IdS2 = 0;
      IdS8 = 0;
      sofiaWR = 0x500;
      
      filename = "/lustre/land/202002_s444/stitched/main0040_0001.lmd";
      outputFilename = "data_s444_online.root";
      
      upexps_dir = ucesb_dir + "/../upexps/";                      // for local computers
      // upexps_dir = "/u/land/fake_cvmfs/upexps";                 // for lxlandana computers
      // upexps_dir = "/u/land/lynx.landexp/202002_s444/upexps/";  // for lxg computers
      ucesb_path = upexps_dir + "/202002_s444/202002_s444 --allow-errors --input-buffer=100Mi";
      
      sofiacaldir = dir + "/sofia/macros/s444/parameters/";
    }
    else if (expId==467){
      NumSofSci = 4; // s467: SECONDARY BEAM EXP, 2 at S2, 1 at S8, 1 at CAVE C
      IdS2 = 2;
      IdS8 = 3;
      sofiaWR = 0xe00;

      if(!RunList.is_open()) std::cerr <<"No run summary found\n";
      int runnumcsv[500], targetpos[500], musicgain[500], junk[500];
      int FRSsetting[500]; // calib:0, ToFCalib:6-8, 40Ca:9, 39Ca:10, 38Ca:11,12, 50Ca:13, ToFWcalib:14
      string dummyline;
      char dumchar;
      double brhocsv[500];
      std::getline (RunList, dummyline);
      Int_t i=0;
      
      while(true){
	RunList>>runnumcsv[i]>>dumchar>>FRSsetting[i]>>dumchar>>brhocsv[i]>>dumchar>>targetpos[i]>>dumchar>>musicgain[i]>>dumchar>>junk[i];
	//std::cout<<runnumcsv[i]<<dumchar<<FRSsetting[i]<<dumchar<<brhocsv[i]<<dumchar<<targetpos[i]<<dumchar<<musicgain[i]<<dumchar<<junk[i]<<std::endl;
	if(runnumcsv[i] == runnum){
	  if(junk[i] == 0){
	    // if(targetpos[i]!=1424) return;
	    brho28 = brhocsv[i];
	    break;
	  } else {
	    std::cout << "Junk run" << std::endl;
	    return;
	  }
	}
	if(i > 400 || !RunList.good()){
	  std::cerr << "No info for run found" <<std::endl;
	  return;
	}
	i++;
      }
      
      filename = dir_rawfile;
      filename.Append(Form("main%04d_*.lmd", runnum));
      //
      if(FRSsetting[i] < 9){
	sofiacalfilename = sofiacaldir + "CalibParam_lowgain.par";
      } else if(musicgain[i] == 0){
	sofiacalfilename = sofiacaldir + "CalibParam_lowgain_FRS" + to_string(FRSsetting[i]) + ".par";
      } else {
	sofiacalfilename = sofiacaldir + "CalibParam_highgain_FRS" + to_string(FRSsetting[i]) + ".par";
      }
      vftxcalfilename = sofiacaldir + "tcal_VFTX.par";
      tofwhitfilename = sofiacaldir + "tofw_hit.par";
      musiccalfilename = sofiacaldir + "music_cal.par";
      //
      auto datime = new TDatime();
      TString str_datime = datime->AsString();
      string month = str_datime(4,3);
      outputFilename = dir_output;
      outputFilename.Append(Form("s467_filltree_Setting%i_%04d_%i%s.root", FRSsetting[i], runnum, datime->GetDay(), month.c_str()));

      std::cout << "LMD FILE: " << filename << std::endl;
      std::cout << "PARAM FILE (VFTX): " << vftxcalfilename << std::endl;
      std::cout << "PARAM FILE (TofW): " << tofwhitfilename << std::endl;
      std::cout << "PARAM FILE (MUSIC CAL): " << musiccalfilename << std::endl;
      std::cout << "PARAM FILE (OTHERS): " << sofiacalfilename << std::endl;
      std::cout << "OUTPUT FILE: " << outputFilename << std::endl;
      std::cout << "Brho28: " << brho28 << std::endl;
      
    }
    else{
      std::cout << "Experiment was not selected" << std::endl;
      gApplication->Terminate();
    }
    // Output file -----------------------------------------
    ucesb_path.ReplaceAll("//", "/");
    sofiacalfilename.ReplaceAll("//", "/");
    vftxcalfilename.ReplaceAll("//", "/");
    tofwhitfilename.ReplaceAll("//", "/");
    musiccalfilename.ReplaceAll("//", "/");
    califamapfilename.ReplaceAll("//", "/");
    califacalfilename.ReplaceAll("//", "/");

    // Online server configuration --------------------------
    Int_t refresh = 100; // Refresh rate for online histograms
    Int_t port = 10000 + runnum; // Port number for the online visualization, example lxgXXXX:8888


    // Create source using ucesb for input ------------------
    EXT_STR_h101 ucesb_struct;

    R3BUcesbSource* source =
        new R3BUcesbSource(filename, ntuple_options, ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
    source->SetMaxEvents(nev);
    source->SetSkipEvents(fSkip_tpat0);//skip tpat=0 events.

    // Definition of reader ---------------------------------
    source->AddReader(new R3BUnpackReader(&ucesb_struct.unpack,offsetof(EXT_STR_h101, unpack)));
    source->AddReader(new R3BTrloiiTpatReader(&ucesb_struct.unpacktpat,offsetof(EXT_STR_h101, unpacktpat)));
    
    R3BFrsReaderNov19* unpackfrs;
    R3BMusicReader* unpackmusic;
    R3BSofSciReader* unpacksci;
    R3BWhiterabbitS2Reader* unpackWRS2;
    R3BWhiterabbitS8Reader* unpackWRS8;
    R3BWhiterabbitMasterReader* unpackWRMaster;
    R3BSofWhiterabbitReader* unpackWRSofia;
    R3BAmsReader* unpackams;
    R3BCalifaFebexReader* unpackcalifa;
    //R3BWhiterabbitCalifaReader* unpackWRCalifa;
    R3BSofMwpcReader* unpackmwpc;
    R3BSofTwimReader* unpacktwim;
    R3BSofTofWReader* unpacktofw;
    R3BSofScalersReader* unpackscalers;
    R3BNeulandTamexReader* unpackneuland;
    R3BWhiterabbitNeulandReader* unpackWRNeuland;

    if (fFrs)
      unpackfrs= new R3BFrsReaderNov19((EXT_STR_h101_FRS*)&ucesb_struct.frs,
					     offsetof(EXT_STR_h101, frs));

    if (fMusic)
        unpackmusic = new R3BMusicReader((EXT_STR_h101_MUSIC_t*)&ucesb_struct.music, offsetof(EXT_STR_h101, music));
    
    if(fFrsSci) {
     unpackWRS2 = new R3BWhiterabbitS2Reader(
            (EXT_STR_h101_WRS2*)&ucesb_struct.wrs2, offsetof(EXT_STR_h101, wrs2), 0x200);
     unpackWRS8 = new R3BWhiterabbitS8Reader(
            (EXT_STR_h101_WRS8*)&ucesb_struct.wrs8, offsetof(EXT_STR_h101, wrs8), 0x800);
    }

    if (fSci)
    {
      unpacksci = new R3BSofSciReader((EXT_STR_h101_SOFSCI_t*)&ucesb_struct.sci, offsetof(EXT_STR_h101, sci),NumSofSci);
      unpackWRMaster = new R3BWhiterabbitMasterReader((EXT_STR_h101_WRMASTER*)&ucesb_struct.wrmaster, offsetof(EXT_STR_h101, wrmaster), 0x300);
      unpackWRSofia = new R3BSofWhiterabbitReader((EXT_STR_h101_WRSOFIA*)&ucesb_struct.wrsofia, offsetof(EXT_STR_h101, wrsofia), sofiaWR);
    }

    if (fCalifa)
    {
        unpackcalifa =
	  new R3BCalifaFebexReader((EXT_STR_h101_CALIFA*)&ucesb_struct.califa, offsetof(EXT_STR_h101, califa));
        //unpackWRCalifa = new R3BWhiterabbitCalifaReader(
        //    (EXT_STR_h101_WRCALIFA*)&ucesb_struct.wrcalifa, offsetof(EXT_STR_h101, wrcalifa), 0xa00, 0xb00);
    }

    if (fMwpc0 || fMwpc1 || fMwpc2 || fMwpc3)
        unpackmwpc = new R3BSofMwpcReader((EXT_STR_h101_SOFMWPC_t*)&ucesb_struct.mwpc, offsetof(EXT_STR_h101, mwpc));
    if (fTwim)
        unpacktwim = new R3BSofTwimReader((EXT_STR_h101_SOFTWIM_t*)&ucesb_struct.twim, offsetof(EXT_STR_h101, twim));
    if (fTofW)
        unpacktofw = new R3BSofTofWReader((EXT_STR_h101_SOFTOFW_t*)&ucesb_struct.tofw, offsetof(EXT_STR_h101, tofw));
    if (fScalers)
        unpackscalers =
            new R3BSofScalersReader((EXT_STR_h101_SOFSCALERS_t*)&ucesb_struct.scalers, offsetof(EXT_STR_h101, scalers));
    // Add readers ------------------------------------------

    if (fFrs)
    {
     unpackfrs->SetOnline(NOTstoremappeddata);
     source->AddReader(unpackfrs);
    }

    if (fMusic)
    {
        unpackmusic->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackmusic);
    }
    if (fSci)
    {
        unpacksci->SetOnline(NOTstoremappeddata);
        source->AddReader(unpacksci);
        unpackWRMaster->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackWRMaster);
        unpackWRSofia->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackWRSofia);
    }

    if(fFrsSci) {
        unpackWRS2->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackWRS2);
        unpackWRS8->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackWRS8);
    }
    if (fCalifa)
    {
        unpackcalifa->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackcalifa);
        //unpackWRCalifa->SetOnline(NOTstoremappeddata);
        //source->AddReader(unpackWRCalifa);
    }
    if (fMwpc0 || fMwpc1 || fMwpc2 || fMwpc3)
    {
        unpackmwpc->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackmwpc);
    }
    if (fTwim)
    {
        unpacktwim->SetOnline(NOTstoremappeddata);
        source->AddReader(unpacktwim);
    }
      if (fTofW)
    {
        unpacktofw->SetOnline(NOTstoremappeddata);
        source->AddReader(unpacktofw);
    }
    if (fScalers)
    {
        unpackscalers->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackscalers);
    }
    // Create online run ------------------------------------
    FairRunOnline* run = new FairRunOnline(source);
    run->SetRunId(fRunId);
    run->SetSink(new FairRootFileSink(outputFilename));
    run->ActivateHttpServer(refresh, port);

    // Runtime data base ------------------------------------
    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    
    FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii
    if (!fCalifa)
    {
        TList* parList1 = new TList();
	parList1->Add(new TObjString(musiccalfilename));
        parList1->Add(new TObjString(sofiacalfilename));
	parList1->Add(new TObjString(vftxcalfilename));
	parList1->Add(new TObjString(tofwhitfilename));
        parIo1->open(parList1, "in");
        rtdb->setFirstInput(parIo1);
        rtdb->print();
    }
    else
    {
        if (!fCal_level_califa)
        { // SOFIA and CALIFA mapping: Ascii files
            TList* parList1 = new TList();
	    parList1->Add(new TObjString(musiccalfilename));
            parList1->Add(new TObjString(sofiacalfilename));
	    parList1->Add(new TObjString(vftxcalfilename));
	    parList1->Add(new TObjString(tofwhitfilename));
            parList1->Add(new TObjString(califamapfilename));
            parIo1->open(parList1);
            rtdb->setFirstInput(parIo1);
            rtdb->print();
        }
        else
        { // SOFIA, CALIFA mapping and CALIFA calibration parameters
            TList* parList1 = new TList();
	    parList1->Add(new TObjString(musiccalfilename));
            parList1->Add(new TObjString(sofiacalfilename));
            parList1->Add(new TObjString(vftxcalfilename));
	    parList1->Add(new TObjString(tofwhitfilename));
            parIo1->open(parList1, "in");
	    rtdb->setFirstInput(parIo1);
            rtdb->print();
            Bool_t kParameterMerged = kFALSE;
            FairParRootFileIo* parIo2 = new FairParRootFileIo(kParameterMerged); // Root file
            TList* parList2 = new TList();
            parList2->Add(new TObjString(califacalfilename));
            parIo2->open(parList2);
            rtdb->setSecondInput(parIo2);
        }
    }
    
    // Add analysis task ------------------------------------
    // TPCs at S2
    if (fFrsTpcs)
    {
      R3BTpcMapped2Cal* TpcMap2Cal = new R3BTpcMapped2Cal();
      TpcMap2Cal->SetOnline(NOTstorecaldata);
      run->AddTask(TpcMap2Cal);
      R3BTpcCal2Hit* TpcCal2Hit = new R3BTpcCal2Hit();
      TpcCal2Hit->SetOnline(NOTstorehitdata);
      run->AddTask(TpcCal2Hit);
    }
    // MWPC0
    if (fMwpc0)
    {
        R3BSofMwpc0Mapped2Cal* MW0Map2Cal = new R3BSofMwpc0Mapped2Cal();
        MW0Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW0Map2Cal);

        R3BSofMwpc0Cal2Hit* MW0Cal2Hit = new R3BSofMwpc0Cal2Hit();
        MW0Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MW0Cal2Hit);
    }

    // MUSIC
    if (fMusic)
    {
        R3BMusicMapped2Cal* MusMap2Cal = new R3BMusicMapped2Cal();
        MusMap2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MusMap2Cal);

        R3BMusicCal2Hit* MusCal2Hit = new R3BMusicCal2Hit();
        MusCal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MusCal2Hit);
    }

    // SCI
    if (fSci)
    {
        // --- Mapped 2 Tcal for SofSci
        R3BSofSciMapped2Tcal* SofSciMap2Tcal = new R3BSofSciMapped2Tcal();
        SofSciMap2Tcal->SetOnline(NOTstorecaldata);
        run->AddTask(SofSciMap2Tcal);

        // --- Tcal 2 SingleTcal for SofSci
        R3BSofSciTcal2SingleTcal* SofSciTcal2STcal = new R3BSofSciTcal2SingleTcal();
        SofSciTcal2STcal->SetOnline(NOTstorecaldata);
        run->AddTask(SofSciTcal2STcal);
        // --- SingleTcal 2 Hit for SofSci
        R3BSofSciSingleTcal2Hit* SofSciSTcal2Hit = new R3BSofSciSingleTcal2Hit();
        SofSciSTcal2Hit->SetOnline(NOTstorehitdata);
        SofSciSTcal2Hit->SetCalParams(675.,-1922.);//ToF calibration at Cave-C
        run->AddTask(SofSciSTcal2Hit);
    }

    // CALIFA
    if (fCalifa && fCal_level_califa)
    {
        // R3BCalifaMapped2CrystalCal ---
        R3BCalifaMapped2CrystalCal* CalifaMap2Cal = new R3BCalifaMapped2CrystalCal();
        CalifaMap2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(CalifaMap2Cal);
        // R3BCalifaCrystalCal2Hit ---
        R3BCalifaCrystalCal2Hit* CalifaCal2Hit = new R3BCalifaCrystalCal2Hit();
        CalifaCal2Hit->SetCrystalThreshold(100.); // 100keV
        CalifaCal2Hit->SetDRThreshold(10000.);    // 10MeV
        CalifaCal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(CalifaCal2Hit);

        //R3BCalifaHitp2p* califaP2p = new R3BCalifaHitp2p();
        //califaP2p->SetOnline(NOTstorecaldata);
        //run->AddTask(califaP2p);
    }

    // MWPC1
    if (fMwpc1)
    {
        R3BSofMwpc1Mapped2Cal* MW1Map2Cal = new R3BSofMwpc1Mapped2Cal();
        MW1Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW1Map2Cal);
    }

    // TWIM
    if (fTwim)
    {
        R3BSofTwimMapped2Cal* TwimMap2Cal = new R3BSofTwimMapped2Cal();
        TwimMap2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(TwimMap2Cal);

        R3BSofTwimCal2Hit* TwimCal2Hit = new R3BSofTwimCal2Hit();
        TwimCal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(TwimCal2Hit);
    }

    // MWPC2
    if (fMwpc2)
    {
        R3BSofMwpc2Mapped2Cal* MW2Map2Cal = new R3BSofMwpc2Mapped2Cal();
        MW2Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW2Map2Cal);

        R3BSofMwpc2Cal2Hit* MW2Cal2Hit = new R3BSofMwpc2Cal2Hit();
        MW2Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MW2Cal2Hit);
    }

    if (fMwpc1 && fMwpc2)
    {
        R3BSofMwpc1Cal2Hit* MW1Cal2Hit = new R3BSofMwpc1Cal2Hit();
        MW1Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MW1Cal2Hit);
    }

    // MWPC3
    if (fMwpc3)
    {
        R3BSofMwpc3Mapped2Cal* MW3Map2Cal = new R3BSofMwpc3Mapped2Cal();
        MW3Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW3Map2Cal);

        R3BSofMwpc3Cal2Hit* MW3Cal2Hit = new R3BSofMwpc3Cal2Hit();
        MW3Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MW3Cal2Hit);
    }

    // ToF-Wall
    if (fTofW)
    {
        // --- Mapped 2 Tcal for SofTofW
        R3BSofTofWMapped2Tcal* SofTofWMap2Tcal = new R3BSofTofWMapped2Tcal();
        SofTofWMap2Tcal->SetOnline(NOTstorecaldata);
        run->AddTask(SofTofWMap2Tcal);

        // --- Tcal 2 SingleTcal for SofTofW
        R3BSofTofWTcal2SingleTcal* SofTofWTcal2STcal = new R3BSofTofWTcal2SingleTcal();
        SofTofWTcal2STcal->SetOnline(NOTstorecaldata);
        run->AddTask(SofTofWTcal2STcal);

        // --- Tcal 2 Hit for SofTofW :
        R3BSofTofWSingleTCal2Hit* SofTofWTcal2Hit = new R3BSofTofWSingleTCal2Hit();
	SofTofWTcal2Hit->SetTofLISE(43.);
        SofTofWTcal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(SofTofWTcal2Hit);
    }
    // Add sofana task ------------------------------------
    if (fSci&&fMusic)
    {
      R3BSofFrsAnalysis* frsana = new R3BSofFrsAnalysis();
      frsana->SetOnline(NOTstorehitdata);
      run -> AddTask(frsana);
    }

    if (fSci&&fMusic&&fTwim&&fMwpc0&&fMwpc1&&fMwpc2&&fMwpc3&&fTofW){
      R3BSofFragmentAnalysis* fragmentana = new R3BSofFragmentAnalysis();
      fragmentana->SetOnline(NOTstorehitdata);
      fragmentana->SetTofWPos(560.);
      run -> AddTask(fragmentana);
    }

    // Add online task ------------------------------------
    if (fFrsTpcs)
    {
       FrsTpcOnlineSpectra* tpconline= new FrsTpcOnlineSpectra();
       run->AddTask(tpconline);
    }

    if (fScalers)
    {
        R3BSofScalersOnlineSpectra* scalersonline = new R3BSofScalersOnlineSpectra();
        run->AddTask(scalersonline);
    }
    /*
    if (fSci&&fMusic&&fTwim){
      R3BSofFrsFillTree* frsfilltree = new R3BSofFrsFillTree();
	/ *frsfilltree->SetNbDetectors(NumSofSci);
	frsfilltree->SetNbChannels(3);
	frsfilltree->SetIdS2(IdS2);
	frsfilltree->SetIdS8(IdS8);* /
        run->AddTask(frsfilltree);
    }
* /
    if (fSci&&fMusic&&fTwim&&fMwpc0&&fMwpc1&&fMwpc2&&fMwpc3&&fTofW){
      R3BSofFrsFragmentTree* frsfragmenttree = new R3BSofFrsFragmentTree();
      frsfragmenttree->SetIdS2(IdS2);
      frsfragmenttree->SetIdS8(IdS8);
      run->AddTask(frsfragmenttree);
    }
    */
    // Initialize -------------------------------------------
    run->Init();
    FairLogger::GetLogger()->SetLogScreenLevel("INFO");
    //FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");

    // Run --------------------------------------------------
    run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);

    // Finish -----------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    std::cout << std::endl << std::endl;
    std::cout << "Macro finished succesfully." << std::endl;
    std::cout << "Output file is " << outputFilename << std::endl;
    std::cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << std::endl << std::endl;
    gApplication->Terminate();
}

