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
    EXT_STR_h101_WRMASTER_t wrmaster;
    EXT_STR_h101_WRSOFIA_t wrsofia;
    EXT_STR_h101_WRS2_t wrs2;
    EXT_STR_h101_WRS8_t wrs8;
    EXT_STR_h101_CALIFA_t califa;
    EXT_STR_h101_SOFTWIM_onion_t twim;
    EXT_STR_h101_SOFTOFW_onion_t tofw;
    EXT_STR_h101_SOFSCALERS_onion_t scalers;
    EXT_STR_h101_raw_nnp_tamex_onion_t raw_nnp; // EXT_STR_h101_raw_nnp_tamex_t raw_nnp;
    EXT_STR_h101_WRNEULAND_t wrneuland;
} EXT_STR_h101;

void filltree(int runnum);

void filltree() { filltree(338); }

void filltree(int runnum)
{
    TStopwatch timer;
    timer.Start();

    // NumSofSci, file names and paths -----------------------------
    Int_t sofiaWR, NumSofSci, IdS2, IdS8;
    TString outputFilename, frs_paramfile, fragment_paramfile, common_paramfile, musiccalfilename;
    Double_t brho28;
    //-------------------------------------------------------------------------
    // added by KB for NeuLAND:
    const Int_t nBarsPerPlane = 50; // number of scintillator bars per plane
    const Int_t nPlanes = 16;       // number of planes (for TCAL calibration)
    const double distanceToTarget = 1520.;
    const Int_t trigger = -1; // 1 - onspill, 2 - offspill. -1 - all
    double timeoffset = 0.;
    if (262 < runnum & runnum < 269)
    {
        timeoffset = 4202. + 6.2;
    }
    else if (276 < runnum & runnum < 304)
    {
        timeoffset = 7200. + 78.0; // check spectra
    }
    else if (303 < runnum & runnum < 359)
    {
        timeoffset = 7200. + 58.0; // check spectra
    }
    else if (367 < runnum & runnum < 381)
    {
        timeoffset = 10610. + 0.3; // check spectra
    }

    // associate parameter file:
#define LENGTH(x) (sizeof x / sizeof *x)

    Int_t calib_file[400] = { 0 };
    Int_t two[] = { 328, 330, 357 };
    Int_t three[] = { 281, 284, 287, 294, 297, 300, 308, 311, 318, 321, 332, 335, 342, 345, 348, 351, 354 };
    Int_t four[] = { 277, 290, 314, 324, 338, 368, 372 };
    Int_t five[] = { 303, 376 };
    Int_t six[] = { 263 };

    for (Int_t j = 0; j < LENGTH(two); j++)
      {
	for (Int_t iruns = 0; iruns < 2; iruns++)
	  calib_file[two[j] + iruns] = two[j];
      }
    for (Int_t j = 0; j < LENGTH(three); j++)
      {
	for (Int_t iruns = 0; iruns < 3; iruns++)
	  calib_file[three[j] + iruns] = three[j];
      }
    for (Int_t j = 0; j < LENGTH(four); j++)
      {
	for (Int_t iruns = 0; iruns < 4; iruns++)
	  calib_file[four[j] + iruns] = four[j];
      }
    for (Int_t j = 0; j < LENGTH(five); j++)
      {
	for (Int_t iruns = 0; iruns < 5; iruns++)
	  calib_file[five[j] + iruns] = five[j];
      }
    for (Int_t j = 0; j < LENGTH(six); j++)
      {
	for (Int_t iruns = 0; iruns < 6; iruns++)
	  calib_file[six[j] + iruns] = six[j];
      }
    const TString syncParFileName = TString::Format("/u/boretzky/s467/params_new_sync_%04d.root", calib_file[runnum]);
    cout << "NeuLAND calibration parameters from:   " << syncParFileName << endl;

    //-------------------------------------------------------------------------

    if (runnum == 0)
    {
        // filename = "--stream=lxlanddaq01:9000";
        cerr << "No online analysis available" << endl;
        return 1;
    }
    else if (expId == 444)
    {                  // not modified
        NumSofSci = 1; // s444: PRIMARY BEAM EXP, 1 SofSci at CAVE C ONLY
        IdS2 = 0;
        IdS8 = 0;
        sofiaWR = 0x500;

        filename = "/lustre/land/202002_s444/stitched/main0040_0001.lmd";
        outputFilename = "data_s444_online.root";

        upexps_dir = ucesb_dir + "/../upexps/"; // for local computers
        // upexps_dir = "/u/land/fake_cvmfs/upexps";                 // for lxlandana computers
        // upexps_dir = "/u/land/lynx.landexp/202002_s444/upexps/";  // for lxg computers
        ucesb_path = upexps_dir + "/202002_s444/202002_s444 --allow-errors --input-buffer=100Mi";

        sofiacaldir = dir + "/sofia/macros/s444/parameters/";
    }
    else if (expId == 467)
    {
        NumSofSci = 4; // s467: SECONDARY BEAM EXP, 2 at S2, 1 at S8, 1 at CAVE C
        IdS2 = 2;
        IdS8 = 3;
        sofiaWR = 0xe00;

        if (!RunList.is_open())
            std::cerr << "No run summary found\n";
        int runnumcsv[500], targetpos[500], musicgain[500], junk[500];
        int FRSsetting[500]; // calib:0, ToFCalib:6-8, 40Ca:9, 39Ca:10, 38Ca:11,12, 50Ca:13, ToFWcalib:14
        string dummyline;
        char dumchar;
        double brhocsv[500];
        std::getline(RunList, dummyline);
        Int_t i = 0;

        while (true)
        {
            RunList >> runnumcsv[i] >> dumchar >> FRSsetting[i] >> dumchar >> brhocsv[i] >> dumchar >> targetpos[i] >>
                dumchar >> musicgain[i] >> dumchar >> junk[i];
            // std::cout<<runnumcsv[i]<<dumchar<<FRSsetting[i]<<dumchar<<brhocsv[i]<<dumchar<<targetpos[i]<<dumchar<<musicgain[i]<<dumchar<<junk[i]<<std::endl;
            if (runnumcsv[i] == runnum)
            {
                if (junk[i] == 0)
                {
                    // if(targetpos[i]!=1424) return;
                    brho28 = brhocsv[i];
                    break;
                }
                else
                {
                    std::cout << "Junk run" << std::endl;
                    return;
                }
            }
            if (i > 400 || !RunList.good())
            {
                std::cerr << "No info for run found" << std::endl;
                return;
            }
            i++;
        }
        filename.ReplaceAll("RUNNUM",Form("%04d", runnum));
        //
        if (FRSsetting[i] < 9)
        { // Default
            frs_paramfile = sofiacaldir + "FRS13.par";
            musiccalfilename = sofiacaldir + "music_highgain.par";
        }
        else
        {
  	    frs_paramfile = sofiacaldir + "FRS" + TString::Itoa(FRSsetting[i],10) + ".par";
            if (FRSsetting[i] < 11)
            {
                fragment_paramfile = sofiacaldir + "FRS13_empty.par";
            }
            else
            {
                switch (targetpos[i])
                {
                    case 539:
                        fragment_paramfile = sofiacaldir + "FRS" + TString::Itoa(FRSsetting[i],10) + "_ch2.par";
                        break;
                    case 362:
                        fragment_paramfile = sofiacaldir + "FRS" + TString::Itoa(FRSsetting[i],10) + "_carbon.par";
                        break;
                    default:
                        fragment_paramfile = sofiacaldir + "FRS" + TString::Itoa(FRSsetting[i],10) + "_empty.par";
                        break;
                }
            }
            if (musicgain[i] == 0)
            {
                musiccalfilename = sofiacaldir + "music_lowgain.par";
            }
            else if (musicgain[i] == 1)
            {
                musiccalfilename = sofiacaldir + "music_highgain.par";
            }
	    else if (musicgain[i] == 122)
	    {
	      musiccalfilename = sofiacaldir + "music_highgain122.par";
	    }
	    else
	    {
	      cerr<< "Nno music parameter file" <<endl;
	    }

        }
        common_paramfile = sofiacaldir + "common.par";
        //
        auto datime = new TDatime();
        TString str_datime = datime->AsString();
        string month = str_datime(4, 3);
        outputFilename = dir_output;
        outputFilename.Append(
            Form("s467_filltree_Setting%i_%04d_%i%s.root", FRSsetting[i], runnum, datime->GetDay(), month.c_str()));

        std::cout << "LMD FILE: " << filename << std::endl;
        std::cout << "PARAM FILE (Common): " << common_paramfile << std::endl;
        std::cout << "PARAM FILE (MUSIC CAL): " << musiccalfilename << std::endl;
        std::cout << "PARAM FILE (FRS): " << frs_paramfile << std::endl;
        std::cout << "PARAM FILE (FRAGMENT): " << fragment_paramfile << std::endl;
        std::cout << "PARAM FILE (CALIFA): " << califacalfilename << std::endl;
        std::cout << "OUTPUT FILE: " << outputFilename << std::endl;
        std::cout << "Brho28: " << brho28 << std::endl;
    }
    else
    {
        std::cout << "Experiment was not selected" << std::endl;
        gApplication->Terminate();
    }
    // Output file -----------------------------------------
    ucesb_path.ReplaceAll("//", "/");
    frs_paramfile.ReplaceAll("//", "/");
    fragment_paramfile.ReplaceAll("//", "/");
    common_paramfile.ReplaceAll("//", "/");
    musiccalfilename.ReplaceAll("//", "/");
    califamapfilename.ReplaceAll("//", "/");
    califacalfilename.ReplaceAll("//", "/");

    // Online server configuration --------------------------
    Int_t refresh = 100;         // Refresh rate for online histograms
    Int_t port = 10000 + runnum; // Port number for the online visualization, example lxgXXXX:8888

    // Create source using ucesb for input ------------------
    EXT_STR_h101 ucesb_struct;

    R3BUcesbSource* source =
        new R3BUcesbSource(filename, ntuple_options, ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
    source->SetMaxEvents(nev);

    // Definition of reader ---------------------------------
    source->AddReader(new R3BUnpackReader(&ucesb_struct.unpack, offsetof(EXT_STR_h101, unpack)));
    auto TrloiiReader = new R3BTrloiiTpatReader(&ucesb_struct.unpacktpat, offsetof(EXT_STR_h101, unpacktpat));
    TrloiiReader->SetTpatRange(min_tpat, max_tpat);
    source->AddReader(TrloiiReader);

    R3BMusicReader* unpackmusic;
    R3BSofSciReader* unpacksci;
    R3BWhiterabbitS2Reader* unpackWRS2;
    R3BWhiterabbitS8Reader* unpackWRS8;
    R3BWhiterabbitMasterReader* unpackWRMaster;
    R3BSofWhiterabbitReader* unpackWRSofia;
    R3BAmsReader* unpackams;
    R3BCalifaFebexReader* unpackcalifa;
    R3BMwpcReader* unpackmwpc;
    R3BTwimReader* unpacktwim;
    R3BSofTofWReader* unpacktofw;
    R3BSofScalersReader* unpackscalers;
    R3BNeulandTamexReader* unpackneuland;
    R3BWhiterabbitNeulandReader* unpackWRNeuland;

    if (fMusic)
        unpackmusic = new R3BMusicReader((EXT_STR_h101_MUSIC_t*)&ucesb_struct.music, offsetof(EXT_STR_h101, music));

    if (fFrsSci)
    {
        unpackWRS2 =
            new R3BWhiterabbitS2Reader((EXT_STR_h101_WRS2*)&ucesb_struct.wrs2, offsetof(EXT_STR_h101, wrs2), 0x200);
        unpackWRS8 =
            new R3BWhiterabbitS8Reader((EXT_STR_h101_WRS8*)&ucesb_struct.wrs8, offsetof(EXT_STR_h101, wrs8), 0x800);
    }

    if (fSci)
    {
        unpacksci =
            new R3BSofSciReader((EXT_STR_h101_SOFSCI_t*)&ucesb_struct.sci, offsetof(EXT_STR_h101, sci), NumSofSci);
        unpackWRMaster = new R3BWhiterabbitMasterReader(
            (EXT_STR_h101_WRMASTER*)&ucesb_struct.wrmaster, offsetof(EXT_STR_h101, wrmaster), 0x300);
        unpackWRSofia = new R3BSofWhiterabbitReader((EXT_STR_h101_WRSOFIA*)&ucesb_struct.wrsofia,
                                                    offsetof(EXT_STR_h101, wrsofia),
                                                    sofiaWR,
                                                    0xf00); // No sofiaWR2 available but tentatively assign as 0xf00.
    }

    if (fCalifa)
    {
        unpackcalifa =
            new R3BCalifaFebexReader((EXT_STR_h101_CALIFA*)&ucesb_struct.califa, offsetof(EXT_STR_h101, califa));
    }

    if (fMwpc0 || fMwpc1 || fMwpc2 || fMwpc3)
        unpackmwpc = new R3BMwpcReader((EXT_STR_h101_SOFMWPC_t*)&ucesb_struct.mwpc, offsetof(EXT_STR_h101, mwpc));
    if (fTwim)
        unpacktwim = new R3BTwimReader((EXT_STR_h101_SOFTWIM_t*)&ucesb_struct.twim, offsetof(EXT_STR_h101, twim));
    if (fTofW)
        unpacktofw = new R3BSofTofWReader((EXT_STR_h101_SOFTOFW_t*)&ucesb_struct.tofw, offsetof(EXT_STR_h101, tofw));
    if (fNeuland)
        unpackneuland = new R3BNeulandTamexReader(&ucesb_struct.raw_nnp, offsetof(EXT_STR_h101, raw_nnp));
    if (fScalers)
        unpackscalers =
            new R3BSofScalersReader((EXT_STR_h101_SOFSCALERS_t*)&ucesb_struct.scalers, offsetof(EXT_STR_h101, scalers));
    // Add readers ------------------------------------------
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

    if (fFrsSci)
    {
        unpackWRS2->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackWRS2);
        unpackWRS8->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackWRS8);
    }
    if (fCalifa)
    {
        unpackcalifa->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackcalifa);
    }
    if (fMwpc0 || fMwpc1 || fMwpc2 || fMwpc3)
    {
        unpackmwpc->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackmwpc);
    }
    if (fTwim)
    {
        unpacktwim->SetOnline(NOTstoremappeddata);
        unpacktwim->SetNumTref(2);
        unpacktwim->SetNumTtrig(2);
        source->AddReader(unpacktwim);
    }
    if (fTofW)
    {
        unpacktofw->SetOnline(NOTstoremappeddata);
        source->AddReader(unpacktofw);
    }
    if (fNeuland)
    {
        unpackneuland->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackneuland);
    }
    if (fScalers)
    {
        unpackscalers->SetOnline(NOTstoremappeddata);
        source->AddReader(unpackscalers);
    }
    // Create online run ------------------------------------
    FairRunOnline* run = new FairRunOnline(source);
    R3BEventHeader* EvntHeader = new R3BEventHeader();
    EvntHeader->SetExpId(expId);
    run->SetEventHeader(EvntHeader);
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
        parList1->Add(new TObjString(frs_paramfile));
        parList1->Add(new TObjString(fragment_paramfile));
        parList1->Add(new TObjString(common_paramfile));
        parIo1->open(parList1, "in");
        rtdb->setFirstInput(parIo1);
    }
    else
    {
        if (!fCal_level_califa)
        { // SOFIA and CALIFA mapping: Ascii files
            TList* parList1 = new TList();
            parList1->Add(new TObjString(musiccalfilename));
            parList1->Add(new TObjString(frs_paramfile));
            parList1->Add(new TObjString(fragment_paramfile));
            parList1->Add(new TObjString(common_paramfile));
            parList1->Add(new TObjString(califamapfilename));
            parIo1->open(parList1);
            rtdb->setFirstInput(parIo1);
        }
        else
        { // SOFIA, CALIFA mapping and CALIFA calibration parameters
            TList* parList1 = new TList();
            parList1->Add(new TObjString(musiccalfilename));
            parList1->Add(new TObjString(frs_paramfile));
            parList1->Add(new TObjString(fragment_paramfile));
            parList1->Add(new TObjString(common_paramfile));
            parIo1->open(parList1, "in");
            rtdb->setFirstInput(parIo1);
            Bool_t kParameterMerged = kFALSE;
            FairParRootFileIo* parIo2 = new FairParRootFileIo(kParameterMerged); // Root file
            //TList* parList2 = new TList();
            //parList2->Add(new TObjString(califacalfilename));
	    //parIo2->open(parList2); // See Class reference. If TList is used, allParams will be created.
	    //
	    parIo2->open(califacalfilename);
            rtdb->setSecondInput(parIo2);
        }
    }
    if (fNeuland)
    {
        // added by KB: root file for NeuLAND - test with one file
        auto parIO = new FairParRootFileIo(false);
        parIO->open(syncParFileName, "in");
        rtdb->setSecondInput(parIO);
        rtdb->addRun(999);
        rtdb->getContainer("LandTCalPar");
        rtdb->setInputVersion(999, (char*)"LandTCalPar", 1, 1);
        rtdb->getContainer("NeulandHitPar");
        rtdb->setInputVersion(999, (char*)"NeulandHitPar", 1, 1);
        cout << "did neuland stuff for rtdb" << endl;
    }
    rtdb->print();
    cout<<"rtdb print end."<<endl;

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
        R3BMwpc0Mapped2Cal* MW0Map2Cal = new R3BMwpc0Mapped2Cal();
        MW0Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW0Map2Cal);

        R3BMwpc0Cal2Hit* MW0Cal2Hit = new R3BMwpc0Cal2Hit();
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
        SofSciSTcal2Hit->SetCalParams(675., -1922.); // ToF calibration at Cave-C
        run->AddTask(SofSciSTcal2Hit);
        //
	/* This is for Offline analysis, R3BFileSource required
        // for CALIFA and Neuland
        R3BEventHeaderPropagator* RunIdTask = new R3BEventHeaderPropagator();
        run->AddTask(RunIdTask);
	*/
        auto sofstart = new R3BSofiaProvideTStart();
        run->AddTask(sofstart);
    }

    // CALIFA
    if (fCalifa && fCal_level_califa)
    {
        // R3BCalifaMapped2CrystalCal ---
        R3BCalifaMapped2CrystalCal* CalifaMap2Cal = new R3BCalifaMapped2CrystalCal();
        CalifaMap2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(CalifaMap2Cal);
        // R3BCalifaCrystalCal2Hit ---
        {
            R3BCalifaCrystalCal2Cluster* CalifaCal2Cluster = new R3BCalifaCrystalCal2Cluster();
            //CalifaCal2Cluster->SetCrystalThreshold(100.); // 100keV
            //CalifaCal2Cluster->SetDRThreshold(10000.);    // 10MeV
            //CalifaCal2Cluster->SetProtonClusterThreshold(10000.);    // 10MeV
            CalifaCal2Cluster->SetOnline(NOTstorehitdata);
            //CalifaCal2Cluster->SetRoundWindowAlg(0.20); // radian
            //CalifaCal2Cluster->SetRoundWindow(0.20);
            CalifaCal2Cluster->SelectGeometryVersion(2020);
            CalifaCal2Cluster->SetRandomization(false);
            //CalifaCal2Cluster->SetTCAName("CalifaHitData_Round020");
            run->AddTask(CalifaCal2Cluster);
        }
	/*
        {
            R3BCalifaCrystalCal2Cluster* CalifaCal2Cluster = new R3BCalifaCrystalCal2Cluster();
            CalifaCal2Cluster->SetCrystalThreshold(100.); // 100keV
            CalifaCal2Cluster->SetDRThreshold(10000.);    // 10MeV
            CalifaCal2Cluster->SetOnline(NOTstorehitdata);
            CalifaCal2Cluster->SetSquareWindowAlg(0.20, 0.20);
            CalifaCal2Cluster->SelectGeometryVersion(2020);
            CalifaCal2Cluster->SetRandomization(false);
            CalifaCal2Cluster->SetTCAName("CalifaHitDataSquare");
            run->AddTask(CalifaCal2Cluster);
        }

        {
            R3BCalifaCrystalCal2Cluster* CalifaCal2Cluster = new R3BCalifaCrystalCal2Cluster();
            CalifaCal2Cluster->SetCrystalThreshold(100.); // 100keV
            CalifaCal2Cluster->SetDRThreshold(10000.);    // 10MeV
            CalifaCal2Cluster->SetOnline(NOTstorehitdata);
            CalifaCal2Cluster->SetConeAlg(0.20);
            CalifaCal2Cluster->SelectGeometryVersion(2020);
            CalifaCal2Cluster->SetRandomization(false);
            CalifaCal2Cluster->SetTCAName("CalifaHitData_Cone020");
            run->AddTask(CalifaCal2Cluster);
        }
	*/
    }

    // MWPC1
    if (fMwpc1)
    {
        R3BMwpc1Mapped2Cal* MW1Map2Cal = new R3BMwpc1Mapped2Cal();
        MW1Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW1Map2Cal);

        R3BMwpc1Cal2Hit* MW1Cal2Hit = new R3BMwpc1Cal2Hit();
        MW1Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MW1Cal2Hit);
    }

    // TWIM
    if (fTwim)
    {
        R3BTwimMapped2Cal* TwimMap2Cal = new R3BTwimMapped2Cal();
        TwimMap2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(TwimMap2Cal);

        R3BTwimCal2Hit* TwimCal2Hit = new R3BTwimCal2Hit();
        TwimCal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(TwimCal2Hit);
    }

    // MWPC2
    if (fMwpc2)
    {
        R3BMwpc2Mapped2Cal* MW2Map2Cal = new R3BMwpc2Mapped2Cal();
        MW2Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW2Map2Cal);

        R3BMwpc2Cal2Hit* MW2Cal2Hit = new R3BMwpc2Cal2Hit();
        MW2Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(MW2Cal2Hit);
    }

    // MWPC3
    if (fMwpc3)
    {
        R3BMwpc3Mapped2Cal* MW3Map2Cal = new R3BMwpc3Mapped2Cal();
        MW3Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(MW3Map2Cal);

        R3BMwpc3Cal2Hit* MW3Cal2Hit = new R3BMwpc3Cal2Hit();
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
        SofTofWTcal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(SofTofWTcal2Hit);
    }

    // NeuLAND
    if (fNeuland)
    {
        auto tcal = new R3BNeulandMapped2Cal();
        tcal->SetTrigger(trigger);
        tcal->SetNofModules(nPlanes, nBarsPerPlane);
        tcal->SetNhitmin(1);
        tcal->EnableWalk(true);
        run->AddTask(tcal);

        auto nlhit = new R3BNeulandCal2Hit();
        nlhit->SetDistanceToTarget(distanceToTarget);
        nlhit->SetGlobalTimeOffset(timeoffset);
        nlhit->SetEnergyCutoff(0.0);
        run->AddTask(nlhit);
    }

    // Add sofana task ------------------------------------
    if (fSci && fMusic)
    {
        R3BSofFrsAnalysis* frsana = new R3BSofFrsAnalysis();
        frsana->SetOnline(NOTstorehitdata);
        run->AddTask(frsana);
    }

    if (fSci && fMusic && fTwim && fMwpc0 && fMwpc1 && fMwpc2 && fMwpc3 && fTofW)
    {
        R3BSofFragmentAnalysis* fragmentana = new R3BSofFragmentAnalysis();
        fragmentana->SetOnline(NOTstorehitdata);
        fragmentana->SetTofWPos(560.);
        run->AddTask(fragmentana);
    }

    // Add online task ------------------------------------
    if (fScalers)
    {
        R3BSofScalersOnlineSpectra* scalersonline = new R3BSofScalersOnlineSpectra();
        run->AddTask(scalersonline);
    }
    /*
      if (fSci&&fMusic&&fTwim&&fMwpc0&&fMwpc1&&fMwpc2&&fMwpc3&&fTofW&&(!fCalifa)){
      R3BSofFrsFragmentTree* frsfragmenttree = new R3BSofFrsFragmentTree();
      frsfragmenttree->SetIdS2(IdS2);
      frsfragmenttree->SetIdS8(IdS8);
      run->AddTask(frsfragmenttree);
      }
    */
    // Initialize -------------------------------------------
    run->Init();
    FairLogger::GetLogger()->SetLogScreenLevel("info");
    // FairLogger::GetLogger()->SetLogScreenLevel("debug");

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
