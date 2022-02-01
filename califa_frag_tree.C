/*
*   This macro uses the R3BMusic infront the target and the TwimMusic behind the target to select charges.
*   In the second step it looks for this charge selection into the CALIFA data, searching for good p2p candidates.
*   A good candidate is exactly two "high energy" hits in CALIFA.
*/

#include "califa_frag_tree.h"



//======== MAIN ===================================================================================================
void califa_frag_tree()
{
  loadtrees();
  
  //Running the analysis Looop (where the magic happens)
  //analysisLoop(1e7);
  analysisLoop();

  //Writing the produced histos to disk. This is up to you to change to your liking
  writeHistos();

  //Closing the used file
  //eventFile->Close();

  //Print out some statistics
  cout<<"More then one proton hits: "<<multiEventCounter<<endl;
  cout<<"There are "<<musicsChargeMatch<<" charges matching the criteria"<<endl;
  cout<<"In this case, CALIFA did not had a good p2p event "<<musicsChargeMatch-califaP2PCandidate<<" times"<<endl;
  for(mult=0; mult<MAXMULT-1; mult++) cout<<"Mult"<<mult<<": "<<mult_dist[mult]<<", ";
  cout<<"Mult"<<MAXMULT-1<<"+: "<<mult_dist[MAXMULT-1]<<endl;
}
//======================================================================================================================


void loadtrees()
{
  // Prepare the input file and get access to the tree
  InputFileName.Append(INPUTDIR);
  InputFileName.Append(FILENAME);
  eventTree = new TChain("evt");
  OutTree = new TTree("evt","Merged hits with CALIFA events");
  
  std::ifstream RunList(RUNLISTNAME, std::ios::in);
  if(!RunList.is_open()) std::cerr <<"No run summary found\n";
  int runnumcsv[500], targetpos[500], musicgain[500], junk[500], FRSsetting[500];
  string dummyline;
  char dumchar;
  double brhocsv[500];
  std::getline (RunList, dummyline);
  Int_t i=0;
  while(true){
    i++;
    if(i > 400 || !RunList.good()){
      //std::cerr << "No info for run found" <<std::endl;
      break;
    }
    RunList>>runnumcsv[i]>>dumchar>>FRSsetting[i]>>dumchar>>brhocsv[i]>>dumchar>>targetpos[i]>>dumchar>>musicgain[i]>>dumchar>>junk[i];
    if(FRSsetting[i]!=FRS) continue;
    if(junk[i]!=0) continue;
    if(targetpos[i]<posmin || targetpos[i]>posmax) continue;
    //cout<<runnumcsv[i]<<" "<<dumchar<<" "<<FRSsetting[i]<<" "<<dumchar<<" "<<brhocsv[i]<<" "<<dumchar<<" "<<targetpos[i]<<" "<<dumchar<<" "<<musicgain[i]<<" "<<dumchar<<" "<<junk[i]<<endl;

    TString tmpfilename(InputFileName);
    //tmpfilename = new TString(InputFileName);
    tmpfilename.ReplaceAll("RUN",Form("%04d", runnumcsv[i]));
    cout <<"Loading: " << tmpfilename << endl;
    eventTree -> Add(tmpfilename);
  }
  if(nEvents==-1)
    {
      nEvents = eventTree->GetEntries();
    }
  cout << "Total events:" << nEvents << endl;
  //eventTree->Print();
  SetTreeBranches();
}




void analysisLoop(int max_events)
{
  if(max_events>0) nEvents = max_events;
  // Main loop thorough the events
  for(Long64_t i = 0; i < nEvents; i++)
    {
      if((i+1)%100000 == 0 || i == nEvents-1)
        {
	  cout<<"\r"<< i+1 <<"/"<< nEvents << setprecision(3)<< " events ("<< (Float_t)i/(Float_t)nEvents*100. << "%) processed. " 
	      <<"Good Events: "<<musicsChargeMatch<<" p2p candidate: "<<califaP2PCandidate<<string(10, ' ')<<flush;
	}
      //Reset the variables to avoid wrong values and to have the findP2PCandidate function working properly
      resetVariables();

      //Select the next event and look for entries
      eventTree->GetEntry(i);
      nCalifaHits = fCalifaHitData->GetEntries();
      //nR3BMusicHits = fR3BMusicHitData->GetEntries();
      nFrsHits = fSofFrsData->GetEntries();
      nFragmentHits = fSofTrackingData->GetEntries();

      //Check if both Musics registered a charge
      if(nFrsHits==0 || nFragmentHits == 0) continue;
      //Prepare the parameters for the charge selection
      //r3bMusicHit = (R3BMusicHitData*)fR3BMusicHitData->At(0);
      FrsHit = (R3BSofFrsData*)fSofFrsData->At(0);
      //r3bMusicCharge = r3bMusicHit->GetZcharge();
      r3bMusicCharge = FrsHit->GetZ();
    
      FragmentHit = (R3BSofTrackingData*)fSofTrackingData->At(0);
      TwimCharge = FragmentHit->GetZ();

      //Make the charge selection
      //if(r3bMusicCharge > R3BMusicCutMin && r3bMusicCharge < R3BMusicCutMax && TwimCharge > twimMusicCutMin && TwimCharge < twimMusicCutMax)
      
      musicsChargeMatch++;
      //If the charges are correct, process the CALIFA data and search for p2p candidates, if we have min. two hits in CALIFA with any energy.
      if(nCalifaHits==0) continue;
      
      //If yes, search for high energy hits, as a proton hit candidate.
      findP2PCandidate();
      if(mult<MAXMULT)
	{
	  mult_dist[mult]++;
	}
      else if(mult>=MAXMULT)
	{
	  mult_dist[MAXMULT-1]++;
	}
      
      //Check if the the search was a success and we have exactly two proton hits
      //if(mult!=2) continue;
      if(mult < 2) continue;
      califaP2PCandidate++;
      //Tranform the angles from rad to degree
      theta[0]=theta[0]*TMath::RadToDeg();
      theta[1]=theta[1]*TMath::RadToDeg();
      phi[0]=phi[0]*TMath::RadToDeg();
      phi[1]=phi[1]*TMath::RadToDeg();
      //Calculate the angles between the two protons and the sum energy
      thetaP2P=theta[0]+theta[1];
      phiP2P=abs(phi[0]-phi[1]);
      sumEnergy=energy[0]+energy[1];

      //Whatever you want to do with the data, do it here! Above, you see which values were stored for each hit.
      //cout<<r3bMusicCharge<<" "<<TwimCharge<<endl;
      if(ftpat!=0) cout<<ftpat<<endl;
      //As an example I fill here histos with my selected events
      histCalifaHighEnTheta->Fill(thetaP2P);
      histCalifaTheta1vsTheta2->Fill(theta[0],theta[1]);
      histCalifaTheta1vsTheta2->Fill(theta[1],theta[0]);
      histCalifaPhi1vsPhi2->Fill(phi[0],phi[1]);
      histCalifaPhi1vsPhi2->Fill(phi[1],phi[0]);
      //
      OutTree->Fill();
    }
  cout<<endl;
}





void findP2PCandidate()
{
  for(Int_t j =0; j<nCalifaHits; j++)
    {
      //Extract the required information for the search
      califaHit = (R3BCalifaHitData*)fCalifaHitData->At(j);
      Double_t tmpProtonPhi = califaHit->GetPhi();
      Double_t tmpProtonTheta = califaHit->GetTheta();
      Double_t tmpProtonEnergy = califaHit->GetEnergy()*dopplerCorrection(tmpProtonTheta,BETA);
       
       
      // Check if the energy of the hit is above the energy threshold for a real proton hit
      if(tmpProtonEnergy > EMINPROTON)
	{
	  if(mult < 2)
	    {
	      energy[mult]= tmpProtonEnergy;
	      theta[mult]= tmpProtonTheta;
	      phi[mult]= tmpProtonPhi;
	    }
	  else
	    {
	      mult==2?multiEventCounter++:1;
	    }
	  mult++;
	}
      //Clear the variables for the next iteration
      califaHit->Reset();
    }
}






void writeHistos()
{
  // Writing the results to the output file and closing it.
  sprintf(outputFile, "%s%s_%s.root", OUTPUTDIR, OUTFILENAME, targetname.Data());
  cout<<"OUTPUT: "<<outputFile<<endl;
  resultFile = new TFile(outputFile,"recreate");
  histCalifaPhi1vsPhi2->Scale(0.5);
  histCalifaTheta1vsTheta2->Scale(0.5);
  histCalifaHighEnTheta->Write();
  histCalifaTheta1vsTheta2->Write();
  histCalifaPhi1vsPhi2->Write();
  //OutTree->Print();
  //OutTree->Scan();
  OutTree->Write();
  cout<<OutTree->GetEntries()<<" events are filled in OutTree."<<endl;
  resultFile->Close();
}


void SetTreeBranches(){
  // Access the data structure to extract the values
  //fEventHeader = new TClonesArray("EventHeader", 5);
  //fEventHeader = new TClonesArray("R3BEvtHeader", 5);
  //eventTree->SetBranchAddress("EventHeader",&fEventHeader);
  //eventTree->SetBranchAddress("R3BEvtHeader",&fEventHeader);
  //eventTree->SetBranchAddress("cbmout",&fEventHeader);
  //eventTree->SetBranchAddress("EventHeader.fTpat",&fTpat);
  //eventTree->SetBranchAddress("EventHeader.fTpat",&ftpat);
  //
  fCalifaHitData = new TClonesArray("R3BCalifaHitData", 5);
  eventTree->SetBranchAddress("CalifaHitData",&fCalifaHitData);
  /*
  fR3BMusicHitData = new TClonesArray("R3BMusicHitData", 5);
  eventTree->SetBranchAddress("MusicHitData",&fR3BMusicHitData);
  */
  fSofFrsData = new TClonesArray("R3BSofFrsData", 5);
  eventTree->SetBranchAddress("SofFrsData",&fSofFrsData);
  //
  fSofTrackingData = new TClonesArray("R3BSofTrackingData", 5);
  eventTree->SetBranchAddress("SofTrackingData",&fSofTrackingData);
  //////
  //////
  //OutTree->Branch("EventHeader",&fEventHeader);
  //OutTree->Branch("fTpat",&ftpat);
  OutTree->Branch("CalifaProtonMult",&mult);
  OutTree->Branch("CalifaTheta",theta, "theta[2]/D");
  OutTree->Branch("CalifaPhi",phi, "phi[2]/D");
  OutTree->Branch("CalifaEdopp",energy, "energy[2]/D");
  OutTree->Branch("CalifaSumE", &sumEnergy);
  OutTree->Branch("CalifaOpeningTheta", &thetaP2P);
  OutTree->Branch("CalifaOpeningPhi", &phiP2P);
  OutTree->Branch("CalifaHitData",&fCalifaHitData);
  OutTree->Branch("SofFrsData",&fSofFrsData);
  OutTree->Branch("SofTrackingData",&fSofTrackingData);
}


void resetVariables()
{
  ftpat=DEFAULTINTEGER;
  mult=DEFAULTINTEGER;
  //
  //CALIFA
  energy[1]=DEFAULTDOUBLE;
  theta[1]=DEFAULTDOUBLE;
  phi[1]=DEFAULTDOUBLE;

  energy[0]=DEFAULTDOUBLE;
  theta[0]=DEFAULTDOUBLE;
  phi[0]=DEFAULTDOUBLE;

  thetaP2P=DEFAULTDOUBLE;
  phiP2P=DEFAULTDOUBLE;
  sumEnergy=DEFAULTDOUBLE;

  //multiEventFlag=0;
  // MUSICs
  TwimCharge=DEFAULTDOUBLE;
  r3bMusicCharge=DEFAULTDOUBLE;
    
  //TClonesArrays
  /*
  if(r3bMusicHit)
    {
      r3bMusicHit->Clear();
    }
    */
  if(FrsHit)
    {
      FrsHit->Clear();
    }
  if(FragmentHit)
    {
      FragmentHit->Clear();
    }
  if(califaHit)
    {
      califaHit->Clear();
    }  
}




Double_t dopplerCorrection(Double_t theta, Double_t BETA)
{
  return (1.-BETA*cos(theta))/sqrt(1.-BETA*BETA);
}
