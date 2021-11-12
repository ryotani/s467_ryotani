TChain* ch;
TMultiDimFit* fit;
Float_t AoQ_S2_Cave, brho_frs, MusicZ, FragZ, TwimTheta, Mw1_X, Mw2_X, Mw3_X, Mw1_Y, Mw2_Y, Mw3_Y, Tofw_Y, FragTof, FragAoQ;
UChar_t Tofw_Paddle;
Long64_t nentry=0, nvalidentry=0;
Double_t* finput();

TChain* loadfiles(TString filename);
Int_t initbranch();
TMultiDimFit* initMDF(Int_t nVars = 4);
void multdimfit();
void multdimfit(TString filename);
Bool_t isgoodevent();

//
void multdimfit(){
  TString filename = "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/fragment/output/mktree_fragment_Nov_empty.root";
  multdimfit(filename);
}

const Int_t nVars = 3;
Double_t input[nVars] = {0.};
Double_t* finput(){
  input[0] = 1000. * TwimTheta;
  input[1] = 0.05 * (Mw1_X + Mw2_X);
  //input[2] = brho_frs;
  //input[2] = 0.5 * (Mw1_Y + Mw2_Y);
  //input[3] = 0.5 * (Mw3_Y - 0.5 * (Mw1_Y + Mw2_Y));
  input[2] = Mw3_X;
  return input;
}

void multdimfit(TString filename){
  ch = loadfiles(filename);
  nentry = ch->GetEntries();
  cout << "Entries:" << nentry <<endl;
  //
  if (initbranch()) return;
  //
  fit = initMDF(nVars);
  //
  //for(Long64_t i=0; i<1000000; i++){
  for(Long64_t i=0; i<nentry; i++){
    ch->GetEntry(i);
    //
    if(!isgoodevent()) continue;
    nvalidentry++;
    //
    //cout<<nvalidentry<<" in "<< i <<" "<<brho_frs <<" "<< MusicZ<<" "<< FragZ<<" "<< TwimTheta<<" "<< Mw1_X<<" "<< Mw2_X<<" "<< Mw3_X<<" "<< Mw1_Y<<" "<< Mw2_Y<<" "<< Mw3_Y<<endl;
    //auto input = finput();
    /*
    cout<<brho_frs;
    for(int j=0; j<nVars; j++) cout<<" "<<input[j];
    cout<<endl;
    */
    fit -> AddRow(finput(), brho_frs, 1);
    //fit -> AddRow(finput(), Mw3_X, 1);
  }
  //
  cout<<nvalidentry<<" entries will be used in "<<nentry<<" ("<<(double)nvalidentry/(double)nentry*100.<<"%)"<<endl<<endl;
  // Find the parameterization -- doing the selection and fit.
  fit->FindParameterization();
  //================================================
  // part 4
  // print the fitting results.
  // P Parameters
  // S Statistics
  // C Coefficients
  // R Result of parameterisation
  // F Result of fit
  // K Correlation Matrix
  // M Pretty print formula

  cout << "#==========Let's print coefficeints========#" << endl;
  fit->Print("c");

  cout << "#==========Let's print formula   ==========#" << endl;
  fit->Print("M");
  
  cout <<endl << "Fit procedure"<<endl;
  fit->Fit("M");
  cout << "#==========Let's print coefficeints========#" << endl;
  fit->Print("c");

  cout << "#==========Let's print formula   ==========#" << endl;
  fit->Print("M");

  fit->Print("R");

  // Write code to file, so that we can reused our parameterization.
  //fit->MakeCode("MDF.C");
  for(Long64_t i=0; i<1000; i++){
    //for(Long64_t i=0; i<nentry; i++){
    ch->GetEntry(i%10000);
    if(!isgoodevent()) continue;
    /*
    input[0] = 1000. * TwimTheta;
    input[1] = 0.05 * (Mw1_X + Mw2_X);
    input[2] = brho_frs;
    */
    Double_t Mw3_X_eval;
    Mw3_X_eval = fit->Eval(finput());
    cout<<Mw3_X<<" Eval: "<<Mw3_X_eval<<" brho_frs: "<<brho_frs<<" "<<Mw3_X_eval - Mw3_X<<" AoQ:"<< AoQ_S2_Cave<<" "<<FragAoQ<<endl;
  }

  
  //delete fit;
  
}

TMultiDimFit* initMDF(Int_t nVars = 4){
  TMultiDimFit *fit = new TMultiDimFit(nVars, TMultiDimFit::kMonomials,"v");

  // set up the max power for each variable.
  // i.e. x1^p_L1, and p_L1 < p_1Max.
  auto mPowers = new Int_t(nVars);
  for(int i=0; i<nVars; i++) mPowers[i] = 6; // Highest order
  fit->SetMaxPowers(mPowers);

  // set up the max number for F_L
  fit->SetMaxFunctions(100000);
  fit->SetMaxStudy(100000);

  // set up the max selected F_L
  fit->SetMaxTerms(1000);

  // set the selection condition for a given F_L, we want:
  // Q_L = ( p_L1/p_1Max + p_L2/p_2Max + ... ) < Q
  fit->SetPowerLimit(1); // the power control limit Q.

  // test1 and test2 to whether we include the Lmax term.
  //fit->SetMinAngle(1);
  //fit->SetMaxAngle(80);

  // def a good fit as:
  // esplion = S / D^2
  // S   = sum over all ( data - simulated value )
  // D^2 = sum over all (data^2)
  //fit->SetMinRelativeError(.01); // esplion = 0.01
  fit->SetMinRelativeError(.01); // esplion = 0.01

  cout << "#==========Let's print Parameters==========#" << endl;
  fit->Print("p");
  
  
  return fit;
}

Bool_t isgoodevent() {
    if(isnan(brho_frs)||isnan(MusicZ)||isnan(FragZ)||isnan(TwimTheta)||isnan(Mw1_X)||isnan(Mw2_X)||isnan(Mw3_X)||isnan(Mw1_Y)||isnan(Mw2_Y)||isnan(Mw3_Y)||isnan(Tofw_Y)||isnan(FragTof)||
       min(MusicZ,FragZ)<10||min(Mw1_X,Mw2_X)<-100||min(Mw1_Y,Mw2_Y)<-100||Mw3_X<-500||Mw3_Y<-500||
       abs(MusicZ-FragZ)>0.5||abs(FragAoQ-AoQ_S2_Cave)<0.02||
       abs(MusicZ-20)>0.5) return false;
    else  return true;
}

Int_t initbranch(){
  ch->SetBranchAddress("AoQ_S2_Cave",&AoQ_S2_Cave);
  ch->SetBranchAddress("Brho_S2_Cave",&brho_frs);
  ch->SetBranchAddress("MusicZ", &MusicZ);
  ch->SetBranchAddress("FragZ", &FragZ);
  ch->SetBranchAddress("TwimTheta", &TwimTheta);
  ch->SetBranchAddress("Mw1_X", &Mw1_X);
  ch->SetBranchAddress("Mw2_X", &Mw2_X);
  ch->SetBranchAddress("Mw3_X", &Mw3_X);
  ch->SetBranchAddress("Mw1_Y", &Mw1_Y);
  ch->SetBranchAddress("Mw2_Y", &Mw2_Y);
  ch->SetBranchAddress("Mw3_Y", &Mw3_Y);
  ch->SetBranchAddress("Tofw_Y", &Tofw_Y);
  ch->SetBranchAddress("FragTof", &FragTof);
  ch->SetBranchAddress("Tofw_Paddle", &Tofw_Paddle);
  ch->SetBranchAddress("FragAoQ", &FragAoQ);
  return 0;
}

TChain* loadfiles(TString filename){
  ch = new TChain("Tree");
  ch -> Add(filename);
  return ch;
}
