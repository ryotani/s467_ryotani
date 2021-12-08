#define NUMPADDLE 28
TProof *p =  TProof::Open("");

TCanvas *c;
TH1F *htime[NUMPADDLE], *hpos[NUMPADDLE];
TH2F *htofhit, *hposhit;
TF1 *func[NUMPADDLE], *funcpos[NUMPADDLE];

//TTree *tree;
//TFile *f;
TChain *ch;

void draw_tofw_hit(int runnum){
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle();
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat();
  TString filename = runnum<1000?Form("/u/taniuchi/s467/rootfiles/rootfiletmp/s467_FRSTree_Setting*%04d_ToFWhitpar.root",  runnum): "/u/taniuchi/s467/rootfiles/rootfiletmp/s467_FRSTree_Setting*_ToFWhitpar.root";
  TString pdfout =   runnum<1000?Form("~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/output/tofw_hit_%04d.pdf", runnum): "~/s467/ana/R3BRoot_ryotani/sofia/macros/s467_ryotani/tofw/output/tofw_hit.pdf";
  
  //TFile *f= TFile::Open(filename, "READ");
  //tree = (TTree*)(f->Get("evt"));
  ch = new TChain("evt");
  ch -> Add(filename);
  ch -> SetProof();
  cout << "Filename: " << filename <<endl;
  cout << "Num Entry = " << ch->GetEntries()<<endl;
  c = new TCanvas("c","c",1200,1000);

  c->Print(pdfout + "[");
  htofhit = new TH2F("htofhit", "ToF in TofW;Paddle ID (1..28);Averaged time (ns)", 30,-0.5,29.5,400,-20,20);
  ch->Draw("TofWHitData.fTime:TofWHitData.fPaddleId>>htofhit", "", "colz");
  c->Print(pdfout);
  hposhit = new TH2F("hposhit", "Position in TofW;X position (mm);Time difference (ns)", 100, -500, 500 ,400,-20,20);
  ch->Draw("TofWHitData.fY:TofWHitData.fX>>hposhit", "", "colz");
  c->Print(pdfout);
  c->Print(pdfout + "]");
}


void draw_tofw_hit(){
  draw_tofw_hit(366);
}
