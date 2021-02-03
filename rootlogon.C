#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#endif

//void AddIncludePath(std::string dir);
//void AddLinkedLibs();
//void LoadModule(std::string dir);

void rootlogon() {
  /*
  gSystem->Load("libGretina");
  gSystem->Load("libS800");
  gSystem->Load("libScaler");
  gSystem->Load("libSettings");
  gSystem->Load("libRunInfo");
  char* c = getenv("TARTSYS");
  if(!c) {
    std::cout << "set environment variable \"TARTSYS\"" << std::endl;
    std::cout << "quit ROOT" << std::endl;      
    gROOT->ProcessLine(".q");
  }
  std::string install_dir(c);

  AddIncludePath(install_dir);
  AddLinkedLibs();
  LoadModule(install_dir);

  std::cout << "\n Manual of ANAROOT \n  http://be.nucl.ap.titech.ac.jp/~kondo/moin/moin.cgi/ANAROOT/Manual\n" << std::endl;
  */
  //Using the RootTools
  //gROOT->ProcessLine(".L ~/ROOT/Macros/Utilities/RootTools.C+");

  //Base Style
  gROOT->SetStyle("Plain");
  // gROOT->SetStyle("Modern");
    // //gROOT->SetStyle("Classic");

  //Force Style
  gStyle->SetHistFillColor(0);
  //gStyle->SetHistFillStyle(3002);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetFuncColor(kRed);
  gStyle->SetNdivisions(210,"XY");
  //gStyle->SetFrameLineWidth(2);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetStatX(0);  
  gStyle->SetStatY(0);  
  //gStyle->SetPalette(0);
  gStyle->SetOptLogz(1);
  gStyle->SetOptTitle(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0); //gStyle->SetOptStat(1111111);
  //gStyle->SetPadBorderMode(1);
  gStyle->SetOptDate(0);

  gStyle->SetLabelFont(43,"XYZ");
  gStyle->SetTitleFont(43,"XYZ");
  gStyle->SetTitleFont(42,"T");
  //gStyle->SetTextFont(43);
  //gStyle->SetStatFont(43);

  gStyle->SetLabelSize(15,"XYZ");
  gStyle->SetTitleSize(15,"XYZ");
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1.5,"Y");
  //gStyle->SetTitleOffset(0,"T");
  //gStyle->SetTitleSize(25,"T");
  //gStyle->SetTextSize(25);
  //gStyle->SetStatSize(25);

  
  TColor *dummyColor = new TColor();
  //dummyColor->SetPalette(54,0);
  //dummyColor->InvertPalette();
  //dummyColor->SetPalette(kInvertedDarkBodyRadiator,0);
  dummyColor->SetPalette(kBlueYellow,0);
  //dummyColor->SetPalette(kDeepSea,0);
  //dummyColor->InvertPalette();
  
  /*
  UInt_t Number = 4;
  Double_t Red[4]   = { 0.50, 0.00, 0.99, 0.99};
  Double_t Green[4] = { 0.99, 0.00, 0.00, 0.99};
  Double_t Blue[4]  = { 0.99, 0.99, 0.00, 0.00};
  Double_t Stops[4] = { 0.00, 0.30, 0.70, 0.99};
  dummyColor->CreateGradientColorTable(Number,Stops,Red,Green,Blue,100);
  */
}
/*
void AddIncludePath(std::string install_dir) {
  std::vector<std::string> include;
  include.push_back("-I"+install_dir+"/include");
  //  include.push_back("`xml2-config --cflags`");

  std::vector<std::string>::iterator it = include.begin();
  while(it != include.end()){
    gSystem->AddIncludePath((*it).c_str());
    std::cout << "add include path : " << *it << std::endl;
    ++it;
  }
}

void AddLinkedLibs()  {
//  std::vector<std::string> include;
//  include.push_back("`xml2-config --libs`");

//  std::vector<std::string>::iterator it = include.begin();
//  while(it != include.end()){
//    gSystem->AddLinkedLibs((*it).c_str());
//    std::cout << "add linked libs : " << *it << std::endl;
//    ++it;
//  }
}

void LoadModule(std::string install_dir) {
  std::vector<std::string> modules;
  modules.push_back("libXMLParser.so");
  modules.push_back(install_dir+"/lib/"+"libanaroot.so"); // load at once
  //modules.push_back(install_dir+"/lib/"+"libananadeko.so"); // load each modules one by one
  //modules.push_back(install_dir+"/lib/"+"libanacore.so");
  //modules.push_back(install_dir+"/lib/"+"libanaeurica.so");
  //modules.push_back(install_dir+"/lib/"+"libanabrips.so");
  //modules.push_back(install_dir+"/lib/"+"libanasamurai.so");
  //modules.push_back(install_dir+"/lib/"+"libanaminos.so");
  //modules.push_back(install_dir+"/lib/"+"libanaanaloop.so");
  //modules.push_back(install_dir+"/lib/"+"libanadali.so");

  std::vector<std::string>::iterator it = modules.begin();
  while(it != modules.end()){
    std::cout << "reading " << *it << std::endl;
    if(gSystem->Load((*it).c_str()) < 0){
      std::cout << "cannnot read in " << *it << std::endl;      
    }
    ++it;
  }
}
*/
