#include "QualityControl/MonitorObject.h"
#include "QualityControl/MonitorObjectCollection.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGrid.h"
#include <string>
#include <iostream>
#include <fstream>

void readQC_loop(){
  
  TGrid::Connect("alien://");

  double phi = 4.;  // choose a point in the problematic phi region
  double threshold = 0.15;  // the threshold below which value we consider the phi-distribution as bad, probably 0.15 is low enough.

  // inputfile
  std::ifstream inputlist("QCfiles.dat");
  std::string fname;
  // outputfile
  std::ofstream outfile;
  outfile.open("problematic_ctfs.dat");
  // loop over inputfile
  while(inputlist>>fname){
    o2::quality_control::core::MonitorObjectCollection *mycoll;
    TFile *infl = TFile::Open(fname.c_str());
    infl->GetObject("ITS/Tracks",mycoll);
    o2::quality_control::core::MonitorObject *mo = (o2::quality_control::core::MonitorObject*)mycoll->At(4);
    TH1F *h1 = (TH1F*)mo->getObject();
    // std::cout<<h1->GetBinContent(2)<<std::endl; // sanity check print the content of the 2nd bin
    if(h1->GetBinContent(h1->FindBin(phi)) < threshold) {
      std::cout<<"***************************** FOUND "<<fname<<std::endl;
      outfile<<fname<<std::endl; // write filepaths of problematic ctf's to outfile
    }
    delete infl;
  }
}

