#include "QualityControl/MonitorObject.h"
#include "QualityControl/MonitorObjectCollection.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGrid.h"
#include <string>
#include <iostream>
#include <fstream>

void readQC(){
  
  TGrid::Connect("alien://");

  std::ifstream inputlist("228filelist.dat");
  std::string fname;
  while(inputlist>>fname){
    o2::quality_control::core::MonitorObjectCollection *mycoll;
    TFile *infl = TFile::Open(fname.c_str());
    infl->GetObject("ITS/Tracks",mycoll);
    o2::quality_control::core::MonitorObject *mo = (o2::quality_control::core::MonitorObject*)mycoll->At(3);
    TH1F *h1 = (TH1F*)mo->getObject();
    std::cout<<h1->GetBinContent(2)<<std::endl;
    if(h1->GetBinContent(h1->FindBin(0.7)) < 0.15) {
      std::cout<<"***************************** FOUND "<<fname<<std::endl;
    }
    delete infl;
  }
}

