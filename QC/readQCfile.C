#include "QualityControl/MonitorObject.h"
#include "QualityControl/MonitorObjectCollection.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGrid.h"
#include <string>
#include <iostream>
#include <fstream>

void readQCfile(){
  
  // alien:///alice/data/2023/LHC23zc/537903/cpass0/1020/QC/001/QC.root
  // TFile *f = TFile::Open("alien:///alice/data/2023/LHC23zc/537903/cpass0/1020/QC/001/QC.root");
  // TFile *f = TFile::Open("alien:///alice/data/2023/LHC23zi/539107/apass1/1820/QC/001/QC.root");
  // TFile *f = TFile::Open("alien:///alice/data/2023/LHC23zc/537903/cpass0/1000/o2_ctf_run00537903_orbit0349958816_tf0002069728_epn059/QC.root");
  // TO BE TESTED!!!
  // for now you can open root in QC env and copy these lines, it should run

  // TGrid::Connect("alien://");

  // std::string filename = "alien:///alice/sim/2023/LHC23d1e/0/520259/001/QC/ITSTracksClusters.root";
  // TFile *f = TFile::Open(filename.c_str());

  // TFile *hitmap = TFile::Open("o2-itsmft-NoiseMap_1685542123734.root");
  o2::itsmft::ChipMappingITS mp; // declare this at the beginning of the macro
  o2::itsmft::NoiseMap * mymap;
  TFile *f = new TFile("o2-itsmft-NoiseMap_1685542123734.root");
  f->GetObject("ccdb_object", mymap);
  for (int i=0; i< mymap->size(); i++){
      if (! mymap->getChip(i).empty()) 
          cout<<mp.getChipNameHW(i)<<endl;; 
  }
  // TBrowser b;

  // std::ifstream inputlist("228filelist.dat");
  // std::string fname;
  // while(inputlist>>fname){
  //   o2::quality_control::core::MonitorObjectCollection *mycoll;
  //   TFile *infl = TFile::Open(fname.c_str());
  //   infl->GetObject("ITS/Tracks",mycoll);
  //   o2::quality_control::core::MonitorObject *mo = (o2::quality_control::core::MonitorObject*)mycoll->At(3);
  //   TH1F *h1 = (TH1F*)mo->getObject();
  //   std::cout<<h1->GetBinContent(2)<<std::endl;
  //   if(h1->GetBinContent(h1->FindBin(0.7)) < 0.15) {
  //     std::cout<<"***************************** FOUND "<<fname<<std::endl;
  //   }
  //   delete infl;
  // }
}

