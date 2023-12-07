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

  // execute this script (without exiting ROOT!!! so no `-q`), then open a TBrowser to inspect the file you just loaded.
  
  // alien:///alice/data/2023/LHC23zc/537903/cpass0/1020/QC/001/QC.root
  // TFile *f = TFile::Open("alien:///alice/data/2023/LHC23zzh/544098/cpass0/2230/QC/001/QC.root");
  // TFile *f = TFile::Open("alien:///alice/data/2023/LHC23zi/539107/apass1/1820/QC/001/QC.root");
  // TFile *f = TFile::Open("alien:///alice/data/2023/LHC23zc/537903/cpass0/1000/o2_ctf_run00537903_orbit0349958816_tf0002069728_epn059/QC.root");
  // TO BE TESTED!!!
  // for now you can open root in QC env and copy these lines, it should run

  TGrid::Connect("alien://");

  std::string filename = "alien:///alice/data/2023/LHC23zzh/544098/cpass0/2230/QC/001/QC.root";
  TFile *f = TFile::Open(filename.c_str());
}

