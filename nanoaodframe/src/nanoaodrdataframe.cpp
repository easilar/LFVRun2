//============================================================================
// Name        : nanoaodrdataframe.cpp
// Author      : Suyong Choi
// Version     :
// Copyright   : suyong@korea.ac.kr, Korea University, Department of Physics
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "NanoAODAnalyzerrdframe.h"
#include "TopLFVAnalyzer.h"
#include "TChain.h"

using namespace std;
using namespace ROOT;

int main(void) {
  TChain c1("outputTree");
  c1.Add("processed/2016data/2016b/00/nanoAOD_3_analyzed.root");

  TopLFVAnalyzer nanoaodrdf(&c1, "testout.root", "2016pre", "", "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt", "", 1);
  nanoaodrdf.setupAnalysis();
  nanoaodrdf.run(false, "Events");

  return EXIT_SUCCESS;
}
