/*
 * NanoAODAnalyzerrdframe.h
 *
 *  Created on: Sep 30, 2018
 *      Author: suyong
 */

#ifndef NANOAODANALYZERRDFRAME_H_
#define NANOAODANALYZERRDFRAME_H_

#include "TTree.h"
#include "TFile.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Vector4D.h"
#include "BTagCalibrationStandalone.h"
#include "WeightCalculatorFromHistogram.h"

#include <string>
#include "json/json.h"

#include "utility.h" // floats, etc are defined here
#include "RNodeTree.h"
#include "JetCorrectorParameters.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectionUncertainty.h"
#include "JetResolution.h"
#include "TauSFTool.h"

#include "correction.h"

using namespace ROOT::RDF;


class NanoAODAnalyzerrdframe {
  using RDF1DHist = RResultPtr<TH1D>;

public:
  NanoAODAnalyzerrdframe(std::string infilename, std::string intreename, std::string outfilename, std::string year="", std::string syst="", std::string jsonfname="", string globaltag="", int nthreads=1);
  NanoAODAnalyzerrdframe(TTree *t, std::string outfilename, std::string year="", std::string syst="", std::string jsonfname="", string globaltag="", int nthreads=1);
  virtual ~NanoAODAnalyzerrdframe();
  void setupAnalysis();

  // object selectors
  // RNode is in namespace ROOT::RDF
  bool readjson();
  void selectElectrons();
  void selectMuons();
  void selectJets(std::vector<std::string> jes_var);
  void skimJets();
  void applyBSFs(std::vector<string> jes_var);
  void selectTaus();
  void selectMET();
  void selectFatJets();
  void matchGenReco();
  void calculateEvWeight();
  void storeEvWeight();
  void topPtReweight();
  virtual void defineMoreVars() = 0; // define higher-level variables from basic ones, you must implement this in your subclassed analysis code

  void addVar(varinfo v);

  template <typename T, typename std::enable_if<!std::is_convertible<T, std::string>::value, int>::type = 0>
  void defineVar(std::string varname, T function,  const RDFDetail::ColumnNames_t &columns = {})
  {
    _rlm = _rlm.Define(varname, function, columns);
  };

  void addVartoStore(std::string varname);
  void addCuts(std::string cut, std::string idx);
  virtual void defineCuts() = 0; // define a series of cuts from defined variables only. you must implement this in your subclassed analysis code
  void add1DHist(TH1DModel histdef, std::string variable, std::string weight, string syst="", string mincutstep="", string maxcutstep="");
  virtual void bookHists() = 0; // book histograms, you must implement this in your subclassed analysis code

  void setupCuts_and_Hists();
  void drawHists(RNode t);
  void run(bool saveAll=true, std::string outtreename="Events");
  void setTree(TTree *t, std::string outfilename);
  void setupTree();

  bool _isSkim = false;
  bool _isHTstitching = false;
  std::string _outfilename;
  std::string _syst;

private:
  ROOT::RDataFrame _rd;
  bool _isData;
  bool _jsonOK;
  std::string _year;
  bool _isRun16pre = false;
  bool _isRun16post = false;
  bool _isRun16 = false;
  bool _isRun17 = false;
  bool _isRun18 = false;
  // you MUST copy syst names from the output of 'python skimcsv.py'
  inline static std::vector<std::string> btag_var = {"central",
                "hfup", "hfdown", "lfup", "lfdown", "hfstats1up", "hfstats1down",
                "hfstats2up", "hfstats2down", "lfstats1up", "lfstats1down",
                "lfstats2up", "lfstats2down", "cferr1up", "cferr1down", "cferr2up", "cferr2down"};
  inline static std::vector<std::string> jes_var_2016 = {"jesAbsoluteup", "jesAbsolutedown",
                "jesAbsolute_2016up", "jesAbsolute_2016down", "jesBBEC1up", "jesBBEC1down",
                "jesBBEC1_2016up", "jesBBEC1_2016down", "jesFlavorQCDup", "jesFlavorQCDdown",
                "jesRelativeBalup", "jesRelativeBaldown", "jesRelativeSample_2016up", "jesRelativeSample_2016down"};
  inline static std::vector<std::string> jes_var_2017 = {"jesAbsoluteup", "jesAbsolutedown",
                "jesAbsolute_2017up", "jesAbsolute_2017down", "jesBBEC1up", "jesBBEC1down",
                "jesBBEC1_2017up", "jesBBEC1_2017down", "jesFlavorQCDup", "jesFlavorQCDdown",
                "jesRelativeBalup", "jesRelativeBaldown", "jesRelativeSample_2017up", "jesRelativeSample_2017down"};
  inline static std::vector<std::string> jes_var_2018 = {"jesAbsoluteup", "jesAbsolutedown",
                "jesAbsolute_2018up", "jesAbsolute_2018down", "jesBBEC1up", "jesBBEC1down",
                "jesBBEC1_2018up", "jesBBEC1_2018down", "jesFlavorQCDup", "jesFlavorQCDdown",
                "jesRelativeBalup", "jesRelativeBaldown", "jesRelativeSample_2018up", "jesRelativeSample_2018down",
                "jesHEMup", "jesHEMdown"};
  floats PDFWeights;
  std::string _jsonfname;
  std::string _globaltag;
  TFile *_inrootfile;
  TFile *_outrootfile;
  std::vector<std::string> _outrootfilenames;
  RNode _rlm;
  std::map<std::string, RDF1DHist> _th1dhistos;
  //bool helper_1DHistCreator(std::string hname, std::string title, const int nbins, const double xlow, const double xhi, std::string rdfvar, std::string evWeight);
  void helper_1DHistCreator(std::string hname, std::string title, const int nbins, const double xlow, const double xhi, std::string rdfvar, std::string evWeight, RNode *anode);
  std::vector<std::string> _originalvars;
  std::vector<std::string> _selections;

  std::vector<hist1dinfo> _hist1dinfovector;
  std::vector<varinfo> _varinfovector;
  std::vector<cutinfo> _cutinfovector;

  std::vector<std::string> _varstostore;
  std::map<std::string, std::vector<std::string>> _varstostorepertree;

  Json::Value jsonroot;

  RNodeTree _rnt;
  RNodeTree *currentnode;
  bool isDefined(string v);

  void setupJetMETCorrection(std::string globaltag, const std::vector<std::string> var = std::vector<std::string>(), std::string jetalgo="AK4PFchs", bool dataMc=false);

};

#endif /* NANOAODANALYZERRDFRAME_H_ */
