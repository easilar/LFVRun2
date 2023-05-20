#include "TauSFTool.h"
#include <iostream> // std::cerr, std::endl
#include <iomanip>
#include <assert.h> // assert

TFile* ensureTFile(const TString filename, bool verbose=false){
  if(verbose)
    std::cout << "Opening " << filename << std::endl;
  TFile* file = new TFile(filename);
  if(!file or file->IsZombie()) {
    std::cerr << std::endl << "ERROR! Failed to open input file = '" << filename << "'!" << std::endl;
    assert(0);
  }
  return file;
}

TH1* extractTH1(const TFile* file, const std::string& histname){
  TH1* hist = dynamic_cast<TH1*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!hist){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  return hist;
}

const TF1* extractTF1(const TFile* file, const std::string& funcname){
  const TF1* function = dynamic_cast<TF1*>((const_cast<TFile*>(file))->Get(funcname.data()));
  if(!function){
    std::cerr << std::endl << "ERROR! Failed to load function = '" << funcname << "' from input file!" << std::endl;
    assert(0);
  }
  return function;
}

std::map<std::string, const TF1*> extractTF1DMandPT(const TFile* file, const std::string& funcname, const std::vector<std::string>& uncerts = {}) {
    const TF1* func = dynamic_cast<TF1*>((const_cast<TFile*>(file))->Get(funcname.data()));
    std::map<std::string, const TF1*> funcs;
    funcs["nom"] = func;
    if (uncerts.size() > 0) {
        for (const auto& u : uncerts) {
            for (const auto& x : {"up", "down"}) {
                std::string syst_funcname;
                if (u.find("syst") != std::string::npos) {
                    syst_funcname = funcname;
                    syst_funcname.replace(syst_funcname.find("fit"), 3, u + "_" + x + "_fit");
                } else {
                    syst_funcname = funcname + "_";
                    syst_funcname += u + "_" + x;
                }
                funcs[u + "_" + x] = dynamic_cast<TF1*>((const_cast<TFile*>(file))->Get(syst_funcname.data()));
            }
        }
    }
    if (!func) {
        std::cerr << std::endl << "ERROR! Failed to load function = '" << funcname << "' from input file!" << std::endl;
        assert(0);
    }
    return funcs;
}

void TauIDSFTool::disabled() const {
  std::cerr << std::endl << "ERROR! Method has been disabled! isVsPT = "<<isVsPT<<", isVsDM = "
            << isVsDM<<", isVsEta = "<<isVsEta<< std::endl;
  assert(0);
}


TauIDSFTool::TauIDSFTool(const std::string& year, const std::string& id, const std::string& wp, const std::string& wp_vsele,  const bool dm, const bool ptdm, const bool embedding): ID(id), WP(wp), WP_VSELE(wp_vsele){

  bool verbose = false;
  //std::string datapath                = Form("%s/src/TauPOG/TauIDSFs/data",getenv("CMSSW_BASE"));
  std::string datapath                = {"data/TauIDSFs"};
  std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco","UL2016_preVFP","UL2016_postVFP","UL2017","UL2018"};
  std::vector<std::string> antiJetIDs = {"MVAoldDM2017v2","DeepTau2017v2p1VSjet"};
  std::vector<std::string> antiEleIDs = {"antiEleMVA6",   "DeepTau2017v2p1VSe"};
  std::vector<std::string> antiMuIDs  = {"antiMu3",       "DeepTau2017v2p1VSmu"};

  if(std::find(years.begin(),years.end(),year)==years.end()){
    std::cerr << std::endl << "ERROR! '"<<year<<"' is not a valid year! Please choose from ";
     std::vector<std::string>::iterator it = years.begin();
    for(it=years.begin(); it!=years.end(); it++){
      if(it!=years.begin()) std::cerr << ", ";
      std::cerr << *it;
    }
    std::cerr << std::endl;
    assert(0);
  }
  
  std::string year_mu = "";
  if(year.find("2016") != std::string::npos){
      year_mu = "2016Legacy";
  }else if(year.find("2017") != std::string::npos){
      year_mu = "2017ReReco";
  }else if(year.find("2018") != std::string::npos){
      year_mu = "2018ReReco";
  }

  if(std::find(antiJetIDs.begin(),antiJetIDs.end(),ID)!=antiJetIDs.end()){
    if (ptdm) {
      DMs    = {0,1,10};
      if (ID.find("oldDM") == std::string::npos)
      {
          DMs.push_back(11);
      }
      std::vector<std::string> allowed_wp={"Loose","Medium","Tight","VTight"};
      std::vector<std::string> allowed_wp_vsele={"VVLoose","Tight"};
      if (ID != "DeepTau2017v2p1VSjet") {
        std::cerr << "Scale factors not available for ID '"+ID+"'!" << std::endl;
        assert(0);
      }
      if (std::find(allowed_wp.begin(), allowed_wp.end(), wp) == allowed_wp.end() ||
      std::find(allowed_wp_vsele.begin(), allowed_wp_vsele.end(), wp_vsele) == allowed_wp_vsele.end()) {
      std::ostringstream msg;
      msg << "Scale factors not available for this combination of WPs! Allowed WPs for VSjet are [";
      for (const auto& x : allowed_wp) {
      msg << x << ", ";
      }
      msg << "]. Allowed WPs for VSele are [";
      for (const auto& x : allowed_wp_vsele) {
      msg << x << ", ";
      }
      msg << "]";
      std::cerr << msg.str() << std::endl;
      assert(0);
      }
      if (embedding) {
        std::cerr << "Scale factors for embedded samples not available in this format! Use either pT-binned or DM-binned SFs." << std::endl;
        assert(0);
      }
      TString filename = Form("%s/TauID_SF_dm_%s_VSjet%s_VSele%s_Mar07.root",datapath.data(),ID.data(),WP.data(),WP_VSELE.data());
      TFile* file = ensureTFile(filename, verbose);
      std::string year_ = year;
      if (year.find("UL") == 0) {
      year_ = year_.substr(2);
      }

      std::vector<std::string> uncerts_dm0={"uncert0","uncert1","syst_alleras","syst_"+year_,"syst_dm0_"+year_};
      std::vector<std::string> uncerts_dm1={"uncert0","uncert1","syst_alleras","syst_"+year_,"syst_dm1_"+year_};
      std::vector<std::string> uncerts_dm10={"uncert0","uncert1","syst_alleras","syst_"+year_,"syst_dm10_"+year_};
      std::vector<std::string> uncerts_dm11={"uncert0","uncert1","syst_alleras","syst_"+year_,"syst_dm11_"+year_};
      funcs_dm0 = extractTF1DMandPT(file,"DM0_"+year_+"_fit", uncerts_dm0);
      funcs_dm1 = extractTF1DMandPT(file,"DM1_"+year_+"_fit", uncerts_dm1);
      funcs_dm10 = extractTF1DMandPT(file,"DM10_"+year_+"_fit", uncerts_dm10);
      funcs_dm11 = extractTF1DMandPT(file,"DM11_"+year_+"_fit", uncerts_dm11);

      isVsDMandPT = true;

    } else if(dm){
      TString filename;
      if (embedding) {
          if (ID.find("oldDM") != std::string::npos)
          {
             std::cerr << "Scale factors for embedded samples are not provided for the MVA IDs." << std::endl;
             assert(0);
          }
          filename = Form("%s/TauID_SF_dm_%s_%s_EMB.root",datapath.data(),ID.data(),year.data());
      }
      else {
          filename = Form("%s/TauID_SF_dm_%s_%s.root",datapath.data(),ID.data(),year.data());
      }
      TFile* file = ensureTFile(filename,verbose);
      hist = extractTH1(file,WP);
      hist->SetDirectory(nullptr);
      file->Close();
      delete file;
      DMs    = {0,1,10};
      if (ID.find("oldDM") == std::string::npos)
      {
          DMs.push_back(11);
      }
      isVsDM = true;
    }else{
      TString filename;
      if (embedding) {
          if (ID.find("oldDM") != std::string::npos)
          {
             std::cerr << "Scale factors for embedded samples are not provided for the MVA IDs." << std::endl;
             assert(0);
          }
          filename = Form("%s/TauID_SF_pt_%s_%s_EMB.root",datapath.data(),ID.data(),year.data());
      }
      else {
          filename = Form("%s/TauID_SF_pt_%s_%s.root",datapath.data(),ID.data(),year.data());
      }
      TFile* file = ensureTFile(filename,verbose);
      func[""]     = extractTF1(file,Form("%s_cent",WP.data()));
      func["Up"]   = extractTF1(file,Form("%s_up",  WP.data()));
      func["Down"] = extractTF1(file,Form("%s_down",WP.data()));
      file->Close();
      delete file;
      isVsPT = true;
    }
  }else if(std::find(antiEleIDs.begin(),antiEleIDs.end(),ID)!=antiEleIDs.end()){
      if (embedding){
          std::cerr << "SF for ID " << ID << " not available for the embedded samples!" << std::endl;
          assert(0);
      }
      TString filename = Form("%s/TauID_SF_eta_%s_%s.root",datapath.data(),ID.data(),year.data());
      TFile* file = ensureTFile(filename,verbose);
      hist = extractTH1(file,WP);
      hist->SetDirectory(nullptr);
      file->Close();
      delete file;
      genmatches = {1,3};
      isVsEta    = true;
  }else if(std::find(antiMuIDs.begin(),antiMuIDs.end(),ID)!=antiMuIDs.end()){
      if (embedding){
          std::cerr << "SF for ID " << ID << " not available for the embedded samples!" << std::endl;
          assert(0);
      }
      TString filename = Form("%s/TauID_SF_eta_%s_%s.root",datapath.data(),ID.data(),year_mu.data());
      TFile* file = ensureTFile(filename,verbose);
      hist = extractTH1(file,WP);
      hist->SetDirectory(nullptr);
      file->Close();
      delete file;
      genmatches = {2,4};
      isVsEta    = true;
  }else{
      std::cerr << "Did not recognize tau ID '" << ID << "'!" << std::endl;
      assert(0);
  }
}



float TauIDSFTool::getSFvsPT(double pt, int genmatch, const std::string& unc){
  if(!isVsPT) disabled();
  if(genmatch==5){
    float SF = static_cast<float>(func[unc]->Eval(pt));
    return SF;
  }
  return 1.0;
}

float TauIDSFTool::getSFvsPT(double pt, const std::string& unc){
  return getSFvsPT(pt,5,unc);
}



float TauIDSFTool::getSFvsDM(double pt, int dm, int genmatch, const std::string& unc) const{
  if(!isVsDM) disabled();
  if(std::find(DMs.begin(),DMs.end(),dm)!=DMs.end() or pt<=40){
    if(genmatch==5){
      Int_t bin = hist->GetXaxis()->FindBin(dm);
      float SF  = static_cast<float>(hist->GetBinContent(bin));
      if(unc=="Up")
        SF += hist->GetBinError(bin);
      else if(unc=="Down")
        SF -= hist->GetBinError(bin);
      return SF;
    }
    return 1.0;
  }
  return 0.0;
}

float TauIDSFTool::getSFvsDM(double pt, int dm, const std::string& unc) const{
  return getSFvsDM(pt,dm,5,unc);
}

float TauIDSFTool::getSFvsDMandPT(double pt, int dm, int genmatch, const std::string& unc) const{
  if(!isVsDMandPT) disabled();
  if(std::find(DMs.begin(),DMs.end(),dm)!=DMs.end()){
    if(genmatch==5){
      // get correct functions depending on DM
      std::map<std::string, const TF1*> funcs = {};
      if (dm==0) funcs = funcs_dm0;
      if (dm>0&&dm<=2) funcs = funcs_dm1;
      if (dm==10) funcs = funcs_dm10;
      if (dm==11) funcs = funcs_dm11;


      float SF;
      if (unc.empty()) {
          SF = funcs["nom"]->Eval(std::max(std::min(pt,140.0),20.0));
      } else {
          SF = funcs[unc]->Eval(std::max(std::min(pt,140.0),20.0));
      }
      return SF;
    }
    return 1.0;
  }
  else {
      return 1.0;
  }
}

float TauIDSFTool::getSFvsDMandPT(double pt, int dm, const std::string& unc) const{
  return getSFvsDMandPT(pt,dm,5,unc);
}

float TauIDSFTool::getSFvsEta(double eta, int genmatch, const std::string& unc) const{
  if(!isVsEta) disabled();
  if(std::find(genmatches.begin(),genmatches.end(),genmatch)!=genmatches.end()){
    Int_t bin = hist->GetXaxis()->FindBin(eta);
    float SF  = static_cast<float>(hist->GetBinContent(bin));
    if(unc=="Up")
      SF += hist->GetBinError(bin);
    else if(unc=="Down")
      SF -= hist->GetBinError(bin);
    return SF;
  }
  return 1.0;
}

TauESTool::TauESTool(const std::string& year, const std::string& id): ID(id){

    bool verbose = false;
    std::string datapath                = {"data/TauIDSFs"};
    std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco","UL2016_preVFP","UL2016_postVFP","UL2017","UL2018"};

    std::string year_highpt = "";
    if(year.find("2016") != std::string::npos){
        year_highpt = "2016Legacy";
    }else if(year.find("2017") != std::string::npos){
        year_highpt = "2017ReReco";
    }else if(year.find("2018") != std::string::npos){
        year_highpt = "2018ReReco";
    }

    if(std::find(years.begin(),years.end(),year)==years.end()){
        std::cerr << std::endl << "ERROR! '"<<year<<"' is not a valid year! Please choose from ";
        std::vector<std::string>::iterator it = years.begin();
        for(it=years.begin(); it!=years.end(); it++){
            if(it!=years.begin()) std::cerr << ", ";
            std::cerr << *it;
        }
        std::cerr << std::endl;
        assert(0);
    }

    TFile* file_lowpt = ensureTFile(Form("%s/TauES_dm_%s_%s.root",datapath.data(),ID.data(),year.data()),verbose);
    TFile* file_highpt = ensureTFile(Form("%s/TauES_dm_%s_%s_ptgt100.root",datapath.data(),ID.data(),year_highpt.data()), verbose);
    hist_lowpt = extractTH1(file_lowpt, "tes");
    hist_highpt = extractTH1(file_highpt, "tes");
    hist_lowpt->SetDirectory(nullptr);
    hist_highpt->SetDirectory(nullptr);
    file_lowpt->Close();
    file_highpt->Close();
    delete file_lowpt;
    delete file_highpt;
    pt_low = 34; // average pT in Z -> tautau measurement (incl. in DM)
    pt_high = 170; // # average pT in W* -> taunu measurement (incl. in DM)
    DMs = {0,1,10};
}

float TauESTool::getTES(double pt, int dm, int genmatch, const std::string& unc){
    if(genmatch == 5 && std::find(DMs.begin(),DMs.end(),dm)!=DMs.end()){
        int bin = hist_lowpt->GetXaxis()->FindBin(dm);
        float tes = static_cast<float>(hist_lowpt->GetBinContent(bin));
        float err = 0.0;
        if(pt > pt_high){
            int bin_high = hist_highpt->GetXaxis()->FindBin(dm);
            err = hist_highpt->GetBinError(bin_high);
        }else if(pt > pt_low){
            int bin_high = hist_highpt->GetXaxis()->FindBin(dm);
            float err_high = hist_highpt->GetBinError(bin_high);
            float err_low = hist_lowpt->GetBinError(bin);
            err = err_low + (err_high-err_low)/(pt_high-pt_low)*(pt-pt_low);
        }else{
            err = hist_lowpt->GetBinError(bin);
        }
        if(unc=="Up"){
            tes += err;
        }else if(unc=="Down"){
            tes -= err;
        }
        return tes;
    }
    return 1.0;
}

TauFESTool::TauFESTool(const std::string&year, const std::string& id): ID(id){

    bool verbose = false;
    std::string datapath                = {"data/TauIDSFs"};
    std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco","UL2016_preVFP","UL2016_postVFP","UL2017","UL2018"};

    std::string year_fes = "";
    if(year.find("2016") != std::string::npos){
        year_fes = "2016Legacy";
    }else if(year.find("2017") != std::string::npos){
        year_fes = "2017ReReco";
    }else if(year.find("2018") != std::string::npos){
        year_fes = "2018ReReco";
    }

    if(std::find(years.begin(),years.end(),year)==years.end()){
    std::cerr << std::endl << "ERROR! '"<<year<<"' is not a valid year! Please choose from ";
    std::vector<std::string>::iterator it = years.begin();
    for(it=years.begin(); it!=years.end(); it++){
        if(it!=years.begin()) std::cerr << ", ";
        std::cerr << *it;
    }
    std::cerr << std::endl;
    assert(0);
    }

    TFile* file = ensureTFile(Form("%s/TauFES_eta-dm_%s_%s.root",datapath.data(),ID.data(),year_fes.data()), verbose);
    graph = (TGraphAsymmErrors*)file->Get("fes");
    DMs = {0,1};
    regions = {"barrel","endcap"};
    int i = 0;
    for(auto reg : regions){
        for(auto dm : DMs){
            float y = graph->GetY()[i];
            float yup = graph->GetErrorYhigh(i);
            float ylow = graph->GetErrorYlow(i);
            vector<float> yout = {y-ylow,y,y+yup};
            FESs.insert(std::make_pair(std::make_pair(reg,dm),yout));
            i += 1;
        }
    }
    file->Close();
    genmatches = {1,3};
}

float TauFESTool::getFES(double eta, int dm, int genmatch, const std::string& unc){
    if(std::find(DMs.begin(),DMs.end(),dm)!=DMs.end()){
        if(std::find(genmatches.begin(),genmatches.end(),genmatch)!=genmatches.end()){
            std::string region = (abs(eta)<1.5) ? "barrel" : "endcap";
            fes = FESs[std::make_pair(region,dm)];
            float fesout = 1.0;
            if(unc=="Up"){
                fesout = fes[2];
            }else if(unc=="down"){
                fesout = fes[0];
            }else if(unc!="All"){
                fesout = fes[1];
            }
            return fesout;
        }
    }
    return 1.0;
}
