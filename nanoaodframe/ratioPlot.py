#! /usr/bin/env python

from ROOT import *
gROOT.SetBatch(True)
import ROOT
import os, sys

from style import *

if len(sys.argv) < 2:
  print("Not enough arguements: Era")
  sys.exit()
era = sys.argv[1]

log = False
each_plot = False
printscale = False

from collections import OrderedDict
datasamples=OrderedDict()
bkgsamples=OrderedDict()
sigsamples=OrderedDict()

def SetData(fname, name, lumi):
  tmp = {}
  f = TFile(fname)
  fname = os.path.basename(fname)[:-5]
  tmp["file"] = f
  tmp["hname"] = [x.GetName() for x in f.GetListOfKeys()]
  try:
    tmp["hname"].remove("LHEPdfWeightSum")
    tmp["hname"].remove("Events")
  except: pass
  tmp["hname"].remove("hcounter")
  tmp["lumi"] = lumi 
  tmp["name"] = name
  datasamples[fname] = tmp
 
def AddBkg(fname, name, color, xsection):
  tmp = {}
  f = TFile(fname)
  fname = os.path.basename(fname)[:-5]
  nevt = 1
  tmp["file"] = f
  tmp["hname"] = [x.GetName() for x in f.GetListOfKeys()]
  h = f.Get("hcounter")
  nevt = h.GetBinContent(2)
  try:
    tmp["hname"].remove("LHEPdfWeightSum")
    tmp["hname"].remove("Events")
  except: pass
  tmp["hname"].remove("hcounter")
  tmp["total"] = nevt
  tmpcol = TColor.GetColor(color)
  tmp["col"] = tmpcol
  tmp["xsection"] = xsection
  tmp["name"] = name
  bkgsamples[fname] = tmp

def AddSig(fname, name, color, xsection):
  tmp = {}
  f = TFile(fname)
  fname = os.path.basename(fname)[:-5]
  nevt = 1
  tmp["file"] = f
  tmp["hname"] = [x.GetName() for x in f.GetListOfKeys()]
  h = f.Get("hcounter")
  nevt = h.GetBinContent(2)
  try:
    tmp["hname"].remove("LHEPdfWeightSum")
    tmp["hname"].remove("Events")
  except: pass
  tmp["hname"].remove("hcounter")
  tmp["total"] = nevt
  tmpcol = TColor.GetColor(color)
  tmp["col"] = tmpcol
  tmp["xsection"] = xsection
  tmp["name"] = name
  sigsamples[fname] = tmp

####Users should provide these information
AddBkg("hist_TTToSemiLeptonic.root", "TTLJ", '#990000', 365.34)
AddBkg("hist_TTTo2L2Nu.root", "TTLL", '#330000', 88.29)
AddBkg("hist_TTToHadronic.root", "TTHad", '#cc0000', 377.96)
AddBkg("hist_ttHTobb.root", "ttV/H", '#ff66ff', 0.2934)
AddBkg("hist_ttHToNonbb.root", "ttV/H", '#ff66ff', 0.2151)
AddBkg("hist_TTWJetsToLNu.root", "ttV/H", '#ff66ff', 0.2043)
AddBkg("hist_TTWJetsToQQ.root", "ttV/H", '#ff66ff', 0.4062)
AddBkg("hist_TTZToLLNuNu.root", "ttV/H", '#ff66ff', 0.2529)
AddBkg("hist_TTZToQQ.root", "ttV/H", '#ff66ff', 0.5297)
AddBkg("hist_ST_s_4f_lepton.root", "ST", '#990099', 3.36)
AddBkg("hist_ST_t_top_4f.root", "ST", '#990099', 136.02)
AddBkg("hist_ST_t_antitop_4f.root", "ST", '#990099', 80.95)
AddBkg("hist_ST_tW_top_5f.root", "ST", '#990099', 35.85)
AddBkg("hist_ST_tW_antitop_5f.root", "ST", '#990099', 35.85)
AddBkg("hist_WW.root", "VV", '#00cccc', 118.7)
AddBkg("hist_WZ.root", "VV", '#00cccc', 47.13)
AddBkg("hist_ZZ.root", "VV", '#00cccc', 16.523)
AddBkg("hist_DYJetsToLL_M10to50.root", "ZJets", '#000099', 18610.0)
AddBkg("hist_DYJetsToLL_M50_amc.root", "ZJets", '#000099', 6077.22)
AddBkg("hist_WJetsToLNu_HT0To100.root", "WJets", '#ff9933', 59373.3)
AddBkg("hist_WJetsToLNu_HT100To200.root", "WJets", '#ff9933', 1335.585)
AddBkg("hist_WJetsToLNu_HT200To400.root", "WJets", '#ff9933', 360.4194) 
AddBkg("hist_WJetsToLNu_HT400To600.root", "WJets", '#ff9933', 49.35019)
AddBkg("hist_WJetsToLNu_HT600To800.root", "WJets", '#ff9933', 13.4960)
AddBkg("hist_WJetsToLNu_HT800To1200.root", "WJets", '#ff9933', 6.612202)
AddBkg("hist_WJetsToLNu_HT1200To2500.root", "WJets", '#ff9933', 1.770228)
AddBkg("hist_WJetsToLNu_HT2500ToInf.root", "WJets", '#ff9933', 0.135072)
#AddBkg("hist_QCD_Pt15To20_MuEnriched.root","QCD", '#d0cfd4', 3819570)
AddBkg("hist_QCD_Pt20To30_MuEnriched.root", "QCD", '#d0cfd4', 2960198.4)
AddBkg("hist_QCD_Pt30To50_MuEnriched.root", "QCD", '#d0cfd4', 1652471.5)
AddBkg("hist_QCD_Pt50To80_MuEnriched.root", "QCD", '#d0cfd4', 437504.1)
AddBkg("hist_QCD_Pt80To120_MuEnriched.root", "QCD", '#d0cfd4', 106033.7)
AddBkg("hist_QCD_Pt120To170_MuEnriched.root", "QCD", '#d0cfd4', 25190.5)
AddBkg("hist_QCD_Pt170To300_MuEnriched.root", "QCD", '#d0cfd4', 8654.5)
AddBkg("hist_QCD_Pt300To470_MuEnriched.root", "QCD", '#d0cfd4', 797.4)
AddBkg("hist_QCD_Pt470To600_MuEnriched.root", "QCD", '#d0cfd4', 79)
AddBkg("hist_QCD_Pt600To800_MuEnriched.root", "QCD", '#d0cfd4', 25.1)
AddBkg("hist_QCD_Pt800To1000_MuEnriched.root", "QCD", '#d0cfd4', 4.7)
AddBkg("hist_QCD_Pt1000_MuEnriched.root", "QCD", '#d0cfd4', 1.6) 

if era == "2016pre":
  SetData("hist_SingleMuon2016pre.root","data", 19502)
  AddSig("hist_ST_LFV_TCMuTau_Scalar.root", "STLFVc_s", '#99CC66', 0.00740)
  AddSig("hist_ST_LFV_TCMuTau_Vector.root", "STLFVc_v", '#99CC66', 0.0368)
  AddSig("hist_ST_LFV_TCMuTau_Tensor.root", "STLFVc_t", '#99CC66', 0.1784)
  AddSig("hist_ST_LFV_TUMuTau_Scalar.root", "STLFVu_s", '#6B8551', 0.0838)
  AddSig("hist_ST_LFV_TUMuTau_Vector.root", "STLFVu_v", '#6B8551', 0.393)
  AddSig("hist_ST_LFV_TUMuTau_Tensor.root", "STLFVu_t", '#6B8551', 1.796)
  AddSig("hist_TT_LFV_TCMuTau_Scalar.root", "TTLFVc_s", '#908DCC', 0.00269)
  AddSig("hist_TT_LFV_TCMuTau_Vector.root", "TTLFVc_v", '#908DCC', 0.0215)
  AddSig("hist_TT_LFV_TCMuTau_Tensor.root", "TTLFVc_t", '#908DCC', 0.1290)
  AddSig("hist_TT_LFV_TUMuTau_Scalar.root", "TTLFVu_s", '#4F4D80', 0.00269)
  AddSig("hist_TT_LFV_TUMuTau_Vector.root", "TTLFVu_v", '#4F4D80', 0.0215)
  AddSig("hist_TT_LFV_TUMuTau_Tensor.root", "TTLFVu_t", '#4F4D80', 0.1290)

elif era == "2016post":
  SetData("hist_SingleMuon2016post.root","data", 16812)
  AddSig("hist_ST_LFV_TCMuTau_Scalar.root", "STLFVc_s", '#99CC66', 0.00740)
  AddSig("hist_ST_LFV_TCMuTau_Vector.root", "STLFVc_v", '#99CC66', 0.0368)
  AddSig("hist_ST_LFV_TCMuTau_Tensor.root", "STLFVc_t", '#99CC66', 0.1784)
  AddSig("hist_ST_LFV_TUMuTau_Scalar.root", "STLFVu_s", '#6B8551', 0.0838)
  AddSig("hist_ST_LFV_TUMuTau_Vector.root", "STLFVu_v", '#6B8551', 0.393)
  AddSig("hist_ST_LFV_TUMuTau_Tensor.root", "STLFVu_t", '#6B8551', 1.796)
  AddSig("hist_TT_LFV_TCMuTau_Scalar.root", "TTLFVc_s", '#908DCC', 0.00269)
  AddSig("hist_TT_LFV_TCMuTau_Vector.root", "TTLFVc_v", '#908DCC', 0.0215)
  AddSig("hist_TT_LFV_TCMuTau_Tensor.root", "TTLFVc_t", '#908DCC', 0.1290)
  AddSig("hist_TT_LFV_TUMuTau_Scalar.root", "TTLFVu_s", '#4F4D80', 0.00269)
  AddSig("hist_TT_LFV_TUMuTau_Vector.root", "TTLFVu_v", '#4F4D80', 0.0215)
  AddSig("hist_TT_LFV_TUMuTau_Tensor.root", "TTLFVu_t", '#4F4D80', 0.1290)

elif era == "2017":
  SetData("hist_SingleMuon2017.root","data", 41480)
  AddSig("hist_ST_LFV_TCMuTau_Scalar.root", "STLFVc_s", '#99CC66', 0.00740)
  AddSig("hist_ST_LFV_TCMuTau_Vector.root", "STLFVc_v", '#99CC66', 0.0368)
  AddSig("hist_ST_LFV_TCMuTau_Tensor.root", "STLFVc_t", '#99CC66', 0.1784)
  AddSig("hist_ST_LFV_TUMuTau_Scalar.root", "STLFVu_s", '#6B8551', 0.0838)
  AddSig("hist_ST_LFV_TUMuTau_Vector.root", "STLFVu_v", '#6B8551', 0.393)
  AddSig("hist_ST_LFV_TUMuTau_Tensor.root", "STLFVu_t", '#6B8551', 1.796)
  AddSig("hist_TT_LFV_TCMuTau_Scalar.root", "TTLFVc_s", '#908DCC', 0.00269)
  AddSig("hist_TT_LFV_TCMuTau_Vector.root", "TTLFVc_v", '#908DCC', 0.0215)
  AddSig("hist_TT_LFV_TCMuTau_Tensor.root", "TTLFVc_t", '#908DCC', 0.1290)
  AddSig("hist_TT_LFV_TUMuTau_Scalar.root", "TTLFVu_s", '#4F4D80', 0.00269)
  AddSig("hist_TT_LFV_TUMuTau_Vector.root", "TTLFVu_v", '#4F4D80', 0.0215)
  AddSig("hist_TT_LFV_TUMuTau_Tensor.root", "TTLFVu_t", '#4F4D80', 0.1290)


elif era == "2018":
  SetData("hist_SingleMuon2018.root", "data", 59832)
  AddSig("hist_ST_LFV_TCMuTau_Scalar.root", "STLFVc_s", '#99CC66', 0.00740)
  AddSig("hist_ST_LFV_TCMuTau_Vector.root", "STLFVc_v", '#99CC66', 0.0368)
  AddSig("hist_ST_LFV_TCMuTau_Tensor.root", "STLFVc_t", '#99CC66', 0.1784)
  AddSig("hist_ST_LFV_TUMuTau_Scalar.root", "STLFVu_s", '#6B8551', 0.0838)
  AddSig("hist_ST_LFV_TUMuTau_Vector.root", "STLFVu_v", '#6B8551', 0.393)
  AddSig("hist_ST_LFV_TUMuTau_Tensor.root", "STLFVu_t", '#6B8551', 1.796)
  AddSig("hist_TT_LFV_TCMuTau_Scalar.root", "TTLFVc_s", '#908DCC', 0.00269)
  AddSig("hist_TT_LFV_TCMuTau_Vector.root", "TTLFVc_v", '#908DCC', 0.0215)
  AddSig("hist_TT_LFV_TCMuTau_Tensor.root", "TTLFVc_t", '#908DCC', 0.1290)
  AddSig("hist_TT_LFV_TUMuTau_Scalar.root", "TTLFVu_s", '#4F4D80', 0.00269)
  AddSig("hist_TT_LFV_TUMuTau_Vector.root", "TTLFVu_v", '#4F4D80', 0.0215)
  AddSig("hist_TT_LFV_TUMuTau_Tensor.root", "TTLFVu_t", '#4F4D80', 0.1290)

else:
  print("Wrong era!")
  sys.exit()

N_hist = len(datasamples[list(datasamples.keys())[0]]["hname"])
N_bkgsamples = len(bkgsamples)
N_sigsamples = len(sigsamples)
if any('DNN' in x for x in datasamples[list(datasamples.keys())[0]]["hname"]):
  rebin = 4
else:
  rebin = 1

fNevt = open("Nevt_ratio.txt",'w')
for fname in bkgsamples.keys():
  fNevt.write(fname + " : generated events : " + str(int(bkgsamples[fname]["total"])) + "\n")
for fname in sigsamples.keys():
  fNevt.write(fname + " : generated events : " + str(int (sigsamples[fname]["total"])) + "\n")

count = 0
for i in range(0, N_hist):
  string_fname = ''
  string_nevt =  ''
  hnames = datasamples[list(datasamples.keys())[0]]["hname"][i]
  print(hnames)

  printHistName = "ncleanbjetspass"
  if len(hnames) < 2: continue
  if any(i in hnames for i in [printHistName, 'DNN']):
    print(hnames)
    string_fname += "%s \n" %hnames
    string_nevt += "%s \n" %hnames


  hs = THStack()
  #l = TLegend(0.30, 0.99 - 0.8 * N_bkgsamples / 20., 0.89, 0.85)
  l = TLegend(0.15,0.71,0.89,0.87)
  l.SetNColumns(4);
  l.SetTextSize(0.05);
  l.SetLineColor(0);
  l.SetFillColor(0);

  h_data = datasamples[list(datasamples.keys())[0]]["file"].Get(datasamples[list(datasamples.keys())[0]]["hname"][i])
  nbins = h_data.GetNbinsX()
  h_data.AddBinContent( nbins, h_data.GetBinContent( nbins+1 ) )  #overflow
  h_data.Rebin(rebin)

  h_sub = h_data.Clone("h_sub")

  ntotalbkg = 0
  k = 0
  for fname in bkgsamples.keys():
    h_tmp = bkgsamples[fname]["file"].Get(datasamples[list(datasamples.keys())[0]]["hname"][i])
    nbins = h_tmp.GetNbinsX()
    h_tmp.AddBinContent( nbins, h_tmp.GetBinContent( nbins+1 ) ) #overflow
    h_tmp.SetFillColor(bkgsamples[fname]["col"])
    ## normalization
    scale = datasamples[list(datasamples.keys())[0]]["lumi"]/(bkgsamples[fname]["total"]/bkgsamples[fname]["xsection"])

    h_tmp.Scale(scale)
    h_tmp.Rebin(rebin)

    ## check if the sample is the same as previous process. 
    if k < N_bkgsamples-1 :
      post_name = list(bkgsamples.keys())[k+1]
      if bkgsamples[fname]["name"] is bkgsamples[post_name]["name"]:
        h_tmp.SetLineColor(bkgsamples[fname]["col"]) 
      else:
        l.AddEntry(h_tmp, bkgsamples[fname]["name"]  ,"F") 
    else: 
      l.AddEntry(h_tmp, bkgsamples[fname]["name"]  ,"F")

    ## print out number of events
    numevt = h_tmp.Integral()
    rawevt = h_tmp.GetEntries()
    ntotalbkg = ntotalbkg + numevt
    if any(i in hnames for i in [printHistName, 'DNN']):
      string_nevt += "%f \n"%(numevt)
      string_fname += "%s :  %s = %f \n"%(fname,bkgsamples[fname]["name"],numevt)
      print(fname, " : ", bkgsamples[fname]["name"], " = ", "{0:.10g}".format(numevt), " scale : " ,"{0:.3g}".format(scale))
    if printscale: print(scale)

    ## Add to Stack
    hs.Add( h_tmp ) #hh_tmp -> add h tmp sig, hs->other
    k = k+1

  h_bkg = hs.GetStack().Last()

  #Sig
  hs_sig = []
  n = 0

  for fname in sigsamples.keys():
    h_sig = sigsamples[fname]["file"].Get(datasamples[list(datasamples.keys())[0]]["hname"][i])
    nbins = h_sig.GetNbinsX()
    h_sig.AddBinContent( nbins, h_sig.GetBinContent( nbins+1 ) ) 
    h_sig.SetLineColor(sigsamples[fname]["col"])

    ## normalization
    scale = datasamples[list(datasamples.keys())[0]]["lumi"]/(sigsamples[fname]["total"]/sigsamples[fname]["xsection"])
    h_sig.Scale(scale)
    h_sig.Rebin(rebin)

    ## check if the sample is the same as previous process. 
    if n < N_sigsamples-1 :
      post_name = list(sigsamples.keys())[n+1]
      if sigsamples[fname]["name"] is sigsamples[post_name]["name"]:
        h_sig.SetLineColor(sigsamples[fname]["col"])
      else:
        l.AddEntry(h_sig, sigsamples[fname]["name"]  ,"F")
    else:
      l.AddEntry(h_sig, sigsamples[fname]["name"]  ,"F")

    ## print out number of events
    numevt = h_sig.Integral()
    rawevt = h_sig.GetEntries()
    if any(i in hnames for i in [printHistName, 'DNN']):
      string_nevt += "%f \n"%(numevt)
      string_fname += "%s :  %s = %f \n"%(fname,sigsamples[fname]["name"],numevt)
      print(fname, " : ", sigsamples[fname]["name"], " = ", "{0:.10g}".format(numevt), " scale : " ,"{0:.3g}".format(scale))
    if printscale: print(scale)

    if h_data.Integral() > 0 and h_sig.Integral() > 0 and h_bkg.Integral() != 0: h_sig.Scale(h_data.Integral()/h_sig.Integral())

    hs_sig.append(h_sig)
    n = n+1



  ndata= h_data.Integral()
  nsub = ndata-ntotalbkg
  if any(i in hnames for i in [printHistName, 'DNN']):
    string_nevt += "%f \n" % ntotalbkg
    string_nevt += "%d \n" % ndata
    string_nevt += "%f \n" % nsub
    string_fname += "ntotal = %f \n" % ntotalbkg
    string_fname += "ndata = %d \n" % ndata
    string_fname += "nsub = %f \n" % nsub
    fNevt.write(string_nevt)
    fNevt.write(string_fname)
    print("ntotal = " , "{0:.10g}".format(ntotalbkg))
    print("ndata = " , "{0:.0f}".format(ndata))
    print("nsub = ", "{0:.10g}".format(nsub))

  def createCanvasPads():
    #creat canvas
    c = TCanvas("c_"+"{}".format(i),"c", 450, 450)

    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0.0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.Draw()
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0.0, 0.0, 1, 0.28)
    pad2.SetBottomMargin(0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()

    return c, pad1, pad2

  def createRatio(mc, data):
    h4 = mc.Clone("h4")
    h3 = data.Clone("h3")
    h3.SetDirectory(0)
    h3.SetLineColor(1)
    h3.SetMarkerStyle(20)
    h3.SetMarkerSize(0.5)
    h3.SetTitle("")
    h3.SetMinimum(0.6)
    h3.SetMaximum(1.4)
    #h3.Sumw2()
    h3.SetStats(0)

    y = h3.GetYaxis()
    y.SetTitle("Data/MC")
    y.SetNdivisions(505)
    y.SetTitleSize(0.11)
    y.SetTitleOffset(0.35)
    y.SetLabelSize(0.1)

    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitleSize(0.11)
    x.SetTitleOffset(1.0)
    x.SetLabelSize(0.1)
    #x.SetLabelSize(0)
    h3.Divide(h4)

    return h3

  def ratioplot():
    h3 = createRatio(h_bkg, h_data)
    c, pad1, pad2 = createCanvasPads()

    #Draw each plot
    pad1.cd()
    if log:
      pad1.SetLogy()
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(0.5)
    max_data = h_data.GetMaximum()
    max_hs = hs.GetMaximum()
    sigmaximums = [i.GetMaximum() for i in hs_sig]
    max_sig = max(sigmaximums)
    #max_sig = max_data
    maxfrac = 0.5
    if(max_hs > max_sig):
      if log :
        if max_data > 100000:
          maxfrac = 1000
        else:
          maxfrac = 100
      if max_hs > max_data :
        h_data.SetMaximum(max_hs+max_hs*maxfrac)
      else:
        h_data.SetMaximum(max_data+max_data*maxfrac)
    else:
      if log:
        maxfrac = 100
        h_data.SetMaximum(max_sig+max_sig*maxfrac)
      else:
        h_data.SetMaximum(max_sig*1.5)
    if log: h_data.SetMinimum(0.5)
    h_data.Draw("p")
    h_data.SetTitle("")
    h_data.GetYaxis().SetTitle("Events")
    h_data.GetYaxis().SetTitleOffset(1.2)
    h_data.GetYaxis().SetTitleSize(0.045)
    h_data.GetXaxis().SetLabelSize(0)#
    #h_data.GetXaxis().SetTitle("")
    h_data.GetXaxis().SetTitleOffset(5.0)
    hs.Draw("histsame")
    for hi in hs_sig:
      hi.SetLineWidth(2)
      hi.Draw("hist same")
    h_data.Draw("psame")
    h_data.Draw("AXIS P SAME")

    l.AddEntry(h_data,"Data","P")
    l.Draw()
    label = TPaveText()
    label.SetX1NDC(gStyle.GetPadLeftMargin())
    label.SetY1NDC(1.0-gStyle.GetPadTopMargin())
    label.SetX2NDC(1.0-gStyle.GetPadRightMargin()+0.03)
    label.SetY2NDC(1.0)
    label.SetTextFont(62)
    if era == '2017':
      label.AddText("Work in Progress        CMS, 41.5 fb^{-1} at #sqrt{s} = 13 TeV")
    elif era == '2018':
      label.AddText("Work in Progress        CMS, 59.7 fb^{-1} at #sqrt{s} = 13 TeV")
    label.SetFillStyle(0)
    label.SetBorderSize(0)
    label.SetTextSize(0.05)
    label.SetTextAlign(32)
    label.Draw("same")
    
    pad2.cd()
    h3.Draw("ep")

    logname = ""
    if log:
      logname = "_log"

    if each_plot:
      c.Print(datasamples[datasamples.keys()]["hname"][i]+logname+".pdf")

    filename = "result_ratio"+logname+".pdf"
    if i == 0 and N_hist > 1:
      c.Print( (filename+"(") )
    elif i > 0 and i == N_hist-1:
      c.Print( (filename+")") ) 
    else:
      c.Print(filename)
    
    c.Clear()

  ratioplot()
  count += 1

fNevt.close()
