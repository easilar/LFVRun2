from math import exp
import numpy as np
import ROOT
from ROOT import TF1, TCanvas, TMath, gStyle, TPaveText

dirpath = "/data1/users/itseyes/LFV/LFVRun2/nanoaodframe/genmatching_debug05010356/2018/"

def DoubleSidedCB2(x, mu, width, a1, p1, a2, p2):
    u   = float(x-mu)/width
    A1  = TMath.Power(p1/TMath.Abs(a1),p1) * exp(-a1*a1/2)
    A2  = TMath.Power(p2/TMath.Abs(a2),p2) * exp(-a2*a2/2)
    B1  = p1/TMath.Abs(a1) - TMath.Abs(a1)
    B2  = p2/TMath.Abs(a2) - TMath.Abs(a2)
    result = 1.
    if (u < -a1): result *= A1*TMath.Power(B1-u, -p1)
    elif (u < a2): result *= exp(-u*u/2)
    else: result *= A2*TMath.Power(B2+u, -p2)
    return result

def DoubleSidedCB(x, par):
    return par[0] * DoubleSidedCB2(x[0], par[1], par[2], par[3], par[4], par[5], par[6])

def combineHist(infiles, histname):
    hists = None
    for infile in infiles:
        finput = ROOT.TFile(dirpath+infile)
        hist = finput.Get(histname)
        if hist != None:
            if hists == None:
                hists = hist.Clone()
                hists.SetDirectory(0)
            else:
                hists.Add(hist)
        finput.Close()
    return hists


def voigt_fit(hist, range_min, range_max):
    fit = TF1("voigt", "[0]*TMath::Voigt(x-[1], [2], [3])", range_min, range_max)
    fit.SetParameters(10000.0, hist.GetMean(), hist.GetRMS(), 3)
    return fit

def gaus_fit(hist, range_min, range_max):
    fit = TF1("gaus", "gaus", range_min, range_max)
    fit.SetParameters(10000.0, hist.GetMean(), hist.GetRMS())
    return fit

def dcb_fit(hist, range_min, range_max):
    fit = TF1("dcb", DoubleSidedCB, range_min, range_max, 7)
    fit.SetParameter(0, 35000)
    fit.SetParameter(1, hist.GetMean())
    fit.SetParameter(2, hist.GetRMS())
    fit.SetParameter(3, 1)
    fit.SetParameter(4, 1)
    fit.SetParameter(5, 1)
    fit.SetParameter(6, 1)
    return fit

def fitting(hist, histname, fit_func="voigt", range_min=80, range_max=280, ch="", outdir="fit_results/"):
    print (hist)
    if fit_func=="voigt":    fit = voigt_fit(hist, range_min, range_max)
    elif fit_func == 'gaus': fit = gaus_fit(hist, range_min, range_max)
    elif fit_func == 'dcb':  fit = dcb_fit(hist, range_min, range_max)
    else:
        print("Invalid fit function. Please choose 'voigt', 'gaus', or 'dcb'.")
        return

    fit.SetParName(0, "Number of signal")
    fit.SetParName(1, "Mean")
    fit.SetParName(2, "Width")

    # Plot and save output histogram
    canvas = ROOT.TCanvas("canvas", histname, 800, 600)
    hist.Draw("PE1")
    gStyle.SetOptStat(0)
    hist.SetTitle(histname)
    hist.Fit(fit_func, "R")
    #gStyle.SetOptFit(111)
    #fit.Draw("same")

    # get fit parameters and create text box for plot
    fitmean = "Mean : %.2f GeV" % fit.GetParameter(1)
    fitwidth = "Width : %.2f GeV" % fit.GetParameter(2)
    chi2 = "Chi2/ndf : %.2f / %d" % (fit.GetChisquare(), fit.GetNDF())
    p_val = "P-value : %f" %fit.GetProb()

    pt1 = ROOT.TPaveText(0.65, 0.4, 0.8, 0.8, "NDC")
    pt1.SetBorderSize(0)
    pt1.SetFillColor(0)
    pt1.SetTextSize(0.04)
    if ch!="": pt1.AddText(ch[:-1])
    pt1.AddText(fit_func)
    pt1.AddText("Entries : %d" % hist.GetEntries())
    pt1.AddText(fitmean)
    pt1.AddText(fitwidth)
    pt1.AddText(chi2)
    pt1.AddText(p_val)
    pt1.Draw("same")

    #gStyle.SetOptFit(111)
    canvas.SaveAs(outdir+ch+histname+"_"+fit_func+".pdf")


TT_infiles = ["hist_TT_LFV_TCMuTau_Scalar.root", "hist_TT_LFV_TCMuTau_Vector.root", "hist_TT_LFV_TCMuTau_Tensor.root",
              "hist_TT_LFV_TUMuTau_Scalar.root", "hist_TT_LFV_TUMuTau_Vector.root", "hist_TT_LFV_TUMuTau_Tensor.root",
]
ST_infiles = ["hist_ST_LFV_TCMuTau_Scalar.root", "hist_ST_LFV_TCMuTau_Vector.root", "hist_ST_LFV_TCMuTau_Tensor.root",
              "hist_ST_LFV_TUMuTau_Scalar.root", "hist_ST_LFV_TUMuTau_Vector.root", "hist_ST_LFV_TUMuTau_Tensor.root",
]

#ST_infiles = ["hist_ST_LFV_TCMuTau_Scalar.root"]

histnames = [("h_GenSMTop_mass_S6", 100., 200.), ("h_GenSMW_mass_S6", 30., 120.),
             ("h_SMTop_mass_S6", 100., 230.), ("h_SMW_mass_S6", 40., 120.)]

for histname, range_min, range_max in histnames:
    ST_hist = combineHist(ST_infiles, histname)
    fitting(ST_hist, histname, "dcb"  , range_min, range_max, "ST_")
    fitting(ST_hist, histname, "voigt", range_min, range_max, "ST_")
    fitting(ST_hist, histname, "gaus" , range_min, range_max, "ST_")
    TT_hist = combineHist(TT_infiles, histname)
    fitting(TT_hist, histname, "voigt", range_min, range_max, "TT_")
    fitting(TT_hist, histname, "gaus" , range_min, range_max, "TT_")
    fitting(TT_hist, histname, "dcb"  , range_min, range_max, "TT_")
    print (ST_hist, TT_hist)
    hist = ST_hist.Clone()
    hist.Add(TT_hist)
    fitting(hist, histname, "voigt", range_min, range_max)
    fitting(hist, histname, "gaus" , range_min, range_max)
    fitting(hist, histname, "dcb"  , range_min, range_max)

