import os, sys, glob
import ROOT
from ROOT import *
import numpy as np
import subprocess

input = sys.argv[1]
year = sys.argv[2]
if year not in ['2016pre', '2016post', '2017', '2018']:
    print('Wrong year, check again')
    sys.exit()

yield_name = 'h_ncleanjetspass'
base_path = './'
nom_path = os.path.join(base_path, input, year)
if not os.path.exists(nom_path):
    print("Folder '{}' does not exists.".format(nom_path))
    sys.exit()
else:
    print("Start postprocessing at '{}'.".format(nom_path))


# Set output folders
out_path = os.path.join(base_path, input, year + '_postprocess')
fig_path = os.path.join(base_path, input, 'figure_' + year)
if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(fig_path):
    os.makedirs(fig_path)


file_list = [i.replace('.root', '') for i in os.listdir(nom_path) if '.root' in i]
data_list = [i[:i.find('201')] for i in os.listdir(nom_path) if '.root' in i and '201' in i and 'jes' not in i]
data_list = list(set(data_list))
#print(data_list)
#print(file_list)


def get_bSFratio(inputf, inputh):
    # ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration
    # rescale histogram by Sum(event weights before applying b weight)/Sum(weights with b weight)
    # This should be done per jet bin - nojet / 3jet

    step = inputh[inputh.rfind('_S')+1:inputh.rfind('_S')+3]

    # This depends on cutflow
    if int(step[-1]) < 4: step = 'S' + step[-1]
    else                : step = 'S4'

    posthist = inputf.Get('h_nevents_' + step)
    prehist = inputf.Get('h_nevents_' + step + '_nobtag')
    if '__btag' in inputh:
        posthist = inputf.Get('h_nevents_' + step + '__' + str(inputh.split('__')[-1]))
    if prehist.Integral() * posthist.Integral() == 0:
        return 1.0
    #print(prehist.Integral() / posthist.Integral())
    return prehist.Integral(0, prehist.GetNbinsX()+1) / posthist.Integral(0, prehist.GetNbinsX()+1)


def rescale(inputh, inputf, sumW, nom_sumW): # rescale up/dn histos

    #Only for ext syst. such as tune and hdamp, not jes/jer/tes
    h = inputf.Get(inputh)
    if not any(i in inputh for i in ['event', 'counter', '_nobtag', 'LHEPdfWeightSum']):
        h.Scale(get_bSFratio(inputf, inputh))
        h.Scale(nom_sumW.GetBinContent(2) / sumW.GetBinContent(2))
        #h.Rebin(nrebin)
        #h = h.Rebin(len(rebin[h.GetName().split('_')[2]])-1, h.GetName(), array.array('d',rebin[h.GetName().split('_')[2]]))
        if yield_name in inputh:
            h1 = h.Clone('h1')
            h1.SetName(inputh.replace(yield_name, yield_name + '_yield'))
            outfile.cd()
            h1.Write()

    outfile.cd()
    h.Write()


def write_envelope(inputh, inputf, syst, nhists, sumW, new_sumW):

    if (inputh + "__" + syst + "0")  in hlists:
        var_list = []
        for x in range(0,nhists):
            h = inputf.Get(inputh + "__" + syst + str(x))
            if any(x in syst for x in ['scale', 'ps']):
                pass
            #elif 'pdf' in syst:
            #  if x == 0: continue
            #  h.Scale(sumW.GetBinContent(2) / new_sumW.GetBinContent(1))
            else: h.Scale(sumW.GetBinContent(2) / new_sumW.GetBinContent(x+1))
            #h.Rebin(nrebin)
            var_list.append(h)

        nominal = inputf.Get(inputh)
        nominal.SetDirectory(ROOT.nullptr)
        #nominal.Rebin(nrebin)
        n_bins = nominal.GetNcells()
        up = nominal.Clone()
        up.SetDirectory(ROOT.nullptr)
        up.Reset()
        dn = nominal.Clone()
        dn.SetDirectory(ROOT.nullptr)
        dn.Reset()

        for i in range(0, n_bins+2):
            minimum = float("inf")
            maximum = float("-inf")

            for v in var_list:
                c = v.GetBinContent(i)
                minimum = min(minimum, c)
                maximum = max(maximum, c)

            up.SetBinContent(i, maximum)
            dn.SetBinContent(i, minimum)

        up.Scale(get_bSFratio(bSFfile, up.GetName()))
        dn.Scale(get_bSFratio(bSFfile, dn.GetName()))
        up.SetName(inputh + "__" + syst + "up")
        dn.SetName(inputh + "__" + syst + "down")
        #We don't draw pdf in full ana due to computing resources

        up.Write()
        dn.Write()

        if yield_name in inputh:
            up_yield = up.Clone('up_yield')
            dn_yield = dn.Clone('dn_yield')
            up_yield.SetName(inputh.replace(yield_name, yield_name + '_yield') + "__" + syst + "up")
            dn_yield.SetName(inputh.replace(yield_name, yield_name + '_yield') + "__" + syst + "down")
            up_yield.Write()
            dn_yield.Write()



# Loop over all files.
for fname in file_list:
    #print(os.path.join(nom_path, fname))
    #if not any(i in fname for i in ['TTTo2L2Nu', 'TTToSemiLeptonic']): continue
    #if not any(i in fname for i in ['hdamp', 'tune']): continue

    #flag for ext. syst with different normalization
    run_on_syst = False
    if any(i in fname for i in ['hdamp', 'tune']):
        run_on_syst = True
        nom_fname = fname.replace('__' + fname.split('__')[1], '')
        nom_file = TFile.Open(os.path.join(nom_path, nom_fname + '.root'), 'READ')
        hcounter_nom = nom_file.Get("hcounter")

    infile = TFile.Open(os.path.join(nom_path, fname + '.root'), 'READ')
    hlists = [ h.GetName() for h in infile.GetListOfKeys() if '_S' in h.GetName() ]
    hlists.append("hcounter")

    # Get ratio for rescaling with b-tagSF.
    if not '__' in fname: bSFfile = infile
    elif '__' in fname and any(i in fname for i in ['hdamp', 'tune', 'jes']):
        bSFfile = infile
    else:
        bSFfname = fname.replace('__' + fname.split('__')[1], '')
        bSFfile = TFile.Open(os.path.join(nom_path, bSFfname + '.root'), 'READ')

    # Collecting Histograms in outfile.
    print("Saving histograms at {}/{}.root".format(out_path, fname))
    outfile = TFile.Open(os.path.join(out_path, fname+'.root'), 'RECREATE')

    nominal_list = []
    isScale = False
    isPS = False
    isPDF = False
    if any('__scale' in i for i in hlists): isScale = True
    if any('__ps' in i for i in hlists): isPS = True
    #if any('__pdf' in i for i in hlists) and 'LFV' in fname: isPDF = True

    for hname in hlists:
        if "__" not in hname: nominal_list.append(hname)
        if run_on_syst: continue
        h = infile.Get(hname)
        if yield_name in hname:
            h1 = h.Clone('h1')
            h1.SetName(hname.replace(yield_name, yield_name + '_yield'))
            if "__" not in h1.GetName():
                nominal_list.append(h1.GetName())
            if '201' not in fname or 'jes' in fname: h1.Scale(get_bSFratio(bSFfile, hname))
            h1.Write()
        if any(i in hname for i in ['event', 'counter', '_nobtag', 'LHEPdfWeightSum']): pass
        elif any(i in hname for i in ['__scale', '__ps', '__pdf']): continue
        elif '201' in fname and 'jes' not in fname: pass
        else:
            ratio = get_bSFratio(bSFfile, hname)
            h.Scale(ratio)
        h.Write()

    hcounter = infile.Get('hcounter')
    #LHEPdfWeightSum = infile.Get('LHEPdfWeightSum')
    nominal_list = list(set(nominal_list))

    for hname2 in nominal_list:

      if isScale: write_envelope(hname2, bSFfile, "scale", 6, hcounter, hcounter)
      if isPS: write_envelope(hname2, bSFfile, "ps", 4, hcounter, hcounter)
      #For PDF: we take 101-102 only for control plots from ttbar
      #if isPDF: write_envelope(hname2, bSFfile, "pdf", 101, hcounter, LHEPdfWeightSum) #sig: 101 / bkg: 103
      if run_on_syst: rescale(hname2, bSFfile, hcounter, hcounter_nom) #placeholder for hdamp and 8tune


    infile.Close()
    outfile.Close()


for dataname in data_list:
    try:
        subprocess.call(['rm', os.path.join(out_path, dataname + year + '.root')])
    except: pass
    subprocess.check_call( ["hadd", "-f", os.path.join(out_path, dataname + year + '.root')] + glob.glob(os.path.join(out_path, dataname) + '201*') )
