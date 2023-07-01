import os, sys, glob
import ROOT
from ROOT import *
import numpy as np
import subprocess
import argparse


base_path = './'
parser = argparse.ArgumentParser()
parser.add_argument('-L', '--label', dest='label', type=str, default="top_lfv_multiClass_jan03")
parser.add_argument('-D', '--discriminator', dest='discriminator', type=str, default="p_st_tt_ob")
parser.add_argument('-A', '--alpha', dest='alpha', type=str, default="0p1")
parser.add_argument('-Y', '--year', dest='year', type=str, default="2017")
args = parser.parse_args()
label = args.label
discriminator = args.discriminator
alpha = args.alpha
year = args.year


if year not in ['2016pre', '2016post', '2017', '2018']:
    print('Wrong year, check again')
    sys.exit()

yield_name = 'h_ncleanjetspass'
base_path = './'
nom_path = base_path+label

if not os.path.exists(nom_path):
    print("Folder '{}' does not exists.".format(nom_path))
    sys.exit()
else:
    print("Start postprocessing at '{}'.".format(nom_path))


out_path = nom_path + '/' + year + '_postprocess/' + discriminator + '/' + alpha + '/'
if not os.path.exists(base_path + label + '/' + year + '_postprocess/' + discriminator + '/' + alpha + '/'):
    os.makedirs(base_path + label + '/' + year + '_postprocess/' + discriminator + '/' + alpha + '/')


# Set output folders

file_list = [i.replace('.root', '') for i in os.listdir(nom_path+'/'+year+'/preds/'+discriminator+'/'+alpha +'/') if '.root' in i]
data_list = [i[:i.find('201')] for i in os.listdir(nom_path+'/'+year+'/preds/'+discriminator+'/'+alpha +'/') if '.root' in i and '201' in i and 'jes' not in i]
data_list = list(set(data_list))
#print(data_list)
#print(file_list)


def get_bSFratio(inputf, inputh):
    # ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration
    # rescale histogram by Sum(event weights before applying b weight)/Sum(weights with b weight)
    # This should be done per jet bin - nojet / 3jet

    step = 'S4'

    posthist = inputf.Get('h_nevents_' + step)
    prehist = inputf.Get('h_nevents_' + step + '_nobtag')
    if '__btag' in inputh:
        posthist = inputf.Get('h_nevents_' + step + '__' + str(inputh.split('__')[-1]))
    if prehist.Integral() * posthist.Integral() == 0:
        return 1.0
    #print(prehist.Integral() / posthist.Integral())
    return prehist.Integral(0, prehist.GetNbinsX()+1) / posthist.Integral(0, prehist.GetNbinsX()+1)


def write_envelope(inputh, inputf, syst, nhists, sumW, new_sumW):
  print("inside Envelope")
  if (inputh + "__" + syst + "0")  in hlists:
    var_list = []
    for x in range(0,nhists):
      h = inputf.Get(inputh + "__" + syst + str(x))
      if any(x in syst for x in ['scale', 'ps']):
        pass
      elif 'pdf' in syst:
        if x == 0: continue
        h.Scale(sumW.GetBinContent(2) / new_sumW.GetBinContent(1))
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
    infile = TFile.Open(os.path.join(nom_path+'/'+year+'/preds/'+discriminator+'/'+alpha +'/', fname + '.root'), 'READ')
    hlists = [ h.GetName() for h in infile.GetListOfKeys() if '_S' in h.GetName() or 'dnn' in h.GetName()]
    hlists.append("hcounter")
    print("HLIST", hlists)

    # Get ratio for rescaling with b-tagSF.
    if not '__' in fname: bSFfile = infile
    elif '__' in fname and any(i in fname for i in ['hdamp', 'tune', 'jes']):
        bSFfile = infile
    else:
        bSFfname = fname.replace('__' + fname.split('__')[1], '')
        bSFfile = TFile.Open(os.path.join(nom_path+'/'+year+'/preds/'+discriminator+'/'+alpha +'/', bSFfname + '.root'), 'READ')

    # Collecting Histograms in outfile.
    print("Saving histograms at {}/{}.root".format(out_path, fname))
    outfile = TFile.Open(os.path.join(out_path, fname+'.root'), 'RECREATE')

    nominal_list = []
    isScale = False
    isPS = False
    isPDF = False
    #if ('__scale' in fname): isScale = True
    #if ('__ps' in fname): isPS = True
    #if ('__pdf' in fname): isPDF = True

    #print(isScale, isPS)

    for hname in hlists:
        print("HISTO : ", hname)
        if "__" not in hname: nominal_list.append(hname)
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
    nominal_list = list(set(nominal_list))

    for hname2 in nominal_list:
      print("HNAME2:" , hname2, isScale, isPS,isPDF)
      if isScale: write_envelope(hname2, bSFfile, "scale", 6, hcounter, hcounter)
      if isPS: write_envelope(hname2, bSFfile, "ps", 4, hcounter, hcounter)
      if isPDF: write_envelope(hname2, bSFfile, "pdf", 103, hcounter, hcounter)
      #  if 'STTH' in files: write_envelope(hname2, hcounter, "pdf", 30, LHEPdfWeightSum)
      #  else:               write_envelope(hname2, hcounter, "pdf", 103, LHEPdfWeightSum)
      #if run_on_syst: rescale([], nom_EventInfo) #placeholder for hdamp and py8tune


    infile.Close()
    outfile.Close()


for dataname in data_list:
    try:
        subprocess.call(['rm', os.path.join(out_path, dataname + year + '.root')])
    except: pass
    subprocess.check_call( ["hadd", "-f", os.path.join(out_path, dataname + year + '.root')] + glob.glob(os.path.join(out_path, dataname) + '201*') )
