import os
import sys
import ROOT
from ROOT import *
import numpy as np
import argparse

base_path = './'
parser = argparse.ArgumentParser()
parser.add_argument('-L', '--label', dest='label', type=str, default="rerun_staug22")
parser.add_argument('-D', '--discriminator', dest='discriminator', type=str, default="p_st_tt_ob")
parser.add_argument('-A', '--alpha', dest='alpha', type=str, default="1p0")
args = parser.parse_args()
label = args.label
discriminator = args.discriminator
alpha = args.alpha
#labels = ['rerun_multi_Multiaug22','rerun_staug22', 'rerun_ttaug22']
# Set Runs
years = ['2016pre', '2016post', '2017','2018']

# Set output folders
for year in years:
	print("Directory created: ", base_path + label + '/' + year + '_postprocess_2/' + discriminator + '/' + alpha + '/')
	if not os.path.exists(base_path + label + '/' + year + '_postprocess_2/' + discriminator + '/' + alpha + '/'):
	    os.makedirs(base_path + label + '/' + year + '_postprocess_2/' + discriminator + '/' + alpha + '/')

# Systematic Sources => All systematics in one file.
systs_toTTTfile = ['tunedown','tuneup','hdampdown','hdampup']
systs_tofile = ['jerdown', 'jerup', 'jesAbsolute_yeardown', 'jesAbsolute_yearup', 'jesAbsolutedown', 'jesAbsoluteup', 'jesBBEC1_yeardown', 'jesBBEC1_yearup', 'jesBBEC1down', 'jesBBEC1up', 'jesFlavorQCDdown', 'jesFlavorQCDup', 'jesRelativeBaldown', 'jesRelativeBalup', 'jesRelativeSample_yeardown', 'jesRelativeSample_yearup', 'jesHEMdown','jesHEMup' , 'tesdown', 'tesup']+systs_toTTTfile

systs_pdf = ['pdfEnvup','pdfEnvdown','pdfup','pdfdown']
syst_otherTheory = ['psup','psdown','scaleup','scaledown']
syst_tauidjets = ['tauidjetHighptextrapdown', 'tauidjetHighptextrapup', 'tauidjetHighptstat_bin1down', 'tauidjetHighptstat_bin1up', 'tauidjetHighptstat_bin2down', 'tauidjetHighptstat_bin2up', 'tauidjetHighptstatdown', 'tauidjetHighptstatup', 'tauidjetHighptsystdown', 'tauidjetHighptsystup', 'tauidjetSystULyeardown', 'tauidjetSystULyearup', 'tauidjetSystallerasdown', 'tauidjetSystallerasup', 'tauidjetSystdm0ULyeardown', 'tauidjetSystdm0ULyearup', 'tauidjetSystdm10ULyeardown', 'tauidjetSystdm10ULyearup', 'tauidjetSystdm11ULyeardown', 'tauidjetSystdm11ULyearup', 'tauidjetSystdm1ULyeardown', 'tauidjetSystdm1ULyearup', 'tauidjetUncert0down', 'tauidjetUncert0up', 'tauidjetUncert1down', 'tauidjetUncert1up']
systs_toweight = ['btagcferr1down', 'btagcferr1up', 'btagcferr2down', 'btagcferr2up', 'btaghfdown', 'btaghfstats1down', 'btaghfstats1up', 'btaghfstats2down', 'btaghfstats2up', 'btaghfup', 'btaglfdown', 'btaglfstats1down', 'btaglfstats1up', 'btaglfstats2down', 'btaglfstats2up', 'btaglfup', 'muiddown', 'muidup', 'muisodown', 'muisoup', 'mutrgdown', 'mutrgup', 'prefiredown','prefireup', 'pudown', 'puup', 'tauideldown', 'tauidelup', 'tauidmudown', 'tauidmuup']+systs_pdf+syst_otherTheory



systs = systs_tofile+systs_toweight+syst_tauidjets+['']

# Produce dictionary for file lists.

nom_path = base_path + label
if not os.path.exists(nom_path):
    print("Folder '{}' does not exists.".format(nom_path))
else:
    print("Start postprocessing for '{}'.".format(nom_path))

file_list = {}
#print(os.listdir(base_path + label + '/' + year + '_postprocess/' + discriminator + '/' + alpha + '/'))
for year in years:
    file_list[year] = [i.replace('.root','') for i in os.listdir(base_path + label + '/' + year + '_postprocess/' + discriminator + '/' + alpha + '/') if '.root' in i and '__' not in i]
#print(file_list)

def collect_systhists(outfile, fname, hlists, syst, syst_, year):
    if 'SingleMuon' not in fname and syst != '':
        try:
          tmpf = TFile.Open(os.path.join(nom_path , year+"_postprocess",discriminator , alpha , fname + '__' + syst.replace("year",year[:4]) + '.root'), 'READ')
        except:
          print("No file: " + os.path.join(base_path + label, year+"_postprocess",discriminator , alpha , fname + '__' + syst.replace("year",year[:4]) + '.root'))
          return
    else: tmpf = TFile.Open(os.path.join(base_path + label, year+"_postprocess",discriminator , alpha , fname + '.root'), 'READ')
    for histname in hlists:
        if syst != '' and 'counter' in histname: continue
        if not "18" in year and "HEM" in syst: continue
        tmpf.cd()
        tmphist = tmpf.Get(histname)
        if syst == '': newtmphist = tmphist.Clone(histname)
        else: newtmphist = tmphist.Clone(histname + '__' + syst_)
        newtmphist.Scale(1.)
        outfile.cd()
        newtmphist.Write()
    tmpf.Close()
    return

# Loop over all files.
for year in years:
    out_path = base_path + label + '/' + year + '_postprocess_2/' + discriminator + '/' + alpha + '/'
    for fname in file_list[year]:
        infile = TFile.Open(os.path.join(base_path+label, year+"_postprocess",discriminator , alpha , fname+'.root'), 'READ')
        hlists = [ h.GetName() for h in infile.GetListOfKeys() if any(i in h.GetName() for i in ['dnn_pred', 'counter']) ]
        ratio = 1
        outfname = fname.replace("_" + year + "_" + fname.split("_")[-1], "")
        if "SingleMuon" in fname:
            #outfname = fname.replace("_" + fname.split("_")[-1], "")
            outfname = fname
        print("Saving histograms at {}/{}.root".format(out_path, outfname))
        # Collecting Histograms in outfile.
        outfile = TFile.Open(os.path.join(out_path, outfname + '.root'), 'RECREATE')
        # Looping over all systematics.
        #for syst, syst_ in systs.items():

        for syst in systs:
            if "year" in syst:
              syst = syst.replace('year', year[:4])
            syst_ = syst
            if '16pre' in year and syst in syst_tauidjets: syst_ = syst_.replace('2016', '2016_preVFP')
            if '16post' in year and syst in syst_tauidjets: syst_ = syst_.replace('2016', '2016_postVFP')
            collect_systhists(outfile, fname, hlists, syst, syst_, year)
        infile.Close()
        outfile.Close()
