from ROOT import *
import ROOT
import sys, os
from subprocess import check_call
from os import listdir, path
import collections
import glob
import multiprocessing
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input', dest='input', type=str, default="test")
args = parser.parse_args()
input = args.input

lumi_dict = {'2016pre': 19502, '2016post': 16812, '2017': 41480, '2018':59832}#16: 36314, run2:137625
file_names = collections.OrderedDict()

if not os.path.exists(input + "/Run2"):
    try: os.makedirs(input + "/Run2")
    except: pass

def store_file(it):
    path = it[0]
    file_name = it[1]
    for f in file_name:
        print(os.path.join(path, f[:f.rfind('_')] + '.root'))
        ftmp = TFile.Open(os.path.join(path, f[:f.rfind('_')] + '.root'), 'READ')

        hist_names = [x.GetName() for x in ftmp.GetListOfKeys()]
        hist_names = list(dict.fromkeys(hist_names)) #remove duplicates from more than one instances 
        hist_names[:] = [item for item in hist_names if item not in ['hcounter']]
        hist_names.sort()

        ntmp = ftmp.Get("hcounter").GetBinContent(2)

        dest_name = path.split('/')[0] + '/Run2/' + f
        dest = TFile.Open(dest_name, 'RECREATE')

        print('Writing scaled histogram to ' + dest_name)

        for hist in hist_names:
            htmp = ftmp.Get(hist)

            if "hcounter" in hist:
                htmp.Write()
            else:
                htmp.Scale(lumi_dict[path.split('/')[1].split('_')[0]]/ntmp)
                dest.cd()
                htmp.Write()
        dest.Write()
        dest.Close()


if __name__ == '__main__':

    for era in ['2016pre', '2016post', '2017', '2018']:
        dir_path = os.path.join(input, era+'_postprocess')
        dirs = os.listdir(dir_path)
        print("POST process path: " , dir_path)
        dirs[:] = [item.replace('.root', '_' + era + '.root') for item in dirs if any(i in item for i in ['LFV']) if '__' not in item]
        print("EDITED DIRS: " , dirs)
        file_names[dir_path] = dirs

    pool = multiprocessing.Pool(20)
    pool.map(store_file, file_names.items())
    pool.close()
    pool.join()

    chs = ['ST_LFV_TCMuTau_Scalar', 'ST_LFV_TCMuTau_Tensor', 'ST_LFV_TCMuTau_Vector',
           'ST_LFV_TUMuTau_Scalar', 'ST_LFV_TUMuTau_Tensor', 'ST_LFV_TUMuTau_Vector',
           'TT_LFV_TCMuTau_Scalar', 'TT_LFV_TCMuTau_Tensor', 'TT_LFV_TCMuTau_Vector',
           'TT_LFV_TUMuTau_Scalar', 'TT_LFV_TUMuTau_Tensor', 'TT_LFV_TUMuTau_Vector']

    for ch in chs:
        print(input + '/Run2/hist_' + ch + '_201*.root')
        check_call(['hadd','-f', input + '/Run2/hist_' + ch + '.root'] +  glob.glob(input + '/Run2/hist_' + ch + '_201*.root'))
