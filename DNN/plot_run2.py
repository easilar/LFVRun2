import sys, os
import re
from subprocess import call
import argparse

#labels = ['rerun_multi_Multiaug22','rerun_staug22','rerun_ttaug22']
parser = argparse.ArgumentParser()
parser.add_argument('-L', '--label', dest='label', type=str, default="rerun_staug22")
parser.add_argument('-D', '--discriminator', dest='discriminator', type=str, default="")
parser.add_argument('-A', '--alpha', dest='alpha', type=str, default="")
args = parser.parse_args()
label = args.label
discriminator = args.discriminator
alpha = args.alpha

config_path = '../plotIt/configs/'
common_syst = 'systematics:\n'
groups = ['GData', 'Gttll', 'Gttlj', 'Gttjj', 'GttV', 'GZJets', 'GWJets', 'GSingleT', 'GVV', 'GQCD',
          'GLFVSTcv', 'GLFVSTuv', 'GLFVTTcv', 'GLFVTTuv',]
#not include prefire and elzvtx which exist only in 2017
#common_syst_list = ['pu', 'muid', 'muiso', 'mutrg', 'elid', 'elreco', 'eltrg',
syst = ['jesAbsolute_year', 'jesBBEC1_year', 'jesRelativeSample_year'] #jesEC2year missing
systs_tofile = ['jer',  'jesAbsolute', 'jesBBEC1', 'jesFlavorQCD', 'jesRelativeBal', 'tes']
systs_toweight = ['btagcferr1', 'btagcferr2', 'btaghf', 'btaghfstats1', 'btaghfstats2', 'btaglf', 'btaglfstats1', 'btaglfstats2',   'muid', 'muiso', 'mutrg', 'prefire','pu', 'tauidel', 'tauidmu']
common_syst_list = systs_tofile+systs_toweight

#common_syst_list = []
#syst = []
years = {'2016pre': 19502, '2016post': 16812, '2017': 41480, '2018':59832}

for sy in syst:
  if 'year' in sy:
    for year in years.keys():
      common_syst_list.append(sy.replace('year', year[:4]))
  else: common_syst_list.append(sy)

reco_str = label

dest_path = reco_str + "/" + discriminator + "/"+alpha+"/"
print(dest_path)
if not os.path.exists(os.path.join(dest_path, 'figures')):
  try: os.makedirs(os.path.join(dest_path, 'figures'))
  except: pass
if not os.path.exists(os.path.join(dest_path, 'figures/qcd')):
  try: os.makedirs(os.path.join(dest_path, 'figures/qcd'))
  except: pass

for item in common_syst_list:
  common_syst += '  - ' + item + '\n'

string_for_files = ''
print("year loop starts")
for year, lumi in years.items():
  #Firstly, merge file list + scale
  with open(config_path + 'files_' + year + '.yml') as f:
    print("yml file opened")
    lines = f.readlines()
    skip_signal = False
    for line in lines:
      #if '#' in line[0]: line = line[1:]
      if skip_signal and 'hist' in line: skip_signal = False
      if 'LFV' in line: skip_signal = True
      if 'hist_QCD' in line: skip_signal = True
      if '#' in line[0]: skip_signal = True
      if 'hist' in line:
        line = line[0] + reco_str + '/' + year + '_postprocess_2/'+ '/' + discriminator + '/' + alpha +'/' + line[1:]
        if not any(i in line for i in ['LFV', 'SingleMuon2']):
          line += '  scale: ' + str(int(lumi)/137570.0) + '\n'
      if not skip_signal: # and not any(i in line for i in ['yields-group']):
        #if 'group' in line and not any(i in line for i in groups): string_for_files += '  group: Gother \n'
        #else: string_for_files += line
        string_for_files += line


  file_syst = ''
  with open(config_path + 'config_' + year + '.yml') as f:
    lines = f.readlines()
    for line in lines:
      if 'type' in line:
        if 'const' in line:
          #file_syst += line[:line.find(':')] + '_' + year + line[line.find(':'):line.find('hist')] + '/' + year + '_postprocess_2/'+'/' + discriminator + '/' + alpha + '/' + line[line.find('hist'):]
          if year == '2016pre':
            file_syst += line
        elif 'shape' in line:
          #file_syst += line[:line.find('hist')] + '/' + year + '_postprocess_2/'+'/' + discriminator + '/' +alpha + '/' + line[line.find('hist'):]
          if year == '2016pre':
            file_syst += line
  print("file_syst", file_syst)

with open(config_path + 'files_Run2.yml', 'w+') as fnew:
  print("""
'{0}/Run2/hist_ST_LFV_TCMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVSTcv'
  cross-section: 0.0368
  generated-events: 12645000
  scale: 100
  group: GLFVSTcv
  order: 1

'{0}/Run2/hist_ST_LFV_TUMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.393
  generated-events: 12785962
  scale: 100
  group: GLFVSTuv
  order: 2

'{0}/Run2/hist_TT_LFV_TCMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVTTcv'
  cross-section: 0.0215
  generated-events: 9725000
  scale: 100
  group: GLFVTTcv
  order: 3

'{0}/Run2/hist_TT_LFV_TUMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.0215
  generated-events: 11286000
  scale: 100
  group: GLFVTTuv
  order: 4
  """.format(dest_path), file=fnew)
  fnew.write(string_for_files)

with open(config_path + 'template_Run2.yml') as f:
  lines = f.readlines()
  with open(config_path + 'config_Run2.yml', 'w+') as f1:
    for line in lines: f1.write(line)
    f1.write(common_syst)
    f1.write(file_syst)
    f1.write("\nplots:\n  include: ['histos_dnn.yml']\n")

print("before call")
print(['../plotIt/plotIt', '-o ' + dest_path + '/figures', config_path + 'config_Run2.yml'])
call(['../plotIt/plotIt', '-o ' + dest_path + '/figures', config_path + 'config_Run2.yml'], shell=False)
print("after call")


#For QCD
string_for_qcd = ''
for year, lumi in years.items():
  with open(config_path + 'files_' + year + '.yml') as f:
    lines = f.readlines()
    skip_signal = True
    for line in lines:
      #if '#' in line[0]: line = line[1:]
      if '#' in line[0]: skip_signal = True
      if skip_signal and 'hist_QCD' in line: skip_signal = False
      if 'hist_QCD' in line:
        line = line[0] + reco_str + '/' + year + '_postprocess_2/' + discriminator + '/' + alpha +'/' + line[1:]
        if not any(i in line for i in ['LFV', 'SingleMuon2']):
          line += '  scale: ' + str(int(lumi)/137570.0) + '\n'
      #if not skip_signal and not any(i in line for i in ['yields-group']): string_for_qcd += line
      if not skip_signal: string_for_qcd += line

with open(config_path + 'files_Run2.yml', 'a') as fnew:
  fnew.write(string_for_qcd)

with open(config_path + 'template_Run2.yml') as f:
  lines = f.readlines()
  with open(config_path + 'config_Run2.yml', 'w+') as f1:
    for line in lines: f1.write(line)
    f1.write(common_syst)
    f1.write(file_syst)
    f1.write("\nplots:\n  include: ['histos_dnn.yml']\n")

call(['../plotIt/plotIt', '-o ' + dest_path + '/figures/qcd', config_path + 'config_Run2.yml'], shell=False)
