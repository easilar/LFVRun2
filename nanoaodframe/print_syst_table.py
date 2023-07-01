import os, shutil, re, sys
from subprocess import call
from collections import OrderedDict

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-I", "--input",  dest="input", type="string", default="", help="Input folder name")
parser.add_option("-Y", "--year",  dest="year", type="string", default="", help="Select 2016pre/post, 2017, or 2018 for years")
(options, args) = parser.parse_args()

year = options.year

config_path = '../plotIt/configs/'
dest_path = os.path.join('./', options.input)
tmp_file_name = 'temp_' + year + '_forSyst.yml'
string_to_add = 'systematics:\n'
plot_to_add = "plots:\n  include: ['histos_yield.yml']\n\n"

if   year == "2016pre" : tauYear = "UL2016_preVFP"
elif year == "2016post": tauYear = "UL2016_postVFP"
else                   : tauYear = "UL" + year

if os.path.exists(config_path + tmp_file_name):
    os.remove(config_path + tmp_file_name)
#shutil.copy2(config_path + 'config_' + year + '.yml', config_path + tmp_file_name)

unc_cat = OrderedDict([
('all', ['xsec', 'pu', 'muid', 'muiso', 'mutrg',
         #'tauidjet', 'tauidel', 'tauidmu', 'tes',
         'tauidjetUncert0', 'tauidjetUncert1', 'tauidjetSystalleras',
         'tauidjetSyst'+tauYear, 'tauidjetSystdm0'+tauYear, 'tauidjetSystdm1'+tauYear,
         'tauidjetSystdm10'+tauYear, 'tauidjetSystdm11'+tauYear,
         'tauidjetHighptstat_bin1', 'tauidjetHighptstat_bin2',
         'tauidjetHighptsyst', 'tauidjetHighptextrap',
         'tauidel', 'tauidmu', 'tes',
         'btaghf', 'btaglf', 'btaghfstats1', 'btaglfstats1',
         'btaghfstats2', 'btaglfstats2', 'btagcferr1', 'btagcferr2',
         'jesAbsolute', 'jesAbsolute_'+year[:4], 'jesBBEC1', 'jesBBEC1_'+year[:4],
         'jesFlavorQCD', 'jesRelativeBal', 'jesRelativeSample_'+year[:4], 'jer',
         'scale', 'ps', 'pdf',
         'tune', 'hdamp',]),
('pu', ['pu']),
('prefire', ['prefire']),
('xsec', ['xsec']),
('muon', ['muid', 'muiso', 'mutrg']),
#('tauid', ['tauidjet', 'tauidel', 'tauidmu']),
('tauid', ['tauidjetUncert0', 'tauidjetUncert1', 'tauidjetSystalleras',
           'tauidjetSyst'+tauYear, 'tauidjetSystdm0'+tauYear, 'tauidjetSystdm1'+tauYear,
           'tauidjetSystdm10'+tauYear, 'tauidjetSystdm11'+tauYear, 'tauidjetHighpt',
           'tauidel', 'tauidmu']),
('tauidjetUncert0', ['tauidjetUncert0']),
('tauidjetUncert1', ['tauidjetUncert1']),
('tauidjetSystalleras', ['tauidjetSystalleras']),
('tauidjetSyst'+tauYear, ['tauidjetSyst'+tauYear]),
('tauidjetSystdm0'+tauYear, ['tauidjetSystdm0'+tauYear]),
('tauidjetSystdm1'+tauYear, ['tauidjetSystdm1'+tauYear]),
('tauidjetSystdm10'+tauYear, ['tauidjetSystdm10'+tauYear]),
('tauidjetSystdm11'+tauYear, ['tauidjetSystdm11'+tauYear]),
('tauidjetHighptstat_bin1', ['tauidjetHighptstat_bin1']),
('tauidjetHighptstat_bin2', ['tauidjetHighptstat_bin2']),
('tauidjetHighptsyst', ['tauidjetHighptsyst']),
('tauidjetHighptextrap', ['tauidjetHighptextrap']),
('tauidjet', ['tauidjetUncert0', 'tauidjetUncert1', 'tauidjetSystalleras',
              'tauidjetSyst'+tauYear, 'tauidjetSystdm0'+tauYear, 'tauidjetSystdm1'+tauYear,
              'tauidjetSystdm10'+tauYear, 'tauidjetSystdm11'+tauYear,
              'tauidjetHighptstat_bin1', 'tauidjetHighptstat_bin2',
              'tauidjetHighptsyst', 'tauidjetHighptextrap']),
('tauidel', ['tauidel']),
('tauidmu', ['tauidmu']),
('tes', ['tes']),
('jesAbsolute', ['jesAbsolute']),
('jesAbsolute_'+year[:4], ['jesAbsolute_'+year[:4]]),
('jesBBEC1', ['jesBBEC1']),
('jesBBEC1_'+year[:4], ['jesBBEC1_'+year[:4]]),
('jesFlavorQCD', ['jesFlavorQCD']),
('jesRelativeBal', ['jesRelativeBal']),
('jesRelativeSample_'+year[:4], ['jesRelativeSample_'+year[:4]]),
('jesAll', ['jesAbsolute', 'jesAbsolute_'+year[:4], 'jesBBEC1', 'jesBBEC1_'+year[:4],
            'jesFlavorQCD', 'jesRelativeBal', 'jesRelativeSample_'+year[:4]]),
('jer', ['jer']),
('scale', ['scale']),
('ps', ['ps']),
('hdamp', ['hdamp']),
('pdf', ['pdf']),
('tune', ['tune']),
('btaghf', ['btaghf']),
('btaglf', ['btaglf']),
('btaghfstats1', ['btaghfstats1']),
('btaglfstats1', ['btaglfstats1']),
('btaghfstats2', ['btaghfstats2']),
('btaglfstats2', ['btaglfstats2']),
('btagcferr1', ['btagcferr1']),
('btagcferr2', ['btagcferr2']),
('bAll', ['btaghf', 'btaglf', 'btaghfstats1', 'btaglfstats1',
          'btaghfstats2', 'btaglfstats2', 'btagcferr1', 'btagcferr2']),
])

if year == "2018":
  unc_cat['jesHEM'] = ['jesHEM']
  unc_cat['jesAll'].append('jesHEM')
  unc_cat['all'].append('jesHEM')

#Create syst yaml everytime, to avoid crash due to different treatments
with open(config_path + "config_" + year + ".yml") as f:
    lines = f.readlines()
    with open(config_path + "syst_all_" + year + ".yml", 'w+') as f1:
        write_lines = False
        for line in lines:
            if write_lines and 'plots' in line: write_lines = False
            if write_lines and line in ['\n', '\r\n']: write_lines = False
            if not write_lines: pass
            else: f1.write(line)
            if 'systematics' in line: write_lines = True

#sys.exit()

for key, value in unc_cat.items():
    print("Print yield tables with source: " + key)
    #prep for systematics
    with open(config_path + 'config_' + year + '.yml', 'r') as f:
        with open(config_path + tmp_file_name, 'w') as f1:
            syst_delete = False
            plot_delete = False
            for line in f:
                if 'luminosity-error:' in line:
                    line = line[:line.find(':')+2]
                    line += "0.0\n"
                    line += "  syst-only: true\n"
                if 'root:' in line:
                    line = line[:line.find(':')+2]
                    line += "'" + os.path.join(dest_path, year + "_postprocess") + "'\n"
                if 'systematics' in line: syst_delete = True
                if 'plots:' in line:
                    syst_delete = False
                    plot_delete = True
                if 'legend' in line: plot_delete = False
                if syst_delete == True: continue
                if plot_delete == True: continue
                f1.write(line)

    with open(config_path + "syst_all_" + year + ".yml") as f:
        lines = f.readlines()
        with open(config_path + tmp_file_name, "a") as f1:
            f1.writelines(plot_to_add)
            for line in lines:
                if any(x in line for x in value): string_to_add += str(line)
            f1.writelines(string_to_add)
    syst_postfix = key
    call(['../plotIt/plotIt', '-o ' + os.path.join(dest_path, 'figure_' + year), config_path + tmp_file_name, '-s'], shell=False)

    #Remove signals in TT specific sources
    with open(os.path.join(dest_path, 'figure_' + year, 'systematics.tex'), 'r') as f:
        with open(os.path.join(dest_path, 'figure_' + year, 'systematics_' + syst_postfix + '.tex'), 'w+') as f1:
            isTTsyst = False
            if key in ['scale', 'ps', 'hdamp', 'tune', 'pdf']: isTTsyst = True
            for line in f:
                isTT = False
                if 'ttbar' in line or 'LFV' in line: isTT = True
                if not isTT and isTTsyst and 'hline' not in line and 'Total' not in line: continue
                if any(str(x) in line[0] for x in list(range(1,10))): line = line[1:]
                f1.write(line)

    string_to_add = 'systematics:\n'

unc_summary = OrderedDict([
('xsec', 'Cross section'), ('pu', 'Pileup'), ('prefire', 'Prefire Reweight'), ('muon', 'Muon SF'),
('tauid', 'Tau ID'), ('tes', 'TES'),
('jesAll', 'JES'), ('jer', 'JER'),
('scale', 'ME scale'), ('ps', 'PS scale'), ('hdamp', 'ME-PS matching'), ('pdf', 'PDF'), ('tune', 'Underlying event'),
('bAll', 'b-tagging shape'), ('all', 'Total sys. unc.'),
])

print("Generating summary table...")
#Gather all results into one summary table
with open("total_syst_template.tex") as f:
    lines = f.readlines()
    with open(os.path.join(dest_path, 'figure_' + year, 'total_syst.tex'), "w") as f1:
        for line in lines:
            #if year != '2017' and 'Prefire' in line: continue
            #if 'Prefire' in line: continue
            for key, value in unc_summary.items():
                if value in line:
                    with open(os.path.join(dest_path, 'figure_' + year, 'systematics_' + key + '.tex'),'r') as f2:
                        lines = f2.read().splitlines()
                        last_line = lines[-1]
                        last_line = last_line[last_line.find("&")+1:]
                        line = line.rstrip('\n') + last_line + '\n'
            f1.write(line)
