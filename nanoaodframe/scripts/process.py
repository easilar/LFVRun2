import os, sys, argparse
from subprocess import call

parser = argparse.ArgumentParser(usage="%prog [options]")
parser.add_argument("-V", "--version", dest="version", type=str, default="", help="Skim version: folder under /data1/common/skimmed_NanoAOD/")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="test", help="Output folder in your working directory")
parser.add_argument("-Y", "--year", dest="year", type=str, default="", help="Select 2016pre, 2016post, 2017, or 2018 runs")
parser.add_argument("-S", "--syst", dest="syst", type=str, default="theory", help="Systematic: 'data' for Data, 'nosyst' for mc without uncertainties. Default is 'theory'. To run without theory unc for TT samples, put 'all'.")
parser.add_argument("-D", "--dataset", dest="dataset", action="store", nargs="+", default=[], help="Put dataset folder name (eg. TTTo2L2Nu) to process specific one.")
parser.add_argument("-F", "--dataOrMC", dest="dataOrMC", type=str, default="", help="data or mc flag, if you want to process data-only or mc-only")
parser.add_argument("--dry", dest="dry", action="store_true", default=False, help="dryrun: not submitting jobs to slurm")
options = parser.parse_args()

year = options.year
workdir = os.getcwd()
indir = os.path.join('/data1/common/skimmed_NanoAOD/', options.version)
tgdir = os.path.join(workdir, options.outdir, year)
logdir = os.path.join(workdir, options.outdir, year, 'log')

os.makedirs(tgdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)

datadir = os.path.join(indir, 'data', year)
data_list = [os.path.join(datadir, s) for s in os.listdir(datadir)]
mcdir = os.path.join(indir, 'mc', year)
mc_list = [os.path.join(mcdir, s) for s in os.listdir(mcdir)]

dataset_list = data_list + mc_list

syst_list = ["", "__tesup", "__tesdown", "__jerup","__jerdown", "__jesAbsoluteup","__jesAbsolutedown",
             "__jesAbsolute_"+year[:4]+"up", "__jesAbsolute_"+year[:4]+"down",
             "__jesBBEC1up", "__jesBBEC1down", "__jesBBEC1_"+year[:4]+"up", "__jesBBEC1_"+year[:4]+"down",
             "__jesFlavorQCDup", "__jesFlavorQCDdown", "__jesRelativeBalup", "__jesRelativeBaldown",
             "__jesRelativeSample_"+year[:4]+"up", "__jesRelativeSample_"+year[:4]+"down"]
if year == "2018": syst_list.extend(["__jesHEMup", "__jesHEMdown"])

# tune and hdamp will appear as an individual dataset.
# thus no need to run in loop, but left here for double check
syst_ext = ["__tuneup", "__tunedown", "__hdampup", "__hdampdown",]

# Default syst for MC = "theory" = all SF + theory
# To run central values only, pass "nosyst"
# To calculate theory uncertainties (PDF, ME, PS, if available), put 'theory'
# Data should use systematic flag "data" not to construct event weights
if options.syst == "nosyst": syst_list = [""]

parameters = [] #order: (tgdir, indir, year, syst)
for ds in dataset_list:

    if len(options.dataset) > 0 and not any(i in ds for i in options.dataset): continue

    for src in syst_list:

        if options.dataOrMC == 'data':
            if 'data' not in ds[5:]: continue
        elif options.dataOrMC == 'mc':
            if 'mc' not in ds[5:]: continue

        dataset_name = ds.split('/')[-1]
        outdir = tgdir
        outfname = 'hist_' + dataset_name + src + '.root'

        ext_syst = False
        if any(i in dataset_name for i in syst_ext): ext_syst = True

        if 'data' in ds[5:]:
            if src == "":
                parameters.append([year, ds, outdir, outfname, "data"])
            else: continue

        else:
            if src == "" and not ext_syst:
                if options.syst == "all":
                    parameters.append((year, ds, outdir, outfname, "all", "", False))
                elif options.syst == "theory":
                    if any(i in dataset_name for i in ["TTTo", "TT_LFV", "ST_LFV"]): #include TTW,TTZ for theory unc?
                        parameters.append([year, ds, outdir, outfname, "theory"])
                    else: parameters.append([year, ds, outdir, outfname, "all"])
                elif options.syst == "nosyst":
                    parameters.append([year, ds, outdir, outfname, "nosyst"])
            elif src == "" and ext_syst:
                parameters.append([year, ds, outdir, outfname, "nosyst"])
            else:
                #os.makedir(os.path.join(tgdir, dataset_name+src)
                if any(i in dataset_name for i in syst_ext): continue
                parameters.append([year, ds, outdir, outfname, src[2:]])


for item in parameters:
    runString = "sbatch -J " + item[0] + '_' + item[3] + " scripts/job_slurm_process.sh " + item[0] + " " + item[1] + " " + item[2] + " " + item[3] + " " + workdir + " " +logdir + " " + item[4]

    print(runString)
    if not options.dry:
        call([runString], shell=True)
