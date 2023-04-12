import os, sys, argparse
from subprocess import call

parser = argparse.ArgumentParser(usage="%prog [options]")
parser.add_argument("-V", "--version", dest="version", type=str, default="", help="Skim version: folder under /data1/common/skimmed_NanoAOD/")
parser.add_argument("-Y", "--year", dest="year", type=str, default="", help="Select 2016pre, 2016post, 2017, or 2018 runs")
parser.add_argument("-D", "--dataset", dest="dataset", action="store", nargs="+", default=[], help="Put dataset folder name (eg. TTTo2L2Nu) to process specific one.")
parser.add_argument("-F", "--dataOrMC", dest="dataOrMC", type=str, default="", help="data or mc flag, if you want to process data-only or mc-only")
parser.add_argument("--dry", dest="dry", action="store_true", default=False, help="dryrun: not submitting jobs to slurm")
options = parser.parse_args()

year = options.year
workdir = os.getcwd()
tgdir = '/data1/common/skimmed_NanoAOD/' + options.version + '/DATAMC/' + year
log = '/data1/common/skimmed_NanoAOD/' + options.version + '/log/' + year

os.makedirs(tgdir.replace('DATAMC', 'data'), exist_ok=True)
os.makedirs(tgdir.replace('DATAMC', 'mc'), exist_ok=True)
os.makedirs(log, exist_ok=True)

for fn in os.listdir("data/dataset/v9UL_" + year):
    if 'json' in fn: continue
    fname = fn.replace('dataset_', '').replace('.txt', '')
    test_list = options.dataset
    if len(test_list) > 0:
        if not any(i in fname for i in test_list): continue
    if len(sys.argv) > 3:
        if options.dataOrMC == 'data':
            if '201' not in fname: continue
        elif options.dataOrMC == 'mc':
            if '201' in fname: continue

    with open(os.path.join("data/dataset/v9UL_" + year, fn), 'r') as f:
        for line in f.readlines():
            if line.startswith('#'): continue

            infile = line.rstrip('\n') 
            lsplit = infile.split('/')
            dirNum = lsplit[-2]
            rootName = lsplit[-1]

            if any('Run201' in l for l in lsplit):
                outputdir = tgdir.replace('DATAMC', 'data')
            else:
                outputdir = tgdir.replace('DATAMC', 'mc')

            os.makedirs(os.path.join(outputdir, fname), exist_ok=True)
            logdir = os.path.join(log, fname)
            os.makedirs(logdir, exist_ok=True)

            runString = "sbatch -J " + fname + " scripts/job_slurm_skim.sh " + year + " " + infile + " " + os.path.join(outputdir, fname) + " " + dirNum + '_' + rootName + " " + workdir + " " + logdir

            print(runString)
            if not options.dry:
                call([runString], shell=True)
