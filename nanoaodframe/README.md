# NanoAOD RDataFrame

With this package, we can do the following to NanoAOD files
- Apply good JSON
- Calculate b-tag weights for MC
- Apply Jet/MET corrections
- Skim events
- Process events with applying selections
- SF and uncertainty calculation

For this, we use the data frame concept to simplify the code,
but it requires to learn new way of using ROOT.
In a data frame concept, the operations on data should be defined in terms of functions,
which means that result of applying a function should produce outputs and no internal states are allowed.
In ROOT, RDataFrame allows us to borrow the concept for analysis.

This package code implements LFV analysis using NanoAOD format.


## I. Introduction

- The advantage of using RDataFrame is that you don't have to
worry about looping, filling of histograms at the right place,
writing out trees. In some sense you can concentrate more on thinking about data analysis.

- Draw back is that you have to think of a new way of doing old things.
This could getting used to. Also complicated algorithms may not
be so easy to implement with RDataFrame

- Purely data frame concept is not ideal for HEP data since
we have to treat real data and simulated data differently.
We must apply different operations depending on the input data, such as scale factors and corrections.
Therefore, some form of global state should be stored. 
In this package, the data frame concept is used together with object oriented concept for this purpose.


## II. Code

- The code consists of a main class NanoAODAnalyzerrdframe. 
There is one header file `.h` and one source file `.cpp` for it.
The class has several methods:
    - Object definitions. `selectElectrons()`, `selectMuons()`, ...
    - Additional derived variables definition. `defineMoreVars()`
    - Histogram definitions. `bookHists()`
    - Selections. `defineCuts()`
    - Plus some utility methods/functions. `gen4vec()`, `readjson()`, `removeOverlaps()`, `helper_1DHistCreator()`, `createHists()`

- **Users should modify**: object definitions, define additional variables, histogram definitions, selections.

- To use the class, you need to pass a pointer to TChain object, output root file name and a JSON file.
  If there is `"genWeight"` branch in the ROOT Tree then the input is assumed to be MC, otherwise data.
  Look at `nanoaoddataframe.cpp` to find how to use within C++.

- Compile
  Make sure to source proper evnv (new!)
  ``` bash
    source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/setup.sh
  ```
  Then compile
  ``` bash
    make clean && make -j 4   # max 6 effective.
    ```
    Need ROOT 6.26.06 or later then for correctionlib,
    this will compile and create a `libnanoaodrdframe.so` shared library that can be loaded in a ROOT session or macro:
    ```c++
    gSystem->Load("libnanoadrdframe.so");
    ```
    or within pyROOT (look in `processnanoaod.py`).


## III. Running over large dataset

### Usage of scripts
In the `scripts/` folder, there are scripts for skimming or processing the NanoAOD files.
All of the options are set in the script files `scripts/*.sh`.

#### Skim
```bash
# CAUTION::::::Create working directory under /data1 to reduce file i/o on disks
# You are at LFVRun2/nanoaodframe/
# First fetch dataset info from google spreadsheet:
# https://docs.google.com/spreadsheets/d/1KNvsRvXi3sgU45T2qOUztFSX3As4elZB325WaPnSkA8/edit#gid=569299692
# The argument will be the name of campaign in 'Campaign Summary' tab
# This will create file lists in data/dataset

python getDatasetInfo.py v8UL_2016pre

# Now, run script submitting slurm job
# arguments: folder, year (2016pre, 2016post, 2017, and 2018)
# You can add data or mc flag at the end of command if needed

python scripts/skim.py -V skim_test -Y 2018

OR

python scripts/skim.py -V skim_test -F mc -Y 2018

# To run specific datasets, add them in test_list of skim.py
# This script uses slurm on htop.
# Please specify which node to EXCLUDE from the batch job, in scripts/job_slurm_skim.sh

#To resubmit knowing that `folder_filename.log`
find log/201*/*/* | xargs grep runtime_error
find log/201*/*/* | xargs grep fault
find log/201*/*/* | xargs grep Traceback
find log/201*/*/* | xargs grep ERROR
find log/201*/*/* | xargs grep error

#Do this ONLY for systematic root file, unless will submit all variations in addition to nominal one
python scripts/skim.py -V skim_test -Y 2018 --dry | grep 270000_221AB515 | sh
```

#### Processing
`scripts/process.py` scripts can automatically run over all ROOT files in an input directory.
``` txt
Usage: python scripts/process.py -V skim_test -O test -Y 2018 -S theory

Options:
  -V, --version         Skim version: folder under /data1/common/skimmed_NanoAOD/
  -O, --outdir          Output folder in your working directory
  -Y YEAR, --year=YEAR  Select 2016pre, 2016post, 2017, or 2018
  -S SYST, --syst=SYST  Systematic: 'data' for Data, 'nosyst' for mc without uncertainties. Default is 'theory'. To run without theory unc for TT samples, put 'all'
  -D, --dataset         Put dataset folder name (eg. -D TTTo2L2Nu,QCD_Pt1000_MuEnriched) to process specific dataset.
  -F, --dataOrMC        Flag to choose Data or MC.
```

In some cases, you may want to submit single file per core using slurm.
`scripts/process.py` will do the job
``` txt
Usage: python scripts/process.py -V skim_test -O test -Y 2018 -S theory (-F data/mc)

For failing job detection:
find log/* | xargs grep runtime_error
find log/* | xargs grep fault
find log/* | xargs grep Traceback
```
Now, apply b SF rescaling, compute uncertainty envelope, etc. Let the `test\2018` is the folder containing histograms.
``` txt
python postprocess.py test 2018
```
Drawing histogram by plotIt
``` txt
cd 2018_postprocess
mkdir ../figure_2018
../../../plotIt/plotIt -o ../figure_2018/ ../../../plotIt/configs/config_18.yml -y -s

#For Run2,
python stack_signals_v2.py -I test #use v2 for now, v1 cannot deal with year based uncertainties
python plot_run2.py -I test
python print_syst_table.py -I test4 -Y 2018
```

### Usage of scripts
In the `scripts/` folder, there are scripts for skimming or processing the NanoAOD files.
All of the options are set in the script files `scripts/*.sh`.
