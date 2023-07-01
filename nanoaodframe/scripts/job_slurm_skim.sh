#!/bin/bash

#SBATCH -J LFV
#SBATCH -p gpu,cpu -x gpu-0-2,compute-0-3
#SBATCH -N 1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --open-mode=append
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb 
#SBATCH --comment python
#SBATCH --time 03:00:00
#SBATCH --hint=compute_bound

year=$1
infile=$2
outpath=$3
outfile=$4
workdir=$5
logdir=$6

cd $workdir
source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/setup.sh
echo "python skimonefile.py -Y $year -I $infile -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log"
python skimonefile.py -Y $year -I $infile -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log
