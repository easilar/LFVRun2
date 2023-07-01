import ROOT
import sys, os
from subprocess import check_call
from os import listdir, path
import collections
import glob
import subprocess

dir_path_st = "/home/ecasilar/lfv_ana_upgrade/DNN/rerun_multi_Multiaug22/"
dirs_st = os.listdir(dir_path_st)

chs = ['ST_LFV_TCMuTau_Scalar', 'ST_LFV_TCMuTau_Tensor', 'ST_LFV_TCMuTau_Vector',
         'ST_LFV_TUMuTau_Scalar', 'ST_LFV_TUMuTau_Tensor', 'ST_LFV_TUMuTau_Vector',
         'TT_LFV_TToCMuTau_Scalar', 'TT_LFV_TToCMuTau_Tensor', 'TT_LFV_TToCMuTau_Vector',
         'TT_LFV_TToUMuTau_Scalar', 'TT_LFV_TToUMuTau_Tensor', 'TT_LFV_TToUMuTau_Vector']

for year in ["16pre_postprocess" ,"16post_postprocess" ,"17_postprocess" ,"18_postprocess"] :
	for disc in ["p_st","p_tt","p_bkg","p_st_tt", "p_max", "p_mean", "p_min"]:
		for ch in chs:
			infilett = dir_path_st +"/"+year + "/" + disc +"/"+"hist_"+ ch  +"_tt.root"
			infilest = dir_path_st +"/"+year + "/" + disc +"/"+"hist_"+ ch  +"_st.root"
			outfile  = dir_path_st +"/"+year + "/" + disc +"/"+"hist_"+ ch  +".root"
			print(['hadd','-f', outfile] + [infilest] + [infilett])
			subprocess.run(['hadd','-f', outfile] + [infilest] + [infilett])
