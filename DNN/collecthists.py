import ROOT
import os
import sys

projname = "aug22"
label = "rerun_multi"

hists = ["data_obs","hstacked_mc_h_dnn_pred",
        #"LFV_STu","LFV_TTc","LFV_STc", "LFV_TTu",
        "LFV_STc_s","LFV_STc_v","LFV_STc_t","LFV_STu_s","LFV_STu_v","LFV_STu_t",
        "LFV_TTc_s","LFV_TTc_v","LFV_TTc_t","LFV_TTu_s","LFV_TTu_v","LFV_TTu_t",
        "TT_di","TT_semi","ST","Wjets","Others"]
systs = {"nom":"nom",
        "puup":"puUp","pudown":"puDown",
        "btagup_hf":"btag_hfUp","btagdown_hf":"btag_hfDown",
        "btagup_lf":"btag_lfUp","btagdown_lf":"btag_lfDown",
        "btagup_hfstats1":"btag_hfstats1Up","btagdown_hfstats1":"btag_hfstats1Down",
        "btagup_lfstats1":"btag_lfstats1Up","btagdown_lfstats1":"btag_lfstats1Down",
        "btagup_hfstats2":"btag_hfstats2Up","btagdown_hfstats2":"btag_hfstats2Down",
        "btagup_lfstats2":"btag_lfstats2Up","btagdown_lfstats2":"btag_lfstats2Down",
        "btagup_cferr1":"btag_cferr1Up","btagdown_cferr1":"btag_cferr1Down",
        "btagup_cferr2":"btag_cferr2Up","btagdown_cferr2":"btag_cferr2Down",
        "up_jesAbsolute":"jesAbsoluteUp","down_jesAbsolute":"jesAbsoluteDown",
        "up_jesAbsolute_year":"jesAbsolute_yearUp","down_jesAbsolute_year":"jesAbsolute_yearDown",
        "up_jesBBEC1":"jesBBEC1Up","down_jesBBEC1":"jesBBEC1Down",
        "up_jesBBEC1_year":"jesBBEC1_yearUp","down_jesBBEC1_year":"jesBBEC1_yearDown",
        "up_jesEC2":"jesEC2Up","down_jesEC2":"jesEC2Down",
        "up_jesEC2_year":"jesEC2_yearUp","down_jesEC2_year":"jesEC2_yearDown",
        "up_jesFlavorQCD":"jesFlavorQCDUp","down_jesFlavorQCD":"jesFlavorQCDDown",
        "up_jesRelativeBal":"jesRelativeBalUp","down_jesRelativeBal":"jesRelativeBalDown",
        "up_jesRelativeSample_year":"jesRelativeSample_yearUp","down_jesRelativeSample_year":"jesRelativeSample_yearDown",
        }

runs_histdict = {"run16APV":"stackhist_19.5.root","run16":"stackhist_16.8.root","run17":"stackhist_41.48.root","run18":"stackhist_59.83.root","run2":"stackhist_137.65.root"}
outfolder = "./pred_"+projname
if "multi" in label:
	lfvprocs = {"Multi": outfolder+"/pred_RUN_merged.root"}
else:
	lfvprocs = {"ST":outfolder+"/pred_RUN_cat1.root","TT":outfolder+"/pred_RUN_cat2.root"}

if not os.path.isdir(outfolder):
    os.makedirs(outfolder)

for run in runs_histdict.keys():

    rname = runs_histdict[run]

    for lfvproc in lfvprocs.keys():

        outfname = lfvprocs[lfvproc].replace("RUN",run) 

        outf = ROOT.TFile(outfname,"RECREATE")
        for key, value in systs.items():
            infilename = label+"_"+lfvproc+projname+"/"+key+"/pred_hists/noblind/"+rname
            inf = ROOT.TFile(infilename)
            outf.cd()
            for h in hists:
                tmphist = inf.Get(h)
                if not tmphist:
                    print("wrong hist",h)
                    continue
                copyhname = h
                if key is not "nom":
                    copyhname = h+"_"+value
                print(copyhname)
                tmphist_copy = tmphist.Clone(copyhname)
                tmphist_copy.Write()
            inf.Close()
        outf.Close()
