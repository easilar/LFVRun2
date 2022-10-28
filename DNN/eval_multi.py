import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"] = "3"
import uproot
import pandas as pd
import awkward as ak
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
import tensorflow.keras.backend as K
from keras.callbacks import EarlyStopping, ModelCheckpoint
import multiprocessing

base_dir = os.getcwd().replace("DNN","") # Upper directory
processed = "aug22"
label = "rerun_multi"
#systs = ["nom",#"jesup","jesdown",
        #"puup","pudown","btagup_jes","btagdown_jes",
        #"btagup_hf","btagdown_hf","btagup_lf","btagdown_lf",
#        ]
systs = [ 'puup', 'pudown', 'btagup_hf', 'btagdown_hf', 'btagup_lf', 'btagdown_lf', 'btagup_hfstats1', 'btagdown_hfstats1', 'btagup_hfstats2', 'btagdown_hfstats2', 'btagup_lfstats1', 'btagdown_lfstats1', 'btagup_lfstats2', 'btagdown_lfstats2', 'btagup_cferr1', 'btagdown_cferr1', 'btagup_cferr2', 'btagdown_cferr2', 'up_jesAbsolute', 'down_jesAbsolute', 'up_jesBBEC1', 'down_jesBBEC1', 'up_jesEC2', 'down_jesEC2', 'up_jesFlavorQCD', 'down_jesFlavorQCD', 'up_jesRelativeBal', 'down_jesRelativeBal', 'up_jesAbsolute_year', 'down_jesAbsolute_year', 'up_jesBBEC1_year', 'down_jesBBEC1_year', 'up_jesEC2_year', 'down_jesEC2_year', 'up_jesRelativeSample_year', 'down_jesRelativeSample_year']
#systs = ["jesup","jesdown",
#systs = ["norm","jesup","jesdown","puup","pudown",]
#systs = ["btagup_jes","btagdown_jes","btagup_hf","btagdown_hf","btagup_lf","btagdown_lf",]

inputvars = ["Sel_muon1pt","Sel_muon1eta",
"Sel_tau1pt","Sel_tau1eta","Sel_tau1mass",
"Sel2_jet1pt","Sel2_jet2pt","Sel2_jet3pt",
"Sel2_jet1eta","Sel2_jet2eta","Sel2_jet3eta",
"Sel2_jet1mass","Sel2_jet2mass","Sel2_jet3mass",
"Sel2_jet1btag","Sel2_jet2btag","Sel2_jet3btag",
"Sys_METpt","Sys_METphi",
"chi2","chi2_SMW_mass","chi2_SMTop_mass",
"chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
"mutau_mass","mutau_dEta","mutau_dPhi","mutau_dR",
]
infile_dict = {"st":("hnevents_cut00000",2),"tt":("hnevents_cut000000",1)}

def run(inputs):
    year = inputs[0]
    syst = inputs[1]
    ch   = inputs[2]

    if 'syst' in syst: syst.replace('year', year)
    p="Multi"
    print("Start "+p+" LFV Evaluation")
    project_dir = "/home/itseyes/github/LFVRun2_ndf_integration/nanoaodframe/aug22_"+ch+"lfv/"+syst+"/"+year+"/"    # MODIFY!!
    path = project_dir
    flist = os.listdir(path)
    flist = [i for i in flist if ".root" in i]
    
    model_dir = label+"_"+p+processed+"/nom/best_model.h5"
    #print(model_dir)
    model = tf.keras.models.load_model(model_dir)
    #model.summary()
    
    eval_dir = label+"_"+p+processed+"/"+syst
    
    weights = ["evWeight"]
    hists_path = eval_dir+"/pred_hists_"+ch+"/"+year+"/"
    if not os.path.isdir(hists_path):
        os.makedirs(hists_path)
    
    for f in flist:
        fdict = f.split("_"+year)[0]
        f_dir = path+f
        outf_dir = hists_path+f.replace(".root","")+"_pred.root"
        infile = uproot.open(f_dir)
        tree = infile["outputTree2"]
        hnocut = infile["hcounter_nocut"]
        hpglepweight = infile["hnevents_pglep_cut0000"]
        hfullweight = infile["hnevents_cut0000"]
        hfinal = infile[infile_dict[ch][0]] 
        if len(tree) == 0:
            print("No events : "+f)
            continue
        print(f)
        pd_data = tree.arrays(inputvars,library="pd")
        pd_weight = tree.arrays(weights,library="np")
        pred_data = np.array(pd_data.filter(items = inputvars))
        pred = model.predict(pred_data,batch_size=128, workers=1, use_multiprocessing=False)
        pred = pred[:,infile_dict[ch][1]].tolist()
        #print(pred)
        pd_weight = pd_weight['evWeight'].tolist()
        binedges = [0,0.1,0.3,0.7,0.9,1.0]
        #dnnhist = np.histogram(pred,20,(0.0,1.0),weights=pd_weight,density=False)
        dnnhist = np.histogram(pred,bins=binedges,weights=pd_weight,density=False)
        dnnhist_entries = np.histogram(pred,bins=binedges,density=False)
        with uproot.recreate(outf_dir) as outf:
            outf["h_dnn_pred"] = dnnhist
            outf["h_dnn_entries"] = dnnhist_entries
            outf["hcounter_nocut"] = hnocut
            outf["hnevents_pglep_cut0000"] = hpglepweight
            outf["hnevents_cut0000"] = hfullweight
            outf["hnevents_final"] = hfinal

    K.clear_session()

if __name__ == '__main__':

    parameters = []
    for year in ["16pre","16post","17","18"]:
        for syst in systs:
            for ch in ['st', 'tt']:
                parameters.append((year, syst, ch))

    pool = multiprocessing.get_context("spawn").Pool(15)
    pool.map(run, parameters)
    pool.close()
    pool.join()