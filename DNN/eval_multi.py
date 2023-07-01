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
processed = "June2023_GoingtoPrep_2"
label = "top_lfv_multiClass"
systs_tofile = ['jerdown', 'jerup', 'jesAbsolute_yeardown', 'jesAbsolute_yearup', 'jesAbsolutedown', 'jesAbsoluteup', 'jesBBEC1_yeardown', 'jesBBEC1_yearup', 'jesBBEC1down', 'jesBBEC1up', 'jesFlavorQCDdown', 'jesFlavorQCDup', 'jesRelativeBaldown', 'jesRelativeBalup', 'jesRelativeSample_yeardown', 'jesRelativeSample_yearup', 'jesHEMdown','jesHEMup','tesdown', 'tesup']
systs_toTTTfile = ['tunedown','tuneup','hdampdown','hdampup']
systs_pdf = ['pdf'+str(i) for i in range(1,101)]  #100 PDF variations
systs_pdf = systs_pdf + ['pdfdown','pdfup'] # Alphas
syst_otherTheory = ['ps0','ps1','ps2','ps3','scale0','scale1','scale2','scale3','scale4','scale5']
syst_tauidjets = ['tauidjetHighptextrapdown', 'tauidjetHighptextrapup', 'tauidjetHighptstat_bin1down', 'tauidjetHighptstat_bin1up', 'tauidjetHighptstat_bin2down', 'tauidjetHighptstat_bin2up', 'tauidjetHighptstatdown', 'tauidjetHighptstatup', 'tauidjetHighptsystdown', 'tauidjetHighptsystup', 'tauidjetSystULyeardown', 'tauidjetSystULyearup', 'tauidjetSystallerasdown', 'tauidjetSystallerasup', 'tauidjetSystdm0ULyeardown', 'tauidjetSystdm0ULyearup', 'tauidjetSystdm10ULyeardown', 'tauidjetSystdm10ULyearup', 'tauidjetSystdm11ULyeardown', 'tauidjetSystdm11ULyearup', 'tauidjetSystdm1ULyeardown', 'tauidjetSystdm1ULyearup', 'tauidjetUncert0down', 'tauidjetUncert0up', 'tauidjetUncert1down', 'tauidjetUncert1up']
systs_toweight = ['btagcferr1down', 'btagcferr1up', 'btagcferr2down', 'btagcferr2up', 'btaghfdown', 'btaghfstats1down', 'btaghfstats1up', 'btaghfstats2down', 'btaghfstats2up', 'btaghfup', 'btaglfdown', 'btaglfstats1down', 'btaglfstats1up', 'btaglfstats2down', 'btaglfstats2up', 'btaglfup', 'muiddown', 'muidup', 'muisodown', 'muisoup', 'mutrgdown', 'mutrgup','prefiredown','prefireup', 'pudown', 'puup', 'tauideldown', 'tauidelup', 'tauidmudown', 'tauidmuup']+syst_tauidjets+systs_pdf+syst_otherTheory

systs = systs_tofile+systs_toTTTfile+systs_toweight+['nom'] 
#convertDict = {'pdfdown':'pdfdown','pdfup':'pdfup','scale1':'scaledown','scale5':'scaleup'}
inputvars = [ "Muon1_pt","Muon1_eta",
        "Tau1_pt","Tau1_mass","Tau1_eta",
        "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagDeepFlavB",
        "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagDeepFlavB",
        "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagDeepFlavB",
        "chi2","chi2_SMW_mass","chi2_SMTop_mass",
        "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
        "mutau_mass","mutau_dEta","mutau_dPhi","mutau_dR",
	"MET_pt"
        ]
#infile_dict = {"st":("hnevents_cut00000",2),"tt":("hnevents_cut000000",1), "bkg":("",0)}
#infile_dict = {"st":("hnevents_cut00000",2),"tt":("hnevents_cut00000",1), "bkg":("hnevents_cut00000",0)}
#infile_dict = {"dummy":("hnevents_cut00000",2)}#,"tt":("hnevents_cut00000",1), "bkg":("hnevents_cut00000",0)}
discriminators = {"p_st" : 2, "p_tt" : 1 , "p_bkg" : 0 , "p_st_tt" : 999, "p_st_tt_ob" : 999 }

def run(inputs):
    year = inputs[0]
    syst = inputs[1]
    discriminator_key = inputs[2]
    alpha = inputs[3]
    print(year, syst , discriminator_key , alpha)
    if 'year' in syst:
        syst_ = syst.replace('year', year[:4])
    else: syst_ = syst
    if '16pre' in year and syst in syst_tauidjets: syst_ = syst_.replace('2016', '2016_preVFP')
    if '16post' in year and syst in syst_tauidjets: syst_ = syst_.replace('2016', '2016_postVFP')

    binedges = [0,0.1,0.3,0.7,0.9,1.0]
    #if "ob" in discriminator_key : binedges = np.linspace(0.0, 200.0, num=21).tolist()
    if "ob" in discriminator_key : binedges = [0,5,10,30,60,100] 

    p="Multi"
    print("Start "+p+" LFV Evaluation")
    project_dir = "/data1/users/minerva1993/work/lfv_production/LFVRun2/nanoaodframe/v9_0608/"+year+"/"
    path = project_dir
    flist = os.listdir(path)
    flist = [i for i in flist if (".root" in i) and not ("__" in i)]
    print(flist)
    
    model_dir = label+"_"+processed+"/nom/best_model.h5"
    model = tf.keras.models.load_model(model_dir)
    
    eval_dir = label+"_"+processed+"/"
    
    hists_path = eval_dir+"/"+year+"/preds/"+discriminator_key+"/"+str(alpha).replace(".","p")+"/"
    if not os.path.isdir(hists_path):
        os.makedirs(hists_path)
    
    syst__ = ""
    if syst == "nom": syst__ = ""
    else: syst__ = "__" + syst_
    if syst in systs_toweight : eventWeight = "eventWeight"+syst__
    else : eventWeight = "eventWeight"
    weights = [eventWeight] 
    if syst in systs_pdf+syst_otherTheory: eventWeight_tt = eventWeight
    #print("BEFORE FOR LOOP :", weights  ) 
    for f in flist:
        #Evaluate DATA only for nom
        if "SingleMuon" in f and not "nom" in syst: continue
        if not "18" in year and "HEM" in syst: continue
        if syst in systs_pdf+syst_otherTheory: 
            print("YESSSSS 1")
            if not "TTT" in f  and "TT_LFV" not in f : 
               eventWeight = "eventWeight"
               weights = ["eventWeight"]
               print("YESSSSS 2")
            if "TTT" in f or "TT_LFV" in f:
               eventWeight = eventWeight_tt
               weights = [eventWeight]
               print(weights, eventWeight)
        outf_dir = hists_path+f.replace(".root",syst__+".root")
        #if syst in convertDict.keys():       
        #    outf_dir = hists_path+f.replace(".root","__"+convertDict[syst]+".root")
        if "TTT" in f:
            if syst in systs_toTTTfile : f_dir = path+f.replace(".root",syst__+".root")
            else : f_dir = path+f 
        if syst in systs_tofile: f_dir = path+f.replace(".root",syst__+".root")
        else: f_dir = path+f
        print("OUT PATH : ", outf_dir) 
        print("FILE PATH:", f_dir)	
        infile = uproot.open(f_dir)
        tree = infile["Events"]
        hcounter = infile["hcounter"]
        if not "SingleMuon" in f :
            h_nevents_S4_nobtag = infile["h_nevents_S4_nobtag"]
            h_nevents_S4 = infile["h_nevents_S4"]	
        if len(tree) == 0:
            print("No events : "+f)
            pred = []
            pd_weight = []
            dnnhist = np.histogram(pred,bins=binedges,weights=pd_weight,density=False)
            dnnhist_entries = np.histogram(pred,bins=binedges,density=False)
        else:
            print("Good file :",f)
            pd_data = tree.arrays(inputvars,library="pd")
            pd_weight = tree.arrays(weights,library="np")
            pred_data = np.array(pd_data.filter(items = inputvars))
            pred = model.predict(pred_data,batch_size=128, workers=1, use_multiprocessing=False)
            if discriminator_key == "p_st_tt_ob"  : 
                 rmzeros = pred[:,0] 
                 rmzeros[rmzeros<=0.00] = 0.001
                 pred[:,0] = rmzeros
                 if alpha <10: pred = (((1-alpha)*pred[:,2]+alpha*pred[:,1])/pred[:,0]).tolist()
                 else: pred = ((pred[:,2]+pred[:,1])/pred[:,0]).tolist()
                 pred = np.array(pred)
                 pred[pred>=90] = 90
                 pred = pred.tolist()
            elif discriminator_key == "p_st_tt"  : pred = (pred[:,2]+pred[:,1]).tolist()
            elif discriminator_key == "p_max"  : pred = (np.amax(pred, axis=1)).tolist()
            elif discriminator_key == "p_mean"  : pred = (np.mean(pred, axis=1)).tolist()
            elif discriminator_key == "p_min"  : pred = (np.min(pred, axis=1)).tolist()
            else  : pred = (pred[:,discriminators[discriminator_key]]).tolist()
            pd_weight = pd_weight[eventWeight].tolist()
            dnnhist = np.histogram(pred,bins=binedges,weights=pd_weight,density=False)
            dnnhist_entries = np.histogram(pred,bins=binedges,density=False)
        with uproot.recreate(outf_dir) as outf:
            outf["h_dnn_pred"] = dnnhist
            outf["h_dnn_entries"] = dnnhist_entries
            outf["hcounter"] = hcounter
            if not "SingleMuon" in f : 
                outf["h_nevents_S4_nobtag"] = h_nevents_S4_nobtag
                outf["h_nevents_S4"] = h_nevents_S4
            #if ch=="bkg": 
            #outf["hnevents_final_st"] = hfinal_st
            ##outf["hnevents_final_tt"] = hfinal_tt
            #else:
            #outf["hnevents_final"] = hfinal

    K.clear_session()

if __name__ == '__main__':

    parameters = []
    #for discriminator in ["p_st","p_tt","p_bkg","p_st_tt_ob"]: #,"p_st_tt","p_st_tt_ob"]:
    for discriminator in ["p_st_tt_ob"]: #,"p_st_tt","p_st_tt_ob"]:
       #for alpha in np.around(np.linspace(0.1, 0.9, num=9),decimals=1).tolist():
       for alpha in [0.1]:
         for year in ["2016pre","2016post","2017","2018"]:
            for syst in systs:
               parameters.append((year, syst, discriminator,alpha))

    pool = multiprocessing.get_context("spawn").Pool(12)
    pool.map(run, parameters)
    pool.close()
    pool.join()
