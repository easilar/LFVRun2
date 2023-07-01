import os
import sys

#os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.utils.np_utils import to_categorical
from utils.plots import *
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import KFold
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV

root_dir = os.getcwd().replace("DNN","") # Upper directory
# MODIFY !!!
processed = "June2023_GoingtoPrep_2"
syst = "nom"
label = "top_lfv_multiClass"
class_names = ["bkg","sigTT", "sigST"]

print("Start multi LFV Training")
epochs = 1000
inputvars_st = [ "Muon1_pt","Muon1_eta",
        "Tau1_pt","Tau1_mass","Tau1_eta", 
        "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagDeepFlavB",
        "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagDeepFlavB",
        "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagDeepFlavB",
        "chi2","chi2_SMW_mass","chi2_SMTop_mass",
        "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
        "mutau_mass","mutau_dEta","mutau_dPhi","mutau_dR",
	"MET_pt"
        ]
#"MET_pt" : helps to the expected limits, do not remove from the input vars.
sbratio = 1 # sig:bkg = 1:1

kfold = KFold(n_splits=10, shuffle=True)

train_outdir = label+"_"+processed+"/"+syst
os.makedirs(train_outdir, exist_ok=True)

siglist_st = ["ST_LFV_TCMuTau_Scalar","ST_LFV_TCMuTau_Vector","ST_LFV_TCMuTau_Tensor","ST_LFV_TUMuTau_Vector","ST_LFV_TUMuTau_Scalar","ST_LFV_TUMuTau_Tensor"]
siglist_tt = ['TT_LFV_TCMuTau_Scalar', 'TT_LFV_TCMuTau_Tensor', 'TT_LFV_TCMuTau_Vector', 'TT_LFV_TUMuTau_Scalar', 'TT_LFV_TUMuTau_Tensor', 'TT_LFV_TUMuTau_Vector']
years = ["2017","2018","2016pre","2016post"]
#project_dir = "/data1/users/itseyes/LFV/processed_LFV/v9test2_theory/"
project_dir = "/data1/users/minerva1993/work/lfv_production/LFVRun2/nanoaodframe/v9_0608/"

df_sig_st_list = []
df_sig_tt_list = []
df_bkg_tt_list = []
#We use all the years together
for year in years:
   project_dir_y = project_dir+"/"+year+"/"

   #Concatinate ST signals
   for sig_tree in siglist_st:
      sig_tree_ = uproot.open(project_dir_y+"hist_"+sig_tree+".root")["Events"]
      df_sig_ = sig_tree_.arrays(inputvars_st,library="pd")
      df_sig_st_list.append(df_sig_)
   
   #Concatinate TT signals
   for sig_tree in siglist_tt:
      sig_tree_ = uproot.open(project_dir_y+"hist_"+sig_tree+".root")["Events"]
      df_sig_ = sig_tree_.arrays(inputvars_st,library="pd")
      df_sig_tt_list.append(df_sig_)

   bkg1_filedir_tt = project_dir_y+"hist_TTTo2L2Nu.root"
   bkg2_filedir_tt = project_dir_y+"hist_TTToSemiLeptonic.root"
   bkg1_tree_tt = uproot.open(bkg1_filedir_tt)["Events"]
   bkg2_tree_tt = uproot.open(bkg2_filedir_tt)["Events"]
   df_bkg1_tt = bkg1_tree_tt.arrays(inputvars_st,library="pd")
   df_bkg2_tt = bkg2_tree_tt.arrays(inputvars_st,library="pd")
   df_bkg_tt_list.append(df_bkg1_tt)
   df_bkg_tt_list.append(df_bkg2_tt)

df_sig_st = pd.concat(df_sig_st_list)
df_sig_tt = pd.concat(df_sig_tt_list)
df_bkg = pd.concat(df_bkg_tt_list)

ntotsig_tt = len(df_sig_st)
ntotsig_st = len(df_sig_tt)
ntotbkg = len(df_bkg)
print(df_bkg.replace(np.nan, 0))
print(df_sig_st.replace(np.nan, 0))
print(df_sig_tt.replace(np.nan, 0))
print("sig tt:", ntotsig_tt)
print("sig st:",ntotsig_st)
print("tot bkg:",ntotbkg)
nsig_st = ntotsig_st
nsig_tt = ntotsig_tt
nbkg = ntotbkg
nsig = min(nsig_st,nsig_tt)
print(nsig)

if nsig >= nbkg:
    if sbratio == 1:
        nsig = nbkg
    elif sbratio == 2:
        nsig = int(nbkg/2)
else:
    if sbratio == 1:
        nbkg = nsig
    elif nsig>=int(nbkg/2):
        if sbratio == 2:
            nsig = int(nbkg/2)
    else:
        "Error::Check the number of events!"
        sys.exit()

print("Take LFV : "+str(nsig)+" events")
print("Take TT  : "+str(nbkg)+" events")
df_sig_st = df_sig_st.sample(n=nsig)
df_sig_tt = df_sig_tt.sample(n=nsig)
df_bkg = df_bkg.sample(n=nbkg)
df_sig_st["category"] = 2
df_sig_tt["category"] = 1
df_bkg["category"] = 0

pd_data = pd.concat([df_sig_tt,df_sig_st,df_bkg])
colnames = pd_data.columns
print(pd_data.head())
print("Col names:",colnames)
pd_sig_st = pd_data[pd_data['category'] == 2]
pd_sig_tt = pd_data[pd_data['category'] == 1]
pd_bkg = pd_data[pd_data['category'] == 0]

print("Plotting corr_matrix total")
plot_corrMatrix(pd_data,train_outdir,"total")
print("Plotting corr_matrix ST")
plot_corrMatrix(pd_sig_st,train_outdir,"sig_st")
print("Plotting corr_matrix TT")
plot_corrMatrix(pd_sig_tt,train_outdir,"sig_tt")
print("Plotting corr_matrix bkg")
plot_corrMatrix(pd_bkg,train_outdir,"bkg")

pd_data = pd_data.sample(frac=1).reset_index(drop=True)

print(pd_data.head())

#    # Min-Max Scale
#    scaler = MinMaxScaler()
#    pd_data = pd.DataFrame(scaler.fit_transform(pd_data), columns=pd_data.columns)

x_total = np.array(pd_data.filter(items = inputvars_st))
y_total = np.array(pd_data.filter(items = ['category']))

y_cat = to_categorical(y_total)
print("Y cat: ", y_cat)
# Splitting between training set and cross-validation set
numbertr = len(y_total)
trainlen = int(0.7*numbertr) # Fraction used for training

from sklearn.model_selection import train_test_split
x_train, x_val, y_train, y_val = train_test_split(x_total, y_cat, test_size=0.3)
print(len(x_train),len(x_val),len(y_train),len(y_val))
'''
x_train = x_total[:trainlen,0::]
y_train = y_total[:trainlen]
print("y_train:" , y_train)
x_val = x_total[trainlen:,0::]
y_val = y_total[trainlen:]
'''

patience_epoch = 30
# Early Stopping with Validation Loss for Best Model
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(train_outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)
print("xtrain shape:",x_train.shape)
###################################################
#                      Model                      #
###################################################

model = tf.keras.models.Sequential()

###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
#model.add(tf.keras.layers.Dense(100, input_shape = (x_train.shape[1],), activation = "relu"))
#drop = 0.2
activation_function='relu'
#activation_function='elu'
weight_initializer = 'random_normal'
#weight_initializer = 'he_uniform'
###############   Hidden Layer 1    ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(50, activation=activation_function))
#model.add(tf.keras.layers.Dropout(drop))

###############   Hidden Layer 2    ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(50, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
#model.add(tf.keras.layers.Dropout(drop))

    
###############    Output Layer     ###############
model.add(tf.keras.layers.Dense(3, activation="softmax"))
batch_size = 1024
model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), loss="categorical_crossentropy", metrics = ["accuracy"])

model.summary()
hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, 
                                validation_data=(x_val,y_val), callbacks=[es, mc])

#tf.keras.models.save_model(model,'model_{epoch:%d}.h5', overwrite=True, include_optimizer=True)
#print(hist.history.keys())


pred_train = model.predict_classes(x_train)
print("pred_train", pred_train)
print("orig train", y_train)
y_train = np.argmax(y_train, axis=1)
from sklearn.metrics import confusion_matrix
#train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[1]]).T, columns=["True", "Pred"])
pred_val = model.predict_classes(x_val)
y_val = np.argmax(y_val, axis=1)
print("pred_val", pred_val)
print("orig train", y_val)
print("conf matrix on train set ")
print(confusion_matrix(y_train, pred_train))
print("conf matrix on val set ")
print(confusion_matrix(y_val, pred_val))
#val_result = pd.DataFrame(np.array([y_val.T[0], pred_val.T[1]]).T, columns=["True", "Pred"])
plot_confusion_matrix(y_val, pred_val, classes=class_names,
                    title='Confusion matrix, without normalization', savename=train_outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=train_outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names,
                    title='Confusion matrix, without normalization', savename=train_outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=train_outdir+"/norm_confusion_matrix_train.pdf")
pred_val = model.predict(x_val)
pred_train = model.predict(x_train)
print(pred_val)
print(pred_val.T)
print(y_val)
print(y_val.T)
val_result = pd.DataFrame(np.array([y_val.T, pred_val.T[1]]).T, columns=["True", "Pred"])
train_result = pd.DataFrame(np.array([y_train.T, pred_train.T[1]]).T, columns=["True", "Pred"])
plot_output_dist(train_result, val_result, sig="tt",savedir=train_outdir)
val_result = pd.DataFrame(np.array([y_val.T, pred_val.T[2]]).T, columns=["True", "Pred"])
train_result = pd.DataFrame(np.array([y_train.T, pred_train.T[2]]).T, columns=["True", "Pred"])
plot_output_dist(train_result, val_result, sig="st",savedir=train_outdir)
