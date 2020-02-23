#!/bin/python3

import os
import shutil
import sys
from datetime import datetime
from copy import deepcopy
import glob
import pandas as pd
import nibabel as nib
import numpy as np
from scipy.stats import ttest_rel, ttest_ind, ttest_1samp
from statsmodels.stats.multitest import multipletests
from funcs import *

def load_nii( imgpath ):
    img4D = nib.load(imgpath, mmap=False).get_data()
    return img4D

def save_nii( data, refnii, filename ):
    out_img = nib.Nifti1Image(data, refnii.affine, refnii.header)
    out_img.to_filename(filename)
    print(str(datetime.now()) + ": File %s saved." % filename)

def save_parcel( parcel_vals, filename, parcellation=MNI_PARCELLATION ):
    parcel_vals = np.array(parcel_vals)
    # MNI parcellationlation
    parcellation = nib.load(MNI_PARCELLATION, mmap=False)
    parcellation.set_data_dtype(np.float)
    p_data = parcellation.get_data()
    out_data = deepcopy(p_data) * 1.0
    for i in range(N_PARCELS):
        p = (i+1)*1.0
        out_val = 1-parcel_vals[i]
        out_data[ p_data==p ] = out_val
    save_nii(out_data, parcellation, filename)

def correct_p (t,p):
    if t < 0:
        # if negative t-value, then relevant context was lower than non-relavant
        p = 1.
    else:
        p = p/2
    return(t,p)

def FDR( x, n=None ):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R
    """
    if n is None:
        n = len(x)
    o = [ i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True) ]
    ro = [ i[0] for i in sorted(enumerate(o), key=lambda v:v[1]) ]
    q = sum([1.0/i for i in range(1,n+1)])
    l = [ q*len(x)/i*x[j] for i,j in zip(reversed(range(1,n+1)),o) ]
    l = [ l[k] if l[k] < 1.0 else 1.0 for k in ro ]
    return l


corrs=['spear', 'reg']
corr_labels = []
tasks=[]
# read in arguments
for arg in sys.argv[1:]:
    print(arg)
    if arg in corrs:
        corr_labels += [arg]
    elif arg in TASKS or arg == 'avg':
        tasks += [arg]

if len(corr_labels) == 0:
    corr_labels = corrs

if len(tasks) == 0:
    tasks = TASKS + ['avg']

task_sects = [[],[]]
for t in tasks:
    task_sects[t in TASKS] = task_sects[t in TASKS] + [t]



in_dir = RSA_DIR
out_dir = SECOND_LEVEL_DIR

for proc in PROCEDURES:
    if proc == SL:
        parc_label = SL
    else:
        parc_label = str(N_PARCELS)
    for corr_label in corr_labels:
        if corr_label == 'reg':
            val_labels = ['R2', 'beta']
        else:
            val_labels = ['r']
        for val_label in val_labels:
            print('Starting analysis of csv files matching: '+parc_label+' '+corr_label+' '+val_label)
            out_fname = "_stat-"+STAT+"_corr-"+ corr_label+"_parc-"+ parc_label+"val-"+val_label+"_tail-1"
            # read in all subject files into single data frame
            for ts in task_sects:
                if len(ts) == 0:
                    continue
                data_fnames = []
                for t in ts:
                    fnames = os.path.join(in_dir, "sub-*", corr_label,
                                             "*task-"+t+"*_stat-"+ STAT+
                                             "_corr-"+ corr_label+
                                             "_parc-"+ parc_label+
                                             "_roi_stats.csv")
                    d_fnames = glob.glob(fnames)
                    data_fnames += d_fnames

                print('Reading in files: ')
                print(data_fnames)
                if len(data_fnames)==0:
                      print('no files found matching: '+fnames)
                      continue
                all_data = pd.concat((pd.read_csv(f) for f in data_fnames))
                # drop row numbers
                #all_data = all_data.iloc[:,1:]
                orig_data = deepcopy(all_data)

                # get predictors
                predictors = pd.unique(all_data['predictor'])
                print("Predictors found: ")
                print(predictors)

                if 'task' in all_data.columns:
                    has_task=True
                    # get tasks
                    tasks = pd.unique(all_data['task'])
                    print("Tasks found: ")
                    print(tasks)
                    if len(tasks)==1:
                        has_task=False
                else:
                    has_task=False

                # get list of ROIs
                rois = pd.unique(all_data['roi'])
                print("Number of ROIs found: ")
                print(len(rois))

                if has_task:
                    # separate tasks into 2 columns
                    all_data = all_data.pivot_table(index=['sub','roi','predictor'],
                                                         columns='task', values=val_label)
                    # reset header
                    all_data.reset_index(inplace=True)
                    all_data.columns.name = ''
                    # calculate average across tasks
                    all_data['avg'] = (all_data['friend'] + all_data['number'])/2
                else:
                    # for regression output, save values under "avg" heading
                    all_data['avg'] = all_data[val_label]

                pred_dfs = []
                diff_dfs = []
                avg_dfs = []
                cnames = ['corr','test','estimate','pred','roi','t','df','p']
                for pred in predictors:
                    print('Running one-way t-tests for predictor: '+pred)
                    # subselect predictor's rows
                    pred_df = all_data[all_data['predictor']==pred]
                    # initialize out dataframes
                    diff_df = pd.DataFrame(None, columns=cnames)
                    avg_df = pd.DataFrame(None, columns=cnames)
                    # iterate through parcels
                    for i,r in enumerate(rois):
                        print("Testing parcel "+ str(r)+ " ("+str(i)+"/"+str(len(rois)) +")")
                        df = pred_df[pred_df['roi']==r]
                        # task-specific tests
                        if has_task:
                            if 'deg' in pred:
                                rel_col = 'number'
                                nonrel_col = 'friend'
                            elif 'dist' in pred:
                                rel_col = 'friend'
                                nonrel_col = 'number'
                            # run two-sample t-test
                            t,p = ttest_rel(df[rel_col], df[nonrel_col])
                            # correct to one-sample
                            t,p = correct_p(t,p)
                            # degrees of freedom
                            t_df = len(df[rel_col]) - 1
                            diff = np.mean(df[rel_col] - df[nonrel_col])
                            diff_df.loc[r]=[corr_label,'diff',diff,pred,r,t,t_df,p]
                        # test where the predictor is encoded across tasks
                        t,p = ttest_1samp(df['avg'], popmean = 0.)
                        # correct p-value to one-sample
                        t,p = correct_p(t,p)
                        # degrees of freedom
                        t_df = len(df['avg']) - 1
                        estimate=np.mean(df['avg'])
                        avg_df.loc[r]=[corr_label,'avg',estimate,pred,r,t,t_df,p]

                    # correct full dataframe for multiple comparisons
                    if has_task:
                        sigs, diff_df['FDR'], a, b = multipletests(diff_df['p'], alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                    sigs, avg_df['FDR'], a, b = multipletests(avg_df['p'], alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                    print(avg_df['FDR'][sigs])

                    # save output
                    fname = out_dir+'parc-'+parc_label+'_pred-'+pred+'_test-%s_stat-'+STAT+'_corr-'+corr_label
                    fname_csv = fname+'.csv'
                    fname_nii = fname+'_p-%s'
                    if has_task:
                        diff_df.to_csv(fname_csv % ('diff'))
                        save_parcel(diff_df['FDR'], fname_nii % ('diff','FDR'))
                        save_parcel(diff_df['p'], fname_nii % ('diff','uncorrected'))
                    avg_df.to_csv(fname_csv % ('avg'))
                    save_parcel(avg_df['FDR'], fname_nii % ('avg', 'FDR'))
                    save_parcel(avg_df['p'], fname_nii % ('avg','uncorrected'))
