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

parc_label = str(N_PARCELS)

if len(sys.argv) > 1:
    corr_labels = sys.argv[1:]
else:
    corr_labels = ['spear', 'reg']

in_dir = RSA_DIR
out_dir = SECOND_LEVEL_DIR

for corr_label in corr_labels:
    print('Starting analysis of '+corr_label+' csv files')
    # read in all subject files into single data frame
    out_fname = "stat-"+STAT+"_corr-"+ corr_label+"_parc-"+ parc_label+"_tail-1"
    fnames = os.path.join(in_dir, "sub-*", corr_label, "parc",
                             "*_stat-"+ STAT+
                             "_corr-"+ corr_label+
                             "_parc-"+ parc_label+
                             "_roi_stats.csv")
    data_fnames = glob.glob(fnames)

    print('Reading in files: ')
    print(data_fnames)
    if len(data_fnames)==0:
          print('no files found matching: '+fnames)
          continue
    all_data = pd.concat((pd.read_csv(f) for f in data_fnames))
    # drop row numbers
    all_data = all_data.iloc[:,1:]
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
    else:
        has_task=False

    # get list of ROIs
    rois = pd.unique(all_data['roi'])
    print("Number of ROIs found: ")
    print(len(rois))
    if len(rois) != N_PARCELS:
      print("WARNING: number of ROIs does not match N_PARCELS:")
      print(str(len(rois)), parc_label)

    if has_task:
        # separate tasks into 2 columns
        all_data = all_data.pivot_table(index=['sub','roi','predictor'],
                                             columns='task', values='r')
        # reset header
        all_data.reset_index(inplace=True)
        all_data.columns.name = ''
        # calculate average across tasks
        all_data['avg'] = (all_data['friend'] + all_data['number'])/2
    else:
        # for regression output, save beta values under "avg" heading
        all_data['avg'] = all_data['beta']

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
                if pred=='deg':
                    rel_col = 'number'
                    nonrel_col = 'friend'
                elif pred=='dist':
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
            diff_df['FDR'] = FDR(np.array(diff_df['p']))
        avg_df['FDR'] = FDR(np.array(avg_df['p']))

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
