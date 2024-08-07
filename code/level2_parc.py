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

corrs=['spear', 'reg']
corr_labels = []
tasks_all = TASKS + ['avg']
tasks=[]
dists=[]
# read in arguments
for arg in sys.argv[1:]:
    print(arg)
    if arg in corrs:
        corr_labels += [arg]
    elif arg in tasks_all:
        tasks += [arg]
    elif arg in ['correlation', 'euclidean', 'mahalanobis']:
        dists += arg

if len(corr_labels) == 0:
    corr_labels = corrs

if len(tasks) == 0:
    tasks = tasks_all

if len(dists) == 0:
    dists=['euclidean', 'correlation']

task_sects = [[],[]]
for t in tasks:
    task_sects[t in TASKS] = task_sects[t in TASKS] + [t]

in_dir = get_thresh_dir(RSA_DIR)

for proc in [PARC]:
    print(proc)
    if proc == SL:
        parc_label = SL
    else:
        parc_label = PARC_LAB
    for corr_label in corr_labels:
        print(corr_label)
        if corr_label == 'reg':
            val_labels = ['R2', 'beta']
        else:
            val_labels = ['r']
        for val_label in val_labels:
            print(val_label)
            for dist in dists:
                print(dist)
                print('Starting analysis of csv files matching: '+parc_label+' '+corr_label+' '+val_label+' '+dist)
                out_dir = get_thresh_dir( os.path.join(SECOND_LEVEL_DIR,'parc-' + PARC_LAB, dist) )
                # read in all subject files into single data frame
                for ts in task_sects:
                    if len(ts) == 0:
                        continue
                    data_fnames = []
                    for t in ts:
                        in_fname_atts = OrderedDict({'task':t, 'space':SPACE,'stat':STAT, 'corr':corr_label, 'parc':parc_label, "val":val_label, "extra":"roi_stats"})
                        in_fname = make_bids_str(in_fname_atts)
                        fnames = os.path.join(in_dir, "sub-*", corr_label, dist,
                                                 "*"+in_fname+".csv")
                        d_fnames = glob.glob(fnames)
                        data_fnames += d_fnames
                    ############################################################
                    out_fname_atts = in_fname_atts
                    if len(data_fnames)==0:
                          print('no files found matching: '+fnames)
                          continue
                    print('Reading in files: ')
                    print(data_fnames)
                    all_data = pd.concat((pd.read_csv(f) for f in data_fnames))
                    # drop row numbers
                    #all_data = all_data.iloc[:,1:]
                    orig_data = deepcopy(all_data)
                    ############################################################
                    # get predictors
                    predictors = pd.unique(all_data['predictor'])
                    print("Predictors found: ")
                    print(predictors)
                    ############################################################
                    if len(pd.unique(all_data['task'])) > 1: #in all_data.columns:
                        has_task=True
                        # get tasks
                        tasks = pd.unique(all_data['task'])
                        print("Tasks found: ")
                        print(tasks)
                        if len(tasks)==1:
                            has_task=False
                    else:
                        has_task=False
                    ############################################################
                    # get list of ROIs
                    rois = pd.unique(all_data['roi'])
                    print("Number of ROIs found: ")
                    print(len(rois))
                    ############################################################
                    if has_task:
                        # separate tasks into 2 columns
                        all_data = all_data.pivot_table(index=['sub','roi','predictor'],
                                                             columns='task', values=val_label)
                        # reset header
                        all_data.reset_index(inplace=True)
                        all_data.columns.name = ''
                        ## calculate average across tasks
                        #all_data['avg'] = (all_data['friend'] + all_data['number'])/2
                    else:
                        # for regression output, save values under "avg" heading
                        all_data['avg'] = all_data[val_label]
                    ############################################################
                    # pred_dfs = []
                    # diff_dfs = []
                    # avg_dfs = []
                    cnames = ['corr','test','estimate','pred','roi','t','df','p']
                    for pred in predictors:
                        print('Running one-way t-tests for predictor: '+pred)
                        # subselect predictor's rows
                        pred_df = all_data[all_data['predictor']==pred]
                        # initialize out dataframes
                        con_df = pd.DataFrame(None, columns=cnames)
                        #diff_df = pd.DataFrame(None, columns=cnames)
                        #avg_df = pd.DataFrame(None, columns=cnames)
                        # iterate through parcels
                        for i,r in enumerate(rois):
                            print("Testing parcel "+ str(r)+ " ("+str(i)+"/"+str(len(rois)) +")")
                            df = pred_df[pred_df['roi']==r]
                            # task-specific tests
                            if has_task:
                                # run relevance (diff) test
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
                                con_df.loc[con_df.shape[0]+r]=[corr_label,'diff',diff,pred,r,t,t_df,p]
                                # run one-sample t-tests within each task
                                for task in tasks:
                                    t,p = ttest_1samp(df[task], popmean = 0.)
                                    # correct p-value to one-sample
                                    t,p = correct_p(t,p)
                                    # degrees of freedom
                                    t_df = len(df[task]) - 1
                                    estimate=np.mean(df[task])
                                    con_df.loc[con_df.shape[0]+r]=[corr_label,task,estimate,pred,r,t,t_df,p]
                            else:
                                # test where the predictor is encoded across tasks
                                t,p = ttest_1samp(df['avg'], popmean = 0.)
                                # correct p-value to one-sample
                                t,p = correct_p(t,p)
                                # degrees of freedom
                                t_df = len(df['avg']) - 1
                                estimate=np.mean(df['avg'])
                                con_df.loc[r]=[corr_label,'avg',estimate,pred,r,t,t_df,p]
                        ############################################################
                        # correct full dataframe for multiple comparisons
                        # if has_task:
                        #     sigs, con_df['FDR'], a, b = multipletests(con_df['p'], alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                        sigs, con_df['FDR'], a, b = multipletests(con_df['p'], alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                        print(con_df['FDR'][sigs])
                        ############################################################
                        # save output
                        out_fname_atts['pred'] = pred
                        #out_fname = os.path.join(out_dir, make_bids_str(in_fname_atts)+'.csv')
                        out_fname_atts['dir'] = 'none'
                        if has_task:
                            for t in np.unique(con_df['test']):
                                df = con_df[con_df['test']==t]
                                out_fname_atts['task'] = t
                                if t=="diff":
                                    out_fname_atts['test'] = 'pairedt'
                                else:
                                    out_fname_atts['test'] = 'singlesamplet'
                                # save csv
                                out_fname_atts['val2'] = 'all'
                                del out_fname_atts['dir']
                                out_fname_atts['correction'] = 'all'
                                out_fname = os.path.join(out_dir, make_bids_str(out_fname_atts)+'.csv')
                                df.to_csv(out_fname)
                                print(out_fname + ' saved.')
                                # save niftis
                                out_fname_atts['val2'] = 'p'
                                out_fname_atts['dir'] = 'rev'
                                out_fname_atts['correction'] = 'none'
                                out_fname = os.path.join(out_dir, make_bids_str(out_fname_atts)+'.nii')
                                save_parcel(df['p'], out_fname)
                                out_fname_atts['correction'] = 'fdr'
                                out_fname = os.path.join(out_dir, make_bids_str(out_fname_atts)+'.nii')
                                save_parcel(df['FDR'], out_fname)
                        else:
                            # save csv
                            del out_fname_atts['dir']
                            out_fname_atts['task'] = 'avg'
                            out_fname_atts['test'] = 'singlesamplet'
                            out_fname_atts['val2'] = 'all'
                            out_fname_atts['correction'] = 'all'
                            out_fname = os.path.join(out_dir, make_bids_str(out_fname_atts)+'.csv')
                            con_df.to_csv(out_fname)
                            # save p map niftis
                            out_fname_atts['val2'] = 'p'
                            out_fname_atts['dir'] = 'rev'
                            out_fname_atts['correction'] = 'none'
                            out_fname = os.path.join(out_dir, make_bids_str(out_fname_atts)+'.nii')
                            save_parcel(con_df['p'], out_fname)
                            out_fname_atts['correction'] = 'fdr'
                            out_fname = os.path.join(out_dir, make_bids_str(out_fname_atts)+'.nii')
                            save_parcel(con_df['FDR'], out_fname)

print("DONE: level2_parc.py")
