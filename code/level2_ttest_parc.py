#!/bin/python

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

N_PARCELS = 800
MNI_PARCELLATION=MNI_DIR+'/tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-'+str(N_PARCELS)+'Parcels17Networks_dseg.nii.gz'

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
    parcellation.set_data_dtype(np.double)
    p_data = parcellation.get_data()
    for i in range(N_PARCELS):
        p = i+1
        p_data[ p_data==p ] = parcel_vals[i]
    save_nii(p_data, parcellation, filename)

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

parc_label = '800'#str(N_PARCELS)
corr_label='spear'

out_fname = "stat-"+STAT+"_corr-"+ corr_label+"_parc-"+ parc_label+"_tail-1"


# read in all subject files into single data frame
out_dir = SECOND_LEVEL_DIR
fnames = os.path.join('/u/project/CCN/cparkins/data/encoding/bids/level1/norecon/rsa',#RSA_DIR,
                         "sub-20*/pymvpa/"+
                         "parc/*_stat-"+ STAT+
                         "_corr-"+ corr_label+
                         "_parc-"+ parc_label+
                         "_roi_stats.csv")
data_fnames = glob.glob(fnames)

print('Reading in files: ')
print(data_fnames)
all_data = pd.concat((pd.read_csv(f) for f in data_fnames))
# drop row numbers
all_data = all_data.iloc[:,1:]
orig_data = deepcopy(all_data)

# get predictors
predictors = pd.unique(all_data['predictor'])
print("Predictors found: ")
print(predictors)

# get tasks
tasks = pd.unique(all_data['task'])
print("Tasks found: ")
print(tasks)

# get list of ROIs
rois = pd.unique(all_data['roi'])
print("Number of ROIs found: ")
print(len(rois))
if len(rois) != N_PARCELS:
  print("WARNING: number of ROIs does not match N_PARCELS:")
  print(str(len(rois)), parc_label)

# separate tasks into 2 columns
all_data_task = all_data.pivot_table(index=['sub','roi','predictor'],
                                     columns='task', values='r')
# reset header
all_data_task.reset_index(inplace=True)
all_data_task.columns.name = ''

# calculate average across tasks
all_data_task['avg'] = (all_data_task['friend'] + all_data_task['number'])/2

# separate into 2 dataframes based on predictor
all_data_deg = all_data_task[all_data_task['predictor']=='deg']
all_data_dist = all_data_task[all_data_task['predictor']=='dist']

# run t-tests within each ROI
all_diff_deg = pd.DataFrame(None, columns=['test','predictor','roi','t','p'])
all_avg_deg = pd.DataFrame(None, columns=['test','predictor','roi','t','p'])
all_diff_dist = pd.DataFrame(None, columns=['test','predictor','roi','t','p'])
all_avg_dist = pd.DataFrame(None, columns=['test','predictor','roi','t','p'])
for i,r in enumerate(rois):
    print("Testing parcel "+ str(r)+ " ("+str(i)+"/"+str(len(rois)) +")")
    # DEGREE
    df = all_data_deg[all_data_deg['roi']==r]
    # run two-sample t-test
    t,p = ttest_rel(df['number'], df['friend'])
    # correct to one-sample
    t,p = correct_p(t,p)
    diff = df['number'] - df['friend']
    all_diff_deg.loc[r]=[diff,'deg',r,t,p]
    # test where degree is encoded across tasks
    t,p = ttest_1samp(df['avg'], popmean = 0.)
    # correct to one-sample
    t,p = correct_p(t,p)
    all_avg_deg.loc[r]=[df['avg'],'deg',r,t,p]
    # DISTANCE
    df = all_data_dist[all_data_dist['roi']==r]
    # run two-sample t-test
    t,p = ttest_rel(df['friend'], df['number'])
    # correct to one-sample
    t,p = correct_p(t,p)
    diff = df['friend'] - df['number']
    all_diff_dist.loc[r]=[diff,'dist',r,t,p]
    # test where degree is encoded across tasks
    t,p = ttest_1samp(df['avg'], popmean = 0.)
    # correct to one-sample
    t,p = correct_p(t,p)
    all_avg_dist.loc[r]=[df['avg'],'dist',r,t,p]

# correct for multiple comparisons
all_diff_deg['FDR'] = FDR(np.array(all_avg_dist['p']))
all_avg_deg['FDR'] = FDR(np.array(all_avg_dist['p']))
all_diff_dist['FDR'] = FDR(np.array(all_avg_dist['p']))
all_avg_dist['FDR'] = FDR(np.array(all_avg_dist['p']))

# save images and csv files
fname = out_dir+'/parc-'+parc_label+'_pred-%s_test-%s_stat-'+STAT+'_corr-'+corr_label
fname_csv = fname+'.csv'
all_diff_deg.to_csv(fname_csv % ('deg','diff'))
save_parcel(all_diff_deg['FDR'], fname % ('deg','diff'))
all_avg_deg.to_csv(fname_csv % ('deg','avg'))
save_parcel(all_avg_deg['FDR'], fname % ('deg','avg'))
all_diff_dist.to_csv(fname_csv % ('dist','diff'))
save_parcel(all_diff_dist['FDR'], fname % ('dist','diff'))
all_avg_dist.to_csv(fname_csv % ('dist','avg'))
save_parcel(all_avg_dist['FDR'], fname % ('dist','avg'))
