
import os
import shutil
import sys
from datetime import datetime
from copy import deepcopy
import glob
import pandas as pd
import csv
from funcs import *

out_fname_all = os.path.join(SECOND_LEVEL_DIR,
    'parc-'+str(N_PARCELS)+'_test-avg_stat-'+STAT+'_corr-spear_desc-fdr.csv')

out_fname_sig = os.path.join(SECOND_LEVEL_DIR,
    'parc-'+str(N_PARCELS)+'_test-avg_stat-'+STAT+'_corr-spear_desc-fdrsig.csv')

fname_reg = os.path.join(SECOND_LEVEL_DIR,
    'parc-'+str(N_PARCELS)+'_pred-*_cat-*_test-avg_stat-'+STAT+'_corr-spear.csv')
fnames = glob.glob(fname_reg)

li=[]
preds = []
cats = []
for f in fnames:
    p = get_bids_att(f, 'pred')
    preds += [p]
    c = get_bids_att(f, 'cat')
    cats += [c]
    df = pd.read_csv(f)
    li.append(df)

dfs_all = pd.concat(li, axis=0, ignore_index=True)
# new data frame with split value columns
new = dfs_all["pred"].str.split("_", n = 1, expand = True)
pred_vals = new[0]
cat_vals = new[1].str.split("-", n=1, expand=True)
# Dropping old Name columns
dfs_all.drop(columns =["pred"], inplace = True)
# making separate first name column from new data frame
dfs_all["pred"]= pred_vals
# making separate last name column from new data frame
dfs_all["cat"]= cat_vals[1]

schaefer_labels_fname = os.path.join(SECOND_LEVEL_DIR, 'Schaefer2018_'+str(N_PARCELS)+'Parcels_'+str(N_NETWORKS)+'Networks_order.txt')
lab_df = pd.read_csv(schaefer_labels_fname, sep='\t', names=["roi", "label", "x", "y", "z", "unknown"])
labels = lab_df["label"].str.split("_",expand=True)
labels.columns = ["n_networks", "hem", "network", "area", "sub_area"]
labels['roi'] = lab_df['roi']
labels = labels.set_index('roi')
# select only ROIs that show at least one result of FDR<.1
dfs_sig = dfs_all.loc[dfs_all["FDR"]<.1]
fdr_tbl_sig = dfs_sig.pivot(index='roi', columns='pred', values='FDR')

lab_sig = labels[labels.index.isin(fdr_tbl_sig.index)]
for c in labels.columns:
    fdr_tbl_sig[c] = labels[c]

fdr_tbl_sig.to_csv(out_fname_sig)

fdr_tbl = dfs_all.pivot(index='roi', columns='pred', values='FDR')
for c in labels.columns:
    fdr_tbl[c] = labels[c]

fdr_tbl.to_csv(out_fname_all)
