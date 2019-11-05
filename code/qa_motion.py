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

def get_term(fname, key):
    fname = os.path.basename(fname).split('.')[0]
    components = fname.split('_')
    key = key.split('-')[0] + '-'
    for i, s in enumerate(components):
        if key in s:
            return s.split('-')[1]
    return -1

print(str(datetime.now()) + ": Begin qa_motion.py")

# read in arguments
# subjects
subs = sys.argv[1:]
if len(subs) == 0:
    subs='all'
subs_str = ""

in_dir = FMRIPREP_DIR
out_dir = FMRIPREP_DIR

# column names for translational motion regressors
motion_columns = ["trans_x", "trans_y", "trans_z"]

# list all directories in this folder (subject directories)
dnames = os.path.join(in_dir,
                         "sub-???")
all_fmriprep_dirs = glob.glob(dnames)
print(all_fmriprep_dirs)
# list of subjects who fail motion QA
fail_subs = []
# combine all subjects' confounds tsv files
for d, dir in enumerate(all_fmriprep_dirs):
    sub = os.path.basename(dir)
    print("Analyzing " + sub)
    if sub in subs or subs == "all":
        if sub in EXCLUDE_SUBS:
            print('Excluding subject '+ sub)
        else:
            subs_str += sub + "_"
            sub_dir = os.path.join(in_dir,sub,'func')
            # find all tsv files for this subject
            confound_fnames = glob.glob(os.path.join(sub_dir, '*desc-confounds_regressors.tsv'))
            # count number of runs this subject fails
            sub_fails = 0
            # construct this subject's dataframe
            for i, f in enumerate(confound_fnames):
                print('Reading in file ' + f)
                # get task
                task = get_term(f, 'task')
                # get run
                run = get_term(f, 'run')
                # read in file
                raw_df = pd.read_csv(f, delimiter = '\t')
                # select specific columns
                sel_df = raw_df[motion_columns]
                # get number of rows
                nrow = sel_df.shape[0]

                # setup columns for distances
                sel_df['dist'] = None
                sel_df['mm2.0'] = 0
                sel_df['mm0.5'] = 0
                # calculate movement
                print("Calculating distances")
                diff = np.zeros((1,nrow))
                for m in motion_columns:
                    col = sel_df[m]
                    d = [0]+[j-i for i, j in zip(col[:-1], col[1:])]
                    d = np.square(d)
                    diff += d
                dist = np.sqrt(diff)
                dist = np.transpose(dist)
                sel_df['dist'] = dist
                sel_df['mm2.0'] = dist > 2
                sel_df['mm0.5'] = dist > .5

                sum20 = sum(sel_df['mm2.0'])
                sum05 = sum(sel_df['mm0.5'])
                fail = sum20>0 or sum05 >= 5
                if fail:
                    print("Subject "+sub+" fails on run "+run)
                sub_fails += int(fail)
                run_dict = {'sub':[sub], 'task':[task], 'run':[run],
                          'dist':[np.mean(sel_df['dist'])], 'mm2.0':[sum20],
                          'mm0.5':[sum05], 'run_fail':[fail]}
                if i == 0:
                    sub_df = pd.DataFrame(run_dict)
                else:
                    sub_df = sub_df.append(pd.DataFrame(run_dict))

        print("Subject " + sub + " failed " + str(sub_fails) + " runs.")
        if sub_fails > 4:
            fail_subs += [sub]
            sub_df['sub_fail'] = True
        else:
            sub_df['sub_fail'] = False
        # add subject's dataframe to the full dataframe
        if 'full_df' in locals():
            full_df = full_df.append(sub_df)
        else:
            full_df = deepcopy(sub_df)
# save output
fname = os.path.join(out_dir, subs_str + "motion_counts.csv")
full_df.to_csv(fname)
print("Saved " + fname)
