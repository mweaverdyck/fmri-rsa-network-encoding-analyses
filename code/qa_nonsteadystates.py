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


print(str(datetime.now()) + ": Begin qa_rm_nonsteadystates.py")

# read in arguments
# subjects
subs = sys.argv[1:]
if len(subs) == 0:
    subs='all'

in_dir = FMRIPREP_DIR
out_dir = DERIVS_DIR

# list all directories in this folder (subject directories)
fnames = os.path.join(FMRIPREP_DIR,
                         "sub-???")
all_fmriprep_dirs = glob.glob(fnames)


# combine all subjects' confounds tsv files
all_dfs = pd.DataFrame(None, columns=['nonsteady_rows', 'task', 'run', 'sub'])
for dir in all_fmriprep_dirs:
    sub = os.path.basename(dir)
    if sub in subs or subs == "all":
      if sub in EXCLUDE_SUBS:
        print(paste('Excluding subject ', sub))
      else:
        raw_sub_dir = os.path.join(BIDS_DIR,sub,'func')
        #events_fnames = glob.glob(os.path.join(raw_sub_dir, '*tsv'))
        # find all tsv files for this subject
        confound_fnames = glob.glob(os.path.join(FMRIPREP_DIR,sub,'func','*tsv'))
        # setup this subject's dataframe
        sub_df = pd.DataFrame(None, columns=['nonsteady_rows', 'task', 'run', 'sub'])
        for f in confound_fnames:
          # get base filename
          confound_basef = f.split('/')[-1]
          t = get_bids_att(confound_basef, 'task')
          r = get_bids_att(confound_basef, 'run')
          #base_f = confound_basef[0:25]
          # read in this file
          confound_df = pd.read_csv(f, delimiter = '\t') #, na.strings = 'n/a')
          confound_df_out = confound_df.iloc[N_TRS_DEL:,:]
          out_fname_components = get_all_bids_atts(confound_basef)
          out_fname_components['rmtr'] = str(N_TRS_DEL)
          if 'extra' in out_fname_components:
              out_fname_components.move_to_end('extra', last=True)
          out_fname = os.path.join(out_dir,
                                    make_bids_str(out_fname_components) +
                                    '.' + confound_basef.split('.')[-1])
          confound_df_out.to_csv(out_fname, sep='\t')
          print(out_fname + " saved")
          #out_fname = paste0(substr(confound_basef, 1, 41), '_rmtr-', as.character(N_TRS_DEL), substr(confound_basef, 42, 56))
          #readr::write_tsv(confound_df_out, path = paste0(out_dir,'/',sub,'/func/', out_fname))
          #print(paste(paste0(out_dir,'/',sub,'/func/', out_fname),"saved"))


          # read in events file
          events_f = glob.glob(os.path.join(raw_sub_dir, sub+'*task-'+t+'*run-'+r+'*tsv'))
          if len(events_f) != 1:
              print("ERROR: "+str(len(events_f))+ " event files found: ")
              print(events_f)
          else:
              events_f = events_f[0]
          print('Reading in file: '+ events_f)
          #events_f = os.path.join(raw_sub_dir,events_fnames+'*')
          events_df = pd.read_csv(events_f, delimiter='\t')

          # # subselect movement parameters
          # confound_df = confound_df.filter(regex='^non_steady_state_outlier',axis=1)
          # nonsteady_rows = []
          # nonsteady_time = 0
          # # if nonsteady state columns found...
          # if confound_df.shape[1] > 0:
          #   warn = False
          #   # iterate through columns
          #   for i in range(confound_df.shape[1]):
          #     r = which(grepl(1, confound_df[,i]))
          #     nonsteady_rows = nonsteady_rows + [r]
          #     # check if consecutive
          #     if (i == 1){
          #       if (r != 1){
          #         print('WARNING: First non-steady state is not first TR!')
          #         print(f)
          #         warn = T
          #       }
          #     } else {
          #       r_prev = nonsteady_rows[i-1]
          #       if (r != r_prev + 1){
          #         print('WARNING: Non-steady states are not consecutive!')
          #         print(f)
          #         warn = T
          #       }
          #     }
          #     # get time
          #     nonsteady_time = max(nonsteady_rows) * TR
          #   }
          #   c_df = data.frame(nonsteady_rows, task, run, sub)
          #   # combine this file's data to subject's dataframe
          #   sub_df = rbind(sub_df, c_df)
          # }
          # # recalculate onsets
          # events_df['onset_corrected_by_sub'] = events_df['onset'] - nonsteady_time

          events_df['onset_corrected'] = events_df['onset'] - (N_TRS_DEL * TR)
          events_df.to_csv(events_f, sep='\t')
          print(events_f + " saved")

        # # append subject's dataframe to all_dfs
        # all_dfs = pd.concat([all_dfs, sub_df])
        
# fname=os.path.join(FMRIPREP_DIR,'nonsteady_outliers.tsv')
# write_tsv(all_dfs, path=fname)
# print(paste(fname,"saved"))

print("Done: qa_nonsteadystates.py")
