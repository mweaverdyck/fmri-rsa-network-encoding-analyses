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
out_dir = DERIVATIVES_DIR

# list all directories in this folder (subject directories)
fnames = os.path.join(FMRIPREP_DIR,
                         "sub-???")
all_fmriprep_dirs = glob.glob(fnames)


# combine all subjects' confounds tsv files
all_dfs = pd.DataFrame(None, columns=['nonsteady_rows', 'task', 'run', 'sub'])
for sub in all_fmriprep_dirs:
    if sub in subs or subs == "all":
      if sub in EXCLUDE_SUBS:
        print(paste('Excluding subject ', sub))
      else:
        raw_sub_dir = os.path.join(BIDS_DIR,sub,'func')
        events_fnames = glob.glob(os.path.join(raw_sub_dir, '*tsv')
        # find all tsv files for this subject
        confound_fnames = glob.glob(os.path.join(FMRIPREP_DIR,sub,'func','*tsv')
        # setup this subject's dataframe
        sub_df=pd.DataFrame(None, columns=['nonsteady_rows', 'task', 'run', 'sub'])
        for f in confound_fnames:
          # get base filename
          confound_basef = f.split('/')[-1]
          base_f = confound_basef[0:25]
          print(paste('Reading in files starting with:', base_f))
          events_f = os.path.join(raw_sub_dir,events_fnames+'*')
          events_csv = pd.read_csv(events_f, delimiter='\t')
          # index for task name in filename
          base_f_comps = base_f.split("-")
          task_i = base_f_comps.index('task')+1
          # task name
          task = base_f_comps[task_i]
          # index for run number in filename
          run_i = base_f_comps.index('run')+1
          # run number
          run = base_f_comps[run_i]

          # read in this file
          confound_csv = pd.read_csv(f, delimiter = '\t') #, na.strings = 'n/a')
          confound_csv_out = confound_csv.iloc[N_TRS_DEL:,:]
          out_fname = confound_basef[0:40]+'_rmtr-' + str(N_TRS_DEL) + confound_basef[40:]
          #out_fname = paste0(substr(confound_basef, 1, 41), '_rmtr-', as.character(N_TRS_DEL), substr(confound_basef, 42, 56))
          readr::write_tsv(confound_csv_out, path = paste0(out_dir,'/',sub,'/func/', out_fname))
          print(paste(paste0(out_dir,'/',sub,'/func/', out_fname),"saved"))

          # get number of non-steady state outliers
          nonsteady = grepl( "non_steady_state_outlier" , names( confound_csv ) )
          # subselect movement parameters
          confound_csv = as.data.frame(confound_csv[, nonsteady])
          nonsteady_rows = c()
          nonsteady_time = 0
          # if nonsteady state columns found...
          if (ncol(confound_csv) > 0){
            warn = F
            # iterate through columns
            for (i in seq(1,ncol(confound_csv))){
              r = which(grepl(1, confound_csv[,i]))
              nonsteady_rows = c(nonsteady_rows, r)
              # check if consecutive
              if (i == 1){
                if (r != 1){
                  print('WARNING: First non-steady state is not first TR!')
                  print(f)
                  warn = T
                }
              } else {
                r_prev = nonsteady_rows[i-1]
                if (r != r_prev + 1){
                  print('WARNING: Non-steady states are not consecutive!')
                  print(f)
                  warn = T
                }
              }
              # get time
              nonsteady_time = max(nonsteady_rows) * TR
            }
            c_df = data.frame(nonsteady_rows, task, run, sub)
            # combine this file's data to subject's dataframe
            sub_df = rbind(sub_df, c_df)
          }
          # recalculate onsets
          events_csv$onset_corrected_by_sub = events_csv$onset - nonsteady_time
          events_csv$onset_corrected = events_csv$onset - (N_TRS_DEL * TR)
          write_tsv(events_csv, path = events_f)
          print(paste(events_f,"saved"))
        }
        # append subject's dataframe to all_dfs
        all_dfs = rbind(all_dfs, sub_df)
      }
    }
  }
}
fname=paste0(FMRIPREP_DIR,'/nonsteady_outliers.tsv')
write_tsv(all_dfs, path=fname)
print(paste(fname,"saved"))
