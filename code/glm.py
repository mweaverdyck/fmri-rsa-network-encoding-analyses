import sys
import os
import shutil
import numpy as np
import pandas as pd
import nistats
import copy
#import nistats.first_level_model
from nistats.first_level_model import * #first_level_models_from_bids
#from nilearn import plotting
from nistats.reporting import plot_design_matrix
import nibabel as nib
import matplotlib.pyplot as plt
from funcs import *


# input arguments
try:
    if len(sys.argv) < 2:
        raise RuntimeError()
    subject_ids = sys.argv[1:]
    if len(subject_ids) == 1 and subject_ids[0] == '--all':
        subject_ids = None
except RuntimeError:
    print('Usage:\n\t\tpython nistats_glm.py <subject_id_1> <subject_id_2> ...'
          '\n\tOR\n\t\tpython nistats_glm.py --all')
    quit()


derivatives_prefix = DERIVS_DIR + 'derivatives_'
space_label = SPACE
tasks = ['friend','number']
nodes = range(N_NODES)

out_dir = GLM_DIR
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# ----------------------- fmriprep 1.2.0 regressors ----------------------------
# col_regressors_fixed = [
# 'a_comp_cor_00', 'a_comp_cor_01',
# 'a_comp_cor_02', 'a_comp_cor_03',
# 'a_comp_cor_04', 'a_comp_cor_05',
# 't_comp_cor_00', 't_comp_cor_01',
# 't_comp_cor_02', 't_comp_cor_03',
# 't_comp_cor_04', 't_comp_cor_05',
# 'trans_x','trans_y','trans_z',
# 'rot_x','rot_y','rot_z']
# # include all columns that start with these prefixes
# col_regressors_prefs_all = []
# # across all 4 runs of the task, include the minimum number of columns that start with these prefixes
# col_regressors_prefs_min = ['cosine', 'non_steady_state_outlier']
# # delete TRs based on these columns

# ----------------------- fmriprep 1.4.0 regressors ----------------------------
col_regressors_fixed = [
'csf', 'white_matter', 'global_signal',
'trans_x','trans_y','trans_z',
'rot_x','rot_y','rot_z']
# include all columns that start with these prefixes
col_regressors_prefs_all = col_regressors_fixed
# across all 4 runs of the task, include the minimum number of columns that start with these prefixes
col_regressors_prefs_min = ['non_steady_state_outlier']
# delete TRs based on these columns



for subj in subject_ids:
    print("Starting subject "+str(subj))
    derivatives_dir_sub = derivatives_prefix + subj + '/'
    if not os.path.exists(derivatives_dir_sub):
        os.makedirs(derivatives_dir_sub)
    if not os.path.exists(derivatives_dir_sub+subj):
        # move derivatives to derivatives folder
        print("moving " +DERIVS_DIR+subj+' to '+derivatives_dir_sub)
        shutil.move(DERIVS_DIR+subj, derivatives_dir_sub)

    # will run on all subjects in derivatives_dir_sub, so all must be done with fmriprep
    for t, task_label in enumerate(tasks):
        print("Starting task "+str(t+1)+" of "+str(len(tasks))+": " + task_label)
        print("Calculating first level model...")
        models, models_run_imgs, models_events, models_confounds = \
            first_level_models_from_bids(
                BIDS_DIR, task_label, space_label
                , derivatives_folder=derivatives_dir_sub
            )

        if sys.version_info[0] == 3:
            # note, need list() around zip in python 3 but not in python 2
            model_and_args = list(zip(models, models_run_imgs, models_events, models_confounds))
        elif sys.version_info[0] == 2:
            model_and_args = zip(models, models_run_imgs, models_events, models_confounds)

        print(len(model_and_args))
        # iterate through subjects
        for midx, (model, imgs, events, confounds) in enumerate(model_and_args):
            print(midx)
            # do not mask
            model.mask=False
            # get subject label
            subnum = model.subject_label
            subid = 'sub-' + subnum
            print("Analyzing subject "+ subid)

            # write subject's output directory
            write_sub_dir = os.path.join(out_dir, subid)
            if not os.path.exists(write_sub_dir):
                os.makedirs(write_sub_dir)

            # select only relevant rows and columns from events files
            print("Setting up events file")
            for r, e in enumerate(events):
                events[r] = e[e['trial_type']=='noncatch']
                events[r] = events[r].filter(['node','onset_corrected','duration'])
                events[r] = events[r].rename(index=str,
                #    columns={"trial_type":"noncatch", "node":"trial_type"})
                    columns={"node":"trial_type", "onset_corrected":"onset"})
                print(events[r].columns)
                #events[r] = events[r].filter(['trial_type','onset_corrected','duration'])
            	#events[r] = events[r].rename(columns = {"onset_corrected": "onset"})

            # subselect regressors from each run's counfounds file
            print("Subselecting relevant regressors")
            # find the minimum prefix columns
            # iterate through this subject's counfounds files
            for r, c in enumerate(confounds):
                # r = run number, c = confounds dataframe
                curr = []
                for p in col_regressors_prefs_min:
                    add_cols = list(c.filter(regex='^'+p,axis=1).columns)
                    print("Adding columns for prefix " + p + ":")
                    print(add_cols)
                    curr = curr + add_cols
                # only include columns that have been found in all confounds files
                if r == 0:
                    col_regressors_min = copy.deepcopy(curr)
                else:
                    col_regressors_min = [x for x in curr if x in col_regressors_min]
            # subselect columns
            confounds_copy=copy.deepcopy(confounds)
            for r, c in enumerate(confounds):
                # find columns that start with prefixes
                col_regressors_all = copy.deepcopy(col_regressors_fixed)
                col_regressors_all = col_regressors_all + col_regressors_min
                for p in col_regressors_prefs_all:
                    add_cols = list(c.filter(regex='^'+p,axis=1).columns)
                    print("Adding columns for prefix " + p + ":")
                    print(add_cols)
                    col_regressors_all = col_regressors_all + add_cols
                print("Selected confounds for run " + str(r+1) + ":")
                print(col_regressors_all)
                confounds_copy[r] = c.filter(col_regressors_all)
                # fill in all NaN with 0
                confounds_copy[r] = confounds_copy[r].fillna(0)

            # fit the GLM
            print("Running GLM")
            #model1=deepcopy(model)
            model.fit(imgs, events, confounds_copy)

            # save design matrices for each run
            print("Saving design matrices")
            for r, mat in enumerate(model.design_matrices_):
                design_matrix = mat
                plot_design_matrix(design_matrix)
                #plt.show()
                filename=write_sub_dir+'/'+subid+'_task-'+task_label+'_run-0'+str(r+1)+'_designmat.png'
                plt.savefig(filename)
                print('Run '+str(r)+' design matrix image saved: '+filename)

            # setup contrast vectors (all zeros)
            n_columns = design_matrix.shape[1]
            con_empty = np.zeros(n_columns)

            # compute each contrast of interest (one per node)
            print("Saving each node's contrast files")
            for node_label in nodes:
                # create contrast vector
                con = copy.deepcopy(con_empty)
                con[node_label] = 1
                # for output options, see
                # https://nistats.github.io/modules/generated/nistats.first_level_model.FirstLevelModel.html
                t_map = model.compute_contrast(con
                                , stat_type='t'
                                , output_type='stat'
                                )
                # save t map
                filename = '%s_task-%s_space-%s_stat-t_node-0%s.nii' % (subid, task_label, space_label, str(node_label))
                t_image_path = os.path.join(write_sub_dir, filename)
                nib.save(t_map, t_image_path)
                print('File ' + write_sub_dir + '/' + filename + ' saved.')


    # move subjects' folders back to fmriprep
    print("Moving subject's data from " + derivatives_dir_sub + " to " + DERIVS_DIR)
    shutil.move(derivatives_dir_sub+subj, DERIVS_DIR)
    os.rmdir(derivatives_dir_sub)
