# searchlight
from nistats.second_level_model import SecondLevelModel
import os
import shutil
import sys
from datetime import datetime
from copy import deepcopy
import glob
import pandas as pd
import nibabel as nib
import numpy as np
from scipy.stats import norm, ttest_rel, ttest_ind, ttest_1samp
from statsmodels.stats.multitest import multipletests
from scipy.ndimage.morphology import binary_dilation
from nilearn.image import math_img
from nilearn.input_data import NiftiMasker
from nistats import second_level_model
from nilearn import plotting
from nilearn.plotting import plot_glass_brain
#from nistats.second_level_model import non_parametric_inference
from funcs import *

def run_sig_tests(data_fnames, mask=None):
    second_level_model = SecondLevelModel(smoothing_fwhm=5.0, mask=mask)
    if isinstance(data_fnames,dict):
        data_fnames_num = sorted(data_fnames['number'])
        data_fnames_fri = sorted(data_fnames['friend'])
        fname_atts = get_all_bids_atts(data_fnames_num[0])
        fname_atts['task'] = 'diff'
        # create paired t-test design matrix
        pair_mat = np.zeros((2*len(data_fnames_num),len(data_fnames_num)+1), dtype=int)
        labs = []
        for i in range(len(data_fnames_num)):
            l = 's'+str(i)
            labs = labs + [l]
            a = [0]*len(data_fnames_num)
            a[i] = 1
            pair_mat[:,i] = a + a
        pair_mat[:, len(data_fnames_num)] = [1] * len(data_fnames_num) + [-1]*len(data_fnames_fri)
        design_matrix = pd.DataFrame(pair_mat,
                                 columns=labs + ['diff'])
        if fname_atts['pred'] == 'deg':
            data_fnames = data_fnames_num + data_fnames_fri
        elif fname_atts['pred'] == 'dist':
            data_fnames = data_fnames_fri + data_fnames_num
        print(data_fnames)
        con_name = 'diff'
    else:
        fname_atts = get_all_bids_atts(data_fnames[0])
        design_matrix = pd.DataFrame([1] * len(data_fnames),
                                 columns=['intercept'])
        con_name = 'intercept'
    # show data files and design matrix
    print("Running significance testing on: ")
    print(data_fnames)
    print("Using design matrix: ")
    print(design_matrix)
    # setup file names
    del fname_atts['sub']
    fname_atts['val2'] = "z"
    fname_atts['correction'] = "none"
    if 'extra' in fname_atts.keys():
        fname_atts.move_to_end('extra')
    # save z map
    second_level_model = second_level_model.fit(data_fnames,
                           design_matrix=design_matrix)
    z_map = second_level_model.compute_contrast(con_name, output_type='z_score')
    out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
    nib.save(z_map, out_fname)
    print(out_fname + ' saved.')
    threshold = 2.88 #3.1  # correponds to  p < .001, uncorrected
    display = plotting.plot_glass_brain(
        z_map, threshold=threshold, colorbar=True, plot_abs=False,
        output_file=out_fname,
        title='z map')
    # save p map
    p_val = second_level_model.compute_contrast(con_name, output_type='p_value')
    fname_atts['val2'] = "p"
    out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
    nib.save(p_val, out_fname)
    print(out_fname + ' saved.')
    # correct for multiple comparisons
    # Correcting the p-values for multiple testing and taking negative logarithm
    n_voxels = np.sum(second_level_model.masker_.mask_img_.get_data())
    neg_log_pval = math_img("-np.log10(np.minimum(1, img * {}))"
            .format(str(n_voxels)),
            img=p_val)
    fname_atts['val2'] = "p"
    fname_atts['correction'] = "parametric"
    out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
    nib.save(neg_log_pval, out_fname)
    print(out_fname + ' saved.')
    # FDR correction
    p_arr = p_val.get_data().flatten()
    sigs, fdr_val, a, b = multipletests(p_arr, alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fname_atts['val2'] = "p"
    fname_atts['correction'] = "fdr"
    out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
    save_nii(data = fdr_val.reshape(p_val.shape), refnii=p_val, filename=out_fname)
    pfdr = deepcopy(1-fdr_val)
    pfdr[fdr_val==0.] = 0.
    pfdr = pfdr.reshape(p_val.shape)
    pfdr_map = nib.Nifti1Image(pfdr, p_val.affine, p_val.header)
    threshold = .95
    display = plotting.plot_glass_brain(
        pfdr_map, threshold=threshold, colorbar=True, plot_abs=False,
        output_file=out_fname,
        title='p map FDR-corrected')
    # neg_log_pvals_permuted_ols_unmasked = \
    #                     non_parametric_inference(second_level_input,
    #                                  design_matrix=design_matrix,
    #                                  model_intercept=True, n_perm=1000,
    #                                  two_sided_test=False,
    #                                  smoothing_fwhm=8.0, n_jobs=1)
    # fname_atts['val2'] = "p"
    # fname_atts['correction'] = "permutation"
    # out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
    # save_nii(neg_log_pvals_permuted_ols_unmasked, ref_img, out_fname)

parc_label = 'sl' + SL_RADIUS
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
    corr_labels = ['spear']

if len(tasks) == 0:
    tasks = TASKS + ['avg']

task_sects = [[],[]]
for t in tasks:
    task_sects[t in TASKS] = task_sects[t in TASKS] + [t]

out_dir = SECOND_LEVEL_DIR

# dilate gray matter mask
tmp = nib.load(MNI_GM_MASK)
gm_mask = load_nii(MNI_GM_MASK)
gm_mask = np.where(gm_mask > 0.5, 1., 0.)
gm_mask_dil = binary_dilation(gm_mask, iterations=int(SL_RADIUS)).astype(gm_mask.dtype)
gm_mask_dil_img = nib.Nifti1Image(gm_mask_dil, tmp.affine, tmp.header)

for corr_label in corr_labels:
    print('starting correlation '+corr_label)
    in_dir = os.path.join(RSA_DIR,'*',corr_label)
    if SPACE=='T1w':
        in_dir = os.path.join(in_dir,'T1w-2-MNI')
    for pred in ['deg', 'dist']:
        print('starting predictor '+pred)
        data_tasks = {}
        for task in tasks:
            print("starting task "+task)
            fnames = os.path.join(in_dir, '*task-'+task+'*space-'+SPACE+'*_parc-'+parc_label+'_*val-r_*pred-'+pred+'*.nii*')
            print(fnames)
            data_fnames = glob.glob(fnames)
            if len(data_fnames) != 0:
                if task in TASKS:
                    data_tasks[task] = data_fnames
                    print(data_tasks)
                run_sig_tests(data_fnames, mask = gm_mask_dil_img)
        if pred in ['deg', 'dist'] and 'number' in data_tasks.keys() and 'friend' in data_tasks.keys():
            run_sig_tests(data_tasks, mask = gm_mask_dil_img)

print(str(datetime.now()) + ": End level2_rsa_sl.py")
