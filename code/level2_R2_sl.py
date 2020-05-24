# second level tests on searchlight results
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
import math
from scipy.stats import kstest, beta, norm, ttest_rel, ttest_ind, ttest_1samp, rankdata
from statsmodels.stats.multitest import multipletests # DON'T USE BECAUSE NOT NEURO
from scipy.ndimage.morphology import binary_dilation
from scipy.special import comb
from nilearn.image import math_img
from nilearn.input_data import NiftiMasker
from nistats import second_level_model
from nilearn import plotting
from nilearn.plotting import plot_glass_brain, plot_stat_map
from nistats import thresholding
from nistats.thresholding import map_threshold
#from nistats.second_level_model import non_parametric_inference
from funcs import *
from rsa_funcs import cfd_soc, cfd_phys, cfd_img

out_dir = os.path.join(SECOND_LEVEL_DIR)#, "no_dilation")
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def load_all(data_fnames):
    for i,d in enumerate(data_fnames):
        if i == 0:
            a = np.expand_dims(load_nii(d), 3)
        else:
            dexp = np.expand_dims(load_nii(d), 3)
            a = np.concatenate([a,dexp], 3)
    return(a)

def run_R2_sig_tests(data_fnames, mask_nii=None, k = None):
    # calculate null distribution
    if k is None:
        if '_pred-soc_' in data_fnames[0]:
            k = len(cfd_soc)
        elif '_pred-phys_' in data_fnames[0]:
            k = len(cfd_phys)
        elif '_pred-img_' in data_fnames[0]:
            k = len(cfd_img)
        elif '_pred-sn_' in data_fnames[0]:
            k = 3 # Change this if wrong
        else:
            print("ERROR: Unknown number of predictors")
    n = comb(int(N_NODES), 2)
    dist = beta((k-1)/2, (n-1)/2)
    # test each voxel against null distribution
    if mask_nii is None:
        mask = load_nii(data_fnames[0])
        mask = np.ones(mask.shape)
    else:
        mask = mask_nii.get_data()
    n_tests = mask.size
    d_map = np.empty(mask.shape)
    p_val = np.empty(mask.shape)
    indices = np.argwhere(mask) # get all indices that are not equal to 0
    # load all data
    all_data = load_all(data_fnames)
    def run_kstest(dv, null_dist, alternative="greater"):
        return kstest(dv, null_dist.cdf, alternative="greater")
    dp_map = np.apply_along_axis(run_kstest, axis = 3, arr = all_data, null_dist = dist)
    d_map = dp_map[:,:,:,0]
    p_val = dp_map[:,:,:,1]
    return (d_map, p_val, all_data)


def run_sig_tests(data_fnames, mask=None):
    cmap = "Wistia"
    second_level_model = SecondLevelModel(smoothing_fwhm=5.0, mask=mask)
    if isinstance(data_fnames,dict):
        test_type = 'z'
        data_fnames_num = sorted(data_fnames['number'])
        data_fnames_fri = sorted(data_fnames['friend'])
        fname_atts = get_all_bids_atts(data_fnames_num[0])
        fname_atts['task'] = 'diff'
        fname_atts['test'] = 'pairedt'
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
            cmap = "winter"
        elif fname_atts['pred'] == 'dist':
            data_fnames = data_fnames_fri + data_fnames_num
            cmap = "cool"
        print(data_fnames)
        con_name = 'diff'
    else:
        fname_atts = get_all_bids_atts(data_fnames[0])
        # if fname_atts['val'] == 'R2':
        #     test_type = 'ks'
        #     fname_atts['test'] = 'ksfit'
        # else:
        test_type = 'z'
        fname_atts['test'] = 'singlesamplet'
        design_matrix = pd.DataFrame([1] * len(data_fnames),
                             columns=['intercept'])
        con_name = 'intercept'
    # setup file names
    del fname_atts['sub']
    fname_atts['correction'] = "none"
    # show data files and design matrix
    print("Running significance testing on: ")
    print(data_fnames)
    if test_type == 'z':
        print("Using design matrix: ")
        print(design_matrix)
        out_stat = "z"
        # save z map
        second_level_model = second_level_model.fit(data_fnames,
                               design_matrix=design_matrix)
        z_map = second_level_model.compute_contrast(con_name, output_type='z_score')
        out_map = z_map
        # save p map
        p_val = second_level_model.compute_contrast(con_name, output_type='p_value')
        refnii = p_val
    elif test_type == 'ks':
        out_stat = "D"
        d_map, p_map, all_data = run_R2_sig_tests(data_fnames, mask_nii=mask)
        refnii = nib.load(data_fnames[0])
        out_map = nib.Nifti1Image(d_map, refnii.affine, refnii.header)
        p_val = nib.Nifti1Image(p_map, refnii.affine, refnii.header)
    # save R2 mean
    if fname_atts['val'] == 'R2':
        all_data = load_all(data_fnames)
        mean_data = np.mean(all_data, axis = 3)
        fname_atts['val2'] = 'mean'
        out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
        save_nii(mean_data, refnii, out_fname)
    # save stat map
    fname_atts['val2'] = out_stat
    out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
    nib.save(out_map, out_fname)
    print(out_fname + ' saved.')
    # save p map
    fname_atts['val2'] = "p"
    out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
    nib.save(p_val, out_fname)
    print(out_fname + ' saved.')
    # FDR correction
    if mask is None:
        test_indices = np.argwhere(np.ones(refnii.get_data().flatten().shape)).flatten()
    else:
        test_indices = np.argwhere(mask.get_data().flatten()).flatten()
    p_arr = p_val.get_data().flatten()
    p_arr2 = np.take(p_arr, test_indices)
    sigs2, fdr_val2, a2, b2 = multipletests(p_arr2, alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    pfdr2 = deepcopy(1-fdr_val2)
    pfdr_arr = np.zeros(p_arr.shape)
    np.put(pfdr_arr, test_indices, pfdr2)
    pfdr_3d2 = pfdr_arr.reshape(p_val.shape)
    pfdr_map = nib.Nifti1Image(pfdr_3d2, refnii.affine, refnii.header)
    fname_atts['val2'] = "p"
    fname_atts['correction'] = "fdr"
    fname_atts['dir'] = "rev"
    out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
    nib.save(pfdr_map, out_fname)
    print(out_fname + ' saved.')
    # plot maps
    thresholded_map2, threshold2 = map_threshold(
        out_map, level=.05, height_control='fdr')#, two_sided=False)
    display = plot_stat_map(out_map, title='raw stat')
    print('The FDR=.05 threshold is %.3g' % threshold2)
    if not math.isinf(threshold2):
        #     thresholded_map2, threshold2 = map_threshold(
        #     out_map, level=.001, height_control=None)#, two_sided=False)
        #     fname_atts['val2'] = out_stat
        #     fname_atts['thresh'] = "p01"
        #     fname_atts['plot'] = "statmap"
        #     fname_atts['correction'] = "none"
        #     out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
        #     display = plot_stat_map(thresholded_map2, title='z map, p < 0.001', threshold = threshold2)
        #     display.savefig(out_fname)
        # else:
        fname_atts['val2'] = out_stat
        fname_atts['thresh'] = "p05"
        fname_atts['plot'] = "statmap"
        fname_atts['correction'] = "fdr"
        out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
        display = plot_stat_map(thresholded_map2, cut_coords=display.cut_coords,
                   title='Thresholded '+out_stat+' map, expected fdr = .05, '+out_stat+' = '+str(threshold2),
                   threshold=threshold2)
        display.savefig(out_fname)
        fname_atts['plot'] = "glassbrain"
        out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
        display = plot_glass_brain(thresholded_map2, cut_coords=display.cut_coords,
                   title='Thresholded '+out_stat+' map, expected fdr = .05, z = '+str(threshold2),
                   threshold=threshold2)
        display.savefig(out_fname)


parc_label = 'sl' + SL_RADIUS
corrs=['spear', 'reg']
preds_all = ['deg', 'dist', 'sn', 'img', 'soc', 'phys']+cfd_soc + cfd_phys
corr_labels = []
tasks_all = TASKS + ['avg']
tasks=[]
preds=[]
# read in arguments
for arg in sys.argv[1:]:
    print(arg)
    if arg in corrs:
        corr_labels += [arg]
    elif arg in tasks_all+['diff']:
        tasks += [arg]
    elif arg in preds_all:
        preds += [arg]

if len(corr_labels) == 0:
    corr_labels = ['spear', 'reg']

if len(tasks) == 0:
    tasks = tasks_all

if len(preds) == 0:
    preds = preds_all

task_sects = [[],[]]
for t in tasks:
    task_sects[t in TASKS] = task_sects[t in TASKS] + [t]


# dilate gray matter mask
tmp = nib.load(MNI_GM_MASK)
gm_mask = load_nii(MNI_GM_MASK)
gm_mask = np.where(gm_mask > 0.5, 1., 0.)
gm_mask_dil = binary_dilation(gm_mask, iterations=int(SL_RADIUS)).astype(gm_mask.dtype)
gm_mask_dil_img = nib.Nifti1Image(gm_mask_dil, tmp.affine, tmp.header)
gm_mask_img = nib.Nifti1Image(gm_mask, tmp.affine, tmp.header)
# choose which mask to use
mask_img_use = gm_mask_dil_img

for corr_label in corr_labels:
    print('starting correlation '+corr_label)
    in_dir = os.path.join(RSA_DIR,'*',corr_label)
    if SPACE=='T1w':
        in_dir = os.path.join(in_dir,'T1w-2-MNI')
    corr_vals = ['r'] if corr_label == 'spear' else ['R2', 'beta']
    for val_label in corr_vals:
        print('starting val_label '+val_label)
        for pred in preds:
            print('starting predictor '+pred)
            data_tasks = {}
            for task in tasks:
                print("starting task "+task)
                fnames = os.path.join(in_dir, '*task-'+task+'*space-'+SPACE+'*_parc-'+parc_label+'_*val-'+val_label+'_*pred-'+pred+'*.nii*')
                data_fnames = glob.glob(fnames)
                print(fnames)
                if len(data_fnames) != 0:
                    if task in TASKS:
                        data_tasks[task] = data_fnames
                        print(data_tasks)
                    run_sig_tests(data_fnames, mask = mask_img_use)
            if pred in ['deg', 'dist'] and 'number' in data_tasks.keys() and 'friend' in data_tasks.keys():
                run_sig_tests(data_tasks, mask = mask_img_use)

print(str(datetime.now()) + ": End level2_rsa_sl.py")
