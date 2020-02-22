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
from scipy.stats import norm, ttest_rel, ttest_ind, ttest_1samp
from statsmodels.stats.multitest import multipletests
from scipy.ndimage.morphology import binary_dilation
from nilearn.image import math_img
from nilearn.input_data import NiftiMasker
from nistats import second_level_model
from nilearn import plotting
from nilearn.plotting import plot_glass_brain, plot_stat_map
from nistats import thresholding
from nistats.thresholding import map_threshold
#from nistats.second_level_model import non_parametric_inference
from funcs import *

def run_sig_tests(data_fnames, out_dir, mask=None):
    cmap = "Wistia"
    second_level_model = SecondLevelModel(smoothing_fwhm=5.0, mask=mask)
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
    out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
    nib.save(z_map, out_fname)
    print(out_fname + ' saved.')
    threshold = 2.88 #3.1  # correponds to  p < .001, uncorrected
    display = plotting.plot_glass_brain(
        z_map, threshold=threshold, colorbar=True, plot_abs=False,
        output_file=out_fname, cmap = cmap,
        title='z map')
    # FDR correction
    p_arr = p_val.get_data().flatten()
    sigs, fdr_val, a, b = multipletests(p_arr, alpha=.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    pfdr = deepcopy(1-fdr_val)
    pfdr[fdr_val==0.] = 0.
    pfdr = pfdr.reshape(p_val.shape)
    pfdr_map = nib.Nifti1Image(pfdr, p_val.affine, p_val.header)
    fname_atts['val2'] = "p"
    fname_atts['correction'] = "fdr"
    fname_atts['dir'] = "rev"
    out_fname = os.path.join(out_dir, make_bids_str(fname_atts))
    nib.save(pfdr_map, out_fname)
    print(out_fname + ' saved.')
    threshold = .95
    display = plotting.plot_glass_brain(
        pfdr_map, threshold=threshold, colorbar=True, plot_abs=False,
        output_file=out_fname, cmap = cmap,
        title='p map FDR-corrected, p < 0.05')

parc_label = 'sl' + SL_RADIUS
corrs=['spear', 'reg']
preds_all = ['sn', 'phys', 'soc']
corr_labels = []
tasks=[]
preds=[]
# read in arguments
for arg in sys.argv[1:]:
    print(arg)
    if arg in corrs:
        corr_labels += [arg]
    elif arg in TASKS or arg == 'avg':
        tasks += [arg]
    elif arg in preds_all:
        preds += [arg]

if len(corr_labels) == 0:
    corr_labels = ['reg']

if len(tasks) == 0:
    tasks = ['avg']

if len(preds) == 0:
    preds = preds_all

# dilate gray matter mask
tmp = nib.load(MNI_GM_MASK)
gm_mask = load_nii(MNI_GM_MASK)
gm_mask = np.where(gm_mask > 0.5, 1., 0.)
gm_mask_dil = binary_dilation(gm_mask, iterations=int(SL_RADIUS)).astype(gm_mask.dtype)
gm_mask_dil_img = nib.Nifti1Image(gm_mask_dil, tmp.affine, tmp.header)

for corr_label in corr_labels:
    print('starting correlation '+corr_label)
    in_dir = os.path.join(RSA_DIR,'*',corr_label)
    out_dir = os.path.join(SECOND_LEVEL_DIR, corr_label)
    if SPACE=='T1w':
        in_dir = os.path.join(in_dir,'T1w-2-MNI')
    for pred in preds:
        print('starting predictor '+pred)
        data_tasks = {}
        for task in tasks:
            print("starting task "+task)
            fnames = os.path.join(in_dir, '*task-'+task+'*space-'+SPACE+'*_parc-'+parc_label+'_*val-R2_*pred-'+pred+'*.nii*')
            data_fnames = glob.glob(fnames)
            print(in_dir)
            if len(data_fnames) != 0:
                if task in TASKS:
                    data_tasks[task] = data_fnames
                    print(data_tasks)
                run_sig_tests(data_fnames, out_dir, mask = gm_mask_dil_img)

print(str(datetime.now()) + ": End level2_rsa_stages.py")
