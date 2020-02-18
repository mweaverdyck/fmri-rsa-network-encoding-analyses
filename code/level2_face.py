from nistats import second_level_model
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
from nilearn import plotting
from nilearn.plotting import plot_glass_brain
#from nistats.second_level_model import non_parametric_inference
from funcs import *

out_dir = SECOND_LEVEL_DIR
mni_mask_dil_img = nib.load(MNI_MASK_DIL)

in_dir = os.path.join(GLM_DIR,'sub-*')
if SPACE=='T1w':
    in_dir = os.path.join(in_dir,'T1w-2-MNI')

# get file names
fnames = os.path.join(in_dir, '*task-??????_space-'+SPACE+'*node-all.nii')
data_fnames = glob.glob(fnames)
print("Running significance testing on: ")
print(data_fnames)

# setup file names
fname_atts = get_all_bids_atts(data_fnames[0])
del fname_atts['sub']
del fname_atts['task']
fname_atts['con'] = "face"
fname_atts['val2'] = "z"
fname_atts['correction'] = "none"
if 'extra' in fname_atts.keys():
    fname_atts.move_to_end('extra')

# setup model
model = SecondLevelModel(mask=mni_mask_dil_img, smoothing_fwhm=5.0)
design_matrix = pd.DataFrame([1] * len(data_fnames),
                         columns=['intercept'])
con_name = 'intercept'
print("Using design matrix: ")
print(design_matrix)
# save z map
model = model.fit(data_fnames,
                       design_matrix=design_matrix)
z_map = model.compute_contrast(con_name, output_type='z_score')
out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
nib.save(z_map, out_fname)
print(out_fname + ' saved.')
threshold = 3.1 #3.1  # correponds to  p < .001, uncorrected
display = plotting.plot_glass_brain(
    z_map, threshold=threshold, colorbar=True, plot_abs=False,
    output_file=out_fname,
    title='z map')
# save p map
p_val = model.compute_contrast(con_name, output_type='p_value')
fname_atts['val2'] = "p"
out_fname = os.path.join(SECOND_LEVEL_DIR, SL, make_bids_str(fname_atts))
nib.save(p_val, out_fname)
print(out_fname + ' saved.')
# correct for multiple comparisons
# Correcting the p-values for multiple testing and taking negative logarithm
n_voxels = np.sum(model.masker_.mask_img_.get_data())
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

# from nistats.second_level_model import non_parametric_inference
# neg_log_pvals_permuted_ols_unmasked = \
#     non_parametric_inference(data_fnames,
#                              design_matrix=design_matrix,
#                              model_intercept=True, n_perm=1000,
#                              two_sided_test=False,
#                              smoothing_fwhm=5.0, n_jobs=1)

print(str(datetime.now()) + ": End level2_rsa_sl.py")
