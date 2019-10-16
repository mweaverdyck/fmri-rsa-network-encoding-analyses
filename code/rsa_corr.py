#!/bin/python
# runs correlation between neural data and predictors
import mvpa2
from mvpa2.suite import *
from mvpa2.measures.rsa import PDistTargetSimilarity
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, pearsonr, norm
from numpy import corrcoef
import json
from collections import OrderedDict
#import statsmodels.api as sm
import os
import shutil
import subprocess
import sys
from datetime import datetime
import copy
import glob
import pandas as pd
import nibabel as nib
from funcs import *

# Fisher transformation function
FisherTransform=FxMapper('features',lambda r: 0.5*np.log((1+r)/(1-r)))

print str(datetime.now()) + ": Begin rsa_corr.py"

# read in arguments
# searchlight or parcellation?
procedure = sys.argv[1]
# which statistic to use
#stat = sys.argv[2]
# subjects
all_subj = sys.argv[2:]

# get project directories
in_dir = RSA_DIR + '%s/'
out_dir = RSA_DIR + '%s/' + procedure + '/' #%(subj)
# filenames
data_fnames = in_dir + '%s_task-%s_space-'+SPACE+'_stat-'+STAT+'_node-4D.nii'#%(subj,task)
parcellation_fname = os.path.basename(MNI_PARCELLATION).split('.')[0]
sub_parc = in_dir + parcellation_fname + '_%s_space-' + SPACE + '_transformed.nii'
sub_mask = in_dir + '*brain_mask*dil-*'
out_fname = out_dir + '%s_task-%s_stat-'+STAT+'_corr-spear_parc-%s_val-r_pred-%s.nii.gz' #% (subj, task, n_parcels, "r" or "b", predictor_name)
csv_fname = out_dir + "%s_stat-"+STAT+"_corr-spear_parc-%s_roi_stats.csv"
atts = RSA_DIR + 'atts.txt' # attributes file

# set variables
deg_label = 'deg'
dist_label = 'dist'
tasks = ['friend','number']
rad = 5 # radius of searchlight in voxels
# subject-general predictors:
deg = array([1, 4, 2, 2, 3, 4, 2, 3, 2, 3])
dist_mat = array([
       [0, 2, 3, 3, 3, 2, 3, 3, 3, 1],
       [2, 0, 1, 1, 1, 2, 3, 3, 3, 1],
       [3, 1, 0, 2, 1, 3, 4, 4, 4, 2],
       [3, 1, 2, 0, 1, 3, 4, 4, 4, 2],
       [3, 1, 1, 1, 0, 3, 4, 4, 4, 2],
       [2, 2, 3, 3, 3, 0, 1, 1, 1, 1],
       [3, 3, 4, 4, 4, 1, 0, 1, 2, 2],
       [3, 3, 4, 4, 4, 1, 1, 0, 1, 2],
       [3, 3, 4, 4, 4, 1, 2, 1, 0, 2],
       [1, 1, 2, 2, 2, 1, 2, 2, 2, 0]])

# get values
preds = [deg_label,dist_label]
sna_dsms = preds[:2]
nsubj = len(all_subj)
n_nodes = len(deg)
attr = SampleAttributes(atts)

# print variables to log file
print str(datetime.now()) + ": Project directory = " + PROJECT_DIR
print str(datetime.now()) + ": Analyzing " + str(nsubj) + " subjects: " + str(all_subj)
print str(datetime.now()) + ": Degrees of all " + str(n_nodes) + " nodes:" + str(deg)

# Make Degree DSM
# make degree matrix by taking difference between degrees
deg_mat = np.zeros([n_nodes,n_nodes], int)
for n in range(n_nodes):
    for n2 in range(n_nodes):
        deg_mat[n,n2] = abs(deg[n] - deg[n2])
# make squareform (upper diagonal elements) for correlating DSMs
deg_sq = squareform(deg_mat)
print str(datetime.now()) + ": Degree DSM: "
print deg_mat

# Make Distance DSM
# make squareform (upper diagonal elements) for correlating DSMs
dist_sq = squareform(dist_mat)
lenconstant = len(dist_sq)
print str(datetime.now()) + ": Distance DSM: "
print dist_mat

# initialize dictionaries/arrays
proc_data = []
ds_dict = {}
dstemp_dict = {}
# dictionary of dsms by subject key
subj_dsms = {}
out_csv_df = []
colnames = ['sub','task','roi','predictor','r']


# run regression for each subject
#subj=all_subj[0]
for subj in all_subj:
    print str(datetime.now()) + ": Analyzing subject %s" % subj

    # Make Mask
    if procedure == 'parc':
        parcellation = sub_parc % (subj, subj)
        print str(datetime.now()) + ": Using parcellation " + parcellation
        parc_ds = fmri_dataset(samples = parcellation)
        roi_list = list(set(parc_ds.samples[0]))
        # remove 0 (i.e., the background)
        roi_list.remove(0.)
        n_parcels = len(roi_list)
    else:
        masks = glob.glob(sub_mask % subj)
        if len(masks) != 1:
            print "WARNING: " + str(len(masks)) + " found. Attempting to use first mask."
        mask = masks[0]
        print str(datetime.now()) + ": Using mask " + mask
        mask_ds = fmri_dataset(samples = mask, mask = mask)

    # make output directories if they don't exist
    if not os.path.exists(out_dir % subj):
        os.makedirs(out_dir % subj)
    #if not os.path.exists(out_dir+subj):
    #    os.makedirs(out_dir+subj)

    # initialize this subject's dictionary
    subj_dsms[subj] = {}
    subj_dsms[subj][deg_label] = deg_sq
    subj_dsms[subj][dist_label] = dist_sq

    ds_dict[subj] = {}
    #task=tasks[0]
    for task in tasks:
        print str(datetime.now()) + ": Creating dataset for task '" + task + "' for subject " + subj

        data_fname = data_fnames % (subj, subj, task)
        print str(datetime.now()) + ": data_fname file = " + data_fname
        if procedure=='parc':
            ds = fmri_dataset(samples = data_fname)
        else:
            ds = fmri_dataset(samples = data_fname, mask = mask)

        # join attributes to dataset
        ds.sa['chunks'] = attr.chunks
        ds.sa['targets'] = attr.targets
        ds_dict[subj][task] = ds
        dset_copy = deepcopy(ds_dict[subj][task])

        if procedure == 'parc':
            # copy for each regression output
            blank_img = copy.deepcopy(parc_ds[0]) #copy.deepcopy(ds[0])
            img_dims = blank_img.samples.shape
            blank_img.samples = np.zeros(img_dims)
            results_dict={}

        measure_dict={} # Create dictionary of measures to use for within subjet searchlights
        # Add measures computing similarity with target DSMs
        for target_name, target_mat in subj_dsms[subj].iteritems():
            measure_dict[target_name] = PDistTargetSimilarity(target_mat,
                                                comparison_metric='spearman',
                                                corrcoef_only=True)
            if procedure == 'parc':
                results_dict[target_name] = copy.deepcopy(blank_img)


        if procedure == 'sl':
            for measurename, measure in measure_dict.iteritems():
                print str(datetime.now()) + ": Creating searchlight of size " + str(rad)
                sl = sphere_searchlight(measure, rad) #, postproc=FisherTransform)
                print str(datetime.now()) + ": Running searchlight."
                sl_map = sl(dset_copy)

                print str(datetime.now()) + ": Saving output images"
                # output image
                nimg_res = map2nifti(sl_map, imghdr = dset_copy.a.imghdr)
                # save images
                fname = out_fname % (subj, subj, task, 'sl', measurename)
                nimg_res.to_filename(fname)
                print str(datetime.now()) + ": File %s saved." % fname

        elif procedure == 'parc':

            print(str(datetime.now()) + ": Starting parcellation "+ str(n_parcels))

            # iterate through each ROI of parcellation and run regression
            for r, parc_roi in enumerate(roi_list):
                print(str(datetime.now()) + ": Running regression on ROI %d of %d..." % (r+1, n_parcels))
                # create mask for this ROI
                parc_mask = copy.deepcopy(parc_ds)
                parc_mask.samples=parc_mask.samples==parc_roi
                # use mask to subselect data_fname from this ROI
                ds_roi = fmri_dataset(samples = data_fname, mask=map2nifti(parc_mask))
                # join attributes to dataset
                ds_roi.sa['chunks'] = attr.chunks
                ds_roi.sa['targets'] = attr.targets

                for measurename, measure in measure_dict.iteritems():
                    # run regression
                    res = measure(ds_roi)
                    r_val = res.samples[0][0]
                    results_dict[measurename].samples[parc_mask.samples] = [r_val]
                    out_csv_df.append([subj, task, parc_roi, measurename, r_val])

            for measurename, measure in measure_dict.iteritems():
                nimg_res = map2nifti(results_dict[measurename], imghdr = ds.a.imghdr)
                # output filename
                fname = out_fname % (subj, subj, task, n_parcels, measurename)
                # save file
                nimg_res.to_filename(fname)
                print str(datetime.now()) + ": File %s saved." % fname

            out_csv_df = pd.DataFrame(out_csv_df, columns = colnames)
            out_csv_df.to_csv(csv_fname % (subj, subj, n_parcels))

print str(datetime.now()) + ": End rsa_corr.py"
