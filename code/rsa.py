#!/bin/python
# runs correlation between neural data and predictors
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, pearsonr, norm
import os
import shutil
import sys
from datetime import datetime
import copy
import glob
import pandas as pd
import nibabel as nib
import numpy as np
from funcs import *

def load_nii( imgpath ):
    img4D = nib.load(imgpath, mmap=False).get_data()
    return img4D

print(str(datetime.now()) + ": Begin rsa.py")

# read in arguments
# searchlight or parcellation?
procedure = sys.argv[1]
# subjects
all_sub = sys.argv[2:]

# get project directories
in_dir = RSA_DIR + '%s/'
out_dir = RSA_DIR + '%s/' + procedure + '/' #%(sub)
# filenames
data_fnames = in_dir + '%s_task-%s_space-'+SPACE+'_stat-'+STAT+'_node-4D.nii'#%(sub,task)
parcellation_fname = os.path.basename(MNI_PARCELLATION).split('.')[0]
sub_parc = in_dir + parcellation_fname + '_%s_space-' + SPACE + '_transformed.nii'
sub_mask = in_dir + '*brain_mask*dil-*'
out_fname = out_dir + '%s_task-%s_stat-'+STAT+'_corr-spear_parc-%s_val-r_pred-%s.nii.gz' #% (sub, task, N_PARCELS, "r" or "b", predictor_name)
csv_fname = out_dir + "%s_stat-"+STAT+"_corr-spear_parc-%s_roi_stats.csv"
# reference image for saving parcellation output
parcellation_template = nib.load(MNI_PARCELLATION, mmap=False)
parcellation_template.set_data_dtype(np.double)

# set variables
deg_label = 'deg'
dist_label = 'dist'
rad = 5 # radius of searchlight in voxels

# subject-general predictors:
deg = np.array([1, 4, 2, 2, 3, 4, 2, 3, 2, 3])
dist_mat = np.array([
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
if N_NODES != len(deg):
    print("WARNING: number nodes entered for degrees ("+ str(len(deg)) +") is different than N_NODES ("+str(N_NODES)+")")

# print variables to log file
print(str(datetime.now()) + ": Project directory = " + PROJECT_DIR)
print(str(datetime.now()) + ": Analyzing " + str(len(all_sub)) + " subjects: " + str(all_sub))
print(str(datetime.now()) + ": Degrees of all " + str(N_NODES) + " nodes:" + str(deg))

# Make Degree RDM
# make degree matrix by taking difference between degrees
deg_mat = np.zeros([N_NODES,N_NODES], int)
for n in range(N_NODES):
    for n2 in range(N_NODES):
        deg_mat[n,n2] = abs(deg[n] - deg[n2])

# make squareform (upper diagonal elements) for correlating RDMs
deg_tri=squareform(deg_mat)
print(str(datetime.now()) + ": Degree RDM: ")
print(deg_mat)

# Make Distance RDM
# make squareform (upper diagonal elements) for correlating RDMs
dist_tri = squareform(dist_mat)
lenconstant = len(dist_tri)
print(str(datetime.now()) + ": Distance RDM: ")
print(dist_mat)

# save RDMs to a dictionary
model_rdms = {deg_label: deg_tri, dist_label: dist_tri}

# dictionary of rdms by subject key
out_csv_df = []
colnames = ['sub','task','roi','predictor','r']

# run regression for each subject
for sub in all_sub:
    print(str(datetime.now()) + ": Analyzing subject %s" % sub)

    # make output directories if they don't exist
    if not os.path.exists(out_dir % sub):
        os.makedirs(out_dir % sub)

    # Make Mask
    if procedure == PARC:
        parcellation = sub_parc % (sub, sub)
        print(str(datetime.now()) + ": Using parcellation " + parcellation)
        parc_data = load_nii(parcellation)
        roi_list = np.unique(parc_data)
        # remove 0 (i.e., the background)
        roi_list = np.delete(roi_list,0)
        # check if number of parcels matches global variable
        if N_PARCELS != len(roi_list):
            print("WARNING: Number of parcels found ("+str(len(roi_list))+") does not equal N_PARCELS ("+str(N_PARCELS)+")")

    for task in TASKS:
        print(str(datetime.now()) + ": Reading in data for task '" + task + "' for subject " + sub)
        # read in 4D image
        # TODO: read in 3D images and combine in script, rather than in wrapper
        data_fname = data_fnames % (sub, sub, task)
        print(str(datetime.now()) + ": data_fname file = " + data_fname)
        sub_data = load_nii(data_fname)

        if procedure == 'sl':
            sub_data_copy = deepcopy(sub_data)#_dict[sub][task])
            for measurename, measure in measure_dict.iteritems():
                print(str(datetime.now()) + ": Creating searchlight of size " + str(rad))
                sl = sphere_searchlight(measure, rad) #, postproc=FisherTransform)
                print(str(datetime.now()) + ": Running searchlight.")
                sl_map = sl(sub_data_copy)

                print(str(datetime.now()) + ": Saving output images")
                # output image
                nimg_res = map2nifti(sl_map, imghdr = sub_data_copy.a.imghdr)
                # save images
                fname = out_fname % (sub, sub, task, 'sl', measurename)
                nimg_res.to_filename(fname)
                print(str(datetime.now()) + ": File %s saved." % fname)

        elif procedure == PARC:
            print(str(datetime.now()) + ": Starting parcellation "+ str(N_PARCELS))
            deg_res_data = parcellation_template.get_data()
            dist_res_data = parcellation_template.get_data()
            out_data_dict = {deg_label: deg_res_data.astype(np.double), dist_label: dist_res_data.astype(np.double)}
            # iterate through each ROI of parcellation and run regression
            for r, parc_roi in enumerate(roi_list):
                print(str(datetime.now()) + ": Running regression on ROI %d (%d/%d)..." % (parc_roi, r+1, N_PARCELS))
                # create mask for this ROI
                roi_mask = parc_data==parc_roi
                roi_mask = roi_mask.astype(int)
                roi_data = np.zeros((N_NODES,sum(roi_mask.flatten())))
                for n in range(N_NODES):
                    sub_data_node = sub_data[:,:,:,n]
                    roi_data_node = sub_data_node * roi_mask
                    roi_data_node = roi_data_node.flatten()
                    roi_data[n:] = roi_data_node[roi_data_node != 0]
                # create neural RDM
                roi_tri = pdist(roi_data, metric='correlation')#1 - np.corrcoef(roi_data, rowvar=False)
                # extract lower triangle, excluding diagonal
                roi_rdm = squareform(roi_tri)
                # correlate with model RDMs
                for model_name in model_rdms:
                    # correlate model and neural RDMs
                    model_tri = model_rdms[model_name]
                    r_val, p = spearmanr(roi_tri, model_tri)
                    # save to dataframe
                    out_csv_df.append([sub, task, parc_roi, model_name, r_val])
                    # update voxels
                    model_data = out_data_dict[model_name]
                    model_data[model_data==parc_roi] = r_val
                    out_data_dict[model_name] = model_data


            for model_name in model_rdms:
                out_img = nib.Nifti1Image(out_data_dict[model_name], parcellation_template.affine, parcellation_template.header)
                # output filename
                fname = out_fname % (sub, sub, task, N_PARCELS, model_name)
                # save file
                out_img.to_filename(fname)
                print(str(datetime.now()) + ": File %s saved." % fname)

            out_csv_df = pd.DataFrame(out_csv_df, columns = colnames)
            out_csv_df.to_csv(csv_fname % (sub, sub, N_PARCELS))

print(str(datetime.now()) + ": End rsa.py")
