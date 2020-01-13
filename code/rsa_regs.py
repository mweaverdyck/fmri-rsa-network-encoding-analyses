#!/bin/python
# runs regression between neural data and various model predictors
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, pearsonr, norm, zscore
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

def isParc():
    return(procedure == PARC)

def isSl():
    return(procedure == SL)

def load_nii( imgpath ):
    img4D = nib.load(imgpath, mmap=False).get_data()
    return img4D

def save_nii( data, refnii, filename ):
    out_img = nib.Nifti1Image(data, refnii.affine, refnii.header)
    out_img.to_filename(filename)
    print(str(datetime.now()) + ": File %s saved." % filename)

def make_RDM( data_array ):
    # find number of nodes
    n = len(data_array)
    if n != N_NODES:
        print("WARNING: number nodes entered ("+ str(n) +") is different than N_NODES ("+str(N_NODES)+"):")
        print(data_array)
    # create empty matrix
    mat = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            mat[i,j] = abs(data_array[i] - data_array[j])
    # make squareform (upper diagonal elements)
    tri = squareform(mat)
    return tri

def get_node_mapping( sub ):
    events_file = glob.glob(BIDS_DIR+sub+"/func/"+sub+"*events.tsv")[0]
    events_df = pd.read_csv(events_file, delimiter='\t')
    node_mapping = {}
    for n in range(N_NODES):
        # find first row that matches this node
        r = events_df.loc[events_df['node']==n].index[0]
        # find corresponding stimulus file name
        stim = events_df.loc[r,'stim_file']
        # remove '.png'
        stim = stim.split('.')[0]
        # save to mapping
        node_mapping[stim] = n
    # order node numbers to match CFD measures
    node_order = []
    cfd_targets = pd.read_csv(CFD_FNAME)['Target']
    for t in cfd_targets:
        node_order.append(node_mapping[t])
    return(node_order)

def get_model_RDM_dict( node_mapping, meas_name_array,
                        df = None, fname=CFD_FNAME,
                        compress=False, out_key='sum' ):
    if df is None:
        df = pd.read_csv(fname)
    # copy chicago face database data frame
    df_sub = copy.deepcopy(df)
    # add to CFD measures data frame and sort
    df_sub['Node'] = node_order
    df_sub.sort_values(by='Node', inplace=True)
    # social features
    rdm_dict = {}
    for i in meas_name_array:
        # extract column
        v = df_sub[i]
        if compress:
            # normalize so that all measures are on the same scale
            df_sub['z'+i] = zscore(v)
        else:
            # make RDMs
            rdm_dict[i] = make_RDM(v)
    if compress:
        # select the normalized names
        norm_meas_name_array = np.core.defchararray.add('z', meas_name_array)
        v2 = df_sub[norm_meas_name_array].sum(axis=1)
        rdm_dict[out_key] = make_RDM(v)
    return(rdm_dict)

def ortho_mat(model_mat):
    model_mat = zscore(model_mat)
    model_mat,R = np.linalg.qr(model_mat)
    return(model_mat)

def run_rsa_reg(neural_v, model_mat):
    # orthogonalize
    model_mat = ortho_mat(model_mat)
    # add column of ones (constant)
    X=np.hstack((model_mat,
                 np.ones((model_mat.shape[0],1))))
    # Convert neural DSM to column vector
    neural_v=neural_v.reshape(-1,1)
    # Compute betas and constant
    betas = np.linalg.lstsq(X, neural_v, rcond=None)[0]
    # Determine model (if interested in R)
    for k in range(len(betas)-1):
        if k==0:
            model = X[:,k].reshape(-1,1)*betas[k]
        else:
            model = model + X[:,k].reshape(-1,1)*betas[k]
    # Get multiple correlation coefficient (R)
    R = pearsonr(neural_v.flatten(),model.flatten())[0]
    out = np.concatenate((betas.flatten()[:-1], np.array([R])))
    return out


print(str(datetime.now()) + ": Begin rsa_regs.py")

# read in arguments
# searchlight or parcellation?
procedure = sys.argv[1]
# subjects
all_sub = sys.argv[2:]

# print variables to log file
print(str(datetime.now()) + ": Project directory = " + PROJECT_DIR)
print(str(datetime.now()) + ": Analyzing " + str(len(all_sub)) + " subjects: " + str(all_sub))

# get project directories
in_dir = RSA_DIR + '%s/'
out_dir = RSA_DIR + '%s/reg/' + procedure + '/' #%(sub)
# filenames
data_fnames = GLM_DIR + '%s/%s_task-%s_space-'+SPACE+'_stat-'+STAT+'_node-%02d.nii'#%(sub, sub,task,node)
parcellation_fname = os.path.basename(MNI_PARCELLATION).split('.')[0]
sub_parc = in_dir + parcellation_fname + '_%s_space-' + SPACE + '_transformed.nii'
sub_mask = in_dir + '*brain_mask*dil-*'
out_fname = out_dir + '%s_stat-'+STAT+'_corr-reg_parc-%s_val-beta_pred-%s.nii.gz' #% (sub, N_PARCELS, "r" or "b", predictor_name)
csv_fname = out_dir + "%s_stat-"+STAT+"_corr-reg_parc-%s_roi_stats.csv"

# set variables
deg_label = 'deg'
dist_label = 'dist'
model_keys = [deg_label, 'soc', 'phys']
rad = 5 # radius of searchlight in voxels


# reference image for saving parcellation output
parcellation_template = nib.load(MNI_PARCELLATION, mmap=False)
parcellation_template.set_data_dtype(np.double)

# chicago face database measures
cfd_soc = ['Dominant', 'Trustworthy']
cfd_phys = ['Unusual', 'Faceshape', 'Heartshapeness', 'Noseshape', 'LipFullness', 'EyeShape', 'EyeSize', 'UpperHeadLength', 'MidfaceLength', 'ChinLength', 'ForeheadHeight', 'CheekboneHeight', 'CheekboneProminence', 'FaceRoundness']

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

# Degree RDM
deg_tri = make_RDM(deg)
print(str(datetime.now()) + ": Degree RDM: ")
print(squareform(deg_tri))

# Distance RDM
dist_tri = squareform(dist_mat)
print(str(datetime.now()) + ": Distance RDM: ")
print(dist_mat)

# social network RDM dictionary
sn_rdms = {deg_label: deg_tri, dist_label: dist_tri}
deg_rdm = {deg_label: deg_tri}

# dictionary of rdms by subject key
out_csv_df = []
colnames = ['sub','roi','predictor','beta']

# run regression for each subject
for sub in all_sub:
    print(str(datetime.now()) + ": Analyzing subject %s" % sub)

    # make output directories if they don't exist
    if not os.path.exists(out_dir % sub):
        os.makedirs(out_dir % sub)

    print(str(datetime.now()) + ": Reading in data for subject " + sub)
    # read in 4D image
    sub_template = nib.load(data_fnames % (sub, sub, TASKS[0], 0), mmap=False)
    sub_dims = sub_template.get_data().shape + (N_NODES,)
    sub_data = np.empty(sub_dims)
    for n in range(N_NODES):
        print(str(datetime.now()) + ": Reading in node " + str(n))
        # average this node's data from both runs
        d1 = load_nii(data_fnames % (sub, sub, TASKS[0], n))
        d2 = load_nii(data_fnames % (sub, sub, TASKS[1], n))
        d_avg = (d1+d2)/2
        # save to fourth dimension
        sub_data[:,:,:,n] = d_avg

    # find subject's node-image mapping
    node_order = get_node_mapping( sub )

    # get model RDMs from CFD measures
    soc_rdms = get_model_RDM_dict(node_mapping=node_order,
                                  meas_name_array=cfd_soc, compress=True,
                                  out_key='soc')
    phys_rdms = get_model_RDM_dict(node_mapping=node_order,
                                   meas_name_array=cfd_phys, compress=True,
                                   out_key='phys')

    # combine all predictor RDM dictionaries
    model_rdms = {**deg_rdm, **soc_rdms, **phys_rdms}

    # turn dictionary into matrix
    for i,k in enumerate(model_keys):
        if i == 0:
            model_rdms_mat = [model_rdms[k]]
        else:
            model_rdms_mat = np.vstack((model_rdms_mat, model_rdms[k]))

    # transpose so that each column corresponds to each measure
    model_rdms_mat = np.transpose(model_rdms_mat)


    if isSl():
        sub_data_copy = deepcopy(sub_data)
        for measurename, measure in measure_dict.iteritems():
            print(str(datetime.now()) + ": Creating searchlight of size " + str(rad))
            sl = sphere_searchlight(measure, rad) #, postproc=FisherTransform)
            print(str(datetime.now()) + ": Running searchlight.")
            sl_map = sl(sub_data_copy)

            print(str(datetime.now()) + ": Saving output images")
            # output image
            nimg_res = map2nifti(sl_map, imghdr = sub_data_copy.a.imghdr)
            # save images
            fname = out_fname % (sub, sub, 'sl', measurename)
            nimg_res.to_filename(fname)
            print(str(datetime.now()) + ": File %s saved." % fname)

    elif isParc():
        # make mask
        parcellation = sub_parc % (sub, sub)
        print(str(datetime.now()) + ": Using parcellation " + parcellation)
        parc_data = load_nii(parcellation)
        roi_list = np.unique(parc_data)
        # remove 0 (i.e., the background)
        roi_list = np.delete(roi_list,0)
        # check if number of parcels matches global variable
        if N_PARCELS != len(roi_list):
            print("WARNING: Number of parcels found ("+str(len(roi_list))+") does not equal N_PARCELS ("+str(N_PARCELS)+")")

        # Run regression on each parcellation
        print(str(datetime.now()) + ": Starting parcellation "+ str(N_PARCELS))
        # get the voxels from parcellation nii
        out_data = parcellation_template.get_data().astype(np.double)
        # create a dictionary of nii's: one per predictor
        out_data_dict = {}
        for i,k in enumerate(model_keys):
            out_data_dict[k] = copy.deepcopy(out_data)
        # iterate through each ROI of parcellation and run regression
        for r, parc_roi in enumerate(roi_list):
            print(str(datetime.now()) + ": Running regression on ROI %d (%d/%d)..." % (parc_roi, r+1, N_PARCELS))
            # create mask for this ROI
            roi_mask = parc_data==parc_roi
            roi_mask = roi_mask.astype(int)
            n_vox = sum(roi_mask.flatten())
            roi_data = np.zeros((N_NODES,n_vox))
            for n in range(N_NODES):
                # get data for specific node based on 4th dimension
                sub_data_node = sub_data[:,:,:,n]
                # mask the data to only include this ROI
                roi_data_node = sub_data_node * roi_mask
                # flatten into vector
                roi_data_node = roi_data_node.flatten()
                # remove 0s and save to roi_data by row
                roi_data[n,:] = roi_data_node[roi_data_node != 0]

            roi_tri = pdist(roi_data, metric='correlation')
            # run regression
            res = run_rsa_reg(neural_v=roi_tri, model_mat=model_rdms_mat)
            for i,k in enumerate(model_keys):
                beta=res[i]
                # save to dataframe
                out_csv_df.append([sub, parc_roi, k, beta])
                # update voxels
                model_data = out_data_dict[k]
                model_data[model_data==parc_roi] = beta
                out_data_dict[k] = model_data


        for k in model_keys:
            fname = out_fname % (sub, sub, N_PARCELS, k)
            save_nii( out_data_dict[k], parcellation_template, fname )

        out_csv_df = pd.DataFrame(out_csv_df, columns = colnames)
        out_csv_df.to_csv(csv_fname % (sub, sub, N_PARCELS))

print(str(datetime.now()) + ": End rsa_regs.py")
