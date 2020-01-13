#!/bin/python
# runs regression between neural data and various model predictors
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, pearsonr, norm, zscore
import os
import shutil
import sys
from datetime import datetime
from copy import deepcopy
import glob
import pandas as pd
import nibabel as nib
import numpy as np
import csv
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
    df_sub = deepcopy(df)
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

def get_roi_data(data, mask):
    roi_mask = mask.astype(int)
    n_vox = sum(roi_mask.flatten())
    roi_data = np.zeros((N_NODES,n_vox))
    for n in range(N_NODES):
        # get data for specific node based on 4th dimension
        data_node = data[:,:,:,n]
        # mask the data to only include this ROI
        roi_data_node = data_node * roi_mask
        # flatten into vector
        roi_data_node = roi_data_node.flatten()
        # remove 0s and save to roi_data by row
        roi_data[n,:] = roi_data_node[roi_data_node != 0]
    return(roi_data)

def get_roi_csv_val(df, roi, val_col='r'):
    # get rows for this roi
    r = df.loc[df['roi']==roi]
    # find corresponding stimulus file name
    preds = np.array(r.loc[:,'predictor'])
    out_dict = {}
    for p in preds:
        v = np.array(r.loc[r['predictor']==p,'r'])[0]
        out_dict[p] = v
    return(out_dict)

def ortho_mat(model_mat):
    model_mat = zscore(model_mat)
    model_mat,R = np.linalg.qr(model_mat)
    return(model_mat)

def rsa_reg(neural_v, model_mat):
    # orthogonalize
    model_mat = ortho_mat(model_mat)
    # add column of ones (constant)
    X=np.hstack((model_mat, np.ones((model_mat.shape[0],1))))
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
    return(out)

def rsa_cor(neural_v, model_mat):
    out = spearmanr(neural_v, model_mat)[0][0,1:]
    return(out)

def run_rsa(roi_data, model_mat, corr='corr'):
    # calculate roi RDM (lower triangle)
    roi_tri = pdist(roi_data, metric='correlation')
    roi_tri[np.isnan(roi_tri)] = 0.
    # run regression
    if corr=='corr':
        res = rsa_cor(neural_v=roi_tri, model_mat=model_rdms_mat)
    elif corr=='reg':
        res = rsa_reg(neural_v=roi_tri, model_mat=model_rdms_mat)
    out_dict = {'result':res}
    return(out_dict)

print(str(datetime.now()) + ": Begin rsa_regs.py")

# read in arguments
# searchlight or parcellation?
procedure = sys.argv[1]
# subjects
all_sub = sys.argv[2:]

# overwrite previous csv files? If False, will read in csv and append
overwrite = False

# print variables to log file
print(str(datetime.now()) + ": Project directory = " + PROJECT_DIR)
print(str(datetime.now()) + ": Analyzing " + str(len(all_sub)) + " subjects: " + str(all_sub))

# get project directories
in_dir = RSA_DIR + '%s/'
out_dir = RSA_DIR + '%s/' + procedure + '/' #%(sub)
# filenames
data_fnames = GLM_DIR + '%s/%s_task-%s_space-'+SPACE+'_stat-'+STAT+'_node-%02d.nii'#%(sub, sub,task,node)
parcellation_fname = os.path.basename(MNI_PARCELLATION).split('.')[0]
sub_parc = in_dir + parcellation_fname + '_%s_space-' + SPACE + '_transformed.nii'
sub_mask_fname = in_dir + '%s_desc-brain_mask_dil-5.nii.gz'
out_fname = out_dir + '%s_task-avg_stat-'+STAT+'_corr-%s_parc-%s_val-%s_pred-%s.nii.gz' #% (sub, corr, N_PARCELS, "r" or "b", predictor_name)
csv_fname = out_dir + "%s_task-avg_stat-"+STAT+"_corr-%s_parc-%s_roi_stats.csv"

corr = 'corr'
corr_label = 'spear' if corr=='corr' else 'reg'
val_label = 'r' if corr=='corr' else 'beta'
# set label for filename
parc_label = SL if isSl() else str(N_PARCELS)

# set variables
deg_label = 'deg'
dist_label = 'dist'
rad = 5 # radius of searchlight in voxels
sphere_dim = 2*rad+1
sphere_template = np.zeros((sphere_dim,sphere_dim,sphere_dim))
center_vox = (rad,rad,rad)
for xs in range(sphere_dim):
    for ys in range(sphere_dim):
        for zs in range(sphere_dim):
            d = np.sqrt((xs-center_vox[0])**2 + (ys-center_vox[1])**2 + (zs-center_vox[2])**2)
            if d <= rad:
                sphere_template[xs,ys,zs] = 1

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


# run regression for each subject
for sub in all_sub:
    print(str(datetime.now()) + ": Analyzing subject %s" % sub)

    # make output directories if they don't exist
    if not os.path.exists(out_dir % sub):
        os.makedirs(out_dir % sub)

    print(str(datetime.now()) + ": Reading in data for subject " + sub)

    # out csv filename
    sub_csv_fname = csv_fname % (sub, sub, corr_label, parc_label)
    if (not overwrite) and os.path.isfile(sub_csv_fname):
        read_bool = True
        # read in csv if already exists
        sub_csv = pd.read_csv(sub_csv_fname)
        # remove row column
        sub_csv = sub_csv.iloc[:,1:]
        # save all completed ROIs except for last one (since may not have been finished)
        completed_rois = np.unique(sub_csv['roi'])[:-1]
        print(completed_rois)
        out_csv_array = sub_csv.values.tolist()
    else:
        read_bool = False
        out_csv_array = []
        completed_rois = []
    wtr = csv.writer(open(sub_csv_fname, 'a'), delimiter=',', lineterminator='\n')

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
                                  meas_name_array=cfd_soc, compress=False,
                                  out_key='soc')
    phys_rdms = get_model_RDM_dict(node_mapping=node_order,
                                   meas_name_array=cfd_phys, compress=False,
                                   out_key='phys')

    # combine all predictor RDM dictionaries
    model_rdms = {**sn_rdms, **soc_rdms}
    model_keys = model_rdms.keys()

    # turn dictionary into matrix
    for i,k in enumerate(model_keys):
        if i == 0:
            model_rdms_mat = [model_rdms[k]]
        else:
            model_rdms_mat = np.vstack((model_rdms_mat, model_rdms[k]))

    # transpose so that each column corresponds to each measure
    model_rdms_mat = np.transpose(model_rdms_mat)


    if isSl():
        # column names for csv file
        colnames = ['sub','x','y','z','roi','predictor',val_label]
        if not read_bool:
            # write out to csv
            wtr.writerow(colnames)
        ref_img = sub_template
        # load subject's whole brain mask
        sub_mask = np.round(load_nii(sub_mask_fname % (sub, sub)))
        # setup parcellation-like image where each voxel that is the center of
        # a sphere is assigned an integer
        parc_data = deepcopy(sub_mask)
        # setup out data
        out_data = deepcopy(sub_data)*0.0
        # setup dictionary of output images
        out_data_dict = {}
        # add a key for each model tested
        for i,k in enumerate(model_keys):
            out_data_dict[k] = deepcopy(out_data)
        # start parcel count
        parc_roi = 0
        nvox = np.sum(sub_mask)
        # go through all voxels in sub_mask
        for index, val in np.ndenumerate(sub_mask):
            # if this voxel has a non-zero value, analyze it
            if val != 0:
                # save integer to parcellation just for the center voxel
                parc_data[index] = parc_roi
                # use next integer for next voxel
                parc_roi += 1
                if parc_roi in completed_rois:
                    print(str(datetime.now()) + ': ROI '+str(parc_roi)+' already saved.')
                    # read in values from dataframe for nii
                    res = get_roi_csv_val(sub_csv, parc_roi, val_label)
                else:
                    perc_done = round((parc_roi / nvox) * 100, 3)
                    print(str(datetime.now()) + ': Analyzing sphere '+str(index)+' -- '+str(perc_done)+'%')
                    # get voxel's coordinates
                    x = index[0]
                    y = index[1]
                    z = index[2]
                    xmin = x-rad if x-rad >= 0 else 0
                    ymin = y-rad if y-rad >= 0 else 0
                    zmin = z-rad if z-rad >= 0 else 0
                    xmax = x+rad if x+rad < sub_mask.shape[0] else sub_mask.shape[0]
                    ymax = y+rad if y+rad < sub_mask.shape[1] else sub_mask.shape[1]
                    zmax = z+rad if z+rad < sub_mask.shape[2] else sub_mask.shape[2]
                    # shape sphere according to above
                    sphere = deepcopy(sphere_template)
                    sphere = sphere[(center_vox[0]-x+xmin):(xmax-x+center_vox[0]),
                                    (center_vox[1]-y+ymin):(ymax-y+center_vox[1]),
                                    (center_vox[2]-z+zmin):(zmax-z+center_vox[2])]
                    mask_roi = deepcopy(sub_mask)*0.0
                    mask_roi[xmin:xmax, ymin:ymax, zmin:zmax] = sphere
                    mask_roi = mask_roi * sub_mask
                    # get the subject's data that is in this sphere
                    roi_data = get_roi_data(sub_data, mask_roi)
                    # run correlations comparing this sphere of data to each model
                    res_dict = run_rsa(roi_data, model_rdms_mat, corr=corr)
                    # get the resulting values
                    res = res_dict['result']
                # for each model, save the result to its image in out_data_dict
                for i,k in enumerate(model_keys):
                    # save to dataframe if not already there
                    if parc_roi not in completed_rois:
                        # array
                        val = res[i]
                        csv_row = [sub, x, y, z, parc_roi, k, val]
                        out_csv_array.append(csv_row)
                        # write out to csv
                        wtr.writerow(csv_row)
                    else:
                        # dictionary
                        val = res[k]
                    # update voxels
                    model_data = out_data_dict[k]
                    model_data[model_data==parc_roi] = val
                    out_data_dict[k] = model_data

    elif isParc():
        # column names for csv file
        colnames = ['sub','roi','predictor',val_label]
        if not read_bool:
            # write out to csv
            wtr.writerow(colnames)
        ref_img = parcellation_template
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
        out_data = ref_img.get_data().astype(np.double)
        # create a dictionary of nii's: one per predictor
        out_data_dict = {}
        for i,k in enumerate(model_keys):
            out_data_dict[k] = deepcopy(out_data)
        # iterate through each ROI of parcellation and run regression
        for r, parc_roi in enumerate(roi_list):
            if parc_roi in completed_rois:
                print(str(datetime.now()) + ': ROI '+str(parc_roi)+' already saved.')
                # read in values from dataframe for nii
                res = get_roi_csv_val(sub_csv, parc_roi, val_label)
            else:
                perc_done = round(((r+1) / len(roi_list)) * 100, 3)
                print(str(datetime.now()) + ': Analyzing ROI '+str(parc_roi)+' -- '+str(perc_done)+'%')
                # create mask for this ROI
                roi_mask = parc_data==parc_roi
                roi_mask = roi_mask.astype(int)
                roi_data = get_roi_data(sub_data, roi_mask)
                res_dict = run_rsa(roi_data, model_rdms_mat, corr=corr)
                res = res_dict['result']
            # for each model, save the result to its image in out_data_dict
            for i,k in enumerate(model_keys):
                # save to dataframe if not already there
                if parc_roi not in completed_rois:
                    val = res[i]
                    csv_row = [sub, parc_roi, k, val]
                    out_csv_array.append(csv_row)
                    # write out to csv
                    wtr.writerow(csv_row)
                else:
                    val = res[k]
                # update voxels
                model_data = out_data_dict[k]
                model_data[model_data==parc_roi] = val
                out_data_dict[k] = model_data

    for k in out_data_dict.keys():
        fname = out_fname % (sub, sub, corr_label, parc_label, val_label, k)
        save_nii( out_data_dict[k], ref_img, fname )


print(str(datetime.now()) + ": End rsa_regs.py")
