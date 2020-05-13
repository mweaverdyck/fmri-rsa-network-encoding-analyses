#!/bin/python
# runs regression between neural data and various model predictors
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, pearsonr, norm, zscore
from scipy.ndimage.morphology import binary_dilation
from brainiak.searchlight.searchlight import Searchlight, Ball
from brainiak.fcma.preprocessing import prepare_searchlight_mvpa_data
from brainiak import io
import os
import shutil
import sys
from datetime import datetime
import time
from copy import deepcopy
import glob
import pandas as pd
import csv
from funcs import *


# DIRECTORIES
in_dir = os.path.join(RSA_DIR, '%s/')
out_dir = get_thresh_dir(RSA_DIR)
out_dir = os.path.join(out_dir, '%s','%s/') #%(sub, corr_label)
print(str(datetime.now()) + ": Project directory = " + PROJECT_DIR)
print(str(datetime.now()) + ": Input directory = " + in_dir)
print(str(datetime.now()) + ": Output directory = " + out_dir)

# FILENAMES
# subject's data images
data_fnames = os.path.join(get_thresh_dir(GLM_DIR), '%s','%s_task-%s_space-'+SPACE+'_stat-'+STAT+'_node-%02d.nii')#%(sub, sub,task,node)
# subject's mask
sub_mask_fname = os.path.join(in_dir, '%s_desc-brain_mask_dil-5.nii.gz')
# parcellation
# standard parcellation filename
parcellation_fname = os.path.basename(MNI_PARCELLATION).split('.')[0]
if SPACE == 'T1w':
    # parcellation in subject's space
    sub_parc = os.path.join(in_dir, parcellation_fname + '_%s_space-' + SPACE + '_transformed.nii')
else:
    sub_parc = MNI_PARCELLATION

# output file names
out_fname = os.path.join(out_dir, '%s_task-%s_space-'+SPACE+'_stat-'+STAT+'_corr-%s_parc-%s_val-%s_pred-%s.nii.gz') #% (sub, task, corr, N_PARCELS, "r" or "b", predictor_name)
csv_fname = os.path.join(out_dir, "%s_task-%s_space-"+SPACE+"_stat-"+STAT+"_corr-%s_parc-%s_val-%s_roi_stats.csv")

# for searchlight
use_mask = False
bcvar = None #range(N_NODES)
sl_rad = int(SL_RADIUS)
max_blk_edge = 5

# chicago face database measures
cfd_phys1 = ['Race', 'Gender', 'Age', 'Luminance_median', 'fWHR']
cfd_soc1 = ['Babyface', 'Feminine', 'Masculine', 'Prototypic', 'Unusual']
cfd_soc2 = ['Attractive', 'Dominant', 'Threatening', 'Trustworthy']
cfd_emot = ['Angry', 'Disgusted', 'Happy', 'Sad', 'Surprised']
cfd_phys = cfd_phys1 + cfd_soc1
cfd_soc = cfd_soc2 + cfd_emot
cfd_img = cfd_phys1 + cfd_soc1 + cfd_soc2 + cfd_emot

# social network measures
deg_lab = 'deg'
dist_lab = 'dist'
dist2_lab = 'dist2'
deg_label = deg_lab + '_cat-sn'
dist_label = dist_lab + '_cat-sn'
dist2_label = dist2_lab + '_cat-sn'
sn_labels = [deg_label, dist_label, dist2_label]


# GENERAL FUNCTIONS
def make_model_RDM( data_array , categorical=False):
    # find number of nodes
    n = len(data_array)
    if n != N_NODES:
        print("WARNING: number nodes entered ("+ str(n) +") is different than N_NODES ("+str(N_NODES)+"):")
        print(data_array)
    if isinstance(data_array[0], str):
        categorical=True
    # create empty matrix
    mat = np.zeros([n,n])
    # find pairwise difference between every element
    for i in range(n):
        for j in range(n):
            if categorical:
                # if data are categorical, enter 0 for same category, 1 for different
                mat[i,j] = int(data_array[i] != data_array[j])
            else:
                # if data are quantitative, enter absolute difference
                mat[i,j] = abs(data_array[i] - data_array[j])
    # make squareform (upper diagonal elements)
    tri = squareform(mat)
    return tri

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
dist2_tri = pdist(dist_mat)

# Degree RDM
deg_tri = make_model_RDM(deg)
print(str(datetime.now()) + ": Degree RDM: ")
print(squareform(deg_tri))

# Distance RDM
dist_tri = squareform(dist_mat)
print(str(datetime.now()) + ": Distance RDM: ")
print(dist_mat)

# Second-Order Distance RDM
print(str(datetime.now()) + ": Second-Order Distance RDM: ")
print(squareform(dist2_tri))

def get_model_RDM_dict( node_mapping, meas_name_array,
                        df = None, fname=CFD_FNAME,
                        compress=False, out_key='sum' ):
    if df is None:
        df = pd.read_csv(fname)
    # copy chicago face database data frame
    df_sub = deepcopy(df)
    # add to CFD measures data frame and sort
    df_sub['Node'] = node_mapping
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
            if out_key != 'sum':
                k = i+'_cat-'+out_key
                print(v)
                # make RDMs
                rdm_dict[k] = make_model_RDM(v)
    if compress:
        # select the normalized names
        norm_meas_name_array = np.core.defchararray.add('z', meas_name_array)
        v2 = df_sub[norm_meas_name_array].sum(axis=1)
        rdm_dict[out_key] = make_model_RDM(v)
    return(rdm_dict)

def ortho_mat(model_mat):
    model_mat = zscore(model_mat)
    model_mat,R = np.linalg.qr(model_mat)
    return(model_mat)

def rsa_reg(neural_v, model_mat, val_label=None, rank_order=True):
    # orthogonalize
    model_mat = ortho_mat(model_mat)
    if rank_order:
        neural_v = np.argsort(neural_v, axis=0)
        model_mat = np.argsort(model_mat, axis=0)
    # add column of ones (constant)
    X=np.hstack((model_mat, np.ones((model_mat.shape[0],1))))
    # Convert neural DSM to column vector
    neural_v=neural_v.reshape(-1,1)
    # Compute betas and constant
    betas, resid = np.linalg.lstsq(X, neural_v, rcond=None)[:2]
    # get R^2 based on residuals
    r2 = 1 - resid / (neural_v.size * neural_v.var())
    # # Determine model (if interested in R)
    # for k in range(len(betas)-1):
    #     if k==0:
    #         model = X[:,k].reshape(-1,1)*betas[k]
    #     else:
    #         model = model + X[:,k].reshape(-1,1)*betas[k]
    # # Get multiple correlation coefficient (R)
    # R = pearsonr(neural_v.flatten(),model.flatten())[0]
    if val_label == 'R2':
        out = np.array(r2)
    else:
        out = np.concatenate((betas.flatten()[:-1], np.array(r2)))
    return(out)

def rsa_cor(neural_v, model_mat):
    out = spearmanr(neural_v, model_mat)
    out = out[0][0,1:]
    return(out)

def run_rsa_roi(roi_data, model_mat, corr='corr', val_label=None):
    # calculate roi RDM (lower triangle)
    roi_tri = pdist(roi_data, metric='correlation')
    roi_tri[np.isnan(roi_tri)] = 0.
    # run regression
    if corr=='corr':
        res = rsa_cor(neural_v=roi_tri, model_mat=model_mat)
    elif corr=='reg':
        res = rsa_reg(neural_v=roi_tri, model_mat=model_mat, val_label=val_label)
    out_dict = {'result':res}
    return(out_dict)

def run_rsa_searchlight(sl, model_mat, corr='corr', val_label=None):
    def calc_rsa(data, sl_mask, myrad, bcvar, model_mat=model_mat, val_label=val_label):
        data4D = data[0]
        labels = NODES
        bolddata_sl = data4D.reshape(sl_mask.shape[0] * sl_mask.shape[1] * sl_mask.shape[2], data[0].shape[3]).T
        res_dict = run_rsa_roi(bolddata_sl, model_mat, corr, val_label=val_label)
        res = res_dict['result']
        return res
    print(str(datetime.now())+ ": Begin Searchlight")
    sl_result = sl.run_searchlight(calc_rsa, pool_size=1)
    print(str(datetime.now())+ ": End Searchlight")
    return sl_result


# STUDY-SPECIFIC
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

def run_rsa_sub(sub, model_rdms, procedure, corr, tasks=TASKS, overwrite=False, val_label=None, pred=None):
    """
    sub: str subject ID
    model_rdms: dictionary with model names as keys and lower triangles of RDMs as elements
    tasks: array of task names: friend, number, avg
    """
    # check type of input arguments

    # subject id
    sub = convert2subid(sub)

    # get labels
    corr_label = 'spear' if corr=='corr' else 'reg'
    if val_label is None:
        val_label = 'r' if corr=='corr' else 'beta'
    elif val_label == "R2" and corr!='reg':
        print("ERROR: cannot calculate R2 with corr label '"+corr+"'. Must be 'reg'")
        exit(1)

    parc_label = SL+str(SL_RADIUS) if isSl(procedure) else str(N_PARCELS)
    # make output directories if they don't exist
    if not os.path.exists(out_dir % (sub, corr_label)):
        os.makedirs(out_dir % (sub, corr_label))

    # tasks to run
    if not isinstance(tasks, list):
        tasks = [tasks]

    # dictionary of model representational dissimilarity matrices (lower triangles)
    if not isinstance(model_rdms, dict):
        if isinstance(model_rdms, list):
            model_rdms={'model':model_rdms}
        else:
            print("ERROR: model_rdms input must be dictionary")
            exit(1)

    # procedure to run
    if procedure not in PROCEDURES_ALL:
        print("ERROR: procedure is not in PROCEDURES_ALL")
        print(PROCEDURES_ALL)
        exit(1)

    # reference image for saving parcellation output
    parcellation_template = nib.load(MNI_PARCELLATION, mmap=False)
    parcellation_template.set_data_dtype(np.double)

    # get model names
    model_keys = [pred] if val_label == "R2" else model_rdms.keys()
    print(str(datetime.now())+ ": Using the following models to get "+val_label+" values")
    print(model_rdms.keys())
    # turn model dictionary into matrix, s.t. column = model RDM lower triangle
    for i,k in enumerate(model_rdms.keys()):
        if i == 0:
            # if first model, setup matrix
            model_rdms_mat = [model_rdms[k]]
        else:
            # add next matrix as row
            model_rdms_mat = np.vstack((model_rdms_mat, model_rdms[k]))

    # transpose so that each column corresponds to each measure
    model_rdms_mat = np.transpose(model_rdms_mat)

    out_tasks = {}
    # iterate through inputted tasks
    for task in tasks:
        print(str(datetime.now()) + ": Task " + task)
        # read in subject's image
        sub_template = nib.load(data_fnames % (sub, sub, TASKS[0], 0), mmap=False)
        sub_dims = sub_template.get_data().shape + (N_NODES,)
        sub_data = np.empty(sub_dims)
        for n in range(N_NODES):
            print(str(datetime.now()) + ": Reading in node " + str(n))
            if task=='avg':
                # average this node's data from both runs
                d1 = load_nii(data_fnames % (sub, sub, TASKS[0], n))
                d2 = load_nii(data_fnames % (sub, sub, TASKS[1], n))
                d = (d1+d2)/2
            else:
                d = load_nii(data_fnames % (sub, sub, task, n))
            # save to fourth dimension
            sub_data[:,:,:,n] = d

        out_data_dict = {}
        if isParc(procedure):
            # out csv filename
            sub_csv_fname = csv_fname % (sub, corr_label, sub, task, corr_label, parc_label, val_label)
            if (not overwrite) and os.path.isfile(sub_csv_fname):
                    read_bool = True
                    # read in csv if already exists
                    sub_csv = pd.read_csv(sub_csv_fname)
                    # remove row column
                    sub_csv = sub_csv.iloc[:,1:]
                    # save all completed ROIs except for last one (since may not have been finished)
                    completed_rois = np.unique(sub_csv['roi'])
                    sub_csv = sub_csv[sub_csv['roi'].isin(completed_rois)]
                    sub_csv.to_csv(sub_csv_fname)
                    out_csv_array = sub_csv.values.tolist()
                    completed_preds = np.unique(sub_csv['predictor'])
            else:
                if os.path.isfile(sub_csv_fname):
                    os.remove(sub_csv_fname)
                    print("Deleted "+sub_csv_fname)
                read_bool = False
                out_csv_array = []
                completed_rois = []
                completed_preds = []
            wtr = csv.writer(open(sub_csv_fname, 'a'), delimiter=',', lineterminator='\n')
            # column names for csv file
            colnames = ['sub','task','roi','predictor',val_label]
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
            for i,k in enumerate(model_keys):
                out_data_dict[k] = deepcopy(out_data)
            # iterate through each ROI of parcellation and run regression
            for r, parc_roi in enumerate(roi_list):
                roi_done = parc_roi in completed_rois and all(mk in completed_preds for mk in model_keys)
                if roi_done:
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
                    res_dict = run_rsa_roi(roi_data, model_rdms_mat, corr=corr, val_label=val_label)
                    res = res_dict['result']
                # for each model, save the result to its image in out_data_dict
                for i,k in enumerate(model_keys):
                    # save to dataframe if not already there
                    if not roi_done:
                        val = res[i]
                        csv_row = [sub, task, parc_roi, k, val]
                        out_csv_array.append(csv_row)
                        # write out to csv
                        wtr.writerow(csv_row)
                    else:
                        val = res[i]
                    # update voxels
                    model_data = out_data_dict[k]
                    model_data[model_data==parc_roi] = val
                    out_data_dict[k] = model_data
        elif isSl(procedure):
            ref_img = sub_template
            # mask
            if use_mask:
                # load all functional masks and make largest mask
                t = task if task in TASKS else "*"
                print(str(datetime.now()) + ": Reading in masks")
                func_mask_names = glob.glob(sub_mask_fname % (sub, sub, t))
                for i,m in enumerate(func_mask_names):
                    print(m)
                    m_data = load_nii(m)
                    if i == 0:
                        m_sum = deepcopy(m_data)
                    else:
                        m_sum += m_data
                whole_brain_mask = np.where(m_sum > 0, 1, 0)
                whole_brain_mask_dil = binary_dilation(whole_brain_mask, iterations=int(SL_RADIUS)).astype(whole_brain_mask.dtype)
                mask = whole_brain_mask_dil
            else:
                mask = deepcopy(d)
                mask.fill(1)

            # Create the searchlight object
            begin_time = time.time()
            sl = Searchlight(sl_rad=sl_rad, max_blk_edge=max_blk_edge, shape=Ball)
            print(str(datetime.now())+ ": Setup searchlight inputs")
            print(str(datetime.now())+ ": Input data shape: " + str(sub_data.shape))

            # Distribute the information to the searchlights (preparing it to run)
            print(str(datetime.now())+ ": Distributing searchlight")
            sl.distribute([sub_data], mask)
            # Data that is needed for all searchlights is sent to all cores via
            # the sl.broadcast function. In this example, we are sending the
            # labels for classification to all searchlights.
            print(str(datetime.now())+ ": Broadcasting bcvar")
            sl.broadcast(bcvar)

            print(str(datetime.now())+ ": Shape of searchlight")
            print(sl.shape)
            # turn model dictionary into matrix, s.t. column = model RDM lower triangle
            sl_result = run_rsa_searchlight(sl, model_rdms_mat, corr=corr, val_label=val_label)
            end_time = time.time()

            # Print outputs
            print(str(datetime.now())+ ": Number of searchlights run: " + str(len(sl_result[mask==1])))
            print(str(datetime.now())+ ': Total searchlight duration (including start up time): %.2f min' % ((end_time - begin_time)/60))
            # separate values
            for i,k in enumerate(model_keys):
                out_data_dict[k] = deepcopy(sl_result)
                for x in range(sl_result.shape[0]):
                    for y in range(sl_result.shape[1]):
                        for z in range(sl_result.shape[2]):
                            val = sl_result[x,y,z]
                            if val is None:
                                out_val = 0.
                            else:
                                out_val = val[i]
                            out_data_dict[k][x,y,z] = out_val

                out_data_dict[k] = out_data_dict[k].astype('double')

            # unindent if saving images in parcellation (in addition to csv)
            # save images
            for k in out_data_dict.keys():
                fname = out_fname % (sub, corr_label, sub, task, corr_label, parc_label, val_label, k)
                save_nii( out_data_dict[k], ref_img, fname )
        # add to output array
        out_tasks[task] = out_data_dict
        out_tasks['ref_img'] = ref_img
    return(out_tasks)

#
# def run_rsa_sub(sub, model_rdms, procedure, corr, tasks=TASKS, overwrite=False, val_label=None):
#     """
#     sub: str subject ID
#     model_rdms: dictionary with model names as keys and lower triangles of RDMs as elements
#     tasks: array of task names: friend, number, avg
#     """
#     # check type of input arguments
#
#     # subject id
#     sub = convert2subid(sub)
#
#     # get labels
#     corr_label = 'spear' if corr=='corr' else 'reg'
#     if val_label is None:
#         val_label = 'r' if corr=='corr' else 'beta'
#
#     parc_label = SL+str(SL_RADIUS) if isSl(procedure) else str(N_PARCELS)
#     # make output directories if they don't exist
#     if not os.path.exists(out_dir % (sub, corr_label)):
#         os.makedirs(out_dir % (sub, corr_label))
#
#     # tasks to run
#     if not isinstance(tasks, list):
#         tasks = [tasks]
#
#     # dictionary of model representational dissimilarity matrices (lower triangles)
#     if not isinstance(model_rdms, dict):
#         if isinstance(model_rdms, list):
#             model_rdms={'model':model_rdms}
#         else:
#             print("ERROR: model_rdms input must be dictionary")
#             exit(1)
#
#     # procedure to run
#     if procedure not in PROCEDURES_ALL:
#         print("ERROR: procedure is not in PROCEDURES_ALL")
#         print(PROCEDURES_ALL)
#         exit(1)
#
#     # reference image for saving parcellation output
#     parcellation_template = nib.load(MNI_PARCELLATION, mmap=False)
#     parcellation_template.set_data_dtype(np.double)
#
#     # get model names
#     model_keys = model_rdms.keys()
#     print(str(datetime.now())+ ": Using models ")
#     print(model_keys)
#     # turn model dictionary into matrix, s.t. column = model RDM lower triangle
#     for i,k in enumerate(model_keys):
#         if i == 0:
#             # if first model, setup matrix
#             model_rdms_mat = [model_rdms[k]]
#         else:
#             # add next matrix as row
#             model_rdms_mat = np.vstack((model_rdms_mat, model_rdms[k]))
#
#     # transpose so that each column corresponds to each measure
#     model_rdms_mat = np.transpose(model_rdms_mat)
#
#     out_tasks = {}
#     # iterate through inputted tasks
#     for task in tasks:
#         # read in subject's image
#         sub_template = nib.load(data_fnames % (sub, sub, TASKS[0], 0), mmap=False)
#         sub_dims = sub_template.get_data().shape + (N_NODES,)
#         sub_data = np.empty(sub_dims)
#         for n in range(N_NODES):
#             print(str(datetime.now()) + ": Reading in node " + str(n))
#             if task=='avg':
#                 # average this node's data from both runs
#                 d1 = load_nii(data_fnames % (sub, sub, TASKS[0], n))
#                 d2 = load_nii(data_fnames % (sub, sub, TASKS[1], n))
#                 d = (d1+d2)/2
#             else:
#                 d = load_nii(data_fnames % (sub, sub, task, n))
#             # save to fourth dimension
#             sub_data[:,:,:,n] = d
#
#         out_data_dict = {}
#         if isSl(procedure):
#             ref_img = sub_template
#             # node labels
#             bcvar = None #range(N_NODES) #None
#             # searchlight radius
#             sl_rad = int(SL_RADIUS)
#             max_blk_edge = 5
#             # mask
#             if use_mask:
#                 # load all functional masks and make largest mask
#                 t = task if task in TASKS else "*"
#                 print(str(datetime.now()) + ": Reading in masks")
#                 func_mask_names = glob.glob(sub_mask_fname % (sub, sub, t))
#                 for i,m in enumerate(func_mask_names):
#                     print(m)
#                     m_data = load_nii(m)
#                     if i == 0:
#                         m_sum = deepcopy(m_data)
#                     else:
#                         m_sum += m_data
#                 whole_brain_mask = np.where(m_sum > 0, 1, 0)
#                 whole_brain_mask_dil = binary_dilation(whole_brain_mask, iterations=int(SL_RADIUS)).astype(whole_brain_mask.dtype)
#                 mask = whole_brain_mask_dil
#             else:
#                 mask = deepcopy(d)
#                 mask.fill(1)
#
#             # Create the searchlight object
#             begin_time = time.time()
#             sl = Searchlight(sl_rad=sl_rad, max_blk_edge=max_blk_edge, shape=Ball)
#             print(str(datetime.now())+ ": Setup searchlight inputs")
#             print(str(datetime.now())+ ": Input data shape: " + str(sub_data.shape))
#
#             # Distribute the information to the searchlights (preparing it to run)
#             print(str(datetime.now())+ ": Distributing searchlight")
#             sl.distribute([sub_data], mask)
#             # Data that is needed for all searchlights is sent to all cores via
#             # the sl.broadcast function. In this example, we are sending the
#             # labels for classification to all searchlights.
#             print(str(datetime.now())+ ": Broadcasting bcvar")
#             sl.broadcast(bcvar)
#
#             print(str(datetime.now())+ ": Shape of searchlight")
#             print(sl.shape)
#             # turn model dictionary into matrix, s.t. column = model RDM lower triangle
#             sl_result = run_rsa_searchlight(sl, model_rdms_mat, val_label=val_label)
#             end_time = time.time()
#
#             # Print outputs
#             print(str(datetime.now())+ ": Number of searchlights run: " + str(len(sl_result[mask==1])))
#             print(str(datetime.now())+ ': Total searchlight duration (including start up time): %.2f' % (end_time - begin_time))
#             # separate values
#             for i,k in enumerate(model_keys):
#                 out_data_dict[k] = deepcopy(sl_result)
#                 for x in range(sl_result.shape[0]):
#                     for y in range(sl_result.shape[1]):
#                         for z in range(sl_result.shape[2]):
#                             val = sl_result[x,y,z]
#                             if val is None:
#                                 out_val = 0.
#                             else:
#                                 out_val = val[i]
#                             out_data_dict[k][x,y,z] = out_val
#
#                 out_data_dict[k] = out_data_dict[k].astype('double')
#
#         elif isParc(procedure):
#             # out csv filename
#             sub_csv_fname = csv_fname % (sub, corr_label, sub, task, corr_label, parc_label)
#             if (not overwrite) and os.path.isfile(sub_csv_fname):
#                     read_bool = True
#                     # read in csv if already exists
#                     sub_csv = pd.read_csv(sub_csv_fname)
#                     # remove row column
#                     sub_csv = sub_csv.iloc[:,1:]
#                     # save all completed ROIs except for last one (since may not have been finished)
#                     completed_rois = np.unique(sub_csv['roi'])[:-1]
#                     sub_csv = sub_csv[sub_csv['roi'].isin(completed_rois)]
#                     sub_csv.to_csv(sub_csv_fname)
#                     out_csv_array = sub_csv.values.tolist()
#             else:
#                 if os.path.isfile(sub_csv_fname):
#                     os.remove(sub_csv_fname)
#                     print("Deleted "+sub_csv_fname)
#                 read_bool = False
#                 out_csv_array = []
#                 completed_rois = []
#             wtr = csv.writer(open(sub_csv_fname, 'a'), delimiter=',', lineterminator='\n')
#             # column names for csv file
#             colnames = ['sub','task','roi','predictor',val_label]
#             if not read_bool:
#                 # write out to csv
#                 wtr.writerow(colnames)
#             ref_img = parcellation_template
#             # make mask
#             parcellation = sub_parc % (sub, sub)
#             print(str(datetime.now()) + ": Using parcellation " + parcellation)
#             parc_data = load_nii(parcellation)
#             roi_list = np.unique(parc_data)
#             # remove 0 (i.e., the background)
#             roi_list = np.delete(roi_list,0)
#             # check if number of parcels matches global variable
#             if N_PARCELS != len(roi_list):
#                 print("WARNING: Number of parcels found ("+str(len(roi_list))+") does not equal N_PARCELS ("+str(N_PARCELS)+")")
#
#             # Run regression on each parcellation
#             print(str(datetime.now()) + ": Starting parcellation "+ str(N_PARCELS))
#             # get the voxels from parcellation nii
#             out_data = ref_img.get_data().astype(np.double)
#             # create a dictionary of nii's: one per predictor
#             for i,k in enumerate(model_keys):
#                 out_data_dict[k] = deepcopy(out_data)
#             # iterate through each ROI of parcellation and run regression
#             for r, parc_roi in enumerate(roi_list):
#                 if parc_roi in completed_rois:
#                     print(str(datetime.now()) + ': ROI '+str(parc_roi)+' already saved.')
#                     # read in values from dataframe for nii
#                     res = get_roi_csv_val(sub_csv, parc_roi, val_label)
#                 else:
#                     perc_done = round(((r+1) / len(roi_list)) * 100, 3)
#                     print(str(datetime.now()) + ': Analyzing ROI '+str(parc_roi)+' -- '+str(perc_done)+'%')
#                     # create mask for this ROI
#                     roi_mask = parc_data==parc_roi
#                     roi_mask = roi_mask.astype(int)
#                     roi_data = get_roi_data(sub_data, roi_mask)
#                     res_dict = run_rsa_roi(roi_data, model_rdms_mat, corr=corr, val_label=val_label)
#                     res = res_dict['result']
#                 # for each model, save the result to its image in out_data_dict
#                 for i,k in enumerate(model_keys):
#                     # save to dataframe if not already there
#                     if parc_roi not in completed_rois:
#                         val = res[i]
#                         csv_row = [sub, task, parc_roi, k, val]
#                         out_csv_array.append(csv_row)
#                         # write out to csv
#                         wtr.writerow(csv_row)
#                     else:
#                         val = res[k]
#                     # update voxels
#                     model_data = out_data_dict[k]
#                     model_data[model_data==parc_roi] = val
#                     out_data_dict[k] = model_data
#         # save images
#         for k in out_data_dict.keys():
#             fname = out_fname % (sub, corr_label, sub, task, corr_label, parc_label, val_label, k)
#             save_nii( out_data_dict[k], ref_img, fname )
#         # add to output array
#         out_tasks[task] = out_data_dict
#         out_tasks['ref_img'] = ref_img
#     return(out_tasks)
#


# def run_second_order_rsa(subs, model_rdms, procedure, corr, tasks=['avg'], overwrite=False):
#     """
#     sub: list of subject IDs (str)
#     model_rdms: dictionary with model names as keys and lower triangles of RDMs as elements
#     tasks: array of task names: friend, number, avg
#     """
#     # check type of input arguments
#
#     # subject id
#     for i,s in enumerate(subs):
#         s = convert2subid(s)
#         subs[i] = s
#
#     # get labels
#     corr_label = 'spear' if corr=='corr' else 'reg'
#     val_label = 'r' if corr=='corr' else 'beta'
#     parc_label = SL+str(SL_RADIUS) if isSl(procedure) else str(N_PARCELS)
#
#
#     out_dir = SECOND_LEVEL_DIR + '%s/' #%(corr_label)
#     out_fname = out_dir + 'order-2_task-%s_stat-'+STAT+'_corr-%s_parc-%s_val-%s_pred-%s.nii.gz' #% (sub, task, corr, N_PARCELS, "r" or "b", predictor_name)
#     csv_fname = out_dir + "order-2_task-%s_stat-"+STAT+"_corr-%s_parc-%s_roi_stats.csv"
#     # make output directories if they don't exist
#     if not os.path.exists(out_dir % (corr_label)):
#         os.makedirs(out_dir % (corr_label))
#
#     # tasks to run
#     if not isinstance(tasks, list):
#         tasks = [tasks]
#
#     # dictionary of model representational dissimilarity matrices (lower triangles)
#     if not isinstance(model_rdms, dict):
#         if isinstance(model_rdms, list):
#             model_rdms={'model':model_rdms}
#         else:
#             print("ERROR: model_rdms input must be dictionary")
#             exit(1)
#
#     # procedure to run
#     if procedure not in PROCEDURES_ALL:
#         print("ERROR: procedure is not in PROCEDURES_ALL")
#         print(PROCEDURES_ALL)
#         exit(1)
#
#     # reference image for saving parcellation output
#     parcellation_template = nib.load(MNI_PARCELLATION, mmap=False)
#     parcellation_template.set_data_dtype(np.double)
#
#     # get model names
#     model_keys = model_rdms.keys()
#     print(str(datetime.now())+ ": Using models ")
#     print(model_keys)
#     # turn model dictionary into matrix, s.t. column = model RDM lower triangle
#     for i,k in enumerate(model_keys):
#         if i == 0:
#             # if first model, setup matrix
#             model_rdms_mat = [model_rdms[k]]
#         else:
#             # add next matrix as row
#             model_rdms_mat = np.vstack((model_rdms_mat, model_rdms[k]))
#     # transpose so that each column corresponds to each measure
#     model_rdms_mat = np.transpose(model_rdms_mat)
#
#     out_tasks = {}
#     # iterate through inputted tasks
#     # setup out data
#     ref_img = parc_template
#     # get the voxels from parcellation nii
#     out_data = ref_img.get_data().astype(np.double)
#     # setup dictionary of output images
#     out_data_dict = {}
#     # out csv filename
#
#     for sub in subs:
#         # read in subject's image
#         sub_template = nib.load(data_fnames % (sub, sub, TASKS[0], 0), mmap=False)
#         sub_dims = sub_template.get_data().shape + (N_NODES,)
#         sub_data = np.empty(sub_dims)
#         for n in range(N_NODES):
#             print(str(datetime.now()) + ": Reading in node " + str(n))
#             if task=='avg':
#                 # average this node's data from both runs
#                 d1 = load_nii(data_fnames % (sub, sub, TASKS[0], n))
#                 d2 = load_nii(data_fnames % (sub, sub, TASKS[1], n))
#                 d = (d1+d2)/2
#             else:
#                 d = load_nii(data_fnames % (sub, sub, task, n))
#             # save to fourth dimension
#             sub_data[:,:,:,n] = d
#
#         # load subject's whole brain mask
#         sub_mask = np.round(load_nii(sub_mask_fname % (sub, sub))) if isSl() else np.round(load_nii(MNI_PARCELLATION))
#         parc_data = deepcopy(sub_mask) if isSl() else load_nii(MNI_PARCELLATION)
#
#         if isSl(procedure):
#             # start parcel count
#             parc_roi = 0
#             # go through all voxels in sub_mask
#             for index, val in np.ndenumerate(parc_data):
#                 # if this voxel has a non-zero value, analyze it
#                 if val != 0:
#                     # save integer to parcellation just for the center voxel
#                     parc_data[index] = parc_roi
#                     # use next integer for next voxel
#                     parc_roi += 1
#         # get list of all ROIs
#         roi_list = np.unique(parc_data)
#         # remove 0 (i.e., the background)
#         roi_list = np.delete(roi_list,0)
#
#         for task in tasks:
#             # add a key for each model tested
#             out_data_dict[task] = deepcopy(out_data)
#             # out csv filename
#             wtr = csv.writer(open(out_csv_fname, 'a'), delimiter=',', lineterminator='\n')
#             # column names for csv file
#             colnames = ['x','y','z','task','roi','r','p']
#             if (not overwrite) and os.path.isfile(out_csv_fname):
#                 read_bool = True
#                 # read in csv if already exists
#                 in_csv = pd.read_csv(out_csv_fname)
#                 # remove row column
#                 in_csv = in_csv.iloc[:,1:]
#                 # save all completed ROIs except for last one (since may not have been finished)
#                 completed_rois = np.unique(in_csv['roi'])[:-1]
#                 in_csv = in_csv[in_csv['roi'].isin(completed_rois)]
#                 in_csv.to_csv(out_csv_fname)
#                 out_csv_array = in_csv.values.tolist()
#             else:
#                 if os.path.isfile(out_csv_fname):
#                     os.remove(out_csv_fname)
#                     print("Deleted "+out_csv_fname)
#                 read_bool = False
#                 out_csv_array = []
#                 completed_rois = []
#                 wtr.writerow(colnames)
#
#             # run correlation in each ROI
#             for r, parc_roi in enumerate(roi_list):
#                 if parc_roi in completed_rois:
#                     print(str(datetime.now()) + ': ROI '+str(parc_roi)+' already saved.')
#                     # read in values from dataframe for nii
#                     res = get_roi_csv_val(in_csv, parc_roi, val_label)
#                 else:
#                     perc_done = round((parc_roi / len(roi_list)) * 100, 3)
#                     print(str(datetime.now()) + ': Analyzing ROI '+str(index)+' -- '+str(perc_done)+'%')
#                     # get mask for this roi
#                     mask_roi = parc_data == parc_roi
#                     mask_roi = roi_mask.astype(int)
#                     index = np.where(mask_roi == 1)
#                     mask_roi = get_sphere_mask(index=index, brain_mask=mask) if isSl() else mask_roi
#                     midi = int(len(index[0])/2)
#                     x = index[0][midi]
#                     y = index[1][midi]
#                     z = index[2][midi]
#                     # get the subject's data that is in this sphere
#                     roi_data = get_roi_data(data, mask_roi)
#                     # run classification
#                     classify_with_nilearn(mask=mask_roi, data = roi_data, labels=labels_dict[lab], classify='svm', groups=subIDs)
#                     # get the resulting values
#                     res = res_dict['result']
#                     # save to dataframe if not already there
#                     if parc_roi not in completed_rois:
#                         # array
#                         csv_row = [x, y, z, parc_roi, k, res]
#                         out_csv_array.append(csv_row)
#                         # write out to csv
#                         wtr.writerow(csv_row)
#                     else:
#                         # dictionary
#                         val = res[k]
#                     # update voxels
#                     model_data = out_data_dict[k]
#                     model_data[model_data==parc_roi] = val
#                     out_data_dict[k] = model_data
#
#
#
#
#             if isSl(procedure):
#                 # column names for csv file
#                 colnames = ['sub','x','y','z','task', 'roi','predictor',val_label]
#                 if not read_bool:
#                     # write out to csv
#                     wtr.writerow(colnames)
#                 ref_img = sub_template
#                 # load subject's whole brain mask
#                 sub_mask = np.round(load_nii(sub_mask_fname % (sub, sub)))
#                 # setup parcellation-like image where each voxel that is the center of
#                 # a sphere is assigned an integer
#                 parc_data = deepcopy(sub_mask)
#                 # setup out data
#                 out_data = deepcopy(sub_data)*0.0
#                 # setup dictionary of output images
#                 out_data_dict = {}
#                 # add a key for each model tested
#                 for i,k in enumerate(model_keys):
#                     out_data_dict[k] = deepcopy(out_data)
#                 # start parcel count
#                 parc_roi = 0
#                 nvox = np.sum(sub_mask)
#                 # go through all voxels in sub_mask
#                 for index, val in np.ndenumerate(sub_mask):
#                     # if this voxel has a non-zero value, analyze it
#                     if val != 0:
#                         # save integer to parcellation just for the center voxel
#                         parc_data[index] = parc_roi
#                         # use next integer for next voxel
#                         parc_roi += 1
#                         if parc_roi in completed_rois:
#                             print(str(datetime.now()) + ': ROI '+str(parc_roi)+' already saved.')
#                             # read in values from dataframe for nii
#                             res = get_roi_csv_val(sub_csv, parc_roi, val_label)
#                         else:
#                             perc_done = round((parc_roi / nvox) * 100, 3)
#                             print(str(datetime.now()) + ': Analyzing sphere '+str(index)+' -- '+str(perc_done)+'%')
#                             # get mask for this roi
#                             mask_roi = get_sphere_mask(index=index, brain_mask=sub_mask)
#                             # get the subject's data that is in this sphere
#                             roi_data = get_roi_data(sub_data, mask_roi)
#                             # run correlations comparing this sphere of data to each model
#                             res_dict = run_rsa_roi(roi_data, model_rdms_mat, corr=corr)
#                             # get the resulting values
#                             res = res_dict['result']
#                         # for each model, save the result to its image in out_data_dict
#                         for i,k in enumerate(model_keys):
#                             # save to dataframe if not already there
#                             if parc_roi not in completed_rois:
#                                 # array
#                                 val = res[i]
#                                 csv_row = [sub, task, x, y, z, parc_roi, k, val]
#                                 out_csv_array.append(csv_row)
#                                 # write out to csv
#                                 wtr.writerow(csv_row)
#                             else:
#                                 # dictionary
#                                 val = res[k]
#                             # update voxels
#                             model_data = out_data_dict[k]
#                             model_data[model_data==parc_roi] = val
#                             out_data_dict[k] = model_data
#
#             elif isParc(procedure):
#                 # column names for csv file
#                 colnames = ['sub','task','roi','predictor',val_label]
#                 if not read_bool:
#                     # write out to csv
#                     wtr.writerow(colnames)
#                 ref_img = parcellation_template
#                 # make mask
#                 parcellation = sub_parc % (sub, sub)
#                 print(str(datetime.now()) + ": Using parcellation " + parcellation)
#                 parc_data = load_nii(parcellation)
#                 roi_list = np.unique(parc_data)
#                 # remove 0 (i.e., the background)
#                 roi_list = np.delete(roi_list,0)
#                 # check if number of parcels matches global variable
#                 if N_PARCELS != len(roi_list):
#                     print("WARNING: Number of parcels found ("+str(len(roi_list))+") does not equal N_PARCELS ("+str(N_PARCELS)+")")
#
#                 # Run regression on each parcellation
#                 print(str(datetime.now()) + ": Starting parcellation "+ str(N_PARCELS))
#                 # get the voxels from parcellation nii
#                 out_data = ref_img.get_data().astype(np.double)
#                 # create a dictionary of nii's: one per predictor
#                 out_data_dict = {}
#                 for i,k in enumerate(model_keys):
#                     out_data_dict[k] = deepcopy(out_data)
#                 # iterate through each ROI of parcellation and run regression
#                 for r, parc_roi in enumerate(roi_list):
#                     if parc_roi in completed_rois:
#                         print(str(datetime.now()) + ': ROI '+str(parc_roi)+' already saved.')
#                         # read in values from dataframe for nii
#                         res = get_roi_csv_val(sub_csv, parc_roi, val_label)
#                     else:
#                         perc_done = round(((r+1) / len(roi_list)) * 100, 3)
#                         print(str(datetime.now()) + ': Analyzing ROI '+str(parc_roi)+' -- '+str(perc_done)+'%')
#                         # create mask for this ROI
#                         roi_mask = parc_data==parc_roi
#                         roi_mask = roi_mask.astype(int)
#                         roi_data = get_roi_data(sub_data, roi_mask)
#                         res_dict = run_rsa_roi(roi_data, model_rdms_mat, corr=corr)
#                         res = res_dict['result']
#                     # for each model, save the result to its image in out_data_dict
#                     for i,k in enumerate(model_keys):
#                         # save to dataframe if not already there
#                         if parc_roi not in completed_rois:
#                             val = res[i]
#                             csv_row = [sub, task, parc_roi, k, val]
#                             out_csv_array.append(csv_row)
#                             # write out to csv
#                             wtr.writerow(csv_row)
#                         else:
#                             val = res[k]
#                         # update voxels
#                         model_data = out_data_dict[k]
#                         model_data[model_data==parc_roi] = val
#                         out_data_dict[k] = model_data
#             # save images
#             for k in out_data_dict.keys():
#                 fname = out_fname % (sub, corr_label, sub, task, corr_label, parc_label, val_label, k)
#                 save_nii( out_data_dict[k], ref_img, fname )
#             # add to output array
#             out_tasks[task] = out_data_dict
#     return(out_tasks)
