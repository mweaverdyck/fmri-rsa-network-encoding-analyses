#!/bin/python
# funcs.py
import os
import sys
from os.path import expanduser
HOME = expanduser("~")
import glob
import pandas as pd
import csv
import nibabel as nib
import numpy as np
import copy
import shutil
from datetime import datetime
from collections import OrderedDict

# keywords
ALL = os.environ['ALL']
NEW = os.environ['NEW']
PARC = os.environ['PARC']
SL = os.environ['SL']
MNI_SPACE = os.environ['MNI_SPACE']
T1_SPACE = os.environ['T1_SPACE']

# arrays that must be manually updated
EXCLUDE_SUBS=["204"]
TASKS=['friend', 'number']
PROCEDURES_ALL=[PARC, SL]
DISTANCES=['correlation', 'euclidean']

# get project directories
PROJECT_DIR = os.environ['PROJECT_DIR'] + '/'
BIDS_DIR = os.environ['BIDS_DIR'] + '/'
CODE_DIR = os.environ['CODE_DIR'] + '/'
DERIVATIVES_DIR = os.environ['DERIVATIVES_DIR'] + '/'
QA_DIR = os.environ['QA_DIR'] + '/'
RECON_DIR = os.environ['RECON_DIR'] + '/'
FMRIPREP_DIR = os.environ['FMRIPREP_DIR'] + '/'
FIRST_LEVEL_DIR = os.environ['FIRST_LEVEL_DIR'] + '/'
SECOND_LEVEL_DIR = os.environ['SECOND_LEVEL_DIR'] + '/'
DERIVS_DIR = os.environ['DERIVS_DIR'] + '/'
GLM_DIR = os.environ['GLM_DIR'] + '/'
RSA_DIR = os.environ['RSA_DIR'] + '/'
MNI_DIR = os.environ['MNI_DIR'] + '/'

# get project variables
EXCLUDE_RUNS_THRESH = float(os.environ['EXCLUDE_RUNS_THRESH'])
SUBID_PREFIX = os.environ['SUBID_PREFIX']
N_NODES = int(os.environ['N_NODES'])
NODES = list(range(N_NODES))
N_TRS = int(os.environ['N_TRS'])
TR = float(os.environ['TR'])
N_TRS_DEL = int(os.environ['N_TRS_DEL'])
N_PARCELS = int(os.environ['N_PARCELS'])
N_NETWORKS = int(os.environ['N_NETWORKS'])
PARC_LAB = os.environ['PARC_LAB']
MNI_PARCELLATION = os.environ['MNI_PARCELLATION']
MNI_GM_MASK = os.environ['MNI_GM_MASK']
MNI_MASK = os.environ['MNI_MASK']
MNI_MASK_DIL = os.environ['MNI_MASK_DIL']
SPACE = os.environ['SPACE']
STAT = os.environ['STAT']
SL_RADIUS = os.environ['SL_RADIUS']
PROCEDURE = os.environ['PROCEDURE']
if PROCEDURE == ALL:
    PROCEDURES = PROCEDURES_ALL
else:
    PROCEDURES = [PROCEDURE]

# filenames
CFD_FNAME = os.environ['CFD_FNAME']
ATTS_FNAME = os.environ['ATTS_FNAME']

# are you on a local computer or Hoffman2
if HOME[:3] == '/u/':
    isLocal=0
else:
    isLocal=1

def convert2subid(s):
    sub=SUBID_PREFIX+str(s).split(SUBID_PREFIX)[-1]
    return sub

def get_bids_att(str, label):
    arr = str.split('_')
    for a in arr:
        if label == a.split('-')[0]:
            return a.split('-')[1]
    return None

def get_all_bids_atts(fname):
    fname = fname.split('/')[-1]
    fname = fname.split('.')[0]
    arr = fname.split('_')
    arr_dict = OrderedDict()
    for a in arr:
        a2 = a.split('-')
        if len(a2) == 2:
            arr_dict[a2[0]] = a2[1]
        else:
            arr_dict['extra'] = a2[0]
    return arr_dict

def make_bids_str(odict):
    out_str = ''
    if 'extra' in odict.keys():
        odict.move_to_end('extra')
    for i,a in enumerate(odict):
        out_str += '_' if i != 0 else ''
        if a == 'extra':
            out_str += odict[a]
        else:
            out_str += a + '-' + odict[a]
    return out_str

def isParc(procedure):
    return(procedure == PARC)

def isSl(procedure):
    return(procedure == SL)

def load_nii( imgpath ):
    img_np_array = nib.load(imgpath, mmap=False).get_data()
    return img_np_array

def save_nii( data, refnii, filename ):
    out_img = nib.Nifti1Image(data, refnii.affine, refnii.header)
    out_img.to_filename(filename)
    print(str(datetime.now()) + ": File %s saved." % filename)

def get_thresh_dir( dir_func ):
    if EXCLUDE_RUNS_THRESH != 0:
        dir_func = os.path.join(dir_func, "excl-%02d" % (EXCLUDE_RUNS_THRESH*100))
    if not os.path.exists(dir_func):
        os.makedirs(dir_func)
    return dir_func

def get_node_mapping( sub ):
    sub = convert2subid(sub)
    print("Subject "+sub+"'s mapping")
    events_file = glob.glob(BIDS_DIR+sub+"/func/"+sub+"*events.tsv")[0]
    events_df = pd.read_csv(events_file, delimiter='\t')
    node_mapping_sub = {}
    for n in range(N_NODES):
        # find first row that matches this node
        r = events_df.loc[events_df['node']==n].index[0]
        # find corresponding stimulus file name
        stim = events_df.loc[r,'stim_file']
        # remove '.png'
        stim = stim.split('.')[0]
        # save to mapping
        node_mapping_sub[stim] = n
    # order node numbers to match CFD measures
    node_order = []
    cfd_df = pd.read_csv(CFD_FNAME)
    cfd_targets = cfd_df['Target']
    for t in cfd_targets:
        study_nodeID = cfd_df.loc[cfd_df['Target']==t]['NodeID_study']
        sub_nodeID = node_mapping_sub[t]
        print(t + ' = study ' + str(study_nodeID) + ' = sub '+str(sub_nodeID))
        node_order.append(node_mapping_sub[t])
    return node_order

def get_stim_from_node(sub, node):
    sub = convert2subid(sub)
    events_file = glob.glob(os.path.join(BIDS_DIR,sub,"func",sub+"*events.tsv"))[0]
    events_df = pd.read_csv(events_file, delimiter='\t')
    # find first row that matches this node
    r = events_df.loc[events_df['node']==int(node)].index[0]
    # find corresponding stimulus file name
    stim = events_df.loc[r,'stim_file']
    # remove '.png'
    stim = stim.split('.')[0]
    return stim

# data-extraction functions (e.g. in searchlight)

def get_sphere_template(radius):
    sphere_dim = 2 * radius + 1
    sphere_template = np.zeros((sphere_dim,sphere_dim,sphere_dim))
    center_vox = (radius,radius,radius)
    for xs in range(sphere_dim):
        for ys in range(sphere_dim):
            for zs in range(sphere_dim):
                d = np.sqrt((xs-center_vox[0])**2 + (ys-center_vox[1])**2 + (zs-center_vox[2])**2)
                if d <= radius:
                    sphere_template[xs,ys,zs] = 1
    return sphere_template

def get_sphere_mask(index, brain_mask, radius=SL_RADIUS):
    try:
        sphere_template
    except NameError:
        print("No sphere template exists yet, creating one:")
        sphere_template = get_sphere_template(radius)
    # get index coordinates
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
    return(mask_roi)

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

def exclude_runs(sid):
    sid = convert2subid(sid)
    events_dir = os.path.join(BIDS_DIR, sid, 'func')
    derivs_dir = os.path.join(DERIVS_DIR, sid, 'func')
    # define exclude folders
    exclude_events_dir = os.path.join(events_dir,"exclude")
    exclude_derivs_dir = os.path.join(derivs_dir, 'exclude')
    # make exclude directories
    if not os.path.exists(exclude_events_dir):
        os.makedirs(exclude_events_dir)
    if not os.path.exists(exclude_derivs_dir):
        os.makedirs(exclude_derivs_dir)
    # move events that are currently in exclude back to func
    move_fnames = glob.glob(os.path.join(exclude_events_dir, "*events.tsv"))
    for f in move_fnames:
        shutil.move(f, events_dir)
    move_fnames = glob.glob(os.path.join(exclude_derivs_dir, "*?.???*"))
    for f in move_fnames:
        shutil.move(f, derivs_dir)
    if EXCLUDE_RUNS_THRESH != 0:
        # get all events files
        events_fnames = glob.glob(os.path.join(events_dir, "*events.tsv"))
        acc_df = pd.DataFrame(events_fnames, columns=["filename"])
        acc_df["accuracy"] = np.nan
        pass_colname = "pass_%02d"%(EXCLUDE_RUNS_THRESH*100)
        acc_df[pass_colname] = True
        for f in events_fnames:
            # read events file
            df = pd.read_csv(f, sep='\t')
            # calculate accuracy
            acc = df.mean(axis=0, skipna=True)['correct']
            print("Accuracy for file "+f+" = "+str(acc))
            acc_df.loc[acc_df['filename'] == f, ['accuracy']] = acc
            # if under threshold, move to exclude
            if acc < EXCLUDE_RUNS_THRESH:
                # save fail to dataframe
                acc_df.loc[acc_df['filename'] == f, [pass_colname]] = False
                # move file to exclude folder
                shutil.move(f, exclude_events_dir)
                # get run
                run_label = get_bids_att(f, 'run')
                nii_fnames = glob.glob(os.path.join(derivs_dir, "*run-"+run_label+"*space-"+SPACE+"_rmtr-6_desc-preproc_bold.nii.gz"))
                regs_fnames = glob.glob(os.path.join(derivs_dir, "*run-"+run_label+"*_rmtr-6_regressors.tsv"))
                if len(nii_fnames) != 1 or len(regs_fnames) != 1:
                    print("ERROR: not exactly 1 nii or regression file found. skipping...")
                    print(os.path.join(derivs_dir, "*run-"+run_label+"*space-"+SPACE+"_rmtr-6_desc-preproc_bold.nii.gz"))
                    print(os.path.join(derivs_dir, "*run-"+run_label+"*_rmtr-6_regressors.tsv"))
                else:
                    print("Moving these files to exclude directories")
                    shutil.move(nii_fnames[0], exclude_derivs_dir)
                    shutil.move(regs_fnames[0], exclude_derivs_dir)
                print(f)
                print(nii_fnames)
                print(regs_fnames)
        # write accuracy dataframe
        f = os.path.join(events_dir,"run_accuracy_exclthresh-%02d.csv"%(EXCLUDE_RUNS_THRESH*100))
        acc_df.to_csv(f, index=False)
        print(f + ' saved.')
