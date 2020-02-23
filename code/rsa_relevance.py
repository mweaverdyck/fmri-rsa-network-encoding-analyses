#!/bin/python
from funcs import *
from rsa_funcs import *

"""
This script runs RSA within each task using social network models as predictors.
Neural data is correlated with each model then the difference between the correlations
is calculated to determine how much more social network aspect was encoded when
it was relevant than when it was not.
"""

print(str(datetime.now()) + ": Begin rsa_relevance.py")

all_sub = []
tasks=[]
procedure=PROCEDURES[0]
# read in arguments
for arg in sys.argv[1:]:
    print(arg)
    if arg in PROCEDURES_ALL:
        procedure = arg
    elif arg in TASKS or arg == 'avg':
        tasks += [arg]
    else:
        all_sub += [arg]

if len(tasks) == 0:
    tasks = ['avg'] + TASKS

# overwrite previous csv files? If False, will read in csv and append
overwrite = True

# print variables to log file
print(str(datetime.now()) + ": Project directory = " + PROJECT_DIR)
print(str(datetime.now()) + ": Analyzing " + str(len(all_sub)) + " subjects: " + str(all_sub))

# set variables
corr = 'corr'
# social network RDM dictionary
sn_rdms = {deg_label: deg_tri, dist_label: dist_tri}


# run regression for each subject
for s in all_sub:
    # get correct form of subID
    sub=convert2subid(s)
    print(str(datetime.now()) + ": Analyzing subject %s" % sub)
    print(str(datetime.now()) + ": Reading in data for subject " + sub)

    # run rsa
    res_dict = run_rsa_sub(sub=sub, model_rdms=sn_rdms, procedure=procedure, corr=corr, overwrite=overwrite)

    # get labels
    corr_label = 'spear' if corr=='corr' else 'reg'
    val_label = 'r' if corr=='corr' else 'beta'
    parc_label = SL+str(SL_RADIUS) if isSl(procedure) else str(N_PARCELS)

    # calculate difference images
    data_rel = res_dict['number'][deg_label] - res_dict['friend'][deg_label]
    fname = out_fname % (sub, corr_label, sub, 'diff', corr_label, parc_label, val_label, deg_label)
    save_nii( data_rel, res_dict['ref_img'], fname )

    data_rel = res_dict['friend'][dist_label] - res_dict['number'][dist_label]
    fname = out_fname % (sub, corr_label, sub, 'diff', corr_label, parc_label, val_label, dist_label)
    save_nii( data_rel, res_dict['ref_img'], fname )

print(str(datetime.now()) + ": End rsa_relevance.py")
