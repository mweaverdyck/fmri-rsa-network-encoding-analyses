#!/bin/python
from funcs import *
from rsa_funcs import *

print(str(datetime.now()) + ": Begin rsa_ind_meas.py")

# read in arguments
preds_all = cfd_phys + cfd_soc + [deg_lab, dist_lab, dist2_lab]
all_sub = []
tasks_all = ['avg'] + TASKS
tasks=[]
preds = []
# read in arguments
for arg in sys.argv[1:]:
    print(arg)
    if arg in PROCEDURES_ALL:
        procedure = arg
    elif arg in tasks_all:
        tasks += [arg]
    elif arg in preds_all:
        preds += [arg]
    else:
        all_sub += [arg]

if len(tasks) == 0:
    tasks = tasks_all

if len(preds) == 0:
    preds = preds_all

# for i, pred in enumerate(preds):
#     if pred == deg_lab:
#         preds[i] = deg_label
#     elif pred == dist_lab:
#         preds[i] = dist_label
#     elif pred == dist2_lab:
#         preds[i] = dist2_label

print(preds)
# overwrite previous csv files? If False, will read in csv and append
overwrite = True

print(str(datetime.now()) + ": Analyzing " + str(len(all_sub)) + " subjects: " + str(all_sub))

# set variables
corr = 'corr' # 'corr' or 'reg'

# social network RDM dictionary
sn_rdms = {deg_label: deg_tri, dist_label: dist_tri}

# run regression for each subject
for s in all_sub:
    # get correct form of subID
    sub=convert2subid(s)
    print(str(datetime.now()) + ": Analyzing subject %s" % sub)
    print(str(datetime.now()) + ": Reading in data for subject " + sub)


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
    all_model_rdms = {**soc_rdms, **phys_rdms, **sn_rdms}
    sel_model_rdms = {}
    for m in all_model_rdms.keys():
        if m in preds or m.split('_')[0] in preds:
            print("Selecting " + m)
            sel_model_rdms[m] = all_model_rdms[m]

    try:
        procedure
    except NameError:
        procedures = PROCEDURES
    else:
        procedures = [procedure]

    # run rsa
    for procedure in procedures:
        for task in tasks:
            if task == 'avg':
                model_rdms = sel_model_rdms
            elif task in TASKS:
                model_rdms = sn_rdms
            run_rsa_sub(sub=sub, model_rdms=model_rdms, procedure=procedure, corr=corr, overwrite=overwrite, tasks=task)

print(str(datetime.now()) + ": End rsa_ind_meas.py")
