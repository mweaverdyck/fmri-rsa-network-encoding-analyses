#!/bin/python
from funcs import *
from rsa_funcs import *

print(str(datetime.now()) + ": Begin rsa_all.py")

# read in arguments

all_sub = []
tasks=[]
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
    tasks = TASKS

# overwrite previous csv files? If False, will read in csv and append
overwrite = True

# print variables to log file
print(str(datetime.now()) + ": Project directory = " + PROJECT_DIR)
print(str(datetime.now()) + ": Analyzing " + str(len(all_sub)) + " subjects: " + str(all_sub))

# set variables
corr = 'corr'
deg_label = 'deg_cat-sn'
dist_label = 'dist_cat-sn'
dist2_label = 'dist2_cat-sn'

# chicago face database measures
all_cols = ['NodeID_study', 'Target', 'Race', 'Gender', 'Age', 'NumberofRaters', 'Female_prop', 'Male_prop', 'Asian_prop', 'Black_prop', 'Latino_prop', 'Multi_prop', 'Other_prop', 'White_prop', 'z', 'Angry', 'Attractive', 'Babyface', 'Disgusted', 'Dominant', 'Feminine', 'Happy', 'Masculine', 'Prototypic', 'Sad', 'Suitability', 'Surprised', 'Threatening', 'Trustworthy', 'Unusual', 'Luminance_median', 'Nose_Width', 'Nose_Length', 'Lip_Thickness', 'Face_Length', 'R_Eye_H', 'L_Eye_H', 'Avg_Eye_Height', 'R_Eye_W', 'L_Eye_W', 'Avg_Eye_Width', 'Face_Width_Cheeks', 'Face_Width_Mouth', 'Forehead', 'Pupil_Top_R', 'Pupil_Top_L', 'Asymmetry_pupil_top', 'Pupil_Lip_R', 'Pupil_Lip_L', 'Asymmetry_pupil_lip', 'BottomLip_Chin', 'Midcheek_Chin_R', 'Midcheek_Chin_L', 'Cheeks_avg', 'Midbrow_Hairline_R', 'Midbrow_Hairline_L', 'Faceshape', 'Heartshapeness', 'Noseshape', 'LipFullness', 'EyeShape', 'EyeSize', 'UpperHeadLength', 'MidfaceLength', 'ChinLength', 'ForeheadHeight', 'CheekboneHeight', 'CheekboneProminence', 'FaceRoundness', 'fWHR']

#cfd_soc = ['Dominant', 'Trustworthy']
#cfd_phys = ['Unusual', 'Faceshape']

cfd_phys = ['Race', 'Gender', 'Age', 'Luminance_median', 'fWHR']
cfd_soc = ['Angry', 'Attractive', 'Babyface', 'Disgusted', 'Dominant', 'Feminine', 'Happy', 'Masculine', 'Prototypic', 'Sad', 'Surprised', 'Threatening', 'Trustworthy', 'Unusual']


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


# social network RDM dictionary
sn_rdms = {deg_label: deg_tri, dist_label: dist_tri} #, dist2_label: dist2_tri}


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
    all_model_rdms = {**sn_rdms, **soc_rdms, **phys_rdms}

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
                model_rdms = all_model_rdms
            elif task in TASKS:
                model_rdms = sn_rdms
            run_rsa_sub(sub=sub, model_rdms=model_rdms, procedure=procedure, corr=corr, overwrite=overwrite, tasks=task)

print(str(datetime.now()) + ": End rsa_all.py")
