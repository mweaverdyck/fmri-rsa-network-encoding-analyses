#!/bin/python
from funcs import *
from rsa_funcs import *

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
out_fname = out_dir + '%s_task-%s_stat-'+STAT+'_corr-%s_parc-%s_val-%s_pred-%s.nii.gz' #% (sub, task, corr, N_PARCELS, "r" or "b", predictor_name)
csv_fname = out_dir + "%s_task-%s_stat-"+STAT+"_corr-%s_parc-%s_roi_stats.csv"

corr = 'corr'
corr_label = 'spear' if corr=='corr' else 'reg'
val_label = 'r' if corr=='corr' else 'beta'
# set label for filename
parc_label = SL if isSl(procedure) else str(N_PARCELS)

# set variables
deg_label = 'deg'
dist_label = 'dist'
rad = 5
sphere_template = get_sphere_template(rad)

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
for s in all_sub:
    # get correct form of subID
    sub=convert2subid(s)
    print(str(datetime.now()) + ": Analyzing subject %s" % sub)

    # make output directories if they don't exist
    if not os.path.exists(out_dir % sub):
        os.makedirs(out_dir % sub)

    print(str(datetime.now()) + ": Reading in data for subject " + sub)

    # combine all predictor RDM dictionaries
    model_rdms = {**sn_rdms}
    
    model_keys = model_rdms.keys()
    print(str(datetime.now())+ ": Using models ")
    print(model_keys)
    # turn dictionary into matrix
    for i,k in enumerate(model_keys):
        if i == 0:
            model_rdms_mat = [model_rdms[k]]
        else:
            model_rdms_mat = np.vstack((model_rdms_mat, model_rdms[k]))
    # transpose so that each column corresponds to each measure
    model_rdms_mat = np.transpose(model_rdms_mat)
    task_data_rel = {}
    for task in TASKS:
        # read in 4D image
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

        # out csv filename
        sub_csv_fname = csv_fname % (sub, sub,task, corr_label, parc_label)
        if (not overwrite) and os.path.isfile(sub_csv_fname):
            read_bool = True
            # read in csv if already exists
            sub_csv = pd.read_csv(sub_csv_fname)
            # remove row column
            sub_csv = sub_csv.iloc[:,1:]
            # save all completed ROIs except for last one (since may not have been finished)
            completed_rois = np.unique(sub_csv['roi'])[:-1]
            sub_csv = sub_csv[sub_csv['roi'].isin(completed_rois)]
            sub_csv.to_csv(sub_csv_fname)
            out_csv_array = sub_csv.values.tolist()
        else:
            read_bool = False
            out_csv_array = []
            completed_rois = []
        wtr = csv.writer(open(sub_csv_fname, 'a'), delimiter=',', lineterminator='\n')

        if isSl(procedure):
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
                        mask_roi = get_sphere_mask(index=index, radius=rad, brain_mask=sub_mask, sphere_template=sphere_template)
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
                            csv_row = [sub, task, x, y, z, parc_roi, k, val]
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

        elif isParc(procedure):
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
                        csv_row = [sub, task, parc_roi, k, val]
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
            fname = out_fname % (sub, sub, task, corr_label, parc_label, val_label, k)
            save_nii( out_data_dict[k], ref_img, fname )

        # create contrasts for relevant vs non-relavant
        if task == 'number':
            task_data_rel[task] = out_data_dict[deg_label] - out_data_dict[dist_label]
            fname = out_fname % (sub, sub, task, corr_label, parc_label, val_label, "rel")
            save_nii( task_data_rel[task], ref_img, fname )
        elif task == 'friend':
            task_data_rel[task] = out_data_dict[dist_label] - out_data_dict[deg_label]
            fname = out_fname % (sub, sub, task, corr_label, parc_label, val_label, "rel")
            save_nii( task_data_rel[task], ref_img, fname )

print(str(datetime.now()) + ": End rsa_regs.py")
