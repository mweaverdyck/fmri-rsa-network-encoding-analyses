#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.rsacorr.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=1:59:00,h_data=8G
# Notify when
#$ -m ae
#
# runs runs searchlight RSA regression based on Nistats GLM output
# break if error raised
set -e
owd=$(pwd)
echo $owd

STAT="t"
procedures=( "parc" ) #parc (parcellation) or sl (searchlight)

# load modules and functions
source funcs
setup_modules $fsl_v $python_v

in_dir=$GLM_DIR
out_dir=${RSA_DIR}
log_dir=${out_dir}/logs
mkdir -p ${out_dir}
mkdir -p ${log_dir}

# get subjects
subs=( "$@" )
convert_sub_args -f subid -c "${subs[@]}" ${in_dir}
subs=( "${SUBS[@]}" )
echo "Running subjects: ${subs[@]}"


# log variables
label='RSA_CORR'
sub_str=''; for sub in "${subs[@]}"; do sub_str=$(echo "${sub_str}_${sub}"); done
log_file="${log_dir}/LOG${sub_str}.log"
log_args="$label $log_file"
echo "Writing to logfile: ${log_file}"

# log start
log_begin $log_args
# log subjects
write_log $log_args "Analyzing subjects ${subs[@]}."

# concatenate t-images from feat folders for each subject, for each block
for sub in "${subs[@]}"; do
    # subject's contrasts directory
    subdir="${in_dir}/${sub}"
    mkdir -p "${out_dir}/$sub"

    write_log $log_args "Subject input directory: $subdir"
    write_log $log_args "Subject output directory: ${out_dir}/$sub"

    # create one concatenated image per task
    for task in "${TASKS[@]}"; do
        # file prefix (must match nii in subject's directory)
        pref="${sub}_task-${task}_space-${SPACE}_stat-${STAT}_node-"
        # output file name
        img_4d="${out_dir}/${sub}/${pref}4D"
        # input file names
        img_pref="${subdir}/${pref}"
        # get list of node images
        node_imgs=()
        for n in ${NODES[@]}; do
            node_imgs+=( "${img_pref}0$n" )
        done
        # create each node's contrast
        if [ ! -f "${img_4d}.nii" ]; then
            write_log $log_args "Concatenating node images for task $task..."
            # merge all nodes by the 4th dimension
            fslmerge -t ${img_4d} ${node_imgs[@]}
            #${img_pref}00 ${img_pref}01 ${img_pref}02 ${img_pref}03 ${img_pref}04 ${img_pref}05 ${img_pref}06 ${img_pref}07 ${img_pref}08 ${img_pref}09
            gunzip "${img_4d}.nii.gz"
            write_log $log_args "File ${img_4d} saved and unzipped."
        else
            write_log $log_args "Subject already agregated. Final file exists: ${img_4d}"
        fi
    done
done

for procedure in "${procedures[@]}"; do
  write_log $log_args "Running procedure: ${procedure}"
  if [[ $procedure == "parc" ]]; then
    all_parcels=( ${N_PARCELS} )
    bash transform_mni-2-T1w.sh ${subs[@]}
    outstat="r"
  elif [[ $procedure == "sl" ]]; then
    all_parcels=( "sl" )
    bash dilate_sub_mask.sh ${subs[@]}
    outstat="beta"
  else
  #if [[ $procedure != "parc" ]] && [[ $procedure != "sl" ]]; then
    echo "ERROR: unrecognized procedure name. Must be 'parc' or 'sl', not: $procedure"
    exit 0
  fi

  # Run RSA
  write_log $log_args "Running RSA correlation"
  python rsa_corr.py ${procedure} $STAT ${subs[@]} | tee -a $log_file
  write_log $log_args "Finished RSA correlation"

  # take difference between output images
  write_log $log_args "Calculate difference images"

  for p in "${all_parcels[@]}"; do
    #img_match="_stat-${STAT}_reg-beta_npred-02_corr-ortho_parc-${p}_pred-"
    img_match="_stat-${STAT}_corr-spear_parc-${p}_val-${outstat}_pred-"
    for sub in "${subs[@]}"; do
        pref="${RSA_DIR}/${sub}/${procedure}/${sub}"
        # distance: friend - number
        img_out=${pref}${img_match}dist_diff-fn
        fslmaths "${pref}_task-friend${img_match}dist" -sub "${pref}_task-number${img_match}dist" ${img_out}
        write_log $log_args "File $img_out saved."
        # degree: number - friend
        img_out=${pref}${img_match}deg_diff-nf
        fslmaths "${pref}_task-number${img_match}deg" -sub "${pref}_task-friend${img_match}deg" ${img_out}
        write_log $log_args "File $img_out saved."
    done
  done
done
log_end $log_args
