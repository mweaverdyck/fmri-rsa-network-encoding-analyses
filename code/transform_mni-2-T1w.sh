#!/usr/bin/env bash

source funcs
setup_modules ants

label='TRANSFORM_SPACES'
in_dir=${GLM_DIR}
out_dir=${RSA_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

# MNI parcellation
mni_fname=$(basename $MNI_PARCELLATION)
mni_fname=${mni_fname%.nii*}
write_log $log_args "MNI image: $MNI_PARCELLATION"

for sub in ${subs[@]}; do
  write_log $log_args "Starting subject ${sub}"
  sub_dir="${FMRIPREP_DIR}/${sub}"

  # transform parcellation for each subject
  sub_ref="${sub_dir}/func/*run-01_*space-${SPACE}_*boldref.nii.gz"
  sub_trans="${sub_dir}/anat/${sub}_from-*_to-${SPACE}_mode-image_xfm.h5"
  out_sub_dir="${out_dir}/${sub}"
  mkdir -p $out_sub_dir
  write_log $log_args "Transforming MNI parcellation to subject's native space using : ${sub_trans}"
  out_fname="${out_sub_dir}/${mni_fname}_${sub}_space-${SPACE}_transformed.nii"
  write_log $log_args "${MNI_PARCELLATION} --> ${out_fname}"
  antsApplyTransforms -d 3 -i ${MNI_PARCELLATION} -o ${out_fname} -r ${sub_ref} -t ${sub_trans} -n NearestNeighbor
done

echo "Done."
