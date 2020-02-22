#!/usr/bin/env bash

source funcs

corr=$1; shift

if [[ $SPACE == 'T1w' ]]; then
  setup_modules ants

  label='TRANSFORM_SPACES'
  in_dir=${RSA_DIR}
  out_dir=${RSA_DIR}

  begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
  log_args_tmp="$LOG_ARGS"
  subs=( "${SUBS[@]}" )

  # MNI parcellation
  write_log $log_args_tmp "MNI image: $MNI_PARCELLATION"
  ref_img="$MNI_PARCELLATION"

  for sub in ${subs[@]}; do
    write_log $log_args_tmp "Starting subject ${sub}"

    sub_trans="${FMRIPREP_DIR}/${sub}/anat/${sub}_from-T1w_to-${MNI_SPACE}_mode-image_xfm.h5"

    sub_dir="${in_dir}/${sub}/${corr}"

    out_sub_dir="${sub_dir}/T1w-2-MNI"
    mkdir -p $out_sub_dir

    write_log $log_args_tmp "Transforming images to MNI space using : ${sub_trans}"
    for f in $(ls -f ${sub_dir}/*space-${SPACE}*parc-${SL}*.nii.gz); do
      echo $f
      in_f=$(basename $f)
      in_f=${in_f%.nii*}
      out_fname="${out_sub_dir}/${in_f}_transformed-${SPACES[1]}.nii"
      echo "$in_f"
      echo $out_fname

      if [[ -f "$out_fname" ]]; then
        write_log $log_args_tmp "$out_fname already exists. Skipping transformation..."
      else
        write_log $log_args_tmp "${f} --> ${out_fname}"
        antsApplyTransforms -d 3 -i ${f} -o ${out_fname} -r ${ref_img} -t ${sub_trans} -n NearestNeighbor
      fi
    done
  done
fi

write_log $log_args_tmp "Done."
