source funcs
setup_modules ${fsl_v}

dil=5

# get subjects
in_dir=${NISTATS_DIR}
convert_sub_args -f subid -c "$@" ${in_dir}
subs=( "${SUBS[@]}" )

for sub in ${subs[@]}; do
  echo "Starting subject ${sub}"
  sub_dir="${FMRIPREP_DIR}/${sub}"

  # copy and dilate T1w mask
  sub_mask="${sub_dir}/anat/*space-${SPACE}_*desc-brain_mask.nii*"
  if [[ ! -f ${sub_mask} ]]; then
    sub_mask="${sub_dir}/anat/${sub}_desc-brain_mask.nii*"
  fi
  echo "Subject mask: ${sub_mask}"
  sub_mask_fname=$(basename ${sub_mask})
  sub_mask_fname=${sub_mask_fname%.nii*}
  out_mask="${out_dir}/${sub_mask_fname}_dil-${dil}.nii"
  echo "Dilating ${sub_mask_fname} with ${dil} kernel --> ${out_mask}"
  fslmaths ${sub_mask} -kernel gauss ${dil} -dilF ${out_mask}
done

echo "Done."
