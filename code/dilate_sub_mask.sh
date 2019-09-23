source funcs
setup_modules ${fsl_v}

dil=5

label='DILATE_MASK'
in_dir=${FMRIPREP_DIR}
out_dir=${RSA_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

for sub in ${subs[@]}; do
  echo "Starting subject ${sub}"
  sub_dir="${in_dir}/${sub}"
  out_sub_dir="${out_dir}/${sub}"
  mkdir -p ${out_sub_dir}

  # copy and dilate T1w mask
  sub_mask="${sub_dir}/anat/*space-${SPACE}_*desc-brain_mask.nii*"
  if [[ ! -f ${sub_mask} ]]; then
    sub_mask="${sub_dir}/anat/${sub}_desc-brain_mask.nii*"
  fi
  echo "Subject mask: ${sub_mask}"
  sub_mask_fname=$(basename ${sub_mask})
  sub_mask_fname=${sub_mask_fname%.nii*}
  out_mask="${out_sub_dir}/${sub_mask_fname}_dil-${dil}.nii"
  if [[ -f ${out_mask} ]] || [[ -f ${out_mask}.gz ]]; then
    echo "${out_mask} already exists. Skipping dilation..."
  else
    echo "Dilating ${sub_mask_fname} with ${dil} kernel --> ${out_mask}"
    fslmaths ${sub_mask} -kernel gauss ${dil} -dilF ${out_mask}
    ref_pattern="${out_sub_dir}/*node-4D*"
    ref_imgs=( $ref_pattern )
    ref_img=${ref_imgs[0]}
    echo "Resampling mask to reference image: ${ref_img}"
    flirt -in ${out_mask} -ref ${ref_img} -out ${out_mask}
  fi
done

echo "Done."
