source funcs
setup_modules ants

# get subjects
in_dir=${NISTATS_DIR}
convert_sub_args -f subid -c "$@" ${in_dir}
subs=( "${SUBS[@]}" )

# MNI parcellation
mni_fname=$(basename $MNI_PARCELLATION)
mni_fname=${mni_fname%.nii*}
echo "MNI image: $MNI_PARCELLATION"

for sub in ${subs[@]}; do
  echo "Starting subject ${sub}"
  sub_dir="${FMRIPREP_DIR}/${sub}"

  # transform parcellation for each subject
  sub_ref="${sub_dir}/func/*run-01_*space-${SPACE}_*boldref.nii.gz"
  sub_trans="${sub_dir}/anat/${sub}_from-*_to-${SPACE}_mode-image_xfm.h5"
  out_dir="${RSA_DIR}/${sub}"
  mkdir -p $out_dir
  out_fname="${out_dir}/${mni_fname}_${sub}_space-${SPACE}_transformed.nii"
  antsApplyTransforms -d 3 -i ${MNI_PARCELLATION} -o ${out_fname} -r ${sub_ref} -t ${sub_trans} -n NearestNeighbor
done

echo "Done."
