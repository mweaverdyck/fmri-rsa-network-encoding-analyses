#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.transformmask.$JOB_ID.log
#$ -j y
#$ -pe shared 1
#$ -l h_rt=0:59:00,h_data=4G
# Notify when
#$ -m ae
#

source funcs
setup_modules ${fsl_v} ants

label='DILATE_MASK'

clust=1 # calculate and save cluster map
social_map=${SOCIAL_MAP}


# find clusters
if [[ $clust -eq 1 ]]; then
  bash get_social_mask_clusters.sh $DILATE
  social_map=${SOCIAL_MAP}_clust_thresh-163
fi

# first transform social map to a subject's T1w space using fmriprep's transformation matrix
# note: it does not matter which subject is used
sub="sub-185"
social_mask=${social_map}_T1
sub_ref="${FMRIPREP_DIR}/${sub}/func/*run-01_*space-T1w_*boldref.nii.gz"
sub_trans="${FMRIPREP_DIR}/${sub}/anat/${sub}_from-*_to-T1w_mode-image_xfm.h5"
antsApplyTransforms -d 3 -i ${social_map}.nii.gz -o ${social_mask}.nii.gz -r ${sub_ref} -t ${sub_trans} -n NearestNeighbor
echo "${social_mask} saved."

# next, transform into the correct MNI space, again using fmriprep's transformation matrix
social_mask2=${social_mask}_MNI
sub_trans="${FMRIPREP_DIR}/sub-185/anat/${sub}_from-T1w_to-${MNI_SPACE}_mode-image_xfm.h5"
antsApplyTransforms -d 3 -i ${social_mask}.nii.gz -o ${social_mask2}.nii.gz -r ${MNI_MASK} -t ${sub_trans} -n NearestNeighbor
echo "${social_mask2} saved."

if [[ $clust -eq 1 ]]; then
  # split by lobes
  # find intersection of atlas and image

  # binarize image to mask atlas
  in_img=${social_mask2}
  img_bin=${social_mask2}_bin
  echo "binarizing image"
  fslmaths ${in_img} -bin ${img_bin}
  echo "--> ${img_bin} saved (will be deleted)."
  echo "masking atlas"
  atlas_img="${MNI_DIR}/MNI_lobes_2mm_hemisphere_split"
  atlas_masked=${atlas_img}_masked
  fslmaths ${atlas_img} -mul ${img_bin} ${atlas_masked}
  echo "--> ${atlas_masked} saved (will be deleted)."

  # find intersection of atlas and image
  out_img=${in_img}_split-atlas
  fslmaths ${atlas_masked} -mul 200 -add $in_img -sub 200 -thr 0 $out_img
  echo "--> ${out_img} saved."

  # remove intermediate files
  echo "removing intermediate files"
  ls ${img_bin}*
  rm ${img_bin}*
  ls ${atlas_masked}*
  rm ${atlas_masked}*

  social_mask3=${out_img}
else
  # threshold map to create a binary mask
  # note: original social map is already thresholded at p < .01, so just need to threshold above 0
  social_mask3=${social_mask2}_bin
  fslmaths ${social_mask2} -thr .01 -bin ${social_mask3}
  echo "${social_mask3} saved."
fi

# finally, dilate mask
social_mask4=${social_mask3}_dil-$DILATE
fslmaths ${social_mask3} -kernel gauss $DILATE -dilF ${social_mask4}
echo "${social_mask4} saved."

# # smooth mask
# social_mask5=${social_mask4}_smooth-${SMOOTH_FWMH}
# fslmaths ${social_mask4} -s ${SMOOTH_SIGMA} ${social_mask5}
# echo "${social_mask5} saved."

# print fsleyes command to view result compared to original map
echo "View results overlaid on the standard template with:"
mni_std="$PROJECT_DIR/anats/2mm_T1.nii.gz"
echo "fsleyes $mni_std -cm greyscale ${social_mask4} -cm blue ${SOCIAL_MAP}.nii.gz -cm red-yellow -dr 0 11 &"

echo "DONE."
