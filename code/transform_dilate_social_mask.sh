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

# first transform social map to a subject's T1w space using fmriprep's transformation matrix
# note: it does not matter which subject is used
sub="sub-185"
social_mask=${SOCIAL_MAP}_T1
sub_ref="${FMRIPREP_DIR}/${sub}/func/*run-01_*space-T1w_*boldref.nii.gz"
sub_trans="${FMRIPREP_DIR}/${sub}/anat/${sub}_from-*_to-T1w_mode-image_xfm.h5"
antsApplyTransforms -d 3 -i ${SOCIAL_MAP}.nii.gz -o ${social_mask}.nii.gz -r ${sub_ref} -t ${sub_trans} -n NearestNeighbor
echo "${social_mask} saved."

# next, transform into the correct MNI space, again using fmriprep's transformation matrix
social_mask2=${social_mask}_MNI
sub_trans="${FMRIPREP_DIR}/sub-185/anat/${sub}_from-T1w_to-${MNI_SPACE}_mode-image_xfm.h5"
antsApplyTransforms -d 3 -i ${social_mask}.nii.gz -o ${social_mask2}.nii.gz -r ${MNI_PARCELLATION} -t ${sub_trans} -n NearestNeighbor
echo "${social_mask2} saved."

# threshold map to create a binary mask
# note: original social map is already thresholded at p < .01, so just need to threshold above 0
social_mask3=${social_mask2}_bin
fslmaths ${social_mask2} -thr .01 -bin ${social_mask3}
echo "${social_mask3} saved."

# finally, dilate mask
fslmaths ${social_mask3} -kernel gauss $DILATE -dilF ${SOCIAL_MASK_DIL}
echo "${SOCIAL_MASK_DIL} saved."

# print fsleyes command to view result compared to original map
echo "View results overlaid on the standard template with:"
mni_std="$PROJECT_DIR/anats/2mm_T1.nii.gz"
echo "fsleyes $mni_std -cm greyscale ${SOCIAL_MASK_DIL} -cm red-yellow ${SOCIAL_MAP}.nii.gz -cm cool &"

echo "DONE."
