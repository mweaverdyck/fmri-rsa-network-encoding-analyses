
source funcs
setup_modules ${fsl_v}

ref_img="${MNI_DIR}/schaefer/tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-200Parcels7Networks_dseg.nii.gz"
lobe_mask="$PROJECT_DIR/anats/lobe_mask_2mm.nii.gz"
hem_mask=${PROJECT_DIR}/anats/hemisphere_mask_2mm
# Transform to new MNI space
hem_mask_mni=${hem_mask}_flirt-sch200_round-T
lobe_mask_mni=${PROJECT_DIR}/anats/lobe_mask_2mm_flirt-sch200_round-T
flirt -in ${hem_mask} -ref ${ref_img} -out ${hem_mask_mni} -interp nearestneighbour
echo "${hem_mask_mni} saved."
flirt -in ${lobe_mask} -ref ${ref_img} -out ${lobe_mask_mni} -interp nearestneighbour
echo "${lobe_mask_mni} saved."

# binarize lobe mask and apply to hemisphere mask
lobe_mask_mni_bin=${lobe_mask_mni}_bin
fslmaths ${lobe_mask_mni} -bin ${lobe_mask_mni_bin}
echo "${lobe_mask_mni_bin} saved."
# apply above binary mask to hemisphere image to remove extra
hem_mask_mni_masked=${hem_mask_mni}_masked
fslmaths ${hem_mask_mni} -mul ${lobe_mask_mni_bin} ${hem_mask_mni_masked}
echo "${hem_mask_mni_masked} saved."
# multiply hem image by 10, add to lobe mask, subtract 10 and threshold it to remove all lower numbers
lobe_mask_mni_hem=${lobe_mask_mni}_hem
fslmaths ${hem_mask_mni_masked} -mul 10 -add ${lobe_mask_mni} -sub 10 -thr 0 ${lobe_mask_mni_hem}
echo "${lobe_mask_mni_hem} saved."

echo "Removing intermediate files:"
ls ${lobe_mask_mni_bin}.nii*
rm ${lobe_mask_mni_bin}.nii*
ls ${hem_mask_mni_masked}.nii*
rm ${hem_mask_mni_masked}.nii*

echo "DONE."
