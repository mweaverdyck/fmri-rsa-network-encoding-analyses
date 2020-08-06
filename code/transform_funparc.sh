
source funcs
setup_modules ${fsl_v}

fun_parc_dir=${MNI_DIR}/functional_parcellation
ref_img="${MNI_DIR}/schaefer/tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-200Parcels7Networks_dseg.nii.gz"
atlas_img="${MNI_DIR}/MNI_lobes_2mm_hemisphere_split.nii.gz"

for k in 50 100 200; do
  echo "-----------------------------------------------------------------------"
  in_img="${fun_parc_dir}/whole_brain_cluster_labels_k-${k}_order-orig"
  out_img=${in_img}_flirt-sch200_round-T
  echo "Transforming ${in_img} to schaefer space -->"
  flirt -in ${in_img} -ref ${ref_img} -out ${out_img} -interp nearestneighbour
  echo "--> ${out_img} saved."

  # binarize image to mask atlas
  in_img=${out_img}
  img_bin=${out_img}_bin
  echo "binarizing image"
  fslmaths ${in_img} -bin ${img_bin}
  echo "--> ${img_bin} saved."
  echo "masking atlas"
  atlas_masked=${atlas_img}_masked
  fslmaths ${atlas_img} -mul ${img_bin} ${atlas_masked}
  echo "--> ${atlas_masked} saved."

  # find intersection of atlas and image
  k1=$(($k + 1))
  out_img=${in_img}_split-atlas
  fslmaths ${atlas_masked} -mul $k1 -add $in_img -sub $k1 -thr 0 $out_img
  echo "--> ${out_img} saved."

  # in_img=${out_img}
  # img100=${out_img}_plus-1000
  # echo "adding 100 to image"
  # fslmaths ${in_img} -add 1000 -thr 1001 ${img100}
  # echo "--> ${img100} saved."
  #
  # k1=$(($k + 1))
  # atlas_img100=${atlas_img}_plus-$k
  # fslmaths ${atlas_img} -add $k -thr $k1 ${atlas_img100}
  #
  # k1001=$(($k1 + 1000))
  # in_img=${img100}
  # out_img=${out_img}_split-atlas
  # echo "adding atlas to image"
  # fslmaths ${in_img} -add ${atlas_img100} -thr $k1001 ${out_img}
  # echo "--> ${out_img} saved."

  echo "Removing intermediate files: "
  ls ${img_bin}.nii*
  rm ${img_bin}.nii*
  ls ${atlas_masked}.nii*
  rm ${atlas_masked}.nii*
  echo "-----------------------------------------------------------------------"
done

#echo "Removing ${atlas_img100}.nii*"
#rm ${atlas_img100}.nii*

echo "DONE."
