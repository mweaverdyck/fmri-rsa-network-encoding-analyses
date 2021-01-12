
source funcs
setup_modules ${fsl_v}

preds=( "deg" "dist" )
dist_meas="correlation"
thresh=95

for pred in "${preds[@]}"; do
  p_img_fname="task-avg_space-T1w_stat-t_corr-spear_parc-sl4_val-r_pred-${pred}_cat-sn_transformed-MNI152NLin2009cAsym_mask-MNI_smooth-6_vsmooth-10_nperm-10000_clust-1.699_vox_p_tstat1"
  p_img="${SECOND_LEVEL_DIR}/fsl_sl4/${dist_meas}/${p_img_fname}"

  echo Using following image to make mask for where $pred is encoded: $p_img

  p_mask="${CODE_DIR}/sig_voxs_task-avg_pred-${pred}_p-uncorr${thresh}_mask"
  # threshold uncorrected p-values at 95%
  fslmaths ${p_img} -thr ".${thresh}" -bin ${p_mask}

  echo ${p_mask} saved.
done
echo "DONE."
