source funcs
setup_modules $fsl_v

thresh=163
#
# if [[ $# -eq 0 ]]; then
#   dils=( 2 3 4 5 )
# else
#   dils=( "$@" )
# fi
#
# for d in "${dils[@]}"; do
  fname="${SOCIAL_MAP}" #_T1_MNI"
  fname_clust="${SOCIAL_MAP}_clust" #_T1_MNI"
  tsv_fname="${fname_clust}_info.tsv"
  cluster --in=${fname} --thresh=1.96 --oindex=${fname_clust} | tee ${tsv_fname}
  export fname_thresh="${SOCIAL_MAP}_clust_thresh-${thresh}" #_T1_MNI"
  fslmaths ${fname_clust} -thr $thresh ${fname_thresh}
  #
  # fname_dil="${fname_clust}_dil-$d"
  # fslmaths ${fname} -kernel gauss $d -dilF ${fname_dil}
# done
