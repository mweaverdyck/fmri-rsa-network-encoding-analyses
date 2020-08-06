
source funcs
setup_modules $fsl_v

p="corrp" #corrp vox_p
thr="0.95" #0.95 0.99

if [[ -z $1 ]]; then
  cd $SECOND_LEVEL_DIR/parc-sch200
else
  cd $SECOND_LEVEL_DIR/parc-$1
fi

fsl_cmd="fsleyes $PROJECT_DIR/anats/2mm_T1.nii.gz -cm greyscale"

for f in $(ls -f */*p_correction-fdr_dir-rev_roi_stats.nii*); do
  minmax=`fslstats $f -R`
  min=`echo "$minmax" | cut -d ' ' -f 1`
  max=`echo "$minmax" | cut -d ' ' -f 2`
  if [[ "$max" > "${thr}" ]]; then
    echo $f $max
    fsl_cmd="$fsl_cmd $f -dr "${thr}" 1.00 -cm red-yellow"
  fi
done

echo "$fsl_cmd &"
`echo $fsl_cmd` &

cd $CODE_DIR
