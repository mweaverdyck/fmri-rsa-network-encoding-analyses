
source funcs
setup_modules $fsl_v

p="corrp" #corrp vox_p
thr="0.95" #0.95 0.99
dist="euclidean"

if [[ -z $1 ]]; then
  cd $SECOND_LEVEL_DIR/fsl_sl4/${dist}
else
  cd $SECOND_LEVEL_DIR/fsl_sl${1}/${dist}
fi

fsl_cmd="fsleyes $PROJECT_DIR/anats/2mm_T1.nii.gz -cm greyscale"

for f in $(ls -f *pred-d*transformed-MNI152NLin2009cAsym*smooth-6*vsmooth-10*nperm-10000_clust-1.699*_${p}_tstat*nii*); do
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
