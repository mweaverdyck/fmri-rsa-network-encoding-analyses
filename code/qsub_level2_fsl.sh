#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/level2fsl.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=11:59:00,h_data=4G
# Notify when
#$ -m ae
#
set -e
source funcs
setup_modules $fsl_v

label='LEVEL2_RSA'
in_dir=`get_thresh_dir ${RSA_DIR}`
out_dir=`get_thresh_dir ${SECOND_LEVEL_DIR}/fsl_sl${SL_RADIUS}`

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid "all"
subs=( ${SUBS[@]} )
log_args="$LOG_ARGS"

write_log $log_args "Analyzing subjects: all"

preds_all=( 'dist' 'deg' ) # 'sn' 'img' )
corrs_all=( 'spear' 'reg' ) # 'spear' 'reg'

preds=()
corrs=()
for arg in "$@"; do
  if [[ " ${preds_all[@]} " =~ " ${arg} " ]]; then
    preds=( ${preds[@]} $arg )
  elif [[ " ${corrs_all[@]} " =~ " ${arg} " ]]; then
    corrs=( ${corrs[@]} $arg )
  fi
done
if [[ ${#preds} -eq 0 ]]; then
  preds=( ${preds_all[@]} )
fi
if [[ ${#corrs} -eq 0 ]]; then
  corrs=( ${corrs_all[@]} )
fi

# file extension
suffix=".nii"

if [[ $SPACE == $T1_SPACE ]]; then
  . transform_T1w-2-mni.sh $ALL
  sdata_dir="/T1w-2-MNI"
else
  sdata_dir=""
fi

# one sample t-tests
for pred in "${preds[@]}"; do
  write_log $log_args "Running tests on predictor: $pred"
  sub_imgs=()
  first_sub=1
  for corr in "${corrs[@]}"; do
    if [[ $corr == "reg" ]]; then
      val="beta"
    else
      val="r"
    fi
    for s in ${subs[@]}; do
      # prefix of filename
      prefix="${s}_"
      # find full file name
      sub_dir=${in_dir}/${s}/${corr}${sdata_dir}
      if [ -d "${sub_dir}" ]; then
        fname=( $(ls -f ${sub_dir}/${prefix}task-avg_space-${SPACE}_*parc-sl${SL_RADIUS}*_val-${val}*pred-${pred}*${suffix}) )
        if [[ ${#fname[@]} -ne 1 ]]; then
          write_log $log_args "WARNING: Subject does not have exactly one image: ${fname[@]} Skipping..."
        else
          write_log $log_args "Reading in subject $s, file: $fname"
          # add to array
          sub_imgs+=( "${fname}" )
          # if first subject, save output file name
          if [[ $first_sub -eq 1 ]]; then
            # get filename only (no path)
            fname=$(basename ${fname})
            # remove subject prefix
            fname=${fname#"${prefix}"}
            # remove file extension suffix
            fname=${fname%"${suffix}"}
            fname_root=$fname
            # add 4D and suffix back on and
            # save full output filename
            out_fname=${out_dir}/${fname_root}_4D
            first_sub=0
          fi
        fi
      fi
    done

    # merge all subject images into one 4D file
    img_4D=${out_fname}${suffix}
    if [[ -f "$img_4D" ]]; then
      write_log $log_args "Found $img_4D, using saved image"
    else
      fslmerge -t ${img_4D} ${sub_imgs[@]}
      write_log $log_args "File ${out_fname}${suffix} saved."
    fi

    # smooth
    img_4D_smooth=${out_fname}_smooth-${SMOOTH_FWMH}${suffix}
    if [[ -f "$img_4D_smooth" ]]; then
      write_log $log_args "Found $img_4D_smooth, using saved image"
    else
      fslmaths ${img_4D} -s ${SMOOTH_SIGMA} ${img_4D_smooth}
    fi

    # run permuatation t-test
    #vox_thresh=1.645 #2.33 # z-score to get p<0.01
    out_t_test=${out_dir}/${fname_root}_mask-${MASK_NAME}_smooth-${SMOOTH_FWMH}_vsmooth-${VAR_SMOOTH}_nperm-${N_PERMS} #_vthresh-${vox_thresh}
    # http://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/Randomise(2f)UserGuide.html
    randomise -i $img_4D_smooth -m ${MASK} -o $out_t_test -1 -v ${VAR_SMOOTH} -T -x -n ${N_PERMS} --uncorrp | tee -a $log_file # -C ${vox_thresh}

    write_log $log_args "Files $out_t_test* saved."
  done
done
log_end $log_args
