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

label='LEVEL2_FACE_FSL'
in_dir=${GLM_DIR}
out_dir=${SECOND_LEVEL_DIR}/face

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid "all"
subs=( ${SUBS[@]} )
log_args="$LOG_ARGS"

write_log $log_args "Analyzing subjects: all"

# one sample t-tests
write_log $log_args "Running tests on predictor: $pred"
sub_imgs=()
first_sub=1
# file extension
suffix=".nii"
for s in ${subs[@]}; do
  # prefix of filename
  prefix="${s}_"
  # find full file name
  sub_dir=${in_dir}/${s}/
  if [ -d "${sub_dir}" ]; then
    fname=( $(ls -f ${sub_dir}/${prefix}task-??????_space-${MNI_SPACE}_stat-t_node-all${suffix}) )
    # if [[ ${#fname[@]} -ne 1 ]]; then
    #   write_log $log_args "WARNING: Subject does not have exactly one image: ${fname[@]} Skipping..."
    # else
    write_log $log_args "Reading in subject $s, file: ${fname[@]}"
    # add to array
    sub_imgs+=( "${fname[@]}" )
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
var_smooth=10
n_perms=5000
#vox_thresh=1.645 #2.33 # z-score to get p<0.01
out_t_test=${out_dir}/${fname_root}_smooth-${SMOOTH_FWMH}_vsmooth-${var_smooth}_nperm-${n_perms} #_vthresh-${vox_thresh}
# http://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/Randomise(2f)UserGuide.html
randomise -i $img_4D_smooth -m ${MNI_MASK_DIL} -o $out_t_test -1 -v ${var_smooth} -T -x -n ${n_perms} --uncorrp | tee -a $log_file # -C ${vox_thresh}

write_log $log_args "Files $out_t_test* saved."
log_end $log_args
