#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/level2fsl.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=11:59:00,h_data=8G
# Notify when
#$ -m ae
#
set -e
source level2_fsl_funcs "$@"

label='LEVEL2_RSA_PAIREDT'

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid "all"
subs=( ${SUBS[@]} )
log_args="$LOG_ARGS"

write_log $log_args "Analyzing subjects: all"

design="design_pairedt${tlab}"
contrast="con_pairedt${tlab}"
group="group_pairedt${tlab}"
text2vest ${design} ${contrast} ${group}

# one sample t-tests
for p in "${preds[@]}"; do
  if [[ $p == 'dist' ]]; then
    task='friend'
    preds_rel=( 'dist' 'deg' )
  elif [[ $p == 'deg' ]]; then
    task='number'
    preds_rel=( 'deg' 'dist' )
  fi
  write_log $log_args "Running tests on task: $task"

  for corr in "${corrs[@]}"; do
    for dist in "${dists[@]}"; do
      sub_imgs=()
      first_sub=1
      for pred in ${preds_rel[@]}; do
        for s in ${subs[@]}; do
          # prefix of filename
          prefix="${s}_"
          # find full file name
          sub_dir=${in_dir}/${s}/${corr}/${dist}${sdata_dir}
          if [ -d "${sub_dir}" ]; then
            if [[ $SPACE == $T1_SPACE ]]; then
              # file extension
              transform_str="transformed-${MNI_SPACE}"
              suffix="pred-${pred}_cat-sn_${transform_str}.nii"
            else
              suffix="pred-${pred}_cat-sn.nii.gz"
            fi
            fname=( $(ls -f ${sub_dir}/${prefix}task-${task}_space-${SPACE}_*parc-sl${SL_RADIUS}*${suffix}) )
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
                fname_root="${fname}${transform_str}"
                # add 4D and suffix back on and
                # save full output filename
                out_fname=${out_dir}/${dist}/${fname_root}_4D
                first_sub=0
              fi
            fi
          fi
        done
      done

      # merge all subject images into one 4D file
      img_4D="${out_fname}.nii"
      if [[ -f "$img_4D" ]]; then
        write_log $log_args "Found $img_4D, using saved image."
      else
        fslmerge -t ${img_4D} ${sub_imgs[@]}
        write_log $log_args "File ${img_4D} saved."
      fi

      # smooth
      img_4D_smooth="${out_fname}_smooth-${SMOOTH_FWMH}.nii"
      if [[ -f "$img_4D_smooth" ]]; then
        write_log ${log_args} "Found ${img_4D_smooth}, using saved image."
      else
        fslmaths ${img_4D} -s ${SMOOTH_SIGMA} ${img_4D_smooth}
        write_log $log_args "File ${img_4D_smooth} saved."
      fi

      # run permuatation t-test
      #vox_thresh=1.645 #2.33 # z-score to get p<0.01
      out_t_test=${out_dir}/${dist}/pred-diff_${fname_root}_mask-${MASK_NAME}_smooth-${SMOOTH_FWMH}_test-pairedt_vsmooth-${VAR_SMOOTH}_nperm-${N_PERMS}_clust-${CLUSTER_THRESH} #_vthresh-${vox_thresh}
      # http://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/Randomise(2f)UserGuide.html
      randomise -i ${img_4D_smooth} -m ${MASK} -o ${out_t_test} -d ${design}.mat -t ${contrast}.con -e ${group}.grp -v ${VAR_SMOOTH} -T -x -n ${N_PERMS} -C ${CLUSTER_THRESH} --uncorrp -D | tee -a $log_file # -C ${vox_thresh}

      write_log $log_args "Files ${out_t_test}* saved."

    done
  done
done
log_end $log_args
