#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.nonsteady.$JOB_ID.log
#$ -j y
#$ -pe shared 2
#$ -l h_rt=11:00:00,h_data=14G
# Notify when
#$ -m ae
#
# break if error raised
#set -e

source funcs
setup_modules $fsl_v $python_v #R

label='NONSTEADYSTATES'
in_dir=${FMRIPREP_DIR}
out_dir="${DERIVS_DIR}"

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

suffix="_desc-preproc_bold"
# remove TRs from beginning of each preprocessed run
for sub in "${subs[@]}"; do
  sub_dir="${in_dir}/${sub}/func"
  sub_outdir="${out_dir}/${sub}/func"
  mkdir -p ${sub_outdir}
  # move all images from exclude directory
  mv ${sub_outdir}/exclude/* ${sub_outdir}/
  in_imgs="$(ls -f ${sub_dir}/${sub}_task-??????_run-??_space-${SPACE}${suffix}.nii.gz)"
  for in_img in ${in_imgs[@]}; do
    fname="$(basename ${in_img} ${suffix}.nii.gz)"
    out_img_bname="${fname}_rmtr-${N_TRS_DEL}${suffix}"
    out_img="${sub_outdir}/${out_img_bname}"
    write_log $log_args "Checking if ${out_img}.nii.gz exists"
    if [[ -f ${out_img}.nii.gz ]]; then
      write_log $log_args "$out_img already exists. Skipping..."
    else
      write_log $log_args "Removing ${N_TRS_DEL} TRs from beginning of volume ${in_img}..."
      fslroi ${in_img} ${out_img} ${N_TRS_DEL} ${N_TRS} | tee -a $log_file
      write_log $log_args "${out_img} saved."
    fi

  done

done

write_log $log_args "Updating events files..."
# calculate new onset times and count number of nonsteady outliers
#Rscript --vanilla --verbose qa_nonsteadystates.R ${subs[@]} | tee -a $log_file
python3 qa_nonsteadystates.py ${subs[@]} | tee -a $log_file
log_end $log_args
