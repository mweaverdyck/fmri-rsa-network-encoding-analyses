#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.rsastages.$JOB_ID.log
#$ -j y
#$ -pe shared 2
#$ -l h_rt=11:59:00,h_data=8G
# Notify when
#$ -m ae
#
# runs runs searchlight RSA regression based on Nistats GLM output
# break if error raised
set -e

# load modules and functions
source funcs
setup_modules $fsl_v $python_v

label='RSA_STAGES'
in_dir=${GLM_DIR}
out_dir=${RSA_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

for procedure in "${PROCEDURES[@]}"; do
  write_log $log_args "Running procedure: ${procedure}"
  if [[ $procedure == $PARC ]]; then
    all_parcels=( ${N_PARCELS} )
    bash transform_mni-2-T1w.sh ${subs[@]}
  fi

  # Run RSA
  write_log $log_args "Running RSA"
  python3 rsa_R2.py ${procedure} ${subs[@]} | tee -a $log_file
  write_log $log_args "Finished RSA"

done
log_end $log_args
