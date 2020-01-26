#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.rsasl.$JOB_ID.log
#$ -j y
#$ -pe shared 2
#$ -l h_rt=1:59:00,h_data=8G
# Notify when
#$ -m ae
#
# runs runs searchlight RSA regression based on Nistats GLM output
# break if error raised
set -e

# load modules and functions
source funcs
setup_modules $python_v

label='RSA_SL'
in_dir=${GLM_DIR}
out_dir=${RSA_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

procedure=$SL
write_log $log_args "Running procedure: ${procedure}"

# Run RSA
write_log $log_args "Running RSA"
python3 rsa_sl.py 'avg' ${TASKS[@]} ${subs[@]} | tee -a $log_file
write_log $log_args "Finished RSA"

log_end $log_args
