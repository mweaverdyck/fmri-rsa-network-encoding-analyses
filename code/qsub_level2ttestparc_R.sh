#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.rsacorr.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=0:59:00,h_data=8G
# Notify when
#$ -m ae
#
# break if error raised
set -e

# load modules and functions
source funcs
setup_modules R

label='LEVEL2_TTEST'
in_dir=${RSA_DIR}
out_dir=${SECOND_LEVEL_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

write_log $log_args "Analyzing subjects: ${sub[@]}"
# subject's contrasts directory
subdir="${in_dir}/${sub}"
mkdir -p "${out_dir}"

Rscript --vanilla --verbose level2_ttest_parc.R

log_end $log_args
