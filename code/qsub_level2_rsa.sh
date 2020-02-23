#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.level2rsa.$JOB_ID.log
#$ -j y
#$ -pe shared 1
#$ -l h_rt=23:59:00,h_data=8G
# Notify when
#$ -m ae
#
# break if error raised
set -e

# load modules and functions
source funcs
setup_modules ${python_v}

label='LEVEL2_RSA'
in_dir=${RSA_DIR}
out_dir=${SECOND_LEVEL_DIR}/${SL}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid "all"
log_args="$LOG_ARGS"

write_log $log_args "Analyzing subjects: ${SUBS[@]}"
mkdir -p "${out_dir}"

for procedure in "${PROCEDURES[@]}"; do
  write_log $log_args "Running procedure: ${procedure}"
  if [[ $procedure -eq $SL ]]; then
    . transform_T1w-2-mni.sh $ALL
  fi
  python3 level2_${procedure}.py $@
done

log_end $log_args
