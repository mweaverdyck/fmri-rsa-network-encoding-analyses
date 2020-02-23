#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.level2parc.$JOB_ID.log
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

SPACE=${MNI_SPACE}

label='LEVEL2_FACE'
in_dir=${GLM_DIR}
out_dir=${SECOND_LEVEL_DIR}/face

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $ALL
log_args="$LOG_ARGS"

write_log $log_args "Analyzing subjects: ${ALL}"
mkdir -p "${out_dir}"

python3 level2_face.py

log_end $log_args
