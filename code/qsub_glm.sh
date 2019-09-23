#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.glm.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=8:00:00,h_data=14G
# Notify when
#$ -m ae
#
# runs GLM in Nistats
# break if error raised
#set -e

# load modules and functions
source funcs
setup_modules $python_v

label='GLM'
in_dir=${DERIVATIVES_DIR}
out_dir=${GLM_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f numid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

for s in "${subs[@]}"; do
    deriv_dir="${in_dir}/derivatives_$s"
    deriv_dir_sub="${deriv_dir}/$s"
    if [[ -d ${deriv_dir_sub} ]]; then
        write_log $log_args "WARNING: Moving subjects directory from derivatives folder"
        mv ${deriv_dir_sub} "${in_dir}/"
        rm -r ${deriv_dir}
    fi
done

write_log $log_args "Counting and removing beginning TRs"
bash qsub_rm_nonsteady_trs.sh "${subs[@]}" | tee -a $log_file
write_log $log_args "Finished removing TRs"

for s in "${subs[@]}"; do
    write_log $log_args "Running Nistats GLM"
    python glm.py "${s}" | tee -a $log_file
    write_log $log_args "Finished GLM"

    deriv_dir="${in_dir}/derivatives_$s/$s"
    if [[ -d ${deriv_dir} ]]; then
        write_log $log_args "WARNING: Moving subjects directory from derivatives folder"
        mv ${deriv_dir} "${in_dir}/"
        rm -r ${deriv_dir}
    fi
done

log_end $log_args
