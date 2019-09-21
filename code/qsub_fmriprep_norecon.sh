#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.fmriprep.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=23:00:00,h_data=23G
# Notify when
#$ -m ae
#
# runs fmriprep on inputted subjects

# break if error raised
set -e

# load modules and functions
source funcs
setup_modules $fsl_v $python_v freesurfer/6.0.0 ants/ants-2.2.0 afni ica-aroma itksnap/3.6.0-RC1.QT4 $fmriprep_v
# module load fsl/5.0.9
# module load python/3.7.2
# module load freesurfer/6.0.0
# module load ants/ants-2.2.0
# module load afni
# module load ica-aroma
# module load itksnap/3.6.0-RC1.QT4
# module load fmriprep/1.3.2
#
# export PATH=/u/project/CCN/apps/c3d/c3d-1.0.0-Linux-x86_64/bin/:$PATH
# export NO_FSL_JOBS=true

# get subjects
convert_sub_args -f numid -c "$@"

# log variables
label='FMRIPREP'
sub_str=''; for s in "${SUBS[@]}"; do sub_str=$(echo "${sub_str}_${s}"); done
log_dir="${BIDS_DIR}/fmriprep/logs"
mkdir -p ${log_dir}
log_file="${log_dir}/LOG${sub_str}.log"
log_args="$label $log_file"

# log start
log_begin $log_args
# log subjects
write_log $log_args "Analyzing subjects [${SUBS[@]}]"

# run fmriprep
fmriprep ${BIDS_DIR} ${BIDS_DIR} --work-dir "${FMRIPREP_DIR}/work" --ignore slicetiming --fs-license-file $FREESURFER_HOME/.license participant --participant-label "${SUBS[@]}" --output-space T1w template MNI152NLin6Asym --fs-no-reconall --nthreads 4 --omp-nthreads 4 --mem-mb 32000 | tee -a ${log_file}

# log end
log_end $log_args
