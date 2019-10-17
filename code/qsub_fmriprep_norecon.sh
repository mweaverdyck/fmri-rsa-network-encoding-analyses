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

label='FMRIPREP'
in_dir=${RECON_DIR}
out_dir=${FMRIPREP_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f numid $@
log_args="${LOG_ARGS}"
subs=( "${SUBS[@]}" )

# run fmriprep
##for version 1.3.2
##fmriprep ${BIDS_DIR} ${BIDS_DIR} --work-dir "${FMRIPREP_DIR}/work" --ignore slicetiming --fs-license-file $FREESURFER_HOME/.license participant --participant-label "${SUBS[@]}" --output-space T1w template MNI152NLin6Asym --fs-no-reconall --nthreads 4 --omp-nthreads 4 --mem-mb 32000 | tee -a ${log_file}
fmriprep ${BIDS_DIR} "${FMRIPREP_DIR}/.." --work-dir "${FMRIPREP_DIR}/work" --ignore slicetiming --fs-license-file $FREESURFER_HOME/.license participant --participant-label "${SUBS[@]}" --output-space T1w template --fs-no-reconall --nthreads 4 --omp-nthreads 4 --mem-mb 32000 | tee -a ${log_file}

# log end
log_end $log_args
