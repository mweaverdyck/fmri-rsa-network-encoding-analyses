#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblogs/joblog.recon.$JOB_ID.log
#$ -j y
#$ -pe shared 4
#$ -l h_rt=16:00:00,h_data=3G
# Notify when
#$ -m ae

# break if error raised
set -e
# load modules and functions
source funcs
setup_modules python/2.7.13 freesurfer/6.0.0

label='RECON'
in_dir=${BIDS_DIR}
out_dir=${RECON_DIR}

begin_script -l ${label} -i ${in_dir} -o ${out_dir} -f subid $@
log_args="$LOG_ARGS"
subs=( "${SUBS[@]}" )

export FREESURFER_HOME="/u/project/CCN/apps/freesurfer/6.0.0"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

write_log ${log_args} "Running recon-all on subjects: ${subs[@]}"

for sub in "${subs[@]}"; do
  write_log ${log_args} "Analyzing subject: $sub"
	#### Create directories and prepare T1 file ####
	if [[ -d ${out_dir}/${sub} ]]; then
		write_log ${log_args} "WARNING: Subject's directory already exists in ${out_dir}. Deleting..."
		rm -rf ${out_dir}/${sub}
	fi

	mkdir -p ${out_dir}/${sub}
	write_log ${log_args} "Created directory ${out_dir}/${sub}"

	write_log ${log_args} "Copying anatomical to ${out_dir}/${sub}"
	cp ${in_dir}/${sub}/anat/*T1w.nii.gz ${out_dir}/${sub}/

	write_log ${log_args} "Unzipping anatomical in ${out_dir}/{$sub}"
	gunzip ${out_dir}/${sub}/*T1w.nii.gz

	mkdir -p ${out_dir}/${sub}/mri/orig
	write_log ${log_args} "Created directory ${out_dir}/${sub}/mri/orig"

	write_log ${log_args} "Converting anatomical ${out_dir}/${sub}/mri/orig/001.mgz"
	mri_convert ${out_dir}/${sub}/*T1w.nii ${out_dir}/${sub}/mri/orig/001.mgz

	#### Run recon-all ####
	write_log ${log_args} "Running recon-all on subject: ${sub}"
	recon-all -subjid ${sub} -all -parallel -openmp 4 -nuintensitycor

done
