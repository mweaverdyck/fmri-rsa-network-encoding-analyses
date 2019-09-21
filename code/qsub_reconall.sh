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
echo "sourcing funcs"
source funcs
echo "loading python and freesurfer"
setup_modules python/2.7.13 freesurfer/6.0.0
export FREESURFER_HOME="/u/project/CCN/apps/freesurfer/6.0.0"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
echo "Setting project directories"
export STUDY_DIR=${BIDS_DIR}
echo "STUDY_DIR is ${STUDY_DIR}"
export SUBJECTS_DIR=${STUDY_DIR}/freesurfer
echo "SUBJECTS_DIR is ${SUBJECTS_DIR}"

# get subjects
convert_sub_args -f subid -c "$@"
subs=( "${SUBS[@]}" )

echo "Running recon-all on subjects: ${subs[@]}"

for sub in "${subs[@]}"; do

	#### Create directories and prepare T1 file ####
	if [[ -d ${SUBJECTS_DIR}/${sub} ]]; then
		echo "WARNING: Subject's directory already exists in ${SUBJECTS_DIR}. Deleting..."
		rm -rf ${SUBJECTS_DIR}/${sub}
	fi

	mkdir -p ${SUBJECTS_DIR}/${sub}
	echo "Created directory ${SUBJECTS_DIR}/${sub}"

	echo "Copying anatomical to ${SUBJECTS_DIR}/${sub}"
	cp ${STUDY_DIR}/${sub}/anat/*T1w.nii.gz ${SUBJECTS_DIR}/${sub}/

	echo "Unzipping anatomical in ${SUBJECTS_DIR}/{$sub}"
	gunzip ${SUBJECTS_DIR}/${sub}/*T1w.nii.gz

	mkdir -p ${SUBJECTS_DIR}/${sub}/mri/orig
	echo "Created directory ${SUBJECTS_DIR}/${sub}/mri/orig"

	echo "Converting anatomical ${SUBJECTS_DIR}/${sub}/mri/orig/001.mgz"
	mri_convert ${SUBJECTS_DIR}/${sub}/*T1w.nii ${SUBJECTS_DIR}/${sub}/mri/orig/001.mgz

	#### Run recon-all ####
	echo "Running recon-all on subject: ${sub}"
	recon-all -subjid ${sub} -all -parallel -openmp 4 -nuintensitycor

done
