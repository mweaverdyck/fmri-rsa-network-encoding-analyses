#!/bin/sh
# copy_dcm.sh sub-###
# copies subject's data from dicom server to hoffman directory
# converts data to bids format using dcm2bids


set -e
#. /u/local/Modules/default/init/modules.sh
source funcs
test_hoffman "copy_dcm.sh"
# get subjects
convert_sub_args -f subid -c "$@"
subs=( "${SUBS[@]}" )

# dicom server directory where subject data is
dcm_dir="ccn,mweaverd,/dicom/data/PARKINSONGROUP"

for sub in "${subs[@]}"; do
	sub_dir="${BIDS_DIR}/${sub}/"
	# if subject's bids directory doesn't exist, make it
	mkdir -p $sub_dir
	# copy data
	setup_subject_2 --download ${dcm_dir} ${BIDS_DIR} ${sub}
	# delete extra directories
	rm -r "${sub_dir}analysis"
	rm -r "${sub_dir}behav"
	rm -r "${sub_dir}dicom"
done

# dicom to bids (not needed for this study)
#dcm2bids -d "${BIDS_DIR}/${sub}" -p "${sub}" -c "${BIDS_DIR}/config_dcm2bids"
# BIDS-validator
