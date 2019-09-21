#!/bin/bash
# run from local mount

source funcs
test_local "copy_behav.sh"
# get subjects
convert_sub_args -f subid -c "$@"
subs=( "${SUBS[@]}" )

for sub in "${subs[@]}"; do

	behav_dir="${HOME}/*Google*Drive/EncodingStudy/encoding-all-tasks/encoding_bids/${sub}/func"
	bids_dir="../bids/${sub}/"
	# if subject's bids directory doesn't exist, make it
	mkdir -p $bids_dir
	# copy files from google drive to hoffman
	echo "Copying ${behav_dir} to ${bids_dir}"
	cp -r ${behav_dir} ${bids_dir}

	# # convert behavioral files to onset files
	# convert_sub_args -f numid -c $sub
	# s="${SUBS}"
	# echo "Running mk_onsets.R on subject $sub"
	# Rscript --vanilla --verbose mk_onsets.R $s

	echo "DONE."
done
