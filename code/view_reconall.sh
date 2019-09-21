#!/bin/bash
# use freeview to look at reconall output. Use No Machine to run
source funcs

setup_modules freesurfer/6.0.0
# get subjects
convert_sub_args -f subid -c "$@"
sub=${SUBS}

subdir="${BIDS_DIR}/freesurfer/$sub"
freeview -v \
${subdir}/mri/T1.mgz \
${subdir}/mri/wm.mgz \
${subdir}/mri/brainmask.mgz \
${subdir}/mri/aseg.mgz:colormap=lut:opacity=0.2 \
-f ${subdir}/surf/lh.white:edgecolor=blue \
${subdir}/surf/lh.pial:edgecolor=red \
${subdir}/surf/rh.white:edgecolor=blue \
${subdir}/surf/rh.pial:edgecolor=red
