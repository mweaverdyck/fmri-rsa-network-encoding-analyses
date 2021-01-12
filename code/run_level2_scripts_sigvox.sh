
# test if contextual relevance modulates encoding of distance and degree in 
# voxels that encode each predictor

qsub -N pair_deg qsub_level2_fsl_pairedt.sh spear deg correlation
qsub -N pair_dist qsub_level2_fsl_pairedt.sh spear dist correlation
