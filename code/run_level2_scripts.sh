
# launch all second level scripts
qsub -N single_deg qsub_level2_fsl.sh spear deg
qsub -N single_dist qsub_level2_fsl.sh spear dist
qsub -N pair_deg qsub_level2_fsl_pairedt.sh spear deg
qsub -N pair_dist qsub_level2_fsl_pairedt.sh spear dist
#qsub -N pairpreds qsub_level2_fsl_pairedt_preds.sh
qsub -N parc qsub_level2_rsa.sh
