
# launch all second level scripts

qsub -N single_deg qsub_level2_fsl.sh spear deg
qsub -N single_dist qsub_level2_fsl.sh spear dist
qsub -N pair_deg qsub_level2_fsl_pairedt.sh spear deg
qsub -N pair_dist qsub_level2_fsl_pairedt.sh spear dist
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh spear deg
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh spear dist

qsub -N single_deg qsub_level2_fsl.sh reg deg
qsub -N single_dist qsub_level2_fsl.sh reg dist
qsub -N pair_deg qsub_level2_fsl_pairedt.sh reg deg
qsub -N pair_dist qsub_level2_fsl_pairedt.sh reg dist
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh reg deg
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh reg dist

#qsub -N parc qsub_level2_rsa.sh
