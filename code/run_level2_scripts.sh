
# launch all second level scripts
# spear, correlation
qsub -N single_deg qsub_level2_fsl.sh spear deg correlation
qsub -N single_dist qsub_level2_fsl.sh spear dist correlation
qsub -N pair_deg qsub_level2_fsl_pairedt.sh spear deg correlation
qsub -N pair_dist qsub_level2_fsl_pairedt.sh spear dist correlation
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh spear deg correlation
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh spear dist correlation

# reg, correlation
qsub -N single_deg qsub_level2_fsl.sh reg deg correlation
qsub -N single_dist qsub_level2_fsl.sh reg dist correlation
qsub -N pair_deg qsub_level2_fsl_pairedt.sh reg deg correlation
qsub -N pair_dist qsub_level2_fsl_pairedt.sh reg dist correlation
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh reg deg correlation
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh reg dist correlation


# spear, euclidean
qsub -N single_deg qsub_level2_fsl.sh spear deg euclidean
qsub -N single_dist qsub_level2_fsl.sh spear dist euclidean
qsub -N pair_deg qsub_level2_fsl_pairedt.sh spear deg euclidean
qsub -N pair_dist qsub_level2_fsl_pairedt.sh spear dist euclidean
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh spear deg euclidean
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh spear dist euclidean

# reg, euclidean
qsub -N single_deg qsub_level2_fsl.sh reg deg euclidean
qsub -N single_dist qsub_level2_fsl.sh reg dist euclidean
qsub -N pair_deg qsub_level2_fsl_pairedt.sh reg deg euclidean
qsub -N pair_dist qsub_level2_fsl_pairedt.sh reg dist euclidean
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh reg deg euclidean
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh reg dist euclidean


# spear, mahalanobis
qsub -N single_deg qsub_level2_fsl.sh spear deg mahalanobis
qsub -N single_dist qsub_level2_fsl.sh spear dist mahalanobis
qsub -N pair_deg qsub_level2_fsl_pairedt.sh spear deg mahalanobis
qsub -N pair_dist qsub_level2_fsl_pairedt.sh spear dist mahalanobis
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh spear deg mahalanobis
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh spear dist mahalanobis

# reg, mahalanobis
qsub -N single_deg qsub_level2_fsl.sh reg deg mahalanobis
qsub -N single_dist qsub_level2_fsl.sh reg dist mahalanobis
qsub -N pair_deg qsub_level2_fsl_pairedt.sh reg deg mahalanobis
qsub -N pair_dist qsub_level2_fsl_pairedt.sh reg dist mahalanobis
qsub -N pairpred_deg qsub_level2_fsl_pairedt_preds.sh reg deg mahalanobis
qsub -N pairpred_dist qsub_level2_fsl_pairedt_preds.sh reg dist mahalanobis

# parc
qsub -N parc qsub_level2_rsa.sh
