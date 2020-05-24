#!/usr/bin/env Rscript

# within regions that encode degree and distance, does context modulate encoding?

# load packages
library(plyr)
library(dplyr)
library(tidyr)
if(!require(pracma)) install.packages("pracma"); require(pracma)
if(!require(oro.nifti)) install.packages("oro.nifti"); require(oro.nifti)


save_parcel <- function(values,file){
  outmap <- parcel
  outmap@datatype <- 16
  outmap@bitpix <- 32
  for (i in 1:n_parcels){
    outmap@.Data[outmap@.Data==i] <- values[i]
  }
  writeNIfTI(outmap,file)
  print(paste(file, "saved."))
}


rsa_dir <- '../bids/derivatives/level1/rsa/'
level2_dir <- '../bids/derivatives/level2/'
preds_sel <- c("dist_cat-sn", "deg_cat-sn")
n_parcels=200
parc_label=as.character(n_parcels)
stat_label='t'
corr_label='spear'
# MNI parcellation
parcel_fname <- paste0('../bids/tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-',parc_label,'Parcels7Networks_dseg.nii.gz')
parcel <- round(readNIfTI(parcel_fname))

main_thresh <- .05
main_col <- "p_uncorr"
diff_thresh <- .05

in_fname_main <- paste0(level2_dir,
                   "stat-",stat_label,
                   "_corr-",corr_label,
                   "_parc-", parc_label,
                   "_tail-1_effect-main.csv")
main_csv <- read.csv(in_fname_main)

in_fname_diff <- paste0(level2_dir,
                        "stat-",stat_label,
                        "_corr-",corr_label,
                        "_parc-", parc_label,
                        "_tail-1_effect-relevant.csv")
diff_csv <- read.csv(in_fname_diff)

out_fname <- paste0(level2_dir,
                    "stat-",stat_label, 
                    "_corr-", corr_label,
                    "_parc-", parc_label,
                    "_tail-1_effect-relevant-in-main")

data_fnames <- dir(rsa_dir,
                   pattern = paste0("\\_stat-", stat_label,
                                    "_corr-", corr_label,
                                    "_parc-", parc_label,
                                    "_val-r",
                                    "_roi_stats.csv$"),
                   recursive=T, full.names = T)
print('Reading in files: ')
print(data_fnames)
all_data <- do.call(rbind,lapply(data_fnames,read.csv))
all_data <- as.data.frame(all_data)
# drop row numbers
all_data <- all_data[,!(names(all_data) %in% c('X'))]
orig_data <- all_data

# separate tasks into 2 columns
all_data_task <- all_data %>% spread(key = task, value = r)

# calculate differences
print("Calculating differences")
all_data_task$diff_rel <- NA
all_data_task[all_data_task$predictor=='dist_cat-sn','diff_rel'] <- all_data_task[all_data_task$predictor=='dist_cat-sn','friend'] - all_data_task[all_data_task$predictor=='dist_cat-sn','number']
all_data_task[all_data_task$predictor=='deg_cat-sn','diff_rel'] <- all_data_task[all_data_task$predictor=='deg_cat-sn','number'] - all_data_task[all_data_task$predictor=='deg_cat-sn','friend']

# get significance levels and corrected p-values
get_sig <- function(p.value, sig_levels=c(0.001, 0.01, 0.05, 0.10),
                    sig_symbols = c('***','**','*','+'), ns_symbol='n.s.',
                    bonf_correction=T,
                    n_tests = 1)
{
  # order significance levels smallest to largest
  sig_levels <- sort(sig_levels)
  sig_symbols <- sig_symbols[order(sig_levels)]
  # correct significance levels for number of tests
  if (bonf_correction) {
    sig_levels <- sig_levels/n_tests
    corr_p.value <- p.value * n_tests
  }
  # iterate through levels and return relevant symbol
  sig_symbol <- ns_symbol
  for (i in seq(length(sig_levels))){
    if (p.value < sig_levels[i]) {
      sig_symbol <- sig_symbols[i]
      break
    }
  }
  return(c(sig_symbol, corr_p.value))
}

ttest_output <- function(data_vector, n_tests1 = 1, n_tests2 = 1, 
                         alternative="greater", 
                         mu = 0, conf.level = 0.95)
{
  tres <- t.test(data_vector, alternative="greater", mu = 0, conf.level = 0.95)
  pval <- tres$p.value[[1]]
  rval <- tres$estimate[[1]] # mean
  tval <- tres$statistic[[1]]
  dfs <- tres$parameter[[1]]
  res <- c(rval, tval, dfs, pval)
  return(res)
}



#cnames <- c('roi', 'predictor', 'r', 't', 'df', 'p_uncorr' ) 
out_df <- NULL
for (p in preds_sel){
  # filter main effect dataframe to only include this predictor
  mdf <- filter(main_csv, predictor==p)
  # get all ROIs where this predictor is below threshold
  sel <- mdf[mdf[,main_col] < main_thresh,"roi"]
  # get these ROIS from relevance effect dataframe
  ddf <- filter(diff_csv, predictor==p, roi %in% sel)
  # get uncorrected p-values
  pvals <- ddf$p_uncorr
  # remove all ajusted p-values
  ddf <- select(ddf, colnames(ddf)[!(colnames(ddf) %in% p.adjust.methods)])
  # recalculate adjustments
  for (pam in p.adjust.methods){
    pvals_adj <- p.adjust(pvals, method=pam)
    ddf[,pam] <- pvals_adj
    if (sum(pvals_adj < diff_thresh) > 0){
      out_pvals <- rep(1, n_parcels)
      out_pvals[sel] <- pvals_adj
      fname <- paste0(out_fname,
                      "_pred-", p,
                      "_pcorr-",pam,
                      "_dir-rev")
      save_parcel(1-out_pvals, fname)
    }
  }
  out_df <- rbind(out_df, ddf)
  
  # run t-test on all ROIs as one
  # extract ROIs and predictor
  df <- filter(all_data_task, roi %in% sel, predictor==p)
  # get mean value across ROIs per subject
  df <- ddply(df, .(sub, predictor), summarize, diff_mean=mean(diff_rel))
  # t-test
  res <- ttest_output(df$diff_mean, n_tests1 = 1, n_tests2 = 1)
  r <- c("mean", p, res, rep(NA, length(p.adjust.methods)))
  # add to output
  out_df <- rbind(out_df, r)
}
write.csv(out_df, file = paste0(out_fname,".csv"), row.names = F)
