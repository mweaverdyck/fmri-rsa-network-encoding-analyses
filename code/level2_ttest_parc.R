#!/usr/bin/env Rscript

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


n_parcels=800
parc_label=as.character(n_parcels)
stat_label='t'
corr_label='spear'

out_fname <- paste0("stat-",stat_label, 
                    "_corr-", corr_label,
                    "_parc-", parc_label,
                    "_tail-1")

# MNI parcellation
parcel_fname <- paste0('../bids/mni_icbm152_nlin_asym_09c/Schaefer/tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-',parc_label,'Parcels17Networks_dseg.nii.gz')
parcel <- round(readNIfTI(parcel_fname))

# read in all subject files into single data frame
rsa_dir <- '../bids/level1/rsa/'
out_dir <- '../bids/level2/'
data_fnames <- dir(rsa_dir,
                   pattern = paste0("\\_stat-", stat_label,
                                    "_corr-", corr_label,
                                    "_parc-", parc_label,
                                    "_roi_stats.csv$"),
                   recursive=T, full.names = T)
print('Reading in files: ')
print(data_fnames)
all_data <- do.call(rbind,lapply(data_fnames,read.csv))
all_data <- as.data.frame(all_data)
# drop row numbers
all_data <- all_data[,!(names(all_data) %in% c('X'))]
orig_data <- all_data

# get predictors
predictors <- unique(all_data$predictor)
print("Predictors found: ")
print(predictors)

# get tasks
tasks <- unique(all_data$task)
print("Tasks found: ")
print(tasks)

# get list of ROIs
rois <- unique(all_data$roi)
print("Number of ROIs found: ")
print(length(rois))
if (length(rois) != n_parcels){
  print("WARNING: number of ROIs does not match n_parcels:")
  print(str(length(rois)), parc_label)
}

# separate tasks into 2 columns
all_data_task <- all_data %>% spread(key = task, value = r)

# calculate average across tasks
all_data_task$avg <- (all_data_task$friend + all_data_task$number)/2

# calculate differences
print("Calculating differences")
all_data_task$diff_rel <- NA
all_data_task[all_data_task$predictor=='dist','diff_rel'] <- all_data_task[all_data_task$predictor=='dist','friend'] - all_data_task[all_data_task$predictor=='dist','number']
all_data_task[all_data_task$predictor=='deg','diff_rel'] <- all_data_task[all_data_task$predictor=='deg','number'] - all_data_task[all_data_task$predictor=='deg','friend']

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
  #sig_pcor1 <- get_sig(pval, n_tests=n_tests1)
  #pval_corr1 <- as.numeric(sig_pcor1[2])
  #sig_symb1 <- sig_pcor1[1]
  #sig_pcor2 <- get_sig(pval, n_tests=n_tests2)
  #pval_corr2 <- as.numeric(sig_pcor2[2])
  #sig_symb2 <- sig_pcor2[1]
  # figure out which direction difference is going in
  #res <- c(rval, tval, dfs, pval, pval_corr1, sig_symb1, pval_corr2, sig_symb2)
  res <- c(rval, tval, dfs, pval)
  return(res)
}


# iterate through each ROI and perform t-test
print("Running one-sample, one-tailed t-tests")
results_diffrel <- c()
cnames_diffrel <- c('roi', 'predictor', 'r', 't', 'df', 'p_uncorr' ) 
                    #, 'p_corr_all','sig_all', 'p_corr_gm', 'sig_gm')
results_avg <- c()
cnames_avg <- c('roi', 'predictor', 'r', 't', 'df', 'p_uncorr' ) 
for (r in unique(all_data_task$roi)) {
  # subselect this ROI's data
  roi_df <- filter(all_data_task, roi==r)
  # run t-test on each predictor
  for (p in predictors){
    df <- filter(roi_df, predictor==p)
    res <- ttest_output(df$diff_rel, n_tests1 = length(rois), n_tests2 = length(gm_rois))
    results_diffrel <- rbind(results_diffrel, c(r, p, res))
    # test average r
    res <- ttest_output(df$avg)
    results_avg <- rbind(results_avg, c(r, p, res))
  }
}
# convert to data frames
results_diffrel <- as.data.frame(results_diffrel)
colnames(results_diffrel) <- cnames_diffrel
results_diffrel$p_uncorr <- as.numeric(as.character(results_diffrel$p_uncorr))
results_avg <- as.data.frame(results_avg)
colnames(results_avg) <- cnames_avg
results_avg$p_uncorr <- as.numeric(as.character(results_avg$p_uncorr))

# calculate corrected p-values
for (pred in predictors){
  pvals_diffrel <- results_diffrel[results_diffrel$predictor == pred,"p_uncorr"]
  pvals_avg <- results_avg[results_avg$predictor == pred, 'p_uncorr']
  for (pam in p.adjust.methods){
    pvals_diffrel_adjust <- p.adjust(pvals_diffrel, method=pam)
    results_diffrel[results_diffrel$predictor == pred,pam] <- pvals_diffrel_adjust
    fname <- paste0(out_dir,
                    out_fname,
                    "_pred-", pred,
                    "_effect-relevant",
                    "_pcorr-",pam)
    save_parcel(1-pvals_diffrel_adjust, fname)
    
    pvals_avg_adjust <- p.adjust(pvals_avg, method=pam)
    results_avg[results_avg$predictor == pred,pam] <- pvals_avg_adjust
    fname <- paste0(out_dir,
                    out_fname,
                    "_pred-", pred,
                    "_effect-main",
                    "_pcorr-",pam)
    save_parcel(1-pvals_avg_adjust, fname)
  }
}

# write out csvs
fname = paste0(out_dir, 
               out_fname,
               "_effect-relevant.csv")
write.csv(results_diffrel, row.names=F,
          file = fname)
print(paste("File saved:", fname))

fname = paste0(out_dir, 
               out_fname,
               "_effect-main.csv")
write.csv(results_avg, row.names=F,
          file = fname)
print(paste("File saved:", fname))


