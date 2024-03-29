#!/usr/bin/env Rscript
# load funcs variables
project_dir <- system("source ~/.bash_profile; echo ${data_dir}", intern=T)
code_dir = paste(project_dir, "code", sep="/")
source(paste(code_dir,"funcs_get_vars.R", sep='/'))

# load packages
library(plyr)
library(dplyr)
library(tidyr)
if(!require(pracma)) install.packages("pracma"); require(pracma)
if(!require(oro.nifti)) install.packages("oro.nifti"); require(oro.nifti)

# MNI parcellation
parcellation <- round(readNIfTI(MNI_PARCELLATION))

corr_label='spear'
val_label=if (corr_label=='spear') 'r' else 'beta'

out_fname <- paste0("stat-",STAT, 
                    "_corr-", corr_label,
                    "_parc-", PARC_LAB,
                    "_tail-1")


# function to save stats into nifti
save_parcel <- function(values,file){
  outmap <- parcellation
  outmap@datatype <- 16
  outmap@bitpix <- 32
  for (i in 1:N_PARCELS){
    outmap@.Data[outmap@.Data==i] <- values[i]
  }
  writeNIfTI(outmap,file)
  print(paste(file, "saved."))
}

# read in all subject files into single data frame
rsa_dir <- paste0(RSA_DIR,'/*/',corr_label)
out_dir <- paste0(SECOND_LEVEL_DIR,'/parc',PARC_LAB,'/')
data_fnames <- dir(Sys.glob(rsa_dir),
                   pattern = paste0("\\_stat-", STAT,
                                    "_corr-", corr_label,
                                    "_parc-", PARC_LAB,
                                    "_val-", val_label,
                                    "_roi_stats.csv$"),
                   recursive=F, full.names = T)
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
if (length(rois) != N_PARCELS){
  print("WARNING: number of ROIs does not match N_PARCELS:")
  print(str(length(rois)), PARC_LAB)
}

# separate tasks into 2 columns
all_data_task <- all_data %>% spread(key = task, value = r)

# # calculate average across tasks
# all_data_task$avg <- (all_data_task$friend + all_data_task$number)/2

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
    # test relavent difference
    if (p %in% c("dist_cat-sn", "deg_cat-sn")){
      res <- ttest_output(df$diff_rel, n_tests1 = length(rois), n_tests2 = length(gm_rois))
      results_diffrel <- rbind(results_diffrel, c(r, p, res))
    }
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
for (pred in c("dist_cat-sn", "deg_cat-sn")){#predictors){
  pvals_diffrel <- results_diffrel[results_diffrel$predictor == pred,"p_uncorr"]
  pvals_avg <- results_avg[results_avg$predictor == pred, 'p_uncorr']
  for (pam in p.adjust.methods){
    if (sum(results_diffrel$predictor == pred) > 0){
      pvals_diffrel_adjust <- p.adjust(pvals_diffrel, method=pam)
      results_diffrel[results_diffrel$predictor == pred,pam] <- pvals_diffrel_adjust
      # only save nifti if there is a significant parcel
      if (sum(pvals_diffrel_adjust < .05) > 0){
        fname <- paste0(out_dir,
                        out_fname,
                        "_pred-", pred,
                        "_effect-relevant",
                        "_pcorr-",pam,
                        "_dir-rev")
        save_parcel(1-pvals_diffrel_adjust, fname)
      }
    }
    
    pvals_avg_adjust <- p.adjust(pvals_avg, method=pam)
    results_avg[results_avg$predictor == pred,pam] <- pvals_avg_adjust
    # only save nifti if there is a significant parcel
    if (sum(pvals_avg_adjust < .05) > 0){
      fname <- paste0(out_dir,
                      out_fname,
                      "_pred-", pred,
                      "_effect-main",
                      "_pcorr-",pam,
                      "_dir-rev")
      save_parcel(1-pvals_avg_adjust, fname)
    }
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


