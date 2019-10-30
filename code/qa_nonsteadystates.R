#library(tidyverse)
library(readr)

exclude_subs=c('sub-204')

subs = commandArgs(trailingOnly=TRUE)
if (length(subs)==0) {
  subs=c('all')
  #stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

tr <- 0.75
rm_n_trs <- 6 # remove 6 TRs

bids_dir <- '../bids'
#setwd(bids_dir)

# path to confounds folder
fmriprep_dir <- paste0(bids_dir, '/prep/fmriprep')
# list all directories in this folder (subject directories)
all_fmriprep_dirs <- list.dirs(path = fmriprep_dir, full.names = FALSE, recursive = FALSE)

out_dir <- paste0(bids_dir, '/level1/derivatives')

# combine all subjects' confounds tsv files
all_dfs <- c()
for (d in all_fmriprep_dirs) {
  # if this is a subjects folder...
  if (startsWith(d,'sub-')){
    # subject ID
    sub <- d
    if (sub %in% subs | subs == "all") {
      if (sub %in% exclude_subs){
        print(paste('Excluding subject ', sub))
      } else {
        raw_sub_dir <- paste0(bids_dir,'/',sub,'/func')
        events_fnames <- list.files(path = raw_sub_dir, pattern='*tsv', full.names = FALSE)
        # find all tsv files for this subject
        confound_fnames <- list.files(path = paste0(fmriprep_dir,'/',d,'/func'),pattern='*tsv', full.names = TRUE)
        # setup this subject's dataframe
        sub_df=c()
        for (f in confound_fnames){
          base_f <- substr(basename(f), 1,26)
          print(paste('Reading in files starting with:', base_f))
          events_f <- paste0(raw_sub_dir,'/',events_fnames[startsWith(events_fnames, base_f)])
          events_csv <- read.csv(events_f, sep='\t')
          # index for task name in filename
          task_i <- regexpr("task-", f)[1]+5
          # task name
          task <- substring(f,task_i,task_i+5)
          # index for run number in filename
          run_i <- regexpr("run-", f)[1]+4
          # run number
          run <- substr(f,run_i,run_i+1)

          # read in this file
          confound_csv <- read.csv(f, sep = '\t') #, na.strings = 'n/a')
          confound_csv_out <- tail(confound_csv, -rm_n_trs)
          confound_basef <- basename(f)
          out_fname <- paste0(substr(confound_basef, 1, 41), '_rmtr-', as.character(rm_n_trs), substr(confound_basef, 42, 56))
          readr::write_tsv(confound_csv_out, path = paste0(out_dir,'/',sub,'/func/', out_fname))
          print(paste(paste0(out_dir,'/',sub,'/func/', out_fname),"saved"))

          # get number of non-steady state outliers
          nonsteady <- grepl( "non_steady_state_outlier" , names( confound_csv ) )
          # subselect movement parameters
          confound_csv <- as.data.frame(confound_csv[, nonsteady])
          nonsteady_rows <- c()
          nonsteady_time <- 0
          # if nonsteady state columns found...
          if (ncol(confound_csv) > 0){
            warn <- F
            # iterate through columns
            for (i in seq(1,ncol(confound_csv))){
              r <- which(grepl(1, confound_csv[,i]))
              nonsteady_rows <- c(nonsteady_rows, r)
              # check if consecutive
              if (i == 1){
                if (r != 1){
                  print('WARNING: First non-steady state is not first TR!')
                  print(f)
                  warn <- T
                }
              } else {
                r_prev <- nonsteady_rows[i-1]
                if (r != r_prev + 1){
                  print('WARNING: Non-steady states are not consecutive!')
                  print(f)
                  warn <- T
                }
              }
              # get time
              nonsteady_time <- max(nonsteady_rows) * tr
            }
            c_df <- data.frame(nonsteady_rows, task, run, sub)
            # combine this file's data to subject's dataframe
            sub_df <- rbind(sub_df, c_df)
          }
          # recalculate onsets
          events_csv$onset_corrected_by_sub <- events_csv$onset - nonsteady_time
          events_csv$onset_corrected <- events_csv$onset - (rm_n_trs * tr)
          write_tsv(events_csv, path = events_f)
          print(paste(events_f,"saved"))
        }
        # append subject's dataframe to all_dfs
        all_dfs <- rbind(all_dfs, sub_df)
      }
    }
  }
}
fname=paste0(fmriprep_dir,'/nonsteady_outliers.tsv')
write_tsv(all_dfs, path=fname)
print(paste(fname,"saved"))
