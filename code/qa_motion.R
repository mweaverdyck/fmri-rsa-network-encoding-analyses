library(tidyverse)
library(dplyr)
# Check motion regressors

exclude_subs=c('sub-204')

# path to confounds folder
fmriprep_dir <- '/Users/miriam/Documents/hoffman2/bids/fmriprep'
setwd(fmriprep_dir)
# list all directories in this folder (subject directories)
all_dirs <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
# combine all subjects' confounds tsv files
all_dfs <- c()
for (d in all_dirs) {
  # if this is a subjects folder...
  if (startsWith(d,'sub-')){
    # subject ID
    sub <- d
    if (sub %in% exclude_subs){
      print(paste('Excluding subject ', sub))
    } else {
      # find all tsv files for this subject
      confound_fnames <- list.files(path = paste0('./',d,'/func'),pattern='*tsv', full.names = TRUE)
      # setup this subject's dataframe
      sub_df=c()
      for (f in confound_fnames){
        print(paste('Reading in file:',f))
        # index for task name in filename
        task_i <- regexpr("task-", f)[1]+5
        # task name
        task <- substring(f,task_i,task_i+5)
        # index for run number in filename
        run_i <- regexpr("run-", f)[1]+4
        # run number
        run <- substr(f,run_i,run_i+1)
        # read in this file
        c <- read.csv(f, sep = '\t', na.strings = 'n/a')
        # subselect movement parameters
        c <- select(c, dvars, framewise_displacement, trans_x, trans_y, trans_z,
                        rot_x, rot_y, rot_z)
        # add task name to dataframe
        c$task <- task
        # add run number to dataframe
        c$run <- run
        # combine this file's data to subject's dataframe
        sub_df <- rbind(sub_df, c)
      }
      # add subject's ID to dataframe
      sub_df$subID <- sub
      # append subject's dataframe to all_dfs
      all_dfs <- rbind(all_dfs, sub_df)
    }
  }
}

# get all subjects
subs <- unique(all_dfs$subID)
# get all tasks
tasks <- unique(all_dfs$task)
# get all runs
runs <- unique(all_dfs$run)

# create new dataframe with one row per run per subject
summ_df <- expand.grid(subs, runs)
colnames(summ_df) <- c('subID', 'run')
# go through each run for each subject, calculate change between frames, 
# determine if these movements are too big (>2 or >.5), and count how many
# of each subject's runs fail this test
all_summ <- data.frame()
for (r in seq(nrow(summ_df))){
  # subject ID
  sub_curr <- summ_df[r,'subID']
  # run number
  run_curr <- summ_df[r,'run']
  # data from this run for this subject
  d <- filter(all_dfs, subID==sub_curr, run==run_curr)
  # task name
  task_curr <- unique(d$task)
  # calculate frame-by-frame changes
  # save the distance from the first frame
  first_diff <- sqrt(d[1,'trans_x']**2 + d[1,'trans_y']**2 + d[1,'trans_z']**2)
  # calculate differences between each frame and the previous frame
  xyz_diff <- data.frame(diff(as.matrix(select(d, trans_x, trans_y, trans_z))))
  # calculate what this difference is in terms of distance
  d$trans_dist <- c(first_diff, sqrt(xyz_diff$trans_x**2 + xyz_diff$trans_y**2 + 
                             xyz_diff$trans_z**2))
  # is it bigger than 2mm?
  d$gt2mm <- as.numeric(d$trans_dist > 2)
  # is it bigger than .5mm?
  d$gt_5mm <- as.numeric(d$trans_dist > .5)
  # does this run fail?
  summ <- ddply(d, .(subID, task, run), summarise, 
                n2=sum(gt2mm), n_5=sum(gt_5mm),
                fail=as.numeric(n2>=1 || n_5 >= 5))
  # add to total summary
  all_summ <- rbind(all_summ, summ)
}

write_csv(all_summ, paste0(fmriprep_dir,'/motion_counts.csv'))
# count how many runs fail per subject per task
out <- ddply(all_summ, .(subID, task), summarize, fails=sum(fail))
write_csv(out, paste0(fmriprep_dir,'/motion_summary.csv'))
# if at least half of the runs (i.e. 2) fail, add to failed subjects list
failed_subs <- filter(out, fails>=2)$subID
# print out failed subjects
print(paste(length(failed_subs), 'failed subjects:', failed_subs))

      