library(dplyr)
library(plyr)
library(readr)
library(tidyr)

parc_dirs=c("parc-2soc1",  "parc-2socclust9", "parc-fun200", "parc-sch200")
dist_dirs=c("correlation", "euclidean")

setwd('~/Documents/hoffman2/bids/derivatives/level2')

all_dat <- data.frame()
for (pd in parc_dirs){
  for (dd in dist_dirs){
    myfiles = list.files(path = paste(pd,dd, sep="/"),pattern="*val2-all_correction-all_roi_stats.csv", full.names=TRUE)
    dat_csv = ldply(myfiles, read_csv)
    dat_csv$parc=pd
    dat_csv$dist=dd
    all_dat <- rbind(all_dat, dat_csv)
  }
}

all_cols <- colnames(all_dat)
stats_cols <- c("estimate", "t", "p")
# dat_avg_wide = filter(dat_csv, corr=="spear", test=='avg', FDR < 0.05)
# dat_avg_wide = select(dat_avg_wide, roi, pred, FDR)
# dat_avg_wide = spread(dat_avg_wide, pred, FDR)

sel_dat <- all_dat %>% 
  select(all_cols[!all_cols %in% stats_cols]) %>% 

dat_wide <- sel_dat %>% 
  filter(FDR < 0.05) %>% 
  spread("dist", "FDR")
View(dat_wide)
write.csv(dat_wide, file = "parc-all_correction-fdr_thresh-05.csv", row.names = F)

dat_wide01 = sel_dat %>% filter(p < 0.01) %>% spread("dist", "p")
View(dat_wide01)
write.csv(dat_wide01, file = "parc-all_correction-none_thresh-01.csv", row.names = F)

# dat_diff_wide = filter(dat_csv, corr=="spear", test=='diff', FDR < 0.05)
# dat_diff_wide = select(dat_diff_wide, roi, pred, FDR)
# dat_diff_wide = spread(dat_diff_wide, pred, FDR)

