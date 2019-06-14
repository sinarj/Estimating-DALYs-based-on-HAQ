##########################################################
# Description: Create HAQ index 2016 for all the most
# detailed locations (age standardized)
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())
library(argparse, lib.loc = "/home/j/temp/vishnn/packages")
library(findpython, lib.loc= '/home/j/temp/vishnn/packages')
library(data.table)
#library(dplyr)
library(broom)
library(parallel)
if (Sys.info()[1] == "Linux"){
  j <- "/home/j/"
  h <- paste0("/homes/", Sys.info()[7], "/haq/haq")
}else if (Sys.info()[1] == "Windows"){
  j <- "J:/"
  h <- "H:/"
}else if (Sys.info()[1] == "Darwin"){
  j <- "/Volumes/snfs/"
  h <- paste0("/Volumes/", Sys.info()[6], "/")
}

##########################################################
## PARSE ARGUMENTS
parser <- ArgumentParser()
parser$add_argument("--prog_dir", help='Where scripts are stored',
                    default=paste0(h, "/risk_standardized_amen_mort"), type="character")
parser$add_argument("--data_dir", help="Space where results are to be loaded from and saved to",
                    default="/share/scratch/projects/hssa/haq/HAQ_2017/haq_US/", type="character")
parser$add_argument("--base_lsid", help="Base location set",
                    default=35, type="integer")
parser$add_argument("--agg_lsids", help="Aggregation location set(s)",
                    default=40, nargs="+", type="integer")
parser$add_argument("--yids", help="Years to run",
                    default=c(seq(1990, 2010, 5), 2017), nargs="+", type="integer")

args <- parser$parse_args()
list2env(args, .GlobalEnv)
rm(args)

##########################################################
## DEFINE FUNCTIONS
source(paste0(prog_dir, "/utilities.R"))
library(dplyr)

#locsdf <- getLocations(base_lsid)
amenabledf_2017 <- data.table()
amenabledf_2017_MI_ratios <- data.table()

amenabledf_2017<-fread(paste0(data_dir,"combined_data/global_data" ,"/allMostDet_ageStd_data_noncancer.csv"))
amenabledf_2017 <- amenabledf_2017 %>% filter(!(cause_id %in% c(429,432,435,441,468,484,487,849)))
#amenabledf_2017 <- amenabledf_2017 %>% filter(age_group_id %in% 2:19)

amenabledf_2017_MI_ratios <-fread(paste0(data_dir,"combined_data/global_data" ,"/allMostDet_ageStd_data_cancer.csv"))
amenabledf_2017_MI_ratios <- amenabledf_2017_MI_ratios %>% filter(cause_id %in% c(429,432,435,441,468,484,487,849))
#amenabledf_2017_MI_ratios <- amenabledf_2017_MI_ratios %>% filter(age_group_id %in% 2:19)

minmaxdf<-fread(paste0(data_dir, "/combined_data/min_max/minmax_noncancer_ageStd.csv"))
# Drop cancers off of min/max, calculate min/max for GBD 2017 MI ratios without addition of log offset
minmaxdf$cause_id<-abs(minmaxdf$cause_id)

# Create new df to calculate min/max for GBD 2017 MI ratios 
minmaxdf_MI <- fread(paste0(data_dir, "/combined_data/min_max/minmax_cancer_ageStd.csv"))
minmaxdf_MI$cause_id<-abs(minmaxdf_MI$cause_id)
# Rbind togetther
minmaxdf<-rbind(minmaxdf,minmaxdf_MI)
minmaxdf <- minmaxdf %>% mutate(rs_min=ifelse(rs_min < 0,0,rs_min))
print("Combining MI Ratios/RSD Deaths then scaling to upper/lower bounds")

# Rbind MI ratios for cancers back onto amenabledf
amenabledf_2017 <-rbind(amenabledf_2017, amenabledf_2017_MI_ratios)

amenabledf_2017 <- amenabledf_2017 %>% select(-c(val))
amenabledf_2017 <- amenabledf_2017[complete.cases(amenabledf_2017),]
amenabledf_2017 <- amenabledf_2017 %>% filter(rsval != Inf)
amenabledf_2017 <- amenabledf_2017 %>% mutate(rsval = rsval + 1e-6)
minmaxdf <- minmaxdf %>% mutate(rs_min = rs_min + 1e-6)
minmaxdf <- minmaxdf %>% mutate(rs_max = rs_max + 1e-6)
# Scale amenabledf to min/max
indexdf <- inner_join(amenabledf_2017, minmaxdf,
                      by = c("cause_id", "draw","age_group_id"))
# indexdf<-indexdf[,c("age_end.x","age_start.x","age_end.y","age_start.y"):=NULL]
indexdf <- indexdf[complete.cases(indexdf),]
#indexdf <- indexdf[, rsval := rsval + 1e-6]
indexdf <- indexdf %>% mutate(rsval = ifelse(rsval < rs_min, rs_min, rsval),
                              rsval = ifelse(rsval > rs_max, rs_max, rsval))

indexdf <- indexdf %>% mutate(log_index_cause = (1 - ((log(rsval) - log(rs_min)) / (log(rs_max) - log(rs_min)))) * 100)
indexdf <- indexdf[complete.cases(indexdf),]
indexdf <- indexdf %>% mutate(log_index_cause=ifelse(log_index_cause < 0,0,log_index_cause),
                              log_index_cause=ifelse(log_index_cause> 100, 100,log_index_cause))

weightsdf <- fread(paste0(data_dir, "/weights/cause_weights_ageStd_2017.csv"))
weightsdf$cause_id <- as.integer(weightsdf$cause_id)

indexdf <- inner_join(indexdf,
                      weightsdf[!is.na(cause_weight), c("cause_id", "cause_weight")],
                      by = c("cause_id"))

indexdf <- indexdf %>% mutate(cause_weight = ifelse(is.na(cause_weight),0,cause_weight))
setDT(indexdf, keep.rownames=FALSE, key=NULL, check.names=FALSE)
aggdf <- indexdf[, .(log_index_cause = sum(log_index_cause * cause_weight),
                     log_index_cause_mean = mean(log_index_cause),
                     log_index_cause_geom_mean = exp(mean(log(log_index_cause + 1)))), by = .(location_id, age_group_id, year_id, draw)]

aggdf$cause_id <- 100
indexdf <- rbind(indexdf[, c("location_id", "year_id",  "age_group_id","cause_id", "draw", "log_index_cause"), with = FALSE],
                 aggdf[, c("location_id", "year_id", "age_group_id","cause_id", "draw", "log_index_cause", "log_index_cause_mean", "log_index_cause_geom_mean"), with = FALSE],
                 fill = TRUE)

indexdf <- indexdf[log_index_cause < 0, log_index_cause := 0][log_index_cause > 100, log_index_cause := 100]

indexdf <- indexdf[, .(index = mean(log_index_cause),
                       index_lval = quantile(log_index_cause, 0.025),
                       index_uval = quantile(log_index_cause, 0.975),
                       #index_efa = mean(log_index_cause_efa),
                       index_mean = mean(log_index_cause_mean),
                       index_geom_mean = mean(log_index_cause_geom_mean)),
                   by = .(location_id, age_group_id,year_id, cause_id)] 

print("Saving compiled scaled components")

#first the 32 amenable causes
index2016_all_ageStd <-indexdf[year_id==2016][cause_id!=100]
write.csv(index2016_all_ageStd, file = paste0(data_dir, "/results/haq_2016_bycause_ageStd_allMostDetLocs.csv"))
#also making a local copy on my drive
write.csv(index2016_all_ageStd, file = paste0("/ihme/homes/arjuns13/notebooks/Documents/Data/haq_2016_bycause_ageStd_allMostDetLocs.csv"))

#now the aggregated cause with cause_id=100
index2016_agg_ageStd <-indexdf[year_id==2016][cause_id==100]
write.csv(index2016_agg_ageStd, file = paste0(data_dir, "/results/haq_2016_by_allcause_ageStd_allMostDetLocs.csv"))
#also making a local copy on my drive
write.csv(index2016_agg_ageStd, file = paste0("/ihme/homes/arjuns13/notebooks/Documents/Data/haq_2016_by_allcause_ageStd_allMostDetLocs.csv"))

