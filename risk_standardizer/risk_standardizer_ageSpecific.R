##########################################################
# Risk standardization values for all the amenable causes
# for the given years and locations. Locations for which 
# was run was all the most detailed locations (for use in
# the regressions) and the countries (for use in the min-max
# band calculations). This script was run for the age-group
# specific version
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())
library(argparse, lib.loc= '/home/j/temp/vishnn/packages')
library(findpython, lib.loc= '/home/j/temp/vishnn/packages')
library(plyr)
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
parser$add_argument("--prog_dir", help="Program directory",
                    default=paste0(h, "/risk_standardized_amen_mort/"), type="character")
parser$add_argument("--data_dir", help="Space where results are to be saved",
                    default="/share/scratch/projects/hssa/haq/HAQ_2017/haq_US/", type="character")
parser$add_argument("--prefix", help="PAF prefix",
                    default='PAF', type="character")
parser$add_argument("--lsid", help="Location set",
                    default=35, type="integer")
parser$add_argument("--lid", help="Location",
                    default=13, type="integer")
parser$add_argument("--yid", help="Year for current job",
                    default=2005, type="integer")
parser$add_argument("--mid", help="Measure",
                    default=1, type="integer")
parser$add_argument("--scale_ceiling", help="Upper limit of PAFs (out of 100)",
                    default=90, type="integer")
args <- parser$parse_args()
list2env(args, .GlobalEnv)
rm(args)

##########################################################
## DEFINE FUNCTIONS
source(paste0(prog_dir, "utilities.R"))

joinGlobal <- function(df, data_dir, lsid, mid, prefix) {
  # Load all years, and use average
  #data_dir = '/share/scratch/projects/hssa/haq/HAQ_2017/cc_90_paf_242_amenable/'
  globdf <- rbindlist(lapply(as.list(c(seq(1990, 2015, 5), 2016, 2017)), function(yid, data_dir, prefix) {
    load(paste0(data_dir, "/draws/inputs/", yid, "/", prefix, "_", mid, "_1.RData"))
    return(inputdf)
  },
  data_dir, prefix))
  if (prefix == "PAF") {
    # Aggregate attributable burden over all years and join
    globdf <- globdf[, .(glob_paf = sum(paf * val) / sum(val), val = sum(val)), by = .(age_group_id, sex_id, measure_id, cause_id, draw)]
    globdf <- globdf[val == 0, glob_paf := 0]
    globdf$val <- NULL
    df <- merge(df,
                globdf,
                by = c("age_group_id", "sex_id", "measure_id", "cause_id", "draw"))
  }
  return(df)
}

makeAgeWeights <- function(cause_list){
  ageweightdf <- getAgeGroups()
  cs_ageweightdf <- data.table()
  for (cid in cause_list$cause_id) {
    end <- cause_list[cause_id == cid]$age_group_years_end
    if(end > 65){
      cause_weights <- data.table(ageweightdf)
      cause_weights$cause_id <- cid
      cause_weights <- cause_weights[age_group_years_start >= 45 & age_group_years_end <= 75]
      cause_weights <- cause_weights[, age_group_weight_value := age_group_weight_value / sum(age_group_weight_value)]
      cs_ageweightdf <- rbind(cs_ageweightdf,
                              cause_weights[, c("cause_id", "age_group_id", "age_group_weight_value"), with = FALSE])
    }
  }
  return(cs_ageweightdf)
}

getWeightShareRatio <- function(df){
  df <- df[,.(location_id,year_id,age_group_id,cause_id,population,age_group_weight_value)]
  df <- df[!(duplicated(df))]
  df <- df[,c('pop_share','average_pop') := .(population/sum(population),mean(population)),by=.(cause_id)]
  df <- df[,age_weight_share_ratio := age_group_weight_value/pop_share]
  df <- df[,c('pop_share','population','age_group_weight_value'):=.(NULL,NULL,NULL)]
  return (df)
}

bad_locs <- c()
riskStandardize13 <- function(prefix, data_dir, yid, lid, mid, lsid, scale_ceiling) {
  # Read in PAFs and CoDCorrect for each most-detailed location
  load(paste0("/share/scratch/projects/hssa/haq/HAQ_2017/vanilla_cc_90_paf_242_amenable/draws/inputs/", yid, "/", prefix, "_", mid, "_", lid, ".RData"))
  
  # Load location hierarchy data frame
  locsdf <- prepLocHierarchy(data_dir, lsid)
  
  # Use deaths to aggregate PAF up to global
  inputdf <- joinGlobal(inputdf, data_dir, lsid, mid, prefix)
  
  if (prefix == "PAF") {
    # Apply scalars
    load(paste0(data_dir, "/draws/inputs/PAF_", mid, "_scalar_", scale_ceiling, ".RData"))
    inputdf <- merge(inputdf,
                     scalardf,
                     by = c("age_group_id", "sex_id", "measure_id", "cause_id"))
    inputdf <- inputdf[, c("paf", "glob_paf") := .(paf * pafscalar, glob_paf * pafscalar)][, pafscalar := NULL]
    
    # Produce risk-standardized deaths
    inputdf <- inputdf[, rsval := val * (1 - paf) * (1 / (1 - glob_paf))]
    inputdf <- inputdf[is.na(rsval), rsval := 0]
    
    # For PAFs of one and for diarrhea, use observed
    inputdf <- inputdf[paf == 1 | cause_id == 302, rsval := val]
  }
  
  # Only keep Nolte/McKee causes and ages
  cause_list <- fread(paste0(prog_dir, "/amenable_cause_list_GBD.csv"))
  causesdf <- prepCauseHierarchy(data_dir)
  
  # Aggregate GBD causes and append amenable
  inputdf <- merge(inputdf,
                   causesdf,
                   by = "cause_id")
  inputdf <- rbindlist(lapply(as.list(c("root", paste0("L", seq(1, max(causesdf$level))))), hierarchyAgg, inputdf, "cause_id", c("rsval", "val")))
  
  # Aggregate to Nolte/McKee list, attach to main data frame
  inputdf <- merge(inputdf,
                   cause_list[, c("cause_id", "age_group_id_start", "age_group_id_end"), with = FALSE],
                   by = "cause_id")
  
  # Aggregate sexes
  inputdf <- inputdf[, bothsex := 3]
  inputdf <- hierarchyAgg("bothsex", inputdf, "sex_id", c("rsval", "val"))
  
  
  #load(paste0(data_dir, "/popsdf_", lsid, ".RData"))
  inputdf <- merge(inputdf,
                   popsdf,
                   by = c("location_id", "year_id", "age_group_id", "sex_id"))
  
  inputdf <- inputdf[, .(rsval = (rsval/population),
                         val = (rsval/population)),
                     by = .(location_id, year_id, age_group_id, sex_id, cause_id, draw)]
  
  print(paste0("running for location ID:", lid))
  # Return data frame
  write.csv(inputdf[, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", "draw", "rsval", "val"), with = FALSE],
            file = paste0("/share/scratch/projects/hssa/haq/HAQ_2017/haq_US/", "std_all_most_det_locs/", yid, "/", lid, ".csv"),
            row.names = FALSE)
  bad_locs = bad_locs
  if(nrow(inputdf)<1 | length(unique(inputdf$rsval)) < 1)
  {
    print(paste0("Did not work for location ID:", lid))
    bad_locs <- c(bad_locs,lid)
  }
}

##########################################################
## RUN PROGRAM
# Risk standardize (should be changed from rsdeathsdf to rsdf, missed that in the update)

#you only need to do this once, and here we do for both sexes 1 and 2, cause we merge on that and then aggregate at the end
popsdf <- getPops(lsid, sids = c(1,2,3))

load(paste0("/share/scratch/projects/hssa/haq/HAQ_2017/vanilla_cc_90_paf_242_amenable", "/locsdf_", lsid, ".RData"))

loc_list = locsdf[locsdf$most_detailed == 1]$location_id

#yids <- c(seq(1990, 2010, 5), 2015, 2016, 2017)
source("/home/j/temp/central_comp/libraries/current/r/get_location_metadata.R")
yids <- c(2016)
peerLocs = c(65,70,73,96)
locsdf <- get_location_metadata(gbd_round_id = 5, location_set_id = 35)
locs = locsdf$location_id[locsdf$parent_id == 102 | locsdf$location_id %in% peerLocs | locsdf$level == 3 ]


locs_to_do <- setdiff(loc_list, locs)

for(yid in c(2016)){
  print(paste0("running for year:", yid))
  for(lid in locs_to_do){
    
    riskStandardize13("PAF", data_dir, yid, lid, mid, lsid, scale_ceiling)
  }
}

##########################################################
## END SESSION
quit("no")

##########################################################