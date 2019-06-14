##########################################################
# Script that prepares all the input RData files to be
# used for all causes, all locations, all years, and for
# both age-standardized and age-specific
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())
library(argparse)
library(data.table)
library(plyr)
if (Sys.info()[1] == "Linux"){
  j <- "/home/j/"
  h <- paste0("/homes/", Sys.info()[7], "/haq")
}else if (Sys.info()[1] == "Windows"){
  j <- "J:/"
  h <- "H:/"
}else if (Sys.info()[1] == "Darwin"){
  j <- "/Volumes/snfs/"
  h <- paste0("/Volumes/", Sys.info()[6], "/")
}

if (Sys.info()[1] == "Linux"){
  j <- "/home/j/"
  h <- paste0("/homes/", Sys.info()[8], "/haq/haq")
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
                    default=paste0(h, "/risk_standardized_amen_mort"), type="character")
parser$add_argument("--data_dir", help="Space where results are to be saved",
                    default="/share/scratch/projects/hssa/haq/HAQ_2017/vanilla_cc_90_paf_242_amenable/", type="character")
parser$add_argument("--cc_vers", help="CoDCorrect version number",
                    default=90, type="integer")
parser$add_argument("--como_vers", help="COMO version number",
                    default=347, type="integer")
parser$add_argument("--gbd_rid", help="GBD round ID",
                    default=5, type="integer")
parser$add_argument("--paf_vers", help="PAF version (needs to be character)",
                    default="242_amenable", type="character")
parser$add_argument("--lsid", help="Location set",
                    default=35, type="integer")
parser$add_argument("--lid", help="Location",
                    default=102, type="integer")
parser$add_argument("--yid", help="Year for current job",
                    default=1990, type="integer")
parser$add_argument("--mid", help="Measure",
                    default=1, type="integer")
parser$add_argument("--drawnum", help="Number of draws",
                    default=1000, type="integer")
args <- parser$parse_args()
list2env(args, .GlobalEnv)
rm(args)

##########################################################
## DEFINE FUNCTIONS
source(paste0("~/haq/haq/risk_standardized_amen_mort/utilities_sinarj.R"))

dataGrabPAF <- function(lid, yid, mid, cids, cc_vers, drawnum) {
  # Collect PAF and death data for given location-year, then combine and return
  pafdf <- loadPAF(lid, yid, mid, cids, drawnum)
  deathsdf <- loadDeaths(lid, yid, mid, cc_vers, cids, drawnum)
  # incidencedf <- loadIncidence(lid, yid, mid = 6, como_vers, cids, drawnum)
  df <- merge(pafdf,
              deathsdf,
              by = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw"),
              all = TRUE)
  # df <- merge(df,
  #             incidencedf,
  #             by = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw"),
  #             all = TRUE)
  df <- df[is.na(val), val := 0][is.na(paf), paf := 0]
  return(df)
}

saveData <- function(inputdf, data_dir, lid, yid, mid, prefix) {
  # Save
  save(inputdf,
       file = paste0(data_dir, "/draws/inputs/", yid, "/", prefix, "_", mid, "_", lid, ".RData"))
}

popScalar <- function(yids, lid, locsdf, popsdf) {
  childpopsdf <- popsdf[location_id %in% locsdf$location_id[locsdf$parent_id == lid & locsdf$location_id != lid] & year_id %in% yids,]
  childpopsdf <- childpopsdf[, .(child_pop_agg = sum(population)), by = .(year_id, age_group_id, sex_id)]
  aggpopsdf <- popsdf[location_id == lid & year_id %in% yids,]
  scaledf <- merge(aggpopsdf,
                   childpopsdf,
                   by = c("year_id", "age_group_id", "sex_id"))
  scaledf <- scaledf[, pop_scalar := population / child_pop_agg][pop_scalar < 1, pop_scalar := 1]
  return(scaledf[, c("location_id", "year_id", "age_group_id", "sex_id", "pop_scalar"), with = FALSE])
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
mostDetailedCollector <- function(lid, yid, mid, cc_vers, data_dir) {
  # Load cause data (only keep amenable causes or their children)
  load(paste0(data_dir, "/causesdf.RData"))
  cause_list <- fread(paste0(prog_dir, "/amenable_cause_list_GBD.csv"))
  cids <- unique(cause_list$cause_id)
  causesdf <- causesdf[most_detailed == 1 | cause_id %in% c(587,297)] #TB and Diabetes mellitus are not most detailed cause
  causesdf <- causesdf[!path_to_top_parent %like% ",297,"][!path_to_top_parent %like% ",587,"]
  causesdf <- causesdf[, amenable := FALSE]
  for (cid in cause_list$cause_id) {
    causesdf <- causesdf[path_to_top_parent %like% cid, amenable := TRUE]
  }
  causesdf <- causesdf[amenable == TRUE]
  
  # PAF
  inputdfPAF <- dataGrabPAF(lid, yid, mid, cids, cc_vers, drawnum)
  
  # Save PAFs
  saveData(inputdfPAF[, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw", "paf", "val"), with = FALSE],
           data_dir, lid, yid, mid, "PAF")
  
  # saveData(inputdfPAF[, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw", "paf", "incidence", "incidence_rate", "val"), with = FALSE], 
  #          data_dir, lid, yid, mid, "PAF")
  # 
  # Summarize PAFs
  inputdfPAF <- inputdfPAF[, .(paf = mean(paf)), by = .(location_id, year_id, age_group_id, sex_id, measure_id, cause_id)]
  saveData(inputdfPAF, data_dir, lid, yid, mid, "md_sums/PAF")
  
}

aggsCollector <- function(lid, yid, mid, data_dir, locsdf, popsdf) {
  scaledf <- popScalar(yid, lid, locsdf, popsdf)
  ## ## ## ## ## ##
  # PAF
  inputdfPAF <- rbindlist(lapply(as.list(locsdf$location_id[locsdf$parent_id == lid & locsdf$location_id != lid]), function(lid, data_dir, yid, mid) 
  {load(paste0(data_dir, "/draws/inputs/", yid, "/PAF_", mid, "_", lid, ".RData"))
    return(inputdf)}, 
  data_dir, yid, mid))
  inputdfPAF <- inputdfPAF[, location_id := lid][, .(paf = sum(paf * val) / sum(val), val = sum(val)), by = .(location_id, year_id, age_group_id, sex_id, measure_id, cause_id, draw)]
  # inputdfPAF <- inputdfPAF[, location_id := lid][, .(paf = sum(paf * val) / sum(val), val = sum(val), incidence = sum(incidence)), by = .(location_id, year_id, age_group_id, sex_id, measure_id, cause_id, draw)]
  inputdfPAF <- inputdfPAF[val == 0, paf := 0]
  if (lsid == 35) {
    inputdfPAF <- merge(inputdfPAF,
                        scaledf,
                        by = c("location_id", "year_id", "age_group_id", "sex_id"))
    inputdfPAF <- inputdfPAF[, val := val * pop_scalar]
  }
  
  # Save
  saveData(inputdfPAF[, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw", "paf", "val"), with = FALSE], 
           data_dir, lid, yid, mid, "PAF")
  
}

##########################################################
## RUN PROGRAM
load(paste0(data_dir, "/locsdf_", lsid, ".RData"))

loc_list = locsdf[locsdf$most_detailed == 1]$location_id

getPops <- function(lsid, yids = c(2016), agids = c(2:20, 30:32, 235), sids = c(3)) {
  # Loads populations for all locations in the specified location set, for a given year
  popsdf <- get_population(location_set_id = lsid, location_id = locsdf$location_id, year_id = yids, age_group_id = agids, sex_id = sids, gbd_round_id = 5)
  return(popsdf[, c("location_id", "year_id", "age_group_id", "sex_id", "population"), with = FALSE])
}

load(paste0("/share/scratch/projects/hssa/haq/HAQ_2017/", "/haq_USpopsdf_35.RData"))
load(paste0("/share/scratch/projects/hssa/haq/HAQ_2017/vanilla_cc_90_paf_242_amenable/", "/popsdf_", lsid, ".RData"))

locs = locsdf$location_id[locsdf$level == 3 ]

yids = c(2015,2016)

popsdf <- getPops(lsid, sids = c(1,2,3))
locs_to_do <- setdiff(locs, loc_list)
yids = c(2016)
#temp run for only one loc
locs_to_do <- c(62)
for (yid in yids) 
{
  for (lid in locs_to_do) 
  {
    print(paste0("Running for location:", lid))
    if (locsdf$most_detailed[locsdf$location_id == lid] == 1) 
    { 
      mostDetailedCollector(lid, yid, mid, cc_vers, data_dir)
    } 
    else 
    {
      ## If aggregate, find most-detailed locations that have already been saved and aggregate them 
      aggsCollector(lid, yid, mid, data_dir, locsdf, popsdf)
    }
  }  
}

##########################################################
## END SESSION
quit("no")

##########################################################