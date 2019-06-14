##########################################################
# Generating the cancer MI ratios for age group specific
# data
##########################################################
## DEFINE ROOT AND LIBRARIES
rm(list=ls())
library(argparse, lib.loc= '/home/j/temp/vishnn/packages')
library(findpython, lib.loc= '/home/j/temp/vishnn/packages')
library(data.table)
library(plyr)
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
                    default="/share/scratch/projects/hssa/haq/HAQ_2017/vanilla_cc_90_paf_242_amenable", type="character")
parser$add_argument("--save_dir", help="Space where results are to be saved",
                    default="/share/scratch/projects/hssa/haq/HAQ_2017/haq_2017_latest_cancer", type="character")

parser$add_argument("--lsid", help="Location set",
                    default=35, type="integer")
parser$add_argument("--lid", help="Location",
                    default=70, type="integer") #523
parser$add_argument("--yid", help="Year for current job",
                    default=1990, type="integer") #2005
parser$add_argument("--mid", help="Measure",
                    default=1, type="integer")
args <- parser$parse_args()
list2env(args, .GlobalEnv)
rm(args)

##########################################################
## DEFINE FUNCTIONS
source(paste0(prog_dir, "/utilities.R"))
source(paste0(j,"temp/central_comp/libraries/2017_archive/r/get_draws.R"))

library(feather)

## This is the original data source, which you'd only need if you don't find the data already stored
## in the temp RData files I created locally below
# start_time = Sys.time()
# new_leukemia <- read_feather("/home/j/temp/emumford/leukemia_draws.feather")
# end_time = Sys.time()
# end_time - start_time

#new_leukemia_feather <- setDT(read_feather("/home/j/temp/emumford/leukemia_draws.feather"))
#save(new_leukemia_feather, file="/ihme/homes/arjuns13/haq/haq/data_files/leukemia_draws_age_std_2016.RData")

# new_leukemia_feather <- new_leukemia_feather[age_group_id %in% c(2:20, 30:32, 235) & year_id %in% c(seq(1990, 2015,5), 2016, 2017) & sex_id %in% c(1,2)]
# new_leukemia_feather <- new_leukemia_feather[, c("modelable_entity_id", "location_set_id") := NULL]
# save(new_leukemia_feather, file="/ihme/homes/arjuns13/haq/haq/data_files/leukemia_draws.RData")
# write.csv2(new_leukemia_feather, "/ihme/homes/arjuns13/haq/haq/data_files/leukemia_draws.csv", row.names = F)

load("/ihme/homes/arjuns13/haq/haq/data_files/leukemia_draws.RData")
new_leukemia_feather <- new_leukemia_feather[year_id %in% c(2016) & sex_id %in% c(1,2)]


#you only need to do this once, and here we do for both sexes 1 and 2, cause we merge on that and then aggregate at the end
popsdf <- getPops(lsid, sids = c(1,2,3))

# Pull raw counts of deaths at draw level, for each cause
causes <- c(849, 429, 432, 435, 441, 468, 484, 487)

load(paste0(data_dir, "/locsdf_", lsid, ".RData"))

loc_list = locsdf[locsdf$most_detailed == 1]$location_id

makeMIratiosNew <- function(data_dir, yid, lid) {
  
  deathsdf<-get_draws(gbd_id_type = rep("cause_id",length(causes)), gbd_id = causes, source = 'codcorrect', measure_id = 1,
                      location_id = lid, year_id = yid, sex_id = c(1,2), gbd_round_id = 5, metric_id = 1,  version_id = 90)
  
  # Make it long 
  deathsdf<-melt(deathsdf, id.vars = c('location_id','year_id','sex_id','age_group_id','cause_id','measure_id','metric_id'), 
                 value.name = 'deaths', 
                 variable.name = 'draw',
                 variable.factor = TRUE)
  deathsdf$draw<-gsub("draw_","", deathsdf$draw)
  deathsdf$draw <- as.numeric(deathsdf$draw)
  
  #Pull incidence 
  incidencedf<-get_draws(gbd_id_type = rep("cause_id",length(causes)), gbd_id = causes, source = 'como', measure_id = 6,
                         location_id = lid, year_id = yid, sex_id = c(1,2), gbd_round_id = 5, metric_id =3)
  
  #Make it long
  incidencedf<-melt(incidencedf, id.vars = c('location_id','year_id','sex_id','age_group_id','cause_id', 'measure_id','metric_id'), 
                    value.name = 'incidence', 
                    variable.name = 'draw',
                    variable.factor = TRUE)
  
  incidencedf$draw<-gsub("draw_","", incidencedf$draw)
  incidencedf$draw <- as.numeric(incidencedf$draw)
  
  # Read in Leukemia incidence estimates
  new_leukemia_t <- new_leukemia_feather[year_id == yid & location_id == lid]
  new_leukemia_t <- melt(data = new_leukemia_t,
                         id.vars = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id"),
                         variable.name = "draw",
                         variable.factor = TRUE,
                         value.name = "val")
  
  new_leukemia_t <- new_leukemia_t[, val := as.numeric(gsub(",", ".", gsub("\\.", "", new_leukemia_t$val)))]
  new_leukemia_t <- new_leukemia_t[, draw := as.numeric(gsub("draw_", "", draw))]
  
  
  setnames(new_leukemia_t, "val", "leukemia_val")
  
  incidencedf <- merge(incidencedf[,c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw", "incidence")],
                       new_leukemia_t, all =  T, by = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw"))
  
  # Replace value with corrected ones, drop column
  
  incidencedf[cause_id == 487, incidence := leukemia_val]
  incidencedf[, "leukemia_val" := NULL]
  
  #Multiply incidence by population to get counts instead of incident rates
  incidencedf <- merge(incidencedf, popsdf, by = c('location_id', 'age_group_id', 'sex_id', 'year_id'))
  incidencedf <- incidencedf[, incidence_count := incidence * population, by = .(location_id, age_group_id, sex_id, draw, cause_id)]
  setnames(incidencedf, "incidence", "incidence_rate")
  
  #Merge dataframes
  incidencedf <- incidencedf[,list(location_id, year_id, sex_id, age_group_id, cause_id, draw, incidence_count)]
  deathsdf <- deathsdf[,list(location_id, year_id, sex_id, age_group_id, cause_id, draw, deaths)]
  inputdf <- merge(incidencedf, deathsdf, by = c("location_id", "year_id", "sex_id", "age_group_id", "cause_id", "draw"))
  inputdf <- inputdf[,deaths:=as.numeric(deaths)]
  #For age groups that there are values for 
  
  # Only keep Nolte/McKee causes and ages
  cause_list <- fread(paste0(prog_dir, "/amenable_cause_list_GBD.csv"))
  causesdf <- prepCauseHierarchy(data_dir)
  
  # Aggregate to Nolte/McKee list, attach to main data frame
  inputdf <- merge(inputdf,
                   cause_list[, c("cause_id", "age_group_id_start", "age_group_id_end"), with = FALSE],
                   by = "cause_id")
  inputdf <- inputdf[age_group_id >= age_group_id_start & age_group_id <= age_group_id_end,]
  
  # Collapse to both sexes
  inputdf<-inputdf[,list(deaths = base::sum(deaths), incidence = base::sum(incidence_count)), by = c('location_id','year_id','age_group_id','cause_id','draw')]
  
  # Aggregate sexes
  inputdf <- inputdf[, sex_id := 3]
  inputdf <- merge(inputdf,
                   popsdf,
                   by = c("location_id", "year_id", "age_group_id", "sex_id"), all.x = T)
  
  inputdf <- inputdf[, .(deaths = (deaths / population),
                         incidence = (incidence / population)), 
                     by = .(location_id, year_id, age_group_id, sex_id, cause_id, draw)]
  
  ## Compute aggregated, all ages both sex MI ratio
  inputdf[,rsval := deaths/incidence]
  inputdf[,val := NA] #Needed to keep shape of dfs uniform
  inputdf$draw<-as.numeric(inputdf$draw)
  inputdf<-inputdf[order(cause_id,draw)]
  
  if(nrow(inputdf)<1 | length(unique(inputdf$rsval)) < 1)
  {
    print(paste0("Did not work for location ID:", lid))
  }
  else
  {
    print(paste0("Done for location ID:", lid))
  }
  # Return data frame
  write.csv(inputdf[, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", "draw", "rsval", "val"), with = FALSE],
            file = paste0("/share/scratch/projects/hssa/haq/HAQ_2017/haq_2017_latest_cancer", "/cancers_age_std_all_most_det_locs/", yid, "/", lid, ".csv"),
            row.names = FALSE)
}

##########################################################
## RUN PROGRAM

source("/home/j/temp/central_comp/libraries/current/r/get_location_metadata.R")
yids <- c(2016)

peerLocs = c(65,70,73,96)
locsdf <- get_location_metadata(gbd_round_id = 5, location_set_id = 35)
locs = locsdf$location_id[locsdf$parent_id == 102 | locsdf$location_id %in% peerLocs | locsdf$level == 3 ]

locs_to_do <- setdiff(loc_list, locs)


loc_list = locsdf[locsdf$most_detailed == 1]$location_id
locs = locsdf$location_id[locsdf$level == 3 ]

locs_to_do <- setdiff(locs, loc_list)

start_time = Sys.time()
# Risk-standardized mortality calculator
for (yid in yids) {
  for (lid in locs_to_do) {
    #if (!file.exists(paste0(save_dir, "/cancers/", yid, "/", lid, ".csv")))
    makeMIratiosNew(data_dir, yid, lid)
  }
}
end_time = Sys.time()
end_time - start_time
