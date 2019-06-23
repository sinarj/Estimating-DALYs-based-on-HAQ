##########################################################
# Author: Ryan Barber
# Date: 03 August 2016
# Description: Define common functions
##########################################################
## DEFINE LIBRARIES
library(data.table)
library(RMySQL)
require(foreign)
require(readstata13)
library(plyr)
library(parallel)

##########################################################
## DEFINE FUNCTIONS
qsub <- function(jobname, shell, code, logloc=NULL, project=NULL, hold=NULL, args=NULL, slots=1) { 
  # Set up number of slots
  if (slots > 1) { 
    slot_string = paste0("-pe multi_slot ", slots, " -l mem_free=", slots*2, "g")
  } 
  # Set up jobs to hold for 
  if (!is.null(hold)) { 
    hold_string <- paste0("-hold_jid \"", hold, "\"")
  } 
  # add log 
  if(!is.null(logloc)){
    log_string <- paste0('-o ',logloc,'/output -e ',logloc,'/errors')
  }
  # Set up project under which jobs are to be submitted
  if (!is.null(project)) { 
    project_string <- paste0("-P proj_", project)
  }  
  # Set up arguments to pass in 
  if (!is.null(args)) { 
    args_string <- paste(args, collapse = " ")
  }  
  # Construct the command 
  sub <- paste("qsub", 
               if (slots>1) slot_string, 
               if (!is.null(hold)) hold_string, 
               if (!is.null(logloc)) log_string, 
               if (!is.null(project)) project_string, 
               "-N ", jobname, 
               shell,
               code,
               if (!is.null(args)) args_string, 
               sep=" ")
  # Submit the command to the system
  if (Sys.info()[1] == "Linux") {
    system(sub) 
  } else {
    cat(paste("\n", sub, "\n\n "))
    flush.console()
  } 
}

# Load relevant GBD shared functions
source(paste0(j, "/temp/central_comp/libraries/current/r/get_draws.R"))
source(paste0(j, "/temp/central_comp/libraries/current/r/get_population.R"))

dbQuery <- function(host_name, query_string) {
  # Pass host name and query string to prevent typing connection info every time, return as data table
  con <- dbConnect(dbDriver("MySQL"), 
                   username="dbview", 
                   password="E3QNSLvQTRJm", 
                   host=paste0("modeling-", host_name, "-db.ihme.washington.edu"))
  results <- dbGetQuery(con, query_string)
  dbDisconnect(con)
  return(data.table(results))
}

getLocations <- function(lsid=35) {
  # Load all locations, along with limited metadata, for a given location set
  locsdf <- dbQuery("gbd", sprintf("
                                   SELECT
                                    location_id,
                                    location_name,
                                    ihme_loc_id,
                                    local_id,
                                    location_type,
                                    region_name,
                                    super_region_name,
                                    parent_id,
                                    path_to_top_parent,
                                    most_detailed,
                                    level,
                                    sort_order
                                   FROM
                                    shared.location_hierarchy_history
                                   WHERE
                                    location_set_version_id = shared.active_location_set_version(%s, 4)
                                   ", lsid))
  return(locsdf)
}

getCauses <- function(csid=3) {
  # Load all locations, along with limited metadata, for a given location set
  causesdf <- dbQuery("gbd", sprintf("
                                     SELECT
                                      cause_id,
                                      cause_name,
                                      path_to_top_parent,
                                      parent_id,
                                      most_detailed,
                                      level,
                                      sort_order
                                     FROM
                                      shared.cause_hierarchy_history
                                     WHERE
                                      cause_set_version_id = shared.active_cause_set_version(%s, 4)
                                     ", csid))
  return(causesdf)
}

getRisks <- function(reisid=1) {
  # Load all locations, along with limited metadata, for a given location set
  reidf <- dbQuery("gbd", sprintf("
                                     SELECT
                                       rei_id,
                                       rei_name,
                                       path_to_top_parent,
                                       most_detailed,
                                       level,
                                       sort_order
                                     FROM
                                       shared.rei_hierarchy_history
                                     WHERE
                                       rei_set_version_id = shared.active_rei_set_version(%s, 4)
                                     ", reisid))
  return(reidf)
}

getPops <- function(lsid, yids = c(seq(1990, 2010, 5), 2016), agids = c(2:20, 30:32, 235), sids = c(1, 2)) {
  # Loads populations for all locations in the specified location set, for a given year
  popsdf <- get_population(location_set_id = lsid, location_id = -1, year_id = yids, age_group_id = agids, sex_id = sids, gbd_round_id = 4)
  return(popsdf[, c("location_id", "year_id", "age_group_id", "sex_id", "population"), with = FALSE])
}

getCovs <- function(cns_call, lsid, yids = c(seq(1990, 2010, 5), 2016)) {  
  # Pull all specified years, rename as predictor for use in epi transition (and distinction from other estimated values)
  covsdf <- dbQuery("cod", sprintf("
                                 SELECT
                                  m.location_id, 
                                  m.year_id,
                                  m.age_group_id, 
                                  m.sex_id, 
                                  m.mean_value AS %s
                                 FROM 
                                  covariate.model m
                                 INNER JOIN
                                  covariate.model_version mv ON m.model_version_id = mv.model_version_id 
                                                             AND mv.model_version_id = %s
                                 INNER JOIN
                                  covariate.data_version dv ON mv.data_version_id = dv.data_version_id 
                                                            AND dv.status = 1
                                 INNER JOIN
                                  shared.covariate c ON dv.covariate_id = c.covariate_id 
                                                     AND c.covariate_name_short = '%s'
                                 INNER JOIN 
                                  (
                                    SELECT 
                                     *  
                                    FROM 
                                     shared.location_hierarchy_history 
                                    WHERE 
                                     location_set_version_id = shared.active_location_set_version(%s, 4)
                                   ) lhh ON m.location_id = lhh.location_id
                                WHERE
                                  m.year_id IN (%s)
                                 ", cns_call, mvid, cns, lsid, paste(yids, collapse=",")))
  return(covsdf)
}

getAgeGroups <- function(agids = c(2:20, 30:32, 235), gbdrid=4) {
  agesdf <- dbQuery("gbd", sprintf("
                                   SELECT 
                                    ag.age_group_id,
                                    ag.age_group_name,
                                    ag.age_group_years_start,
                                    ag.age_group_years_end,
                                    agw.age_group_weight_value
                                   FROM 
                                    shared.age_group ag
                                   INNER JOIN
                                    shared.age_group_weight agw ON ag.age_group_id = agw.age_group_id AND agw.gbd_round_id = %s
                                   WHERE 
                                    ag.age_group_id IN (%s)
                                   ", gbdrid, paste(agids, collapse = ",")))
  return(agesdf)
}

getAgeGroupNames <- function(agids) {
  agnmdf <- dbQuery("gbd", sprintf("
                                   SELECT
                                   age_group_id,
                                   age_group_name
                                   FROM
                                   shared.age_group
                                   WHERE
                                   age_group_id IN (%s)
                                   ", paste0(agids, collapse = ",")))
  return(agnmdf)
}

subsetDraws <- function(df, drawnum, indexvars) {
  df <- df[draw %in% sort(sample(seq(0,999), drawnum, replace = FALSE)), ]
  df <- df[, draw := seq.int(draw) - 1, by = indexvars]
  return(df)
}

loadDeaths <- function(lid, yid, mid, cc_vers, cids, drawnum) {
  # Read in CoDCorrect draws and format (takes BEST, daly_vers is not applied...)
  deathsdf <- get_draws(gbd_id_type = rep("cause_id", length(cids)), gbd_id = cids, 
                        source = "codcorrect", version_id = cc_vers,
                        location_id = lid, year_id = yid, age_group_id = c(2:20, 30:32, 235), sex_id = c(1, 2), 
                        measure_id = mid, metric_id = 1)
  deathsdf <- deathsdf[, metric_id := 1]
  deathsdf <- deathsdf[, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "metric_id", "cause_id", paste0("draw_", seq(0, 999))), with = FALSE]
  deathsdf <- melt(data = deathsdf,
                   id.vars = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "metric_id", "cause_id"),
                   variable.name = "draw",
                   variable.factor = TRUE,
                   value.name = "val")
  deathsdf$draw <- as.character(deathsdf$draw)
  deathsdf <- deathsdf[, draw := gsub("draw_", "", draw)]
  deathsdf$draw <- as.numeric(deathsdf$draw)
  if (drawnum < 1000) deathsdf <- subsetDraws(deathsdf, drawnum, indexvars = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "metric_id", "cause_id"))
  return(deathsdf)
}

loadPAF <- function(lid, yid, mid, paf_vers, cids, drawnum) {
  # Read in PAF .dta's and format
  pafdf <- data.table(read.dta13(paste0("/ihme/centralcomp/pafs/", paf_vers, "/", lid, "_", yid, ".dta")))
  if (mid %in% c(1, 4)) paftype <- "yll"
  else if (mid == 3) paftype <- "yld"
  pafdf <- pafdf[age_group_id %in% c(2:20, 30:32, 235) & sex_id %in% c(1,2) & cause_id %in% cids & rei_id == 169, c("location_id", "year_id", "age_group_id", "sex_id", "cause_id", paste0("paf_", paftype, "_", seq(0,999))), with=FALSE]
  pafdf <- melt(data = pafdf,
                id.vars = c("location_id", "year_id", "age_group_id", "sex_id", "cause_id"),
                variable.name = "draw",
                variable.factor = TRUE,
                value.name = "paf")
  pafdf$draw <- as.character(pafdf$draw)
  pafdf <- pafdf[, draw := gsub(paste0("paf_", paftype, "_"), "", draw)]
  pafdf$draw <- as.numeric(pafdf$draw)
  pafdf$measure_id <- mid
  if (drawnum < 1000) pafdf <- subsetDraws(pafdf, drawnum, indexvars = c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id"))
  return(pafdf)
}

prepLocHierarchy <- function(data_dir, lsid) {
  # Create data frame with parsed path to parent
  load(paste0(data_dir, "/locsdf_", lsid, ".RData"))
  parents <- read.table(text = locsdf$path_to_top_parent, 
                        sep = ",", 
                        fill = TRUE,
                        colClasses = "integer", 
                        col.names = c("root", paste0("L", seq(1, max(locsdf$level)))))
  locsdf <- cbind(locsdf[, c("location_id", "most_detailed", "level"), with=FALSE], parents)
  return(locsdf)
}

prepCauseHierarchy <- function(data_dir) {
  # Create data frame with parsed path to parent
  load(paste0(data_dir, "/causesdf.RData"))
  # Fixing the error in causedf.RData. Cannot save this to dir due to permissions issue.
  causesdf[82,]$path_to_top_parent = "294,295,962,366,995"
  parents <- read.table(text = causesdf$path_to_top_parent, 
                        sep = ",", 
                        fill = TRUE,
                        colClasses = "integer", 
                        col.names = c("root", paste0("L", seq(1, max(causesdf$level)))))
  causesdf <- cbind(causesdf[, c("cause_id", "most_detailed", "level"), with=FALSE], parents)
  return(causesdf)
}

hierarchyAgg <- function(levelvar, leveldf, idvar, aggvars, manual = FALSE) {
  # Aggregate by level (manual option for simple one- or two-step aggregations)
  aggleveldf <- leveldf[!is.na(leveldf[[levelvar]]),]
  aggleveldf[[idvar]] <- aggleveldf[[levelvar]]
  aggleveldf <- aggleveldf[, lapply(.SD, sum), by = .(location_id, year_id, age_group_id, sex_id, measure_id, cause_id, draw), .SDcols = aggvars]
  if (manual) {
    aggleveldf <- rbind(leveldf[, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw", aggvars), with = FALSE], 
                        aggleveldf[, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "draw", aggvars), with = FALSE])
  }
  return(aggleveldf)
}

compileRSM <- function(lids, yids = c(seq(1990, 2010, 5), 2016), agids = 27, sids = 3, mids = 1, data_dir) {
  # Load all specified risk-standardized mortality summaries
  ncores <- round(detectCores()/2, 0)
  args <- expand.grid(location_id = lids, year_id = yids)
  rsdf <- rbindlist(mcmapply(function(data_dir, yid, lid, agids, sids, mids) { 
                                                            load(paste0(data_dir, "/summaries/", yid, "/", lid, ".RData")) 
                                                            return(rsdeathsdf[age_group_id %in% agids & sex_id %in% sids & measure_id %in% mids, c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "cause_id", "rs_mval", "mval", "attrib_mval", "paf_mval"), with = FALSE])
                                                          }, 
                                      yid = args$year_id,
                                      lid = args$location_id,
                                      MoreArgs = list(data_dir, agids, sids, mids),
                                      SIMPLIFY = FALSE,
                                      mc.cores = ncores))
  return(rsdf)
}

##########################################################