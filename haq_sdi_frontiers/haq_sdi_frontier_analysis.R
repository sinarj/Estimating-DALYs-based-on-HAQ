# This code generates the HAQ-SDI frontiers for both the age-group specific, as 
# well as the age-standardized versions for the aggregated causes (cause_id = 100)
# and the list of all most detailed locations, and the United States, which is the
# location we use for the final results

library(data.table)
library(micEcon, lib.loc = '/ihme/homes/arjuns13/r_packages')
library(lmtest, lib.loc = '/ihme/homes/arjuns13/r_packages')
library(frontier, lib.loc = '/ihme/homes/arjuns13/r_packages')

# read in data for age-group specific
#mydat <- fread('/ihme/homes/arjuns13/notebooks/Documents/Data/haq_sdi_perAge_aggregatedCauses_allLocsIncUS.csv')
# read in data for age standardized
mydat <- fread('/ihme/homes/arjuns13/notebooks/Documents/Data/haq_sdi_ageStd_aggregatedCauses_allLocsIncUS.csv')


require(data.table)
require(ggplot2)
require(frontier)

sfa2_log <- sfa(data = mydat, formula = ln_haq ~ 1 + logit_sdi)
summary(sfa2_log, truncNorm=F)

fit_the_stuff <- function(age) {
  subset <- copy(mydat[age_group_id == age])
  sfa2_log <- sfa(data = subset, formula = ln_haq ~ 1 + logit_sdi)
  subset[, fit2:= fitted(sfa2_log, asInData = T)]
  return(subset)
}

age_group_ids = unique(mydat[, age_group_id])
# mydat[, fit2:= fitted(sfa2_log, asInData = T)]

new_mydat <- lapply(X = age_group_ids, FUN = fit_the_stuff)
new_mydat <- rbindlist(new_mydat)

#write.csv(new_mydat, file="/ihme/homes/arjuns13/notebooks/Documents/Data/haq_sdi_frontier_estimates_perAge_aggregatedCauses_allLocsIncUS.csv")
write.csv(new_mydat, file="/ihme/homes/arjuns13/notebooks/Documents/Data/haq_sdi_frontier_estimates_ageStd_aggregatedCauses_allLocsIncUS.csv")

plot2 <- ggplot(new_mydat) + 
  geom_point(aes(x = logit_sdi, y = ln_haq, group = location_id, color = location_id), alpha=.5) +
  geom_line(aes(x = logit_sdi, y = fit2)) +
  theme(legend.position = 'none') +
  xlab("Logit SDI") + ylab("Log HAQ") +
  ggtitle("SFA of Log HAQ against logit SDI, with a cross section structure")

plot2

plot2_lev <- ggplot(new_mydat) + 
  geom_point(aes(x = logit_sdi, y = 10**ln_haq, group = location_id, color = location_id), alpha=.5) +
  geom_line(aes(x = logit_sdi, y = 10**(fit2))) +
  theme(legend.position = 'none') +
  xlab("Logit SDI") + ylab("HAQ") +
  ggtitle("SFA of HAQ against logit SDI")

plot2_lev

