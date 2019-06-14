# Estimating-DALYs-based-on-HAQ

### Credits: 
This reseach project was completed under the guidance and support of my research advisor Prof. Marcia Weaver.

The project builds on top of the code of my IHME (Institute for Health Metrics and Evaluation) colleagues: Vishnu Nandakumar, Kate  Rosettie and Jamal Yearwood. 

## Overview: To estimate the burden of disease averted by personal health care in the United States, where burden is measured by disability-adjusted-life-years (DALYs), and personal health care is measured by the healthcare quality and access (HAQ) index

What follows next is the list of steps (with references to the code) that were followed to achieve this.

## Exploratory Data Analysis (in the folder “pre_analysis”)

•	 “DALY_regression_with_Visualizations”: This notebook runs the regressions using the Data Rich locations and produces a number of interesting visualizations

•	“DALY_regression_with_Highest_Expenditure_data”: Load and analyze the data for the causes with highest expenditure in the US (Joe's list)

•	“DALY_regression_locations_Joe_CODEM”: Loading the locations based on the list from Joe's team intersected with CODEM


## Steps in the data generation/retrieval process:

### Creation of the input RData files is done in the script “all_input_data_files.R” in the folder “generate_input_RData_files”.
This is the script that prepares all the input RData files to be used for all causes, all locations, all years, and for both age-standardized and age-specific.

### Using the RData files to create the Risk Standardized ratios (in csv files) for all amenable causes separately:

The script named “risk_standardizer_ageStandardized.R” in the folder “risk_standardizer”:
Risk standardization values for all the amenable causes for the given years and locations. Locations for which was run was all the most detailed locations (for use in the regressions) and the countries (for use in the min-max band calculations). This script was run for the age-std version

The script named “risk_standardizer_ageSpecific.R” in the folder “risk_standardizer”:
Very similar to the previously described file but for age-group specific values.

### Cancer-MI-Ratios: for cancer (total of 8 causes) the MI ratios were created separately in the following scripts, both in the folder “cancer_mi_ratios”:

“cancer_mi_ratios_ageSpecific.R”:  Generating the cancer MI ratios for age group specific
data
“cancer_mi_ratios_ageStandardized.R”: Generating the cancer MI ratios for age standardized data

### Using the country specific csv files for cancer as well as non-cancer to create the min-max bands
For age-specific: notebook titled “Global_minmax_ageSpecific” in folder “create_min_max”
For age-standardized: notebook titled “Global_minmax_ageStandardized” in folder “create_min_max”

### Create csv files that has all the input HAQ data for all of the most detailed locations. The scripts are in the folder “create_haq_input_allLocs”:

“HAQ_sourceData_allMostDetLocs_ageSpecific”: Notebook to create the HAQ input data for all the most detailed locations for age specific version

“HAQ_sourceData_allMostDetLocs_ageStandardized”: Notebook to create the HAQ input data for all the most detailed locations for age standardized version

“HAQ_sourceData_onlyUS_ageSpecific”: Notebook to create the HAQ input data for only US for the age specific version. This is being done separately because US is the to be used for the final predictions but US itself is not a "most detailed location"

“HAQ_sourceData_USonly_ageStandardized”: # Notebook to create the HAQ input data for only US for the age standardized version. This is being done separately because US is the to be used for the final predictions but US itself is not a "most detailed location"


### Cause Weights Preprocessing files (in the folder “preprocess_cause_weights”):

“PreWeights_processing_ageSpecific”: This notebook takes as input the excel worksheet which has all the new cause weights and pre-processes it in a way that can then be directly consumed by the notebook "Weights_processing_ageSpecific" which will then prepare it further for the HAQ-index-generation code. This notebook deals with the age specific version.

“Weights_processing_ageSpecific”: This notebook takes as input the pre-processed data prepared via the notebook "PreWeightsProcessing_ageSpecific" and further processes it in a way that can then be directly consumed by the HAQ-index-generation code. This notebook deals with the age specific version.

“Weights_processing_ageStandardized”: This notebook takes as input the excel worksheet which has all the new cause weights and pre-processes it in a way that can then be directly consumed by the HAQ-index-generation code. This notebook deals with the age standardized version.

### Using the min max bands as well as csv files for cancer and non-cancer from all most detailed locations, and the processed weights, to create the HAQ indices (see scripts in the folder create_haq_index)

There are different versions for both age-group specific and age standardized; also there are separate files for the US, as initially the regressions only depended on the data from the most detailed locations. The US is not a most detailed location, but the final predictions were to be created for the US.

### Generate data files: these combine the HAQ and SDI values which is then used as the input into the HAQ-Frontier analysis. All the files are present in the folder “haq_sdi_frontiers” and are as follows:

“generate_data_for_ageSpecific_causeSpecific_allMostDetLocs_HAQ_frontier”: This generates the HAQ-SDI data to be used for the HAQ-Frontier Analysis. Current notebook is for all most detailed locations and age specific data.

“generate_data_for_ageSpecific_causeSpecific_USonly_HAQ_frontier”: This generates the HAQ-SDI data to be used for the HAQ-Frontier Analysis. Current notebook is for all only the US and age specific data.

“generate_data_for_ageStandardized_causeSpecific_allLocsIncUS_HAQ_frontier”: This generates the HAQ-SDI data to be used for the HAQ-Frontier Analysis. Current notebook is for all most detailed locations and the US and age standardized data.

“generate_data_for_ageStandardized_causeSpecific_allMostDetLocs_HAQ_frontier”: This generates the HAQ-SDI data to be used for the HAQ-Frontier Analysis. Current notebook is for all most detailed locations and age standardized data.


### Generating the HAQ over frontier values for both age group specific as well as age-standardized versions (see script “haq_sdi_frontier_analysis.R”):
This code generates the HAQ-SDI frontiers for both the age-group specific, as well as the age-standardized versions for the aggregated causes (cause_id = 100) and the list of all most detailed locations, and the United States, which is the location we use for the final results


### Using the above along with the logit_SDI values to run regressions for every group of cause-ageGroup-sex

The notebook titled “DALYs_regression_HAQ_LogitSDI_ageSpecific” in the folder “regression_analysis” is the main file for for running all the regressions using OLS for every group of cause-ageGroup-sex. Based on the previous analysis, we use the age groups up to 74 years and logit SDI as the predictor rather than SDI, along with the HAQ indices obtained for GBD round 5

Similar to the one above, is the notebook titled “DALYs_regression_HAQ_LogitSDI_ageStandardized” in the same folder which runs for all combinations of cause and sex, since the age groups are now standardized.

## Final Tables:

Based on the regressions the final tables are prepared using the 2016 US data for HAQ over frontier and the US population (which is age group specific or standardized based on whichever of the two versions we’re looking at). The files are “PrepareFinalTables_usingUS_HAQoverFrontier_ageSpecific.ipynb" and “PrepareFinalTables_usingUS_HAQoverFrontier_ageStandardized.ipynb” and they’re both in the folder “regression_analysis”

## Regression comparison:

•	The notebook titled “Compare_SDI_and_LogitSDI_DALYs_regression” in the folder “regression_analysis” was written to compare how our regression results fared when we used SDI as compared to logit_SDI as one of the predictors in our regressions

## More from the creation of data/sheets for NPC:

•	The notebook titled “candidate_causes_criteria_list_fillUp” in the folder “creating_worksheets_NPC/” was created to code the completion of the CandidateCauses_Criteria sheet, which we were to send across to NPC

### Note that for all the R scripts, you’d need some of the functions from a utilities script, called “utilities_sinarj.R”. Make sure to update the path from where you save this script on your machine.
