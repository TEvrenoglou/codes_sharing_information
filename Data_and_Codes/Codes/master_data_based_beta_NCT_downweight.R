# load the required packages
library(readxl)   
library(rjags) 
library(R2jags)
library(tidyverse)
library(curl)
library(dplyr)
library(devtools)
library(netmeta)
library(stringr)
library(remotes)

# install nmajags from github
#remotes::install_github("guido-s/nmajags")
library(nmajags)

# load functions for the Codes file

source("./Codes/helpers.R") ## read the helpers file

source("./Codes/pairwiseMA_data_expert_based_priors.R") ## pairwise meta-analysis model to construct the priors for the location parameter using the data based approach

source("./Codes/NMA_model_data_based_beta.R") ## NMA model for the first stage

source("./Codes/NMA_model_with_informative_priors.R") ## NMA model for the second stage

'%!in%' = function(x, y)  ! ('%in%'(x, y))
 
# set seed for reproducibility
set.seed(1995)

# read the data
all.data =
  read_excel("./Data/data.xlsx",
             na = c("", " ")) 

# create different datasets for the two populations
ca.data = subset(all.data, all.data$ChildAdolesc == 1)

gp.data = subset(all.data, all.data$GeneralPopulation == 1)

#### Find the common basic comparisons between the two networks. Those are all common pairwise comparisons "Drug vs Placebo"

### bring the data for GP in a wide format to enable locating the pairwise comparisons in the network
pair_gp = pairwise(
  data = gp.data,
  mean = OverallEfficM,
  n = OverallEfficN,
  sd = OverallEfficSD,
  sm = "SMD",
  treat = Drug,
  studlab = Study_No
)

pair_gp$comparison = paste(pair_gp$Drug1, "-", pair_gp$Drug2, sep = "")

### bring the data for CA in a wide format to enable locating the pairwise comparisons in the network
pair_ca = pairwise(
  data = ca.data,
  mean = OverallEfficM,
  n = OverallEfficN,
  sd = OverallEfficSD,
  sm = "SMD",
  treat = Drug,
  studlab = Study_No
)

pair_ca$comparison = paste(pair_ca$Drug1, "-", pair_ca$Drug2, sep = "")

### Find all basic comparisons in the network of GP
A = pair_gp$comparison[grepl("-Placebo", pair_gp$comparison)]

A = unique(A)

### Find all basic comparisons in the network of CA
B = pair_ca$comparison[grepl("-Placebo", pair_ca$comparison)]

B = unique(B)

### identify the common basic comparisons between the network of GP and CA
common = intersect(A, B)

# BEGIN CONSTRUCTING PRIORS FOR THE LOCATION PARAMETER

### prepare lists to store the results
dat.gp = list()

dat.ca = list()

model = list()

res = list()

shifts = list()

#### Run the meta-analysis model and construct the priors for the location parameter

for (i in 1:length(common)) {
  dat.gp[[i]] = pair_gp %>%
    filter(comparison == common[i]) ### restrict GP to contain only the comparison i
  
  dat.ca[[i]] = pair_ca %>%
    filter(comparison == common[i]) ### restrict CA to contain only the comparison i
  
  N1 = length(unique(dat.gp[[i]]$studlab))
  
  N2 = length(unique(dat.ca[[i]]$studlab))
  
  ### prepare data in a form suitable for JAGS
  data_list = list(
    y1 = dat.gp[[i]]$TE,
    v1 = dat.gp[[i]]$seTE * dat.gp[[i]]$seTE,
    N1 = N1,
    N2 = N2,
    y2 = dat.ca[[i]]$TE,
    v2 = dat.ca[[i]]$seTE * dat.ca[[i]]$seTE
  )
  
  ## run the model for comparison i
  model[[i]] = jags(
    data = data_list,
    parameters.to.save = save_ma_data,
    model.file = modfile_meta,
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 10000
  )
  
  #### Store the results
  res[[i]] = as.data.frame(model[[i]]$BUGSoutput$summary)
  
  shifts[[i]] = res[[i]][1, ]
  
  shifts[[i]]$comparison = common[i]
  
  shifts[[i]]$prec = 1 / (shifts[[i]]$sd ^ 2)
}

shift_data = bind_rows(shifts, .id = "column_label")

row.names(shift_data) = NULL

shift_data[, 1] = NULL

shift_data = shift_data %>%
  arrange(comparison)

beta.results = data.frame(
  comparison = c(
    "Aripiprazole vs Placebo",
    "Asenapine vs Placebo",
    "Clozapine vs Placebo",
    "Fluphenazine vs Placebo",
    "Haloperidol vs Placebo",
    "Loxapine vs Placebo",
    "Lurasidone vs Placebo",
    "Molindone vs Placebo",
    "Olanzapine vs Placebo",
    "Paliperidone vs Placebo",
    "Placebo vs Placebo",
    "Quetiapine vs Placebo",
    "Risperidone vs Placebo",
    "Trifluoperazine vs Placebo",
    "Ziprasidone vs Placebo"
  ),
  mean = rep(NA, 15),
  sd = rep(NA, 15)
)

# dataset with all the priors constructed for each basic comparison
# an noninfromative prior is assigned if a comparison exists for GP but not for CA

beta.results[c(1, 2, 7, 9, 10, 12, 13, 15), 2:3] = shift_data[, 1:2]

beta.results$mean = ifelse(is.na(beta.results$mean), 0, beta.results$mean)

beta.results$sd = ifelse(is.na(beta.results$sd), 100, beta.results$sd)

drugs = intersect(unique(ca.data$Drug), unique(gp.data$Drug))

drugs = sort(drugs)

row.names(beta.results) = drugs

# create a dataframe which contains the constructed informative priors for beta and noninformative priors for the non-common comparisons
st_gp = sort(unique(gp.data$Drug))

beta.resultsf = matrix(NA, ncol = 2, nrow = length(st_gp))

beta.resultsf = as.data.frame(beta.resultsf)

names(beta.resultsf) = c("mean", "sd")

row.names(beta.resultsf) = st_gp

### Locate each treatment's row at the two dataframes
Lox_betas = which(row.names(beta.results) == "Loxapine")
Lox_betaf = which(row.names(beta.resultsf) == "Loxapine")

Trif_betas = which(row.names(beta.results) == "Trifluoperazine")
Trif_betaf = which(row.names(beta.resultsf) == "Trifluoperazine")

Flup_betas = which(row.names(beta.results) == "Fluphenazine")
Flup_betaf = which(row.names(beta.resultsf) == "Fluphenazine")

Halop_betas = which(row.names(beta.results) == "Haloperidol")
Halop_betaf = which(row.names(beta.resultsf) == "Haloperidol")

Arip_betas = which(row.names(beta.results) == "Aripiprazole")
Arip_betaf = which(row.names(beta.resultsf) == "Aripiprazole")

Plac_betas = which(row.names(beta.results) == "Placebo")
Plac_betaf = which(row.names(beta.resultsf) == "Placebo")

Quet_betas = which(row.names(beta.results) == "Quetiapine")
Quet_betaf = which(row.names(beta.resultsf) == "Quetiapine")

Zip_betas = which(row.names(beta.results) == "Ziprasidone")
Zip_betaf = which(row.names(beta.resultsf) == "Ziprasidone")

Ase_betas = which(row.names(beta.results) == "Asenapine")
Ase_betaf = which(row.names(beta.resultsf) == "Asenapine")

Lur_betas = which(row.names(beta.results) == "Lurasidone")
Lur_betaf = which(row.names(beta.resultsf) == "Lurasidone")

Ris_betas = which(row.names(beta.results) == "Risperidone")
Ris_betaf = which(row.names(beta.resultsf) == "Risperidone")

Ola_betas = which(row.names(beta.results) == "Olanzapine")
Ola_betaf = which(row.names(beta.resultsf) == "Olanzapine")

Clo_betas = which(row.names(beta.results) == "Clozapine")
Clo_betaf = which(row.names(beta.resultsf) == "Clozapine")

Pal_betas = which(row.names(beta.results) == "Paliperidone")
Pal_betaf = which(row.names(beta.resultsf) == "Paliperidone")

Mol_betas = which(row.names(beta.results) == "Molindone")
Mol_betaf = which(row.names(beta.resultsf) == "Molindone")

sf = c(
  Lox_betaf,
  Trif_betaf,
  Flup_betaf,
  Halop_betaf,
  Arip_betaf,
  Plac_betaf,
  Quet_betaf,
  Zip_betaf,
  Ase_betaf,
  Lur_betaf,
  Ris_betaf,
  Ola_betaf,
  Clo_betaf,
  Pal_betaf,
  Mol_betaf
)

s = c(
  Lox_betas,
  Trif_betas,
  Flup_betas,
  Halop_betas,
  Arip_betas,
  Plac_betas,
  Quet_betas,
  Zip_betas,
  Ase_betas,
  Lur_betas,
  Ris_betas,
  Ola_betas,
  Clo_betas,
  Pal_betas,
  Mol_betas
)

beta.resultsf$mean[sf] = beta.results$mean[s]

beta.resultsf$sd[sf] = beta.results$sd[s]

## add a non-informative prior (N(0,100^2)) for treatments without information 
beta.resultsf$mean = ifelse(is.na(beta.resultsf$mean), 0, beta.resultsf$mean)

beta.resultsf$sd = ifelse(is.na(beta.resultsf$sd), 100, beta.resultsf$sd)

# END OF PRIOR CONSTRUCTION FOR THE LOCATION PARAMETER
# the vector beta.results contain all the priors

# use the ind_drugs() function to create an indicator parameter which is 1 when a study is of evaluates treatments in Ta-Tc and 0 otherwise
ind_NCT = ind_drugs(all.data)

# BEGIN THE 1st STAGE OF THE APPROACH: ANALYZE GP data and construct informative priors

# prepare data for JAGS

runGP.data = long2jags(
  Study_No,
  Drug,
  mean = OverallEfficM,
  sd = OverallEfficSD,
  n = OverallEfficN,
  data = gp.data
)

# store only data who are useful for the analysis
runGP.data = list(
  ns = runGP.data$k,
  nt = runGP.data$n,
  na = runGP.data$n.k,
  t = runGP.data$T,
  y = runGP.data$Y,
  prec = runGP.data$Pr,
  pooled.sd = runGP.data$pooled.sd,
  ref = 24, ## treatment coded as 24 is the Placebo, you can verify this using: levels(as.factor(gp.data$Drug))
  b = beta.resultsf[, 1:2],
  ind = ind_NCT
)

# run the NMA model for the first stage of the method
run.GP = jags(
  data = runGP.data,
  inits = NULL,
  parameters.to.save = c("SMD.ref", "predSMD.ref", "tau"),
  n.chains = 2,
  n.iter = 100000,
  n.burnin = 10000,
  DIC = F,
  model.file = modfile_gp
)

## Store the results of the NMA model
GP.results = run.GP$BUGSoutput$summary[, c(1, 2, 3, 5, 7)]

GP.results = as.data.frame(run.GP$BUGSoutput$summary)

# add one column in the data.frame and use it as an indicator to find the useful lines in the data
GP.results$ind = as.character(rownames(GP.results))

# save predictive distribution results for all estimates
pred_all_GP = GP.results %>%
  filter(grepl("pred", ind))

# save heterogeneity results for GP
tau_GP = GP.results %>%
  filter(ind == "tau")

# keep only the useful columns
tau_GP = tau_GP %>%
  dplyr::select("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")

# create a vector with the names of all common basic comparisons between GP and CA
common_comparisonsGP = c(
  "Aripiprazole vs Placebo",
  "Asenapine vs Placebo",
  "Clozapine vs Placebo",
  "Fluphenazine vs Placebo",
  "Haloperidol vs Placebo",
  "Loxapine vs Placebo",
  "Lurasidone vs Placebo",
  "Molindone vs Placebo",
  "Olanzapine vs Placebo",
  "Paliperidone vs Placebo",
  "Quetiapine vs Placebo",
  "Risperidone vs Placebo",
  "Trifluoperazine vs Placebo",
  "Ziprasidone vs Placebo"
)

# store the results in terms of predictive distributions for common comparisons between GP and CA
basic_pred_commonGP_CA = pred_all_GP %>%
  filter(
    ind %in% c(
      "predSMD.ref[2]",
      "predSMD.ref[3]",
      "predSMD.ref[8]",
      "predSMD.ref[10]",
      "predSMD.ref[11]",
      "predSMD.ref[14]",
      "predSMD.ref[16]",
      "predSMD.ref[17]",
      "predSMD.ref[18]",
      "predSMD.ref[19]",
      "predSMD.ref[25]",
      "predSMD.ref[26]",
      "predSMD.ref[31]",
      "predSMD.ref[32]"
    )
  )

# add comparison name in the dataframe
basic_pred_commonGP_CA$comparisons = common_comparisonsGP

# add precision calculated as 1/var(SMD)
basic_pred_commonGP_CA$prec = (1 / basic_pred_commonGP_CA$sd) ^ 2

# dplyr::select only the useful columns
basic_pred_commonGP_CA = basic_pred_commonGP_CA %>%
  dplyr::select("mean",
         "sd",
         "2.5%",
         "50%",
         "97.5%",
         "Rhat",
         "ind",
         "comparisons",
         "prec")

# re-arrange the columns
basic_pred_commonGP_CA = basic_pred_commonGP_CA[, c(8, 1, 2, 3, 4, 5, 9, 6, 7)]

# give informative names to the columns
names(basic_pred_commonGP_CA) = c(
  "comparison",
  "mean",
  "sd",
  "lower",
  "median",
  "upper",
  "precision",
  "Rhat",
  "indicator"
)

# END OF THE 1st STAGE OF THE APPROACH
# The vector basic_pred_commonGP_CA contains the extrapolated SMD's which will be used as informative priors

# BEGIN THE 2nd STAGE OF THE APPROACH: ANALYZE CA data using the informative priors constructed at the 1st STAGE 

# bring ca.data in a JAGS format
runCA.data = long2jags(
  Study_No,
  Drug,
  mean = OverallEfficM,
  sd = OverallEfficSD,
  n = OverallEfficN,
  data = ca.data
)

#create a list with the useful columns and incorporate extrapolated SMD's in the data for the CA analysis
runCA.data = list(
  ns = runCA.data$k,
  nt = runCA.data$n,
  na = runCA.data$n.k,
  t = runCA.data$T,
  y = runCA.data$Y,
  prec = runCA.data$Pr,
  pooled.sd = runCA.data$pooled.sd,
  ref = 11, ## treatment coded as 11 is the Placebo, you can verify this using: levels(as.factor(ca.data$Drug))
  prior.d = basic_pred_commonGP_CA[2:8]
)

# run the model
run.CA = jags(
  data = runCA.data,
  inits = NULL,
  parameters.to.save = c(
    "SMD.ref",
    "SMD",
    "tau",
    "SUCRA"
  ),
  n.chains = 2,
  n.iter = 50000,
  n.burnin = 10000,
  DIC = F,
  model.file = modfile_ca
)

#check convergence using traceplots

run.CA_coda = jags.model(
  'CA.model.txt',
  data = runCA.data,
  n.chains = 2,
  n.adapt = 10000
)

# create coda samples
samples = coda.samples(run.CA_coda,
                        c("SMD.ref", "tau"),
                        50000)

# Appendix 1, Figure 10
plot(samples)
################################## 

# STORE THE FINAL RESULTS

# create a dataframe with the results
CA.results = as.data.frame(run.CA$BUGSoutput$summary)

# add one column in the data.frame and use it as an indicator to find the useful lines in the data
CA.results$ind = as.character(rownames(CA.results))

# results for all basic comparisons and heterogeneity parameter in CA
basic_comp_CA = CA.results %>%
  filter(grepl("ref|tau", ind))

##### heterogeneity estimate
tau_CA = CA.results %>%
  filter(ind == "tau")

# keep only the useful columns
tau_CA = tau_CA %>%
  dplyr::select("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")

# exclude Placebo and tau from basic comparion data
basic_comp_CA = basic_comp_CA %>%
  filter(ind %!in% c("SMD.ref[11]", "tau"))

# create a vector with all the comparison names
comparisonsCA = c(
  "Aripiprazole vs Placebo",
  "Asenapine vs Placebo",
  "Clozapine vs Placebo",
  "Fluphenazine vs Placebo",
  "Haloperidol vs Placebo",
  "Loxapine vs Placebo",
  "Lurasidone vs Placebo",
  "Molindone vs Placebo",
  "Olanzapine vs Placebo",
  "Paliperidone vs Placebo",
  "Quetiapine vs Placebo",
  "Risperidone vs Placebo",
  "Trifluoperazine vs Placebo",
  "Ziprasidone vs Placebo"
)


# add the comparison names in the dataframe
basic_comp_CA$comparisons = comparisonsCA

# keep only the useful columns
basic_comp_CA = basic_comp_CA %>%
  dplyr::select("mean",
         "sd",
         "2.5%",
         "50%",
         "97.5%",
         "comparisons")

# re-arrrange the dataframe
basic_comp_CA = basic_comp_CA[, c(6, 1, 2, 3, 4, 5)]
# delete rownames
row.names(basic_comp_CA) = NULL

# add informative column names in the dataframe with the results
names(basic_comp_CA) = c("comparison",
                         "mean",
                         "sd",
                         "lower",
                         "median",
                         "upper"
                         )

# add column to specify the type of the approach used
basic_comp_CA$model_type = rep("Data based beta-NCT DW", nrow(basic_comp_CA))
 
# store SUCRAs
SucrasCA = CA.results %>%
  filter(grepl("SUCRA", ind))

# keep only the useful columns
SucrasCA = SucrasCA %>%
  dplyr::select("mean", "sd", "2.5%", "97.5%")

# create a vector with the drug names
treats = c(
  "Aripiprazole",
  "Asenapine",
  "Clozapine",
  "Fluphenazine",
  "Haloperidol",
  "Loxapine",
  "Lurasidone",
  "Molindone",
  "Olanzapine",
  "Paliperidone",
  "Placebo",
  "Quetiapine",
  "Risperidone",
  "Trifluoperazine",
  "Ziprasidone"
)

# add the drug names in the SUCRA dataframe
SucrasCA$drug = treats

# arrange the dataframe in the order from the 1st to the last ranked drug
SucrasCA = SucrasCA %>%
  arrange(desc(mean))

# delete rownames
row.names(SucrasCA) = NULL

# re-arrange the dataframe
SucrasCA = SucrasCA[, c(5, 1)]

# add column to specify the type of the approach used
SucrasCA$model_type = rep("Data based beta-NCT DW", nrow(SucrasCA))

# store league table using the sourced "get_results" function
A = get_results(run.CA, parameter = "SMD")

# bring it in the class of dataframes
league_table = as.data.frame(A$leaguetable)

# column names after the name of the drugs
names(league_table) = treats

# row names after the name of the drugs
row.names(league_table) = treats

## Save the results in the "Results" folder
## These results are readily available in the "Intermediate results" folder

# These results are part of the "basic_comparisons_all.xlsx" file in the "Intermediate results" folder
write.csv(basic_comp_CA,
          "./Results/data_based_beta_NCT_DW.csv",
          row.names = F)

# These results are part of the "SUCRAS.xlsx" file in the "Intermediate results" folder
write.csv(SucrasCA,
          "./Results/SUCRA_data_based_beta_NCT_DW.csv",
          row.names = F)

write.csv(league_table,
          "./Results/Supporting_table_13.csv",
          row.names = F)
