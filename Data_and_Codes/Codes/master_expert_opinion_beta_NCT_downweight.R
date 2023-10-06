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

source("Codes\\helpers.R") ## read the helpers file

source("Codes\\pairwiseMA_data_expert_based_priors.R") ## pairwise meta-analysis model to construct the priors for the location parameter using the expert opinion based approach

source("Codes\\NMA_model_expert_based_beta.R") ## NMA model for the first stage

source("Codes\\NMA_model_with_informative_priors.R") ## NMA model for the second stage

'%!in%' = function(x, y)   ! ('%in%'(x, y))

# set seed for reproducibility
set.seed(1995)

# load data
all.data =
  read_excel(
    "Data\\data.xlsx",
    na = c("", " ")
  )

## dataset with expert opinion
expert.data =
  read_excel(
    "Data\\expert_opinion.xlsx",
    na = c("", " ")
  )

## dataset with the pooled results from the GP dataset
GP.pooled.data =
  read_excel(
    "Data\\GP_pooled_data.xlsx",
    na = c("", " ")
  )


# create different datasets for the two populations
ca.data = subset(all.data, all.data$ChildAdolesc == 1)

gp.data = subset(all.data, all.data$GeneralPopulation == 1)

# BEGIN CONSTRUCTING PRIORS FOR THE LOCATION PARAMETER

#### meta-analysis model to create the priors for the location parameter based on expert's opinion

# bring expert opinion data into a JAGS format
expert.data.jags = long2jags(
  studlab = study,
  treat = drug,
  mean = mean,
  sd = sd,
  n = n,
  data = expert.data
)

expert.data.jags$pooled.sd = expert.data$sd_pooled

expert.data.jags$T = matrix(expert.data.jags$T, ncol = 1)

col2 = matrix(rep(NA, nrow(expert.data.jags$T)), ncol = 1)

col3 = matrix(rep(NA, nrow(expert.data.jags$T)), ncol = 1)

expert.data.jags$T = cbind(expert.data.jags$T, col2, col3)

# create the final list with the data
shift.data_expert = list(
  nCA = expert.data.jags$k,
  nt = expert.data.jags$n,
  naCA = expert.data.jags$n.k,
  tCA = expert.data.jags$T,
  yCA = expert.data.jags$Y,
  yGP = GP.pooled.data$mean,
  precCA = expert.data.jags$Pr,
  pooled.sdCA = rep(17.56161, length(expert.data.jags$pooled.sd)),
  # 17.56161 is the median sd_pooled across the studies in the database
  pooled.sdGP = rep(17.56161, length(GP.pooled.data$sd_pooled)),
  gamma = expert.data$gamma
)

run.shifts_expert = jags(
  data = shift.data_expert,
  inits = NULL,
  parameters.to.save = save_ma_expert,
  n.chains = 2,
  n.iter = 50000,
  n.burnin = 10000,
  DIC = F,
  model.file = modfile_expert
)

# results obtained by the meta-analysis model
beta.results_expert = run.shifts_expert$BUGSoutput$summary[, c(1, 2, 3, 5, 7)]

st=sort(unique(ca.data$Drug))

beta.results_expert=as.data.frame(beta.results_expert)

row.names(beta.results_expert)=st

# create a dataframe which contains the constructed informative priors for beta and non-informative priors for the non-common comparisons
st_gp = sort(unique(gp.data$Drug))

beta.resultsf = matrix(NA, ncol = 2, nrow = length(st_gp))

beta.resultsf = as.data.frame(beta.resultsf)

names(beta.resultsf) = c("mean", "sd")

row.names(beta.resultsf) = st_gp

### Locate each treatment's row at the two dataframes
Lox_betas = which(row.names(beta.results_expert) == "Loxapine")
Lox_betaf = which(row.names(beta.resultsf) == "Loxapine")

Trif_betas = which(row.names(beta.results_expert) == "Trifluoperazine")
Trif_betaf = which(row.names(beta.resultsf) == "Trifluoperazine")

Flup_betas = which(row.names(beta.results_expert) == "Fluphenazine")
Flup_betaf = which(row.names(beta.resultsf) == "Fluphenazine")

Halop_betas = which(row.names(beta.results_expert) == "Haloperidol")
Halop_betaf = which(row.names(beta.resultsf) == "Haloperidol")

Arip_betas = which(row.names(beta.results_expert) == "Aripiprazole")
Arip_betaf = which(row.names(beta.resultsf) == "Aripiprazole")

Plac_betas = which(row.names(beta.results_expert) == "Placebo")
Plac_betaf = which(row.names(beta.resultsf) == "Placebo")

Quet_betas = which(row.names(beta.results_expert) == "Quetiapine")
Quet_betaf = which(row.names(beta.resultsf) == "Quetiapine")

Zip_betas = which(row.names(beta.results_expert) == "Ziprasidone")
Zip_betaf = which(row.names(beta.resultsf) == "Ziprasidone")

Ase_betas = which(row.names(beta.results_expert) == "Asenapine")
Ase_betaf = which(row.names(beta.resultsf) == "Asenapine")

Lur_betas = which(row.names(beta.results_expert) == "Lurasidone")
Lur_betaf = which(row.names(beta.resultsf) == "Lurasidone")

Ris_betas = which(row.names(beta.results_expert) == "Risperidone")
Ris_betaf = which(row.names(beta.resultsf) == "Risperidone")

Ola_betas = which(row.names(beta.results_expert) == "Olanzapine")
Ola_betaf = which(row.names(beta.resultsf) == "Olanzapine")

Clo_betas = which(row.names(beta.results_expert) == "Clozapine")
Clo_betaf = which(row.names(beta.resultsf) == "Clozapine")

Pal_betas = which(row.names(beta.results_expert) == "Paliperidone")
Pal_betaf = which(row.names(beta.resultsf) == "Paliperidone")

Mol_betas = which(row.names(beta.results_expert) == "Molindone")
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

## add means and sd's in the two dataframes
beta.resultsf$mean[sf] = beta.results_expert$mean[s]

beta.resultsf$sd[sf] = beta.results_expert$sd[s]

## add a non-informative prior (N(0,100^2)) for treatments without information 
beta.resultsf$mean = ifelse(is.na(beta.resultsf$mean), 0, beta.resultsf$mean)

beta.resultsf$sd = ifelse(is.na(beta.resultsf$sd), 100, beta.resultsf$sd)

# END OF PRIOR CONSTRUCTION FOR THE LOCATION PARAMETER
# the vector beta.results_expert contain drug specific priors

# create an indicator parameter which is 1 when a study is of evaluates treatments in Ta-Tc and 0 otherwise
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
  ref = 24,
  b = beta.resultsf,
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
  model.file = modfile_GP_expert_downweight
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
  select("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")

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

# select only the useful columns
basic_pred_commonGP_CA = basic_pred_commonGP_CA %>%
  select("mean",
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
# The vector basic_pred_commonGP_CA contains the extrapolated SMD's which are used as informative priors

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
  ref = 11,
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

# Appendix 1, Figure 13
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
  select("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")

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
  select("mean",
         "sd",
         "2.5%",
         "50%",
         "97.5%",
         "comparisons")

# re-arrrange the dataframe
basic_comp_CA = basic_comp_CA[, c(6, 1, 2, 3, 4, 5)]

# delete rownames
row.names(basic_comp_CA) = NULL

# add column to specify the type of the approach used
basic_comp_CA$model_type = rep("Expert opinion beta-NCT DW", nrow(basic_comp_CA))

# store SUCRAs
SucrasCA = CA.results %>%
  filter(grepl("SUCRA", ind))

# keep only the useful columns
SucrasCA = SucrasCA %>%
  select("mean", "sd", "2.5%", "97.5%")

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
SucrasCA$model_type = rep("Expert opinion beta-NCT DW", nrow(SucrasCA))
  
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
          "Results\\expert_based_beta_NCT_DW.csv",
          row.names = F)

# These results are part of the "SUCRAS.xlsx" file in the "Intermediate results" folder
write.csv(SucrasCA,
          "Results\\SUCRA_expert_based_beta_NCT_DW.csv",
          row.names = F)

write.csv(league_table,
          "Results\\Supporting_table_16.csv",
          row.names = F)
