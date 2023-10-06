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

source("Codes\\standard_NMA_model.R") ## standard NMA model

'%!in%' = function(x, y)   ! ('%in%'(x, y))

# set seed for reproducibility
set.seed(1995)

# load data
all.data =
  read_excel("Data\\data.xlsx",
             na = c("", " "))

# restrict data to contain only Children-Adolescent data
all.data = subset(all.data, all.data$ChildAdolesc == 1)

runCA.data = long2jags(
  Study_No,
  Drug,
  mean = OverallEfficM,
  sd = OverallEfficSD,
  n = OverallEfficN,
  data = all.data
)

runCA.data = list(
  ns = runCA.data$k,
  nt = runCA.data$n,
  na = runCA.data$n.k,
  t = runCA.data$T,
  y = runCA.data$Y,
  prec = runCA.data$Pr,
  pooled.sd = runCA.data$pooled.sd,
  ref = 11
)

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
  model.file = modfile_standardNMA
)

CA.results = as.data.frame(run.CA$BUGSoutput$summary)

CA.results$ind = as.character(rownames(CA.results))

####### Manipulate the results and create suitable datasets

#### Basic comparisons CA
basic_comp_CA = CA.results %>%
  filter(grepl("ref|tau", ind))

##### Heterogeneity parameter
tau_CA = CA.results %>%
  filter(ind == "tau")

tau_CA = tau_CA %>%
  select("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")

#### Exclude Placebo and tau from basic comparion data
basic_comp_CA = basic_comp_CA %>%
  filter(ind %!in% c("SMD.ref[11]", "tau"))

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


#### NMA estimates results
basic_comp_CA$comparisons = comparisonsCA

basic_comp_CA = basic_comp_CA %>%
  select("mean",
         "sd",
         "2.5%",
         "50%",
         "97.5%",
         "comparisons")

basic_comp_CA = basic_comp_CA[, c(6, 1, 2, 3, 4, 5)]

row.names(basic_comp_CA) = NULL

# add column to specify the type of the approach used
basic_comp_CA$model_type = rep("NMA with non-informative priors", nrow(basic_comp_CA))

################ SUCRAS for CA ############

SucrasCA = CA.results %>%
  filter(grepl("SUCRA", ind))

SucrasCA = SucrasCA %>%
  select("mean", "sd", "2.5%", "97.5%")

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

SucrasCA$drug = treats

SucrasCA = SucrasCA %>%
  arrange(desc(mean))

row.names(SucrasCA) = NULL

SucrasCA = SucrasCA[, c(5, 1)]

# add column to specify the type of the approach used
SucrasCA$model_type = rep("NMA with non-informative priors", nrow(SucrasCA))
  
#### League table ###

A = get_results(run.CA, parameter = "SMD")

league_table = as.data.frame(A$leaguetable)

names(league_table) = levels(as.factor(unique(all.data$Drug)))

row.names(league_table) = levels(as.factor(unique(all.data$Drug)))

league_table = league_table[treats, treats]

## Save the results in the "Results" file
## These results are readily available in the "Intermediate Results" file 

write.csv(basic_comp_CA,
          "Results\\nma_non_informative_priors.csv",
          row.names = F)

write.csv(SucrasCA,
          "Results\\SUCRA_nma_non_informative_priors.csv",
          row.names = F)

write.csv(league_table,
          "Results\\Supporting_table_17.csv",
          row.names = F)
