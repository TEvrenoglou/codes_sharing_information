library(readxl) 
library(rjags)
library(R2jags)
library(tidyverse)
library(curl)
library(dplyr)
library(devtools)
library(nmajags)


source("./Codes/NMA_model_node_splitting.R")

source("./Codes/helpers.R")

'%!in%' <- function(x,y)!('%in%'(x,y))

# read the data
all.data =
  read_excel("./Data/data.xlsx",
             na = c("", " ")) 
priors =
  read_excel("./Intermediate results/priors.xlsx",
             na = c("", " "))


### transform the data
sub.data=subset(all.data,(all.data$GeneralPopulation==1 | all.data$ChildAdolesc==1) & !is.na(all.data$OverallEfficM) & !is.na(all.data$OverallEfficSD) & !is.na(all.data$OverallEfficN))

ca.data=subset(sub.data,sub.data$ChildAdolesc==1)

### Define pairwise comparions for which consistency checks are needed.
treat1 <- c("Aripiprazole","Olanzapine","Paliperidone","Quetiapine","Risperidone")

treat2 <- rep("Placebo",length(treat1))


#### Run consistency checks for the naive model (Appendix 1: Table 2)

naive_NMA <- run_consistency_check(data = all.data,
                                   data.priors = priors,
                                   beta = "none",
                                   downweight = "none",
                                   treat1 = treat1,
                                   treat2 = treat2,
                                   ref=24,
                                   posterior.samples=F,
                                   seed = 1995,
                                   expand=T)

### get the results (Appendix 1: Table 2)
naive_NMA$results

#### Run consistency checks for the Data based beta-No DW model (Appendix 1: Table 3)

data_based_noDW <- run_consistency_check(data = ca.data,
                      data.priors = priors,
                      beta = "data",
                      downweight = "no",
                      treat1 = treat1,
                      treat2 = treat2,
                      ref=11,
                      posterior.samples=F,
                      seed = 1995)

### get the results (Appendix 1: Table 3)
data_based_noDW$results

#### Run consistency checks for the Data based beta-RoB DW model (Appendix 1: Table 4)

data_based_RoBDW <- run_consistency_check(data = ca.data,
                                         data.priors = priors,
                                         beta = "data",
                                         downweight = "rob",
                                         treat1 = treat1,
                                         treat2 = treat2,
                                         ref=11,
                                         posterior.samples=F,
                                         seed = 1995)
### get the results (Appendix 1: Table 4)
data_based_RoBDW$results

#### Run consistency checks for the Data based beta-NCT DW model (Appendix 1: Table 5)

data_based_NCTDW <- run_consistency_check(data = ca.data,
                                          data.priors = priors,
                                          beta = "data",
                                          downweight = "nct",
                                          treat1 = treat1,
                                          treat2 = treat2,
                                          ref=11,
                                          posterior.samples=F,
                                          seed = 1995)
### get the results (Appendix 1: Table 5)
data_based_NCTDW$results

#### Run consistency checks for the expert based beta-no DW model (Appendix 1: Table 6)

expert_based_noDW <- run_consistency_check(data = ca.data,
                                          data.priors = priors,
                                          beta = "expert",
                                          downweight = "no",
                                          treat1 = treat1,
                                          treat2 = treat2,
                                          ref=11,
                                          posterior.samples=F,
                                          seed = 1995)

### get the results (Appendix 1: Table 6)
expert_based_noDW$results

#### Run consistency checks for the expert based beta-RoB DW model (Appendix 1: Table 7)

expert_based_RoBDW <- run_consistency_check(data = ca.data,
                                           data.priors = priors,
                                           beta = "expert",
                                           downweight = "rob",
                                           treat1 = treat1,
                                           treat2 = treat2,
                                           ref=11,
                                           posterior.samples=F,
                                           seed = 1995)

### get the results (Appendix 1: Table 7)
expert_based_RoBDW$results

#### Run consistency checks for the expert based beta-NCT DW model (Appendix 1: Table 8)

expert_based_NCTDW <- run_consistency_check(data = ca.data,
                                            data.priors = priors,
                                            beta = "expert",
                                            downweight = "nct",
                                            treat1 = treat1,
                                            treat2 = treat2,
                                            ref=11,
                                            posterior.samples=F,
                                            seed = 1995)

### get the results (Appendix 1: Table 8)
expert_based_NCTDW$results

#### Run consistency checks for standard NMA model for CA with non-informative priors (Appendix 1: Table 9)

NMA_non_informative <- run_consistency_check(data = ca.data,
                                            data.priors = priors,
                                            beta = "none",
                                            downweight = "none",
                                            treat1 = treat1,
                                            treat2 = treat2,
                                            ref=11,
                                            posterior.samples=F,
                                            seed = 1995)

### get the results (Appendix 1: Table 9)
NMA_non_informative$results

