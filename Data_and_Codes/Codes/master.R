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
library(ggplot2) ### *** Please, make sure you run version 3.3.5, as also shown in the README file
library(ggsci)
library(mvtnorm)
library(igraph)
library(gplots)

# install nmajags from github
#remotes::install_github("guido-s/nmajags")
library(nmajags)

# load functions for the Codes file

source("./Codes/helpers.R") ## read the helpers file

source("./Codes/pairwiseMA_data_expert_based_priors.R") ## pairwise meta-analysis model to construct the priors for the location parameter using the data based approach

source("./Codes/NMA_model_data_based_beta.R") ## NMA model for the first stage for the case of data-based beta

source("./Codes/NMA_model_expert_based_beta.R") ## NMA model for the first stage for the case of expert-based beta

source("./Codes/NMA_model_with_informative_priors.R") ## NMA model for the second stage

source("./Codes/standard_NMA_model.R") ## standard NMA model

source("./Codes/NMA_model_node_splitting.R") ## node splitting NMA model for consistency checks

# set seed for reproducibility
set.seed(1995)

# read the data
all.data =
  read_excel("./Data/data.xlsx",
             na = c("", " ")) 

## dataset with expert opinion
expert.data =
  read_excel(
    "./Data/expert_opinion.xlsx",
    na = c("", " ")
  )

## dataset with the pooled results from the GP dataset
GP.pooled.data =
  read_excel(
    "./Data/GP_pooled_data.xlsx",
    na = c("", " ")
  )

## dataset with the direct evidence coming from the CA dataset
direct_evidence =
  read_excel(
    "./Data/direct_evidence.xlsx",
    na = c("", " ")
  )

## priors to be used for the consistency checks 
priors =
  read_excel("./Data/priors.xlsx",
             na = c("", " "))

## densities to be used for enabling reproduction of Figures 2-6 in the Appendix 1
posterior_samples =
  read_excel("./Data/posterior_samples.xlsx",
             na = c("", " "))


############## RUN THE ANALYSIS AND REPRODUCE THE RESULTS #################################


###### Analysis using the naive pooling approach #####
res_naive_pooling <- run_simple_NMA(all.data = all.data,
                                    model_type = "naive_pooling")

### export Appendix_1_Table_10

write.csv(res_naive_pooling$league_table,"./Results/Appendix_1_Table_10.csv")

###### Analysis using model with data based beta and no downweight #####

dat_prep_data_no_dw <- data_prep_data_based(all.data = all.data,nct = F,seed = 1995)

# Fit the model
res_data_no_dw <- run.data.based(all.data = all.data,
                       gp.data = dat_prep_data_no_dw$gp.data,
                       ca.data = dat_prep_data_no_dw$ca.data,
                       beta = dat_prep_data_no_dw$beta.results,
                       downweight = "no"
)

### export Appendix_1_Table_11

write.csv(res_data_no_dw$league_table,"./Results/Appendix_1_Table_11.csv")

### save all traceplots in Appendix_1_Figure_8 in a pdf file

Appendix_1_Figure_8 <- "./Results/Appendix_1_Figure_8.pdf"

pdf(file = Appendix_1_Figure_8)

plot(res_data_no_dw$samples)

dev.off()


###### Analysis using model with data based beta and downweight according to risk of bias (rob) #####

dat_prep_data_rob_dw <- data_prep_data_based(all.data = all.data,nct = F,seed = 1995)

# Fit the model
res_data_rob_dw <- run.data.based(all.data = all.data,
                            gp.data = dat_prep_data_rob_dw$gp.data,
                            ca.data = dat_prep_data_rob_dw$ca.data,
                            beta = dat_prep_data_rob_dw$beta.results,
                            downweight = "rob"
)

### export Appendix_1_Table_12

write.csv(res_data_rob_dw$league_table,"./Results/Appendix_1_Table_12.csv")

### save all traceplots in Appendix_1_Figure_9 in a pdf file

Appendix_1_Figure_9 <- "./Results/Appendix_1_Figure_9.pdf"

pdf(file = Appendix_1_Figure_9)

plot(res_data_rob_dw$samples)

dev.off()


###### Analysis using model with data based beta and downweight according to the non-common-treatment (nct) criterion #####

dat_prep_data_nct_dw <- data_prep_data_based(all.data = all.data,nct = T,seed = 1995)

# Fit the model
res_data_nct_dw <- run.data.based(all.data = all.data,
                            gp.data = dat_prep_data_nct_dw$gp.data,
                            ca.data = dat_prep_data_nct_dw$ca.data,
                            beta = dat_prep_data_nct_dw$beta.results,
                            downweight = "nct"
)

### export Appendix_1_Table_13

write.csv(res_data_nct_dw$league_table,"./Results/Appendix_1_Table_13.csv")

### save all traceplots in Appendix_1_Figure_10 in a pdf file

Appendix_1_Figure_10 <- "./Results/Appendix_1_Figure_10.pdf"

pdf(file = Appendix_1_Figure_10)

plot(res_data_nct_dw$samples)

dev.off()


###### Analysis using model with expert opinion based beta and and no downweight ##### #####

dat_prerp_expert_no_dw <- data_prep_expert_based(all.data = all.data,
                                                  expert.opinion.data = expert.data,
                                                  pooled.data = GP.pooled.data,
                                                  nct = F,
                                                  seed = 1995 
)

# Fit the model
res_expert_no_dw <- run.expert.based(all.data = all.data,
                                      gp.data = dat_prerp_expert_no_dw$gp.data,
                                      ca.data = dat_prerp_expert_no_dw$ca.data,
                                      beta = dat_prerp_expert_no_dw$beta.results,
                                      downweight = "no"
)

### export Appendix_1_Table_14

write.csv(res_expert_no_dw$league_table,"./Results/Appendix_1_Table_14.csv")


### save all traceplots in Appendix_1_Figure_11 in a pdf file

Appendix_1_Figure_11 <- "./Results/Appendix_1_Figure_11.pdf"

pdf(file = Appendix_1_Figure_11)

plot(res_expert_no_dw$samples)

dev.off()

###### Analysis using model with expert opinion based beta and downweight according to the risk of bias (rob) ##### #####

dat_prerp_expert_rob_dw <- data_prep_expert_based(all.data = all.data,
                                                  expert.opinion.data = expert.data,
                                                  pooled.data = GP.pooled.data,
                                                  nct = F,
                                                  seed = 1995 
)

# Fit the model
res_expert_rob_dw <- run.expert.based(all.data = all.data,
                                      gp.data = dat_prerp_expert_rob_dw$gp.data,
                                      ca.data = dat_prerp_expert_rob_dw$ca.data,
                                      beta = dat_prerp_expert_rob_dw$beta.results,
                                      downweight = "rob"
)

### export Appendix_1_Table_15

write.csv(res_expert_rob_dw$league_table,"./Results/Appendix_1_Table_15.csv")


### save all traceplots in Appendix_1_Figure_12 in a pdf file

Appendix_1_Figure_12 <- "./Results/Appendix_1_Figure_12.pdf"

pdf(file = Appendix_1_Figure_12)

plot(res_expert_rob_dw$samples)

dev.off()


###### Analysis using model with expert opinion based beta and downweight according to the non-common-treatment (nct) criterion ##### #####

dat_prerp_expert_nct_dw <- data_prep_expert_based(all.data = all.data,
                                                 expert.opinion.data = expert.data,
                                                 pooled.data = GP.pooled.data,
                                                 nct = T,
                                                 seed = 1995 
)

# Fit the model
res_expert_nct_dw <- run.expert.based(all.data = all.data,
                                  gp.data = dat_prerp_expert_nct_dw$gp.data,
                                  ca.data = dat_prerp_expert_nct_dw$ca.data,
                                  beta = dat_prerp_expert_nct_dw$beta.results,
                                  downweight = "nct"
)

### export Appendix_1_Table_16

write.csv(res_expert_nct_dw$league_table,"./Results/Appendix_1_Table_16.csv")


### save all traceplots in Appendix_1_Figure_13 in a pdf file

Appendix_1_Figure_13 <- "./Results/Appendix_1_Figure_13.pdf"

pdf(file = Appendix_1_Figure_13)

plot(res_expert_nct_dw$samples)

dev.off()


###### Analysis using the standard NMA model with non-informative priors #####
res_non_informative_priors <- run_simple_NMA(all.data = all.data,
                                    model_type = "standard_NMA")

### export Appendix_1_Table_10

write.csv(res_non_informative_priors$league_table,"./Results/Appendix_1_Table_17.csv")


### reproduce Figure 2 of main manuscript


### gather all results together

data_basic_all <- rbind.data.frame(res_naive_pooling$basic_comp_CA,
                                   res_data_no_dw$basic_comp_CA,
                                   res_data_rob_dw$basic_comp_CA,
                                   res_data_nct_dw$basic_comp_CA,
                                   res_expert_no_dw$basic_comp_CA,
                                   res_expert_rob_dw$basic_comp_CA,
                                   res_expert_nct_dw$basic_comp_CA,
                                   res_non_informative_priors$basic_comp_CA,
                                   direct_evidence
                                   ) 


### main paper, Figure 2
Figure_2 <- forest_plot(data_basic_all)

ggsave(filename = "./Results/Figure_2_main_paper.pdf",plot = Figure_2,width = 15,height = 12)


#### Reproduce Appendix 1 Table 1 and Figure 1

Appendix_1_Figure_1 <- "./Results/Appendix_1_Figure_1.pdf"

pdf(file = Appendix_1_Figure_1,width = 10)

run_community(all.data = all.data, show="plot")

dev.off()

Appendix_1_Table_1 <- run_community(all.data = all.data, show="table")

write.csv(Appendix_1_Table_1,"./Results/Appendix_1_Table_1.csv",row.names = F)

#### Reproduce Appendix 1 Tables 2-9

## prepare data 

sub.data=subset(all.data,(all.data$GeneralPopulation==1 | all.data$ChildAdolesc==1) & !is.na(all.data$OverallEfficM) & !is.na(all.data$OverallEfficSD) & !is.na(all.data$OverallEfficN))

ca.data=subset(sub.data,sub.data$ChildAdolesc==1)

### Define pairwise comparisons for which consistency checks are needed
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

Appendix_1_Table_2 <- naive_NMA$results

write.csv(Appendix_1_Table_2, "./Results/Appendix_1_Table_2.csv")

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
Appendix_1_Table_3 <- data_based_noDW$results

write.csv(Appendix_1_Table_3, "./Results/Appendix_1_Table_3.csv")

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
Appendix_1_Table_4 <- data_based_RoBDW$results

write.csv(Appendix_1_Table_4, "./Results/Appendix_1_Table_4.csv")

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
Appendix_1_Table_5 <- data_based_NCTDW$results

write.csv(Appendix_1_Table_5, "./Results/Appendix_1_Table_5.csv")

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
Appendix_1_Table_6 <- expert_based_noDW$results

write.csv(Appendix_1_Table_6, "./Results/Appendix_1_Table_6.csv")

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
Appendix_1_Table_7 <- expert_based_RoBDW$results

write.csv(Appendix_1_Table_7, "./Results/Appendix_1_Table_7.csv")

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
Appendix_1_Table_8 <- expert_based_NCTDW$results

write.csv(Appendix_1_Table_8, "./Results/Appendix_1_Table_8.csv")

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
Appendix_1_Table_9 <- NMA_non_informative$results

write.csv(Appendix_1_Table_9, "./Results/Appendix_1_Table_9.csv")


#### Reproduce Appendix 1 Figures 2-7 

## Appendix 1, Figure 2
Appendix_1_Figure_2 <- density_plot(data = posterior_samples,treat="Aripiprazole")

ggsave(filename = "./Results/Appendix_1_Figure_2.pdf",plot = Appendix_1_Figure_2,width = 10,height = 7)

## Appendix 1, Figure 3
Appendix_1_Figure_3 <- density_plot(data = posterior_samples,treat="Olanzapine")

ggsave(filename = "./Results/Appendix_1_Figure_3.pdf",plot = Appendix_1_Figure_3,width = 10,height = 7)

## Appendix 1, Figure 4
Appendix_1_Figure_4 <-density_plot(data = posterior_samples,treat="Paliperidone")

ggsave(filename = "./Results/Appendix_1_Figure_4.pdf",plot = Appendix_1_Figure_4,width = 10,height = 7)

## Appendix 1, Figure 5
Appendix_1_Figure_5 <- density_plot(data = posterior_samples,treat="Quetiapine")

ggsave(filename = "./Results/Appendix_1_Figure_5.pdf",plot = Appendix_1_Figure_5,width = 10,height = 7)

## Appendix 1, Figure 6
Appendix_1_Figure_6 <- density_plot(data = posterior_samples,treat="Risperidone")

ggsave(filename = "./Results/Appendix_1_Figure_6.pdf",plot = Appendix_1_Figure_6,width = 10,height = 7)

## Appendix 1, Figure 7

### gather the SUCRA results across models

sucras <- rbind.data.frame(res_naive_pooling$SucrasCA,
                           res_data_no_dw$SucrasCA,
                           res_data_rob_dw$SucrasCA,
                           res_data_nct_dw$SucrasCA,
                           res_expert_no_dw$SucrasCA,
                           res_expert_rob_dw$SucrasCA,
                           res_expert_nct_dw$SucrasCA,
                           res_non_informative_priors$SucrasCA
                           )

## Appendix 1, Figure 7
Appendix_1_Figure_7 <- heatmap(data = sucras)

ggsave(filename = "./Results/Appendix_1_Figure_7.pdf",plot = Appendix_1_Figure_7,width = 14,height = 7)
