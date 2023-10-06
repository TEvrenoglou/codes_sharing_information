library(readxl)
library(tidyverse)
library(netmeta)

## read the dataset
all.data <- read_excel("Data\\data.xlsx")

# code with the community detection algorithm. Code obtained by Law et al. (https://doi.org/10.1186/s12874-019-0689-9)
source("Codes\\helpers.R")

'%!in%' <- function(x, y)   ! ('%in%'(x, y))

# create different datasets for the two populations
ca.data = subset(all.data, all.data$ChildAdolesc == 1)

gp.data = subset(all.data, all.data$GeneralPopulation == 1)

#identify the common drugs for CA and GP
drugs = intersect(unique(ca.data$Drug), unique(gp.data$Drug))

# reduce the data to contain only the common drugs between GP and CA
drug.data = subset(all.data, all.data$Drug %in% drugs)

drug.data$OverallEfficM = as.numeric(drug.data$OverallEfficM)

drug.data = as.data.frame(drug.data)

### Remove single arm studies
s = as.data.frame(table(drug.data$Study_No))

s1 = which(s$Freq == 1)

single_armID = s$Var1[s1]

drug.data = subset(drug.data, Study_No %!in% single_armID)

ca.data = subset(drug.data, ChildAdolesc == 1)

gp.data = subset(drug.data, GeneralPopulation == 1)


## Children-adolescents

## bring the Children-adolescents data in a suitable form for the netmeta package
p_ca <- pairwise(data=ca.data,
              studlab = Study_No,
              treat = Drug,
              mean = OverallEfficM,
              sd = OverallEfficSD,
              n = OverallEfficN,
              sm="SMD"
              )

## run NMA using netmeta
mod_ca <- netmeta(p_ca,reference.group = "Placebo")

# variance-covariannce matrix for all estimates
cov_ca <- mod_ca$Cov.fixed

# get the names of the rows in the variance covariance matrix that contain "Placebo" (reference group) 
names_ca <- grep('Placebo', row.names(mod_ca$Cov.random), value=TRUE)

# find the rows in the variance covariance matrix which contain the reference group
rows_ca <- which(row.names(cov_ca) %in% names_ca)

# find the columns in the variance covariance matrix which contain the reference group
cols_ca <- which(colnames(cov_ca) %in% names_ca)

# reduce the variance covariance matrix to contain only the basic comparisons
cov_f_ca <- cov_ca[rows_ca,cols_ca]

# extract all the treatment effect estimates
ests_ca <- mod_ca$TE.random

# find the comparisons in the form "treatment vs Placebo"
ests_ca <- ests_ca[,which(colnames(ests_ca)=="Placebo")]

# remove the redundant "Placebo vs Placebo" comparison 
ests_ca <- ests_ca[-which(names(ests_ca)=="Placebo")]

# create an unnamed vector
ests_f_ca <- unname(ests_ca)

# Find the 20%, 40%, 60% and 80% quantiles of the standard errors in the Children-Adolescent network
# Appendix 1: Table 1, second row
round(quantile(sqrt(diag(cov_ca)),probs = 0.2),digits = 2)

round(quantile(sqrt(diag(cov_ca)),probs = 0.4),digits = 2)

round(quantile(sqrt(diag(cov_ca)),probs = 0.6),digits = 2)

round(quantile(sqrt(diag(cov_ca)),probs = 0.8),digits = 2)

# Run the community detection algorithm. 
cd(ests=ests_f_ca, 
   Cov=cov_f_ca, 
   method = "community", 
   weighted=FALSE, 
   seed=1,
   quants = 0.2)



## Analysis for the General patient network

## bring the General patients data in a suitable form for the netmeta package

p_gp <- pairwise(data=gp.data,
                 studlab = Study_No,
                 treat = Drug,
                 mean = OverallEfficM,
                 sd = OverallEfficSD,
                 n = OverallEfficN,
                 sm="SMD")

## run NMA using netmeta

mod_gp <- netmeta(p_gp,reference.group = "Placebo")

# variance-covariannce matrix for all estimates
cov_gp <- mod_gp$Cov.fixed

# get the names of the rows in the variance covariance matrix that contain "Placebo" (reference group) 
names_gp <- grep('Placebo', row.names(mod_gp$Cov.random), value=TRUE)

# find the rows in the variance covariance matrix which contain the reference group
rows_gp <- which(row.names(cov_gp) %in% names_gp)

# find the columns in the variance covariance matrix which contain the reference group
cols_gp <- which(colnames(cov_gp) %in% names_gp)

# reduce the variance covariance matrix to contain only the basic comparisons
cov_f_gp <- cov_gp[rows_gp,cols_gp]

# extract all the treatment effect estimates
ests_gp <- mod_gp$TE.random

# find the comparisons in the form "treatment vs Placebo"
ests_gp <- ests_gp[,which(colnames(ests_gp)=="Placebo")]

# remove the redundant "Placebo vs Placebo" comparison 
ests_gp <- ests_gp[-which(names(ests_gp)=="Placebo")]

# create an unnamed vector
ests_f_gp <- unname(ests_gp)

# match the letters with the treatments to create the graph legend
name_letter <- c("Placebo",names(ests_gp))

name_letter <- cbind.data.frame(name_letter,LETTERS[1:length(name_letter)])

names(name_letter) <- c("treat","letter")

name_letter$legend <- paste(name_letter$letter,name_letter$treat,sep = "=")

# Run the community detection algorithm. 
cd(ests=ests_f_gp, 
   Cov=cov_f_gp, 
   method = "community", 
   weighted=FALSE, 
   seed=1,
   quants = 0.8)

# add the legend in the graph
legend("topleft",legend=name_letter$legend,cex = 1,bg="lightblue")

# Find the 20%, 40%, 60% and 80% quantiles of the standard errors in the General patients network
# Appendix 1: Table 1, first row
round(quantile(sqrt(diag(cov_gp)),probs = 0.2),digits = 2)

round(quantile(sqrt(diag(cov_gp)),probs = 0.4),digits = 2)

round(quantile(sqrt(diag(cov_gp)),probs = 0.6),digits = 2)

round(quantile(sqrt(diag(cov_gp)),probs = 0.8),digits = 2)

