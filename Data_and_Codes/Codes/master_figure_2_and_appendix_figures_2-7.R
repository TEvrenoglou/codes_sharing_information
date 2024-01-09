library(tidyverse)    
library(ggplot2)
library(ggsci)
library(readxl)

## intermediate results with the posterior samples to reproduce Figure 2 of the main paper and Figures 2-7 of the Appendix 1

basic_comparisons <- read_excel("./Intermediate Results/basic_comparisons_all.xlsx")

posterior_samples <- read_excel("./Intermediate Results/posterior_samples.xlsx")

sucras <- read_excel("./Intermediate Results/SUCRAS.xlsx")

## read the functions to reproduce the graphs
source("./Codes/helpers.R")

### main paper, Figure 2
forest_plot(basic_comparisons)

## Appendix 1, Figure 2
density_plot(data = posterior_samples,treat="Aripiprazole")

## Appendix 1, Figure 3
density_plot(data = posterior_samples,treat="Olanzapine")

## Appendix 1, Figure 4
density_plot(data = posterior_samples,treat="Paliperidone")

## Appendix 1, Figure 5
density_plot(data = posterior_samples,treat="Quetiapine")

## Appendix 1, Figure 6
density_plot(data = posterior_samples,treat="Risperidone")

## Appendix 1, Figure 7
heatmap(data = sucras)


