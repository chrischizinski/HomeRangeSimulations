# options(device=NULL)

library(circular)
library(ggplot2)
library(dplyr)
library(adehabitatHR)
library(ade4)
require(rgeos)
require(maptools)
library(tidyverse)
library(parallel)
library(foreach)
library(doSNOW)

# set.seed(55)

source("/R/squareworldfunctions_102116.R")

############
# Parameters
Norgs<-1
Nsteps<-120
by_samples<-c(1:4,6,12,1:4,6,12)
iter<-1
ecotype=c("desert", "sage", "prairie")
unin <- "m"
unout <- "m2"
do.print = FALSE
#################

#Generate world
myworld <- make_world2 (1200, 1200)
# 

RNGkind("L'Ecuyer-CMRG")

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores, outfile="") 
registerDoSNOW(cl)
iterations <- 1000
# pb <- txtProgressBar(min = 1, max = iterations, style = 3)
result <- foreach(q=1:iterations,
                  .packages = c("circular","adehabitatHR","ade4","rgeos","maptools", "tidyverse"),
                  .combine = rbind) %dopar% hr_simulations(q,myworld, do.print = TRUE)

stopCluster(cl)

result2<-do.call("rbind",result)

# save(result2,file="sim_results_data_frame.RData")
