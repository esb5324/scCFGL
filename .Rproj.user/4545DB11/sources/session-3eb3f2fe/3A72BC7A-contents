source("R/functions_networks.R")
source("R/CFGL_functions.R")
source("R/scLink_functions.R")
library(parallel)
library(glmnet)
library(scLink)
# choose one of the two blocks of code or write your own
load("data/simu_data_1.RData") # run make_simulation_data.R first
mynorm <- simu_data
mytxt <- "1" # ex: If running S1 of the simulations, can use mytxt = "1" to indicate this

## load("data/tc_norm.RData) # run time_course.Rmd first
## mynorm <- readRDS("data/tc_norm.rds)
## mytxt <- "tc" # to indicate that you're using the time-course (tc) data

# this function runs the desired methods to get estimated networks and saves results in an RData file
# this can take a while depending on the number of penalty parameters (lambdas) and how many methods you are running
# NOTE: the chosen lambda values are for illustrative purposes only
# Arguments
## The first four arguments are for the "lambda1" parameter for the respective method (lower means more edges); if any are NA, this method will not be run
## The ..._lda2 arguments are for the "lambda2" parameters for CFGL, CFGL-r, and CFGL-rp: higher means more shared edges
## scLink_lda: length of list should be number of conditions because each vector corresponds to a condition - here, c(0.09) is for condition 1, c(0.08) is for condition 2; however, if you don't want to run scLink, just use list(c(NA), c(NA)) or list(c(NA), c(NA), c(NA)) depending on the number of conditions you have
## SM: Set to NA, and the function with calculate the Estimated W. Or, you can input a screening matrix (ex: scrmats.true if using simulated data)
## n_cores: used for parallelization, not necessary if only running 1 lambda1 and 1 lambda2
# data_mat: list of 2 or 3 expression matrices, make sure the number of genes (columns) is the same for all conditions

run_nets(SM=NA, n_cores=20, data_mat=mynorm,
         scLink_lda=list(c(0.5), c(0.5)), CFGL_lda1=c(0.5), CFGLr_lda1=c(0.5), CFGLrp_lda1=c(0.5),
         CFGL_lda2=0.01, CFGLr_lda2=0.01, CFGLrp_lda2=0.01, txt="1")
