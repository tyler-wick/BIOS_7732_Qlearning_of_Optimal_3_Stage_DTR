# BIOS 7732 Group Project
# test_2stage: test simulations and model fitting code for two-stage model
# author: karen kanaster

# note: tested on scenario 6 with using qLearn::qLearn() function
# output is reasonable but run-time is 50 minutes for 10 datasets!

# using hm_2 with n=300:
# N= 100, B1= 500,  B2=100: Time difference of 5.707698 mins
# N=1000, B1= 500,  B2=100: Time difference of 55.66289 mins
# N=1000, B1=1000,  B2=200: Time difference of 2.232852 hours

library(future.apply)
library(dplyr)
library(qLearn)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("mod_2stage.R")
source("boot_2stage.R")
source("sim_2stage.R")

# simulation and bootstrap settings
N <- 1000  # simulated data sets
n <-  300  # dataset sample size
B <- 1000  # bootstrap rep dim 1
C <-  200  # bootstrap rep dim 2
sc  <-  9  # sim scenario number

# optimal bootstraps (Booth & Hall 1994)
bootnum <- get_bootnum(BC=B*C)
B1 <- bootnum$B1
B2 <- bootnum$B2

# loop: simulate N datasets for single scenario
plan(multisession, workers=8)
t1 <- Sys.time()
out <- future_lapply(1:N, function(i) {
  x <- sim_data(sim_par[sc,], n=n)
  double_boot(x, B1=B1, B2=B2)
}, future.seed=628*sc)
Sys.time() - t1

# true stage 1 parameters
true_psi <- list(psi_10 = sim_par[sc, "psi10"],
                 psi_11 = sim_par[sc, "psi11"])

true_psi %>% unlist

# estimated stage 1 parameters
sapply(out, function(x) x$psi) %>% rowMeans

# single bootstrap confidence intervals
ci1_10 <- sapply(out, function(x) x$ci1[,1])
ci1_11 <- sapply(out, function(x) x$ci1[,2])

# double bootstrap confidence intervals
ci2_10 <- sapply(out, function(x) x$ci2[,1])
ci2_11 <- sapply(out, function(x) x$ci2[,2])

# single bootstrap coverage
mean(ci1_10[1,] < true_psi$psi_10 & true_psi$psi_10 < ci1_10[2,])
mean(ci1_11[1,] < true_psi$psi_11 & true_psi$psi_11 < ci1_11[2,])

# double bootstrap coverage
mean(ci2_10[1,] < true_psi$psi_10 & true_psi$psi_10 < ci2_10[2,])
mean(ci2_11[1,] < true_psi$psi_11 & true_psi$psi_11 < ci2_11[2,])

# out9 <- out
# save(out1, out2, out3, out4, out5, out6, out7, out8, out9, file="out_2stage.Rdata")