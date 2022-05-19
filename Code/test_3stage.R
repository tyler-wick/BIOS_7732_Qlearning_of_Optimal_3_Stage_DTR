# BIOS 7732 Group Project
# test_3stage: test simulations and model fitting code for three-stage model
# author: karen kanaster

# using hm_3 with n=300:
# N= 100, B1= 500, B2=100: Time difference of 9.931976 mins
# N=1000, B1=1000, B2=200: Time difference of 4.877481 hrs

library(dplyr)
library(future.apply)
library(qLearn)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("mod_3stage.R")
source("boot_3stage.R")
source("sim_3stage.R")

# simulation and bootstrap settings
N <- 1000  # simulated data sets
n <-  300  # dataset sample size
B <- 1000  # bootstrap rep dim 1           
C <-  200  # bootstrap rep dim 2
sc  <- 15  # sim scenario number

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

# true stage 1 and stage 2 parameters
true_psi <- list(psi_20 = sim_par[sc, "psi20"],
                 psi_21 = sim_par[sc, "psi21"],
                 psi_22 = sim_par[sc, "psi22"],
                 psi_10 = sim_par[sc, "psi10"],
                 psi_11 = sim_par[sc, "psi11"])

true_psi %>% unlist

# estimated stage 1 and stage 2 parameters
sapply(out, function(x) x$psi) %>% rowMeans

# single bootstrap confidence intervals
ci1_20 <- sapply(out, function(x) x$ci1[,1])
ci1_21 <- sapply(out, function(x) x$ci1[,2])
ci1_22 <- sapply(out, function(x) x$ci1[,3])
ci1_10 <- sapply(out, function(x) x$ci1[,4])
ci1_11 <- sapply(out, function(x) x$ci1[,5])

# double bootstrap confidence intervals
ci2_20 <- sapply(out, function(x) x$ci2[,1])
ci2_21 <- sapply(out, function(x) x$ci2[,2])
ci2_22 <- sapply(out, function(x) x$ci2[,3])
ci2_10 <- sapply(out, function(x) x$ci2[,4])
ci2_11 <- sapply(out, function(x) x$ci2[,5])

# single bootstrap coverage
mean(ci1_20[1,] < true_psi$psi_20 & true_psi$psi_20 < ci1_20[2,])
mean(ci1_21[1,] < true_psi$psi_21 & true_psi$psi_21 < ci1_21[2,])
mean(ci1_22[1,] < true_psi$psi_22 & true_psi$psi_22 < ci1_22[2,])
mean(ci1_10[1,] < true_psi$psi_10 & true_psi$psi_10 < ci1_10[2,])
mean(ci1_11[1,] < true_psi$psi_11 & true_psi$psi_11 < ci1_11[2,])

# double bootstrap coverage
mean(ci2_20[1,] < true_psi$psi_20 & true_psi$psi_20 < ci2_20[2,])
mean(ci2_21[1,] < true_psi$psi_21 & true_psi$psi_21 < ci2_21[2,])
mean(ci2_22[1,] < true_psi$psi_22 & true_psi$psi_22 < ci2_22[2,])
mean(ci2_10[1,] < true_psi$psi_10 & true_psi$psi_10 < ci2_10[2,])
mean(ci2_11[1,] < true_psi$psi_11 & true_psi$psi_11 < ci2_11[2,])

out15 <- out
save(out10, out11, out12, out13, out14, out15, file="out_combo.Rdata")