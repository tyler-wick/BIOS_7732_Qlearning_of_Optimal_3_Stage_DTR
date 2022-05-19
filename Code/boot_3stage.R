# BIOS 7732 Group Project
# boot_2stage: functions to perform double bootstrap for three-stage model
# sourced by: test_3stage.R
# author: karen kanaster

# bootstrap sampling function used by double_boot() function
sample_df <- function(x) x[sample.int(n=nrow(x), replace=TRUE),]

fit_model <- function(x) {
  # separate dataset into design matrices and outcome
  mod <- hm_3(X1 = x[, 1:4], 
              O2 = x[, 5], 
              X2 = x[, 6:8], 
              X3 = x[, 9:11], 
              Y3 = x[, 12])
  
  return(c(psi20=mod$Q2[6],   # A2
           psi21=mod$Q2[7],   # O2A2
           psi22=mod$Q2[8],   # A1A2
           psi10=mod$Q1[3],   # A1
           psi11=mod$Q1[4]))  # O1A1
}

# get optimal bootstrap number per level (Booth & Hall 1994)
get_bootnum <- function(BC=500*100, alpha=0.05) {
  gamma <- (1/2*(1-alpha)^-2*alpha*(5/4-alpha))^(1/3)
  B2 <- 1/gamma*BC^(1/3)
  B2 <- round(B2/50)*50  
  B1 <- round(BC/B2, 0)
  return(list(B1=B1, B2=B2))
}

# get single and double bootstrap confidence intervals
# inputs are data set x and model fitting function f()
# returns ci1 and ci2: single and double bootstrap cis
double_boot <- function(x, B1=500, B2=100, alpha=0.05) {
  coef0 <- fit_model(x)
  coefs <- names(coef0)
  probs <- c(alpha/2, 1-alpha/2)
  
  boot1 <- lapply(1:B1, function(b) sample_df(x))
  coef1 <- sapply(boot1, fit_model)

  ci1 <- sapply(coefs, function(x) {
    quantile(coef1[x, ], probs)
  })
  
  u <- sapply(boot1, function(x) {
    replicate(B2, {
      boot2 <- sample_df(x)
      coef2 <- fit_model(boot2)
      coef2 <= coef0
    }) %>% rowMeans
  })
  
  ci2 <- sapply(coefs, function(x) {
    q <- quantile(u[x, ], probs)
    quantile(coef1[x, ], q)
  })
  
  return(list(psi=coef0, ci1=ci1, ci2=ci2))
}
