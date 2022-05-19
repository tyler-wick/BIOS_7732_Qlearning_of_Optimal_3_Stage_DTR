# BIOS 7732 Group Project
# boot_2stage: functions to perform double bootstrap for two-stage model
# sourced by: test_2stage.R
# author: karen kanaster

# bootstrap sampling function used by double_boot() function
sample_df <- function(x) x[sample.int(n=nrow(x), replace=TRUE),]

fit_model <- function(x) {
  # separate dataset into design matrices and outcome
  mod <- hm_2(X1 = x[, 1:4], 
              X2 = x[, 5:7], 
              Y2 = x[, 8])
  
  return(c(psi10=mod[3],   # A1
           psi11=mod[4]))  # O1A1
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
