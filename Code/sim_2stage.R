# BIOS 7732 Group Project
# sim_2stage: functions to simulate dataset for one of 9 two-stage scenarios
# sourced by: test_2stage.R
# author: karen kanaster

# two-stage parameters from Chakraborty Table 8.3
gammas <- rbind(c(0.00, 0.00,  0.00, 0.00, 0.00, 0.00, 0.00),
                c(0.00, 0.00,  0.00, 0.00, 0.01, 0.00, 0.00),
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.50),
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.49),
                c(0.00, 0.00, -0.50, 0.00, 1.00, 0.50, 0.50),
                c(0.00, 0.00, -0.50, 0.00, 0.25, 0.50, 0.50),
                c(0.00, 0.00, -0.25, 0.00, 0.75, 0.50, 0.50),
                c(0.00, 0.00,  0.00, 0.00, 0.25, 0.00, 0.25),
                c(0.00, 0.00,  0.00, 0.00, 0.25, 0.00, 0.24))

deltas <- rbind(c(0.5, 0.5), 
                c(0.5, 0.5),
                c(0.5, 0.5),
                c(0.5, 0.5),
                c(1.0, 0.0),
                c(0.1, 0.1),
                c(0.1, 0.1),
                c(0.0, 0.0),
                c(0.0, 0.0))

sim_par <- data.frame(cbind(gammas, deltas))
colnames(sim_par) <- c(paste0("gamma", 0:6), paste0("delta", 1:2))

# four possible (O2, A1) combinations from Chakraborty Table 8.2
cells <- as.data.frame(rbind(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)))
colnames(cells) <- c("O2", "A1")

get_q <- function(O2, A1, delta1, delta2, psi) {
  # P(O1 = 1) = P(O1 = -1) = P(A1 = 1) = P(A1 = -1) = 1/2
  # P(O2 = 1) = plogis(x) = exp(x)/(1 + exp(x))
  # P(O2 = -1) = 1 - P(O2 = 1) = 1/(1 + exp(x)) 
  #            = exp(-x)/(1 + exp(-x)) = plogis(-x)
  # P(O2) = plogis(O2*x) where x = delta1*O1 + delta2*A1
  
  O1 <- c(1,-1)
  q_a <- plogis(O2*(delta1*O1[1] + delta2*A1))
  q_b <- plogis(O2*(delta1*O1[2] + delta2*A1))
  
  q10 <- 1/4 * (q_a + q_b)  # cell probability averaged over O1
  q11 <- 1/4 * (q_a - q_b)  # q_prime from Chakraborty page 156
  
  if (psi == 10) return (q10) else if (psi == 11) return (q11)
}

# construct table of parameter values for 9 simulation scenarios
# checked scenarios 1-6 against "Inference for non-regular parameters"
sim_par <- sim_par %>%
  mutate(type = c("NR", "NNR", "NR", "NNR", "NR", "R", "R", "NR", "NNR"),
         q10 = get_q(O2=1, A1= 1, delta1, delta2, psi=10),
         q20 = get_q(O2=1, A1=-1, delta1, delta2, psi=10),
         q30 = q20,
         q40 = q10,
         f1 = gamma4 + gamma5*cells$O2[1] + gamma6*cells$A1[1],
         f2 = gamma4 + gamma5*cells$O2[2] + gamma6*cells$A1[2],
         f3 = gamma4 + gamma5*cells$O2[3] + gamma6*cells$A1[3],
         f4 = gamma4 + gamma5*cells$O2[4] + gamma6*cells$A1[4],
         p = 1/4 * ((f1 == 0) + (f2 == 0) + (f3 == 0) + (f4 == 0)),
         E = q10*f1   + q20*f2   + q30*f3   + q40*f4,
         E2 = q10*f1^2 + q20*f2^2 + q30*f3^2 + q40*f4^2,
         V = E2 - E^2,
         phi = E / sqrt(V),
         q11 = get_q(O2=1, A1= 1, delta1, delta2, psi=11),
         q21 = get_q(O2=1, A1=-1, delta1, delta2, psi=11),
         q31 = q11,
         q41 = q21,
         # see psi_2stage for derivation of formulas for true stage 1 parameters
         psi10 = gamma2 + q10*abs(f1) - q20*abs(f2) + q30*abs(f3) - q40*abs(f4),
         psi11 = gamma3 + q11*abs(f1) - q21*abs(f2) - q31*abs(f3) + q41*abs(f4))

# function to construct one simulated dataset
# input is one row of parameters from sim_par
sim_data <- function(par, n=300) {
  p_A1 <- 0.5
  p_A2 <- 0.5
  p_O1 <- 0.5
  
data.frame(e = rnorm(n),
           int = 1, 
           O1 = rbinom(n, size=1, prob=p_O1)*2-1,         # covariate 1
           A1 = rbinom(n, size=1, prob=p_A1)*2-1,         # treatment 1
           A2 = rbinom(n, size=1, prob=p_A2)*2-1) %>%     # treatment 2
    mutate(p_O2 = plogis(par$delta1*O1 + par$delta2*A1),
           O2 = rbinom(n, size=1, prob=p_O2)*2-1,         # covariate 2
           O1A1 = O1*A1,
           O2A2 = O2*A2,
           A1A2 = A1*A2,
           S1 = par$gamma0 + par$gamma1*O1 + par$gamma2*A1 + par$gamma3*O1A1, 
           Y2 = S1 + par$gamma4*A2 + par$gamma5*O2A2 + par$gamma6*A1A2 + e) %>%
  select(int, O1, A1, O1A1, A2, O2A2, A1A2, Y2) %>% as.matrix
}

