# BIOS 7732 Group Project
# sim_3stage: functions to simulate dataset for one of 9 three-stage scenarios
# sourced by: test_3stage.R
# author: karen kanaster

# two-stage parameters from Chakraborty Table 8.3
gammas <- rbind(c(0.00, 0.00,  0.00, 0.00, 0.00, 0.00, 0.00),  # 1
                c(0.00, 0.00,  0.00, 0.00, 0.01, 0.00, 0.00),  # 2
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.50),  # 3
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.49),  # 4
                c(0.00, 0.00, -0.50, 0.00, 1.00, 0.50, 0.50),  # 5
                c(0.00, 0.00, -0.50, 0.00, 0.25, 0.50, 0.50),  # 6
                c(0.00, 0.00, -0.25, 0.00, 0.75, 0.50, 0.50),  # A
                c(0.00, 0.00,  0.00, 0.00, 0.25, 0.00, 0.25),  # B
                c(0.00, 0.00,  0.00, 0.00, 0.25, 0.00, 0.24))  # C

deltas <- rbind(c(0.5, 0.5),  # 1
                c(0.5, 0.5),  # 2
                c(0.5, 0.5),  # 3
                c(0.5, 0.5),  # 4
                c(1.0, 0.0),  # 5
                c(0.1, 0.1),  # 6
                c(0.1, 0.1),  # A
                c(0.0, 0.0),  # B
                c(0.0, 0.0))  # C

sim_par <- data.frame(cbind(gammas, deltas))
colnames(sim_par) <- c(paste0("gamma", 0:6), paste0("delta", 1:2))

# extension of two-stage parameters from Chakraborty Table 8.3
# A2: (gamma4, gamma5, gamma6) -> A3: (gamma7, gamma8, gamma9)
sim_par <- sim_par %>%
  mutate(gamma7 = gamma4, 
         gamma8 = gamma5, 
         gamma9 = gamma6,
         delta3 = delta1,
         delta4 = delta2)

# combinations of scenarios 3-NR, 4-NNR, and 6-R
gammas <- rbind(c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.50, 0.50, 0.00, 0.49),  # 3-4
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.50, 0.25, 0.50, 0.50),  # 3-6
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.49, 0.50, 0.00, 0.50),  # 4-3
                c(0.00, 0.00, -0.50, 0.00, 0.50, 0.00, 0.49, 0.25, 0.50, 0.50),  # 4-6
                c(0.00, 0.00, -0.50, 0.00, 0.25, 0.50, 0.50, 0.50, 0.00, 0.50),  # 6-3
                c(0.00, 0.00, -0.50, 0.00, 0.25, 0.50, 0.50, 0.50, 0.00, 0.49))  # 6-4

deltas <- rbind(c(0.5, 0.5, 0.5, 0.5),  # 3-4
                c(0.5, 0.5, 0.1, 0.1),  # 3-6
                c(0.5, 0.5, 0.5, 0.5),  # 4-3
                c(0.5, 0.5, 0.1, 0.1),  # 4-6
                c(0.1, 0.1, 0.5, 0.5),  # 6-3
                c(0.1, 0.1, 0.5, 0.5))  # 6-4

combos <- data.frame(cbind(gammas, deltas))
colnames(combos) <- c(paste0("gamma", 0:9), paste0("delta", 1:4))

sim_par <- rbind(sim_par, combos)

# four possible (O2, A1) combinations from Chakraborty Table 8.2
cells <- rbind(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
cells <- as.data.frame(cbind(cells, cells))
colnames(cells) <- c("O3", "A2", "O2", "A1")
cells$O1 <- cells$O2

# return q and q_prime for stage 1 true psi components
get_q1 <- function(O2, A1, delta1, delta2, psi_sub) {
  # P(O1 = 1) = P(O1 = -1) = P(A1 = 1) = P(A1 = -1) = 1/2
  # P(O2 = 1) = plogis(x) = exp(x)/(1 + exp(x))
  # P(O2 = -1) = 1 - P(O2 = 1) = 1/(1 + exp(x)) 
  #            = exp(-x)/(1 + exp(-x)) = plogis(-x)
  # P(O2) = plogis(O2*x) where x = delta1*O1 + delta2*A1
  
  O1 <- c(1,-1)
  q_a <- plogis(O2*(delta1*O1[1] + delta2*A1))
  q_b <- plogis(O2*(delta1*O1[2] + delta2*A1))
  
  q_0 <- 1/4 * (q_a + q_b)  # cell probability averaged over O1
  q_1 <- 1/4 * (q_a - q_b)  # q_prime from Chakraborty page 156

  # psi_10 = treatment main effect, psi_11 = treatment interaction
  if (psi_sub == 10) {
    return (q_0) 
  } else if (psi_sub == 11) {
    return (abs(q_1))
  }
}

# return q and q' for stage 2 true psi components
get_q2 <- function(O3, A2, delta3, delta4, psi_sub) {
  get_q1(O2=O3, A1=A2, delta1=delta3, delta2=delta4, psi_sub=psi_sub-10)
}

# construct table of parameter values for 9 simulation scenarios
sim_par <- sim_par %>%
  mutate(scen = c(1:6, LETTERS[1:3], "3-4", "3-6", "4-3", "4-6", "6-3", "6-4"),
         type = c("NR", "NNR", "NR", "NNR", "NR", "R", "R", "NR", "NNR",
                  "NR-NNR", "NR-R", "NNR-NR", "NNR-R", "R-NR", "R-NNR"),
         q20_1 = get_q2(O3=cells$O3[1], A2=cells$A2[1], delta3, delta4, psi_sub=20),
         q20_2 = get_q2(O3=cells$O3[2], A2=cells$A2[2], delta3, delta4, psi_sub=20),
         q20_3 = get_q2(O3=cells$O3[3], A2=cells$A2[3], delta3, delta4, psi_sub=20),
         q20_4 = get_q2(O3=cells$O3[4], A2=cells$A2[4], delta3, delta4, psi_sub=20),
         q21_1 = get_q2(O3=cells$O3[1], A2=cells$A2[1], delta3, delta4, psi_sub=21),
         q21_2 = get_q2(O3=cells$O3[2], A2=cells$A2[2], delta3, delta4, psi_sub=21),
         q21_3 = get_q2(O3=cells$O3[3], A2=cells$A2[3], delta3, delta4, psi_sub=21),
         q21_4 = get_q2(O3=cells$O3[4], A2=cells$A2[4], delta3, delta4, psi_sub=21),
         f21 = gamma7 + gamma8*cells$O3[1] + gamma9*cells$A2[1],
         f22 = gamma7 + gamma8*cells$O3[2] + gamma9*cells$A2[2],
         f23 = gamma7 + gamma8*cells$O3[3] + gamma9*cells$A2[3],
         f24 = gamma7 + gamma8*cells$O3[4] + gamma9*cells$A2[4],
         beta24 = q21_1*(abs(f21) + abs(f22) - abs(f23) - abs(f24)),
         psi20 = q20_1*abs(f21) - q20_2*abs(f22) + q20_3*abs(f23) - q20_4*abs(f24),
         psi21 = q21_1*abs(f21) - q21_2*abs(f22) - q21_3*abs(f23) + q21_4*abs(f24),
         gamma4 = gamma4 - psi20,
         gamma5 = gamma5 - psi21,
         psi20  = gamma4 + psi20,
         psi21  = gamma5 + psi21,
         psi22  = gamma6,
         p_f2 = 1/4 * ((f21 == 0) + (f22 == 0) + (f23 == 0) + (f24 == 0)),
         E_f2 = q20_1*f21 + q20_2*f22 + q20_3*f23 + q20_4*f24,
         E2_f2 = q20_1*f21^2 + q20_2*f22^2 + q20_3*f23^2 + q20_4*f24^2,
         V_f2 = E2_f2 - E_f2^2,
         phi_f2 = E_f2 / sqrt(V_f2),
         q10_1 = get_q1(O2=cells$O2[1], A1=cells$A1[1], delta1, delta2, psi_sub=10),
         q10_2 = get_q1(O2=cells$O2[2], A1=cells$A1[2], delta1, delta2, psi_sub=10),
         q10_3 = get_q1(O2=cells$O2[3], A1=cells$A1[3], delta1, delta2, psi_sub=10),
         q10_4 = get_q1(O2=cells$O2[4], A1=cells$A1[4], delta1, delta2, psi_sub=10),
         q11_1 = get_q1(O2=cells$O2[1], A1=cells$A1[1], delta1, delta2, psi_sub=11),
         q11_2 = get_q1(O2=cells$O2[2], A1=cells$A1[2], delta1, delta2, psi_sub=11),
         q11_3 = get_q1(O2=cells$O2[3], A1=cells$A1[3], delta1, delta2, psi_sub=11),
         q11_4 = get_q1(O2=cells$O2[4], A1=cells$A1[4], delta1, delta2, psi_sub=11),
         f11 = psi20 + psi21*cells$O2[1] + psi22*cells$A1[1],
         f12 = psi20 + psi21*cells$O2[2] + psi22*cells$A1[2],
         f13 = psi20 + psi21*cells$O2[3] + psi22*cells$A1[3],
         f14 = psi20 + psi21*cells$O2[4] + psi22*cells$A1[4],
         psi10 = 2*q10_1*beta24 - 2*q10_2*beta24,
         gamma2 = gamma2 - psi10, 
         psi10 = psi10 + gamma2 + q10_1*abs(f11) - q10_2*abs(f12) + q10_3*abs(f13) - q10_4*abs(f14),
         psi11 = gamma3 + q11_1*abs(f11) - q11_2*abs(f12) - q11_3*abs(f13) + q11_4*abs(f14),
         p_f1 = 1/4 * ((f11 == 0) + (f12 == 0) + (f13 == 0) + (f14 == 0)),
         E_f1 = q10_1*f11 + q10_2*f12 + q10_3*f13 + q10_4*f14,
         E2_f1 = q10_1*f11^2 + q10_2*f12^2 + q10_3*f13^2 + q10_4*f14^2,
         V_f1 = E2_f1 - E_f1^2,
         phi_f1 = E_f1 / sqrt(V_f1))

# function to construct one simulated dataset
# input is one row of parameters from sim_par
sim_data <- function(par, n=300) {
  p_A1 <- 0.5
  p_A2 <- 0.5
  p_A3 <- 0.5
  p_O1 <- 0.5
  
data.frame(e = rnorm(n),
           int = 1, 
           A1 = rbinom(n, size=1, prob=p_A1)*2-1,         # treatment 1
           A2 = rbinom(n, size=1, prob=p_A2)*2-1,         # treatment 2
           A3 = rbinom(n, size=1, prob=p_A3)*2-1,         # treatment 3
           O1 = rbinom(n, size=1, prob=p_O1)*2-1) %>%     # covariate 1  
    mutate(p_O2 = plogis(par$delta1*O1 + par$delta2*A1),
           O2 = rbinom(n, size=1, prob=p_O2)*2-1,         # covariate 2
           p_O3 = plogis(par$delta3*O2 + par$delta4*A2),
           O3 = rbinom(n, size=1, prob=p_O3)*2-1,         # covariate 3
           O1A1 = O1*A1,
           O2A2 = O2*A2,
           A1A2 = A1*A2,
           O3A3 = O3*A3,
           A2A3 = A2*A3,
           S1 = par$gamma0 + par$gamma1*O1 + par$gamma2*A1 + par$gamma3*O1A1, 
           S2 = par$gamma4*A2 + par$gamma5*O2A2 + par$gamma6*A1A2,
           Y3 = S1 + S2 + par$gamma7*A3 + par$gamma8*O3A3 + par$gamma9*A2A3 + e) %>%
  select(int, O1, A1, O1A1, O2, A2, O2A2, A1A2, A3, O3A3, A2A3, Y3) %>% as.matrix
}