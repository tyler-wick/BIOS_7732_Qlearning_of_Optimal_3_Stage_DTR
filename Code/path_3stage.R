# BIOS 7732 Group Project
# path_3stage: example of optimal treatment paths for three-stage model
# author: karen kanaster

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("sim_3stage.R")

# construct all combinations of response and previous treatment
all_paths <- expand.grid(O1 = c(-1,1), A1 = c(-1,1), O2 = c(-1,1), A2 = c(-1,1), O3 = c(-1,1))

# compute optimal decision rules d1, d2, d3 for each combination
sim_paths <- sim_par %>% 
  select(scen, type, gamma7, gamma8, gamma9, contains("psi")) %>%
  slice(rep(1:n(), each=nrow(all_paths))) %>%
  bind_cols(all_paths %>% slice(rep(1:n(), times=15))) %>%
  mutate(f3 = gamma7 + gamma8*O3 + gamma9*A2,
         f2 = psi20 + psi21*O2 + psi22*A1,
         f1 = psi10 + psi11*O1,
         d3 = sign(f3),
         d2 = sign(f2),
         d1 = sign(f1)) %>%
  subset(A1 == d1 & A2 == d2 & d3 != 0)

# extract paths for scenario 6
sim_paths_6 <- sim_paths %>% 
  subset(scen == "6")

# stage 1 optimal decision rules
sim_paths_6 %>% 
  select(psi10, psi11, O1, f1, d1) %>%
  arrange(-O1) %>%
  unique

# stage 2 optimal decision rules
sim_paths_6 %>% 
  select(psi20, psi21, O2, psi22, A1, f2, d2) %>%
  arrange(-A1, -O2) %>%
  unique

# stage 3 optimal decision rules
sim_paths_6 %>% 
  select(gamma7, gamma8, O3, gamma9, A2, f3, d3) %>%
  arrange(-A2, -O3) %>%
  unique
