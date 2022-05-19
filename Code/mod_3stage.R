# BIOS 7732 Group Project
# function to fit three-stage DTR model with qlearning
# sourced by: test_3stage.R
# author: xinyi yang

# Double bootstrap run-time for one dataset:
# hm_2:  Time difference of 8.301749 secs
# hm_3:  Time difference of 14.31958 secs

# function to fit three-stage model with hard-max estimation
hm_3 = function(X1, O2, X2, X3, Y3) {
  
  p1 = ncol(X1)  # int, O1, A1, O1A1
  p2 = ncol(X2)  # A2, O2A2, A1A2
  p3 = ncol(X3)  # A3, O3A3, A2A3

  p <- p1 + p2 + p3
  X <- cbind(X1, X2, X3)
  Q3 = .lm.fit(X, Y3)$coefficients

  X12 <- cbind(X1, X2)
  beta3 = Q3[1:(p-p3)]
  psi3 = Q3[(p-p3+1):p]

  y2hat = X12 %*% beta3 + abs(X3 %*% psi3)

  p <- p1 + 1 + p2
  X12 <- cbind(X1, O2, X2)  # add O2 KK
  Q2 = .lm.fit(X12, y2hat)$coefficients
  
  H20 = cbind(X1, O2)       # add O2 KK
  beta2 = Q2[1:(p-p2)]
  psi2  = Q2[(p-p2+1):p]

  y1hat = H20 %*% beta2 + abs(X2 %*% psi2)
  Q1 = .lm.fit(X1, y1hat)$coefficients
  
  return(list(Q2 = Q2, Q1 = Q1))  
  # Q2: "1", "O1", "A1", "O1*A1", "O2", A2", "O2*A2", "A1*A2"
  # Q1: "1", "O1", "A1", "O1*A1"
}
