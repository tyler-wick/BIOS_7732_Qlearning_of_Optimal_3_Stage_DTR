# BIOS 7732 Group Project
# function to fit two-stage DTR model with qlearning
# sourced by: test_2stage.R
# author: xinyi yang

# Double bootstrap run-time for one dataset:
# Time difference of 21.99293 secs 
# Time difference of 16.57393 secs (.lm.fit)
# Time difference of 8.301749 secs (-cbind)

# function to fit two-stage model with hard-max estimation
hm_2 = function(X1, X2, Y2) {
  p1 <- ncol(X1)  
  p2 <- ncol(X2)
  
  Q2 <- .lm.fit(cbind(X1, X2), Y2)$coefficients
  
  beta <- Q2[1:p1]
  psi  <- Q2[(p1+1):(p1+p2)]
  y1hat <- X1 %*% beta + abs(X2 %*% psi)
  
  Q1 <- .lm.fit(X1, y1hat)$coefficients
  
  return(Q1)  # returns stage 1 coefficients
}
