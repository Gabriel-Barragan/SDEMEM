log_prior <- function(theta){
  
  # Logaritmo de la densidad a priori
  mu_phi_1 <- theta[1]
  log_tau_1 <- theta[2]
  
  mu_phi_2 <- theta[3]
  log_tau_2 <- theta[4]
  
  log_tau_sig <- theta[5]
  log_tau_xi <- theta[6]
  
  Log_prior <- sum(logNorm(mu_phi_1, 190, 15),
                     logGamma(log_tau_1, 2, 0.1),
                     logNorm(mu_phi_2, 350, 20),
                     logGamma(log_tau_2, 5, 0.1),
                     logGamma(log_tau_sig, 5, 10),
                     logGamma(log_tau_xi, 5, .5) 
  )
  
  return(Log_prior)
}
