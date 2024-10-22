log_eta <- function(phi_1, phi_2, theta) {
  mu_phi_1 <- theta[1]
  log_tau_phi_1 <- theta[2]
  mu_phi_2 <- theta[3]
  log_tau_phi_2 <- theta[4]

  log_p_eta_theta <- sum(logNorm(phi_1, mu_phi_1, exp(-0.5*log_tau_phi_1))+ logNorm(phi_2, mu_phi_2, exp(-0.5*log_tau_phi_2))
  )
  return(log_p_eta_theta)
}
