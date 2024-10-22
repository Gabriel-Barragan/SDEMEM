log_gradiente <- function(phi_1, phi_2, theta) {
  # parametros (theta)
  mu_phi_1 <- theta[1]
  log_tau_phi_1 <- theta[2]
  mu_phi_2 <- theta[3]
  log_tau_phi_2 <- theta[4]

  M <- length(phi_1)
  
  # gradiente (theta)
  g_mu_phi_1 <- exp(log_tau_phi_1)*(sum(phi_1) - M*mu_phi_1) - mu_phi_1/(100)^2
  g_tau_phi_1 <- 0.5*M - 0.5*exp(log_tau_phi_1)*sum((phi_1 - mu_phi_1)**2) + 1 - 0.01*exp(log_tau_phi_1)
  
  g_mu_phi_2 <- exp(log_tau_phi_2)*(sum(phi_2) - M*mu_phi_2) - mu_phi_2/(100)^2
  g_tau_phi_2 <- 0.5*M - 0.5*exp(log_tau_phi_2)*sum((phi_2 - mu_phi_2)**2) + 1 - 0.01*exp(log_tau_phi_2)

  return(c(g_mu_phi_1, g_tau_phi_1, g_mu_phi_2, g_tau_phi_2))  
}
