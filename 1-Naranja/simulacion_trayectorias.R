simulacion_trayectorias <- function(Y,
                                    particulas, 
                                    theta,
                                    phi_1, 
                                    phi_2,
                                    D,
                                    N,
                                    entre_tiempo,
                                    FP,
                                    bmarray){
  
  # FP: Puente de difusion modificado
  PDM <- 2 
  
  # FP: Puente de difusion residual
  PDR <- 3 
  
  # Discretizacion
  sim_dt <- entre_tiempo/D
  
  # Almacernar los estados simulados
  X <- matrix(0, D+1, N)
  X[1,] <- particulas # Guardar los estados iniciales
  
  # Parametros
  log_tau_sig <- theta[5]
  log_tau_xi <- theta[6]
  sigma_2 <- exp(-log_tau_xi)
  
  # Almacenar log pesos simulados
  log_pesos_sim <- rep(0, N)

  # Residual Bridge (RB) - resolver EDO drift
  if (FP == PDR){
    # Resolver numericamente para eta_t
    phi_1 <- unname(phi_1)
    phi_2 <- unname(phi_2)
    parms <- c(phi_1 = phi_1, phi_2 = phi_2)
    times <- seq(0, entre_tiempo, sim_dt)
    out <- ode(X[1,], times, model, parms)
    Xd <- out[,-1]
  }
  
  for (k in 0:(D-1)) {

    # Funcion drift
    ak <- 1/(phi_1*phi_2)*X[k+1,]*(phi_1 - X[k+1,])
    
    # Funcion de difusion al cuadrado
    vk <- exp(-log_tau_sig)*X[k+1,]
   
    # Euler Maruyama (EM) Bridge
    mu_EM <-  X[k+1,] + ak*sim_dt
    mu <-  mu_EM
    
    sigma_EM <-  sqrt(vk*sim_dt)
    sigma  <-  sigma_EM
    
    # Puente de Difusion Modificado (PDM) Bridge
    if (FP == PDM){
      dk <- (entre_tiempo - k*sim_dt)
      K <- vk*dk + sigma_2
      
      mu_MDB <- (ak*sigma_2 + vk*(Y - X[k+1,]))/K 
      psi_MDB <- (vk*sigma_2 + vk^2*(dk - sim_dt))/K
      
      # psi_MDB <- exp( log(vk)+log(sigma_2 + vk*(dk - dt)) - log(K) )
      
      mu <- X[k+1,] + mu_MDB*sim_dt
      sigma <- sqrt(psi_MDB*sim_dt)
    
    # Puente de Difusion Residual (PDR) Bridge  
    } else if(FP == PDR){
      r <- X[k+1,] - Xd[k+1,]
      
      # Delta k
      dk <- (entre_tiempo - k*sim_dt)
      
      chord <- (Xd[k+2,] - Xd[k+1,])/sim_dt
      
      K  <-  vk*dk + sigma_2
      C <-  vk*(Y - Xd[nrow(Xd),] - r - (ak - chord)*dk)
      
      mu_RB  <-  ak + C/K
      psi_RB <-  (vk*sigma_2 + vk^2*(dk - sim_dt))/K 
      
      mu <-  X[k+1,] + mu_RB*sim_dt
      sigma <-  sqrt(psi_RB*sim_dt)
    }
    
    # Simular nuevos estados
    X_nuevo <-  mu + sigma*bmarray[, k+1]
    
    if (any(X_nuevo <= 0)) {
      indices <- X_nuevo <= 0
      X_nuevo[indices] <- X[k+1, indices]
    }
    X[k+2,] <- X_nuevo 
    
    # Logaritmos de pesos
    log_pesos_sim <- log_pesos_sim + logNorm(X[k+2,], mu_EM, sigma_EM) - logNorm(X[k+2,], mu, sigma) 
  }
  
  J <- is.na(X[D+1,] + log_pesos_sim) | is.infinite(X[D+1,] + log_pesos_sim)
  if (any(J)) {
   X[D+1,J] <- particulas[J];
   log_pesos_sim[J] <- 0;
  }
  
  
  return(list(log_pesos_sim = log_pesos_sim, particulas_sim = X[D+1,]))
}
