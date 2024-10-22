filtro_particulas <- function(datos,
                              M, 
                              theta,
                              x0,
                              FP,
                              D, 
                              N,
                              ea,
                              bmarray,
                              umat) {
  
  # Definir la funcion de observacion
  log_tau_xi <- theta[6]
  logNorm_t <- function(x,m) {
    return(0.5*(log_tau_xi - log(2*pi)) - 0.5*exp(log_tau_xi)*(x-m)^2) 
  }
  
  # Numero de observaciones
  T_obs <- dim(datos)[1]
  
  # Almacenar el logaritmo de la verosimilitud de cada individuo
  log_ver_fp <- rep(0, M)
  
  # Calcular el log verosimilitud para cada individuo
  for (m in 1:M) {
    
    # Efectos aleatorios
    phi_1 <- ea[m,1]
    phi_2 <- ea[m,2] 

    # Tiempos de observacion
    tiempos <- datos[,1,m]
    
    # Observaciones
    observaciones <- datos[,2,m]
    
    # Inicializar las particulas
    particulas <- rep(x0, N)
    
    # Logaritmo de los pesos
    log_pesos <- logNorm_t(observaciones[1], particulas)
    
    # Logaritmo de los pesos normalizados
    log_pesos_norm <- rep(-log(N), N)
    
    # Actualizar el logaritmo de la verosimilitud de cada individuo
    log_ver_fp[m] <- log_ver_fp[m] + logSumExp(log_pesos) - log(N)
    
    # Ejecutar el filtro de particulas 
    for (t in 2:T_obs) {
      # Remuestreo
      
      # Numero aleatorio para remuestreo
      u_tm <- umat[t-1,m]
      
      # Numero efectivo de particulas 
      n_efectivo <- exp(-logSumExp(2*log_pesos_norm))
      
      if (n_efectivo < 0.5*N){
        
        # Remuestreo sistematico
        I <- order(particulas)
        particulas <- sort(particulas)
        log_pesos_norm <- log_pesos_norm[I]
        I <- sysresamp2(exp(log_pesos_norm), N, u_tm)
        particulas <- particulas[I]
        log_pesos_norm <- rep(-log(N), N)
      }
      
      # Propagar las particulas
      entre_tiempo <- tiempos[t] - tiempos[t-1]
      PD <- simulacion_trayectorias(observaciones[t], 
                                    particulas,
                                    theta,
                                    phi_1, 
                                    phi_2,
                                    D, 
                                    N,
                                    entre_tiempo,
                                    FP,
                                    bmarray[,,t-1,m])
      particulas <- PD$particulas_sim
      
      # Ponderacion y normalizacion
      log_pesos <- logNorm_t(observaciones[t], particulas) + PD$log_pesos_sim

      # Extraer el mayor log peso (estabilidad numerica)
      max_log_pesos <- max(log_pesos)
      log_pesos <- log_pesos - max_log_pesos
      
      # Normalizar los pesos
      log_pesos_norm <- log_pesos - logSumExp(log_pesos)
      
      # Calcular el logaritmo de la verosimilitud
      log_ver_fp[m] <- log_ver_fp[m] + logSumExp(log_pesos) + max_log_pesos - log(N)
    }
  }
  return(log_ver_fp)  
}
