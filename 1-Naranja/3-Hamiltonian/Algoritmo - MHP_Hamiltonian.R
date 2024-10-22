MHP_Hamiltonian <- function(datos, 
                            M,
                            theta_0,
                            x0,
                            Sigma, 
                            iteraciones=1100,
                            periodo_quemado=100,
                            FP,
                            D,
                            N,
                            eps,
                            L,
                            phi1_m, phi2_m,
                            sig_phi1, sig_phi2,
                            rho=0.999){
  
  # Numero de parametros
  npar <- length(theta_0)
  
  # Inicializar la cadena
  cadena <- array(0, dim = c(iteraciones+periodo_quemado, npar))
  
  # Guardar valor inicial en la cadena y crear una copia para cada bloque
  theta <- cadena[1, ] <- theta_0

  # Bloques de parametros
  B1 <- 1:4 # mu_phi_1, sigma_phi_1, mu_phi_2, sigma_phi_2
  B2 <- 5:6 # sigma, xi
  
  # Tasa de aceptacion de cada bloque
  tasa_theta <-  c(0, 0)
  
  # Tasa de aceptacion de los efectos aleatorios
  tasa_ea <- rep(0, M)
  
  # Descomposicion de Cholesky de la matriz de covarianza
  S <- chol(Sigma)
  
  # Matriz inversa de la covarianza (Hamiltoniano)
  S_inv <- solve(Sigma)
  S_inv_sqrt <- chol(S_inv) 
  
  # Efectos aleatorios
  # Almacenar los efectos aleatorios phi 1
  phi1 <- array(0, dim = c(iteraciones+periodo_quemado, M))
  phi1[1,] <- phi1_m
  phi_1 <- phi1[1,]
  
  # Almacenar los efectos aleatorios phi 2  
  phi2 <- array(0, dim = c(iteraciones+periodo_quemado, M))
  phi2[1,] <- phi2_m
  phi_2 <- phi2[1,]

  # Numero de observaciones de cada individuo
  n <- dim(datos)[1] 
  
  # Numeros aleatorios para el filtro de particula para cada individuo
  n_FP <- N*D*(n-1)*M
  
  # Matriz de numeros aleatorios para el filtro de particulas N[0,1]
  bmarray=array(rnorm(n_FP,0,1),dim=c(N,D,n-1,M)) 
  
  # Matriz de numeros aleatorios para remuestreo U[0,1]
  umat=matrix(runif((n-1)*M), nrow=n-1, ncol=M)
  
  # Definir la funcion del logaritmo de verosimilitud
  fun_log_verosimilitud <- function(theta, ea, bmarray, umat) {
    filtro_particulas(datos = datos,
                      M = M, 
                      theta,
                      x0 = x0,
                      FP = FP, 
                      D = D,
                      N = N, 
                      ea,
                      bmarray,
                      umat)
  }
  
  # Calcular el logaritmo de la verosimilitud (inicial)
  log_ver <- fun_log_verosimilitud(theta, cbind(phi_1, phi_2), bmarray, umat)
  
  ################################# MCMC ######################################
  
  # Tiempo inicial
  start.time <- Sys.time()
  
  for (i in 2:(iteraciones+periodo_quemado)) {
    
    # Avance de iteraciones
    if (i %% 100 == 0){
      cat(
        sprintf("############################################################\n")
      )
      cat(
        sprintf("Ejecutando ... Iteracion %d de %d completada --> Avance %s. \n \n",
                i, iteraciones+periodo_quemado,
                paste(round(i/(iteraciones+periodo_quemado)*100), '%', sep='') )
      )
    }
   
    # Actualizar phi 1 y phi 2
    
    # Proponer nuevos efectos aleatorios
    # phi 1
    nuevo_phi1 <- phi_1 + sig_phi1*rnorm(M, 0, 1)
    
    # phi 2
    nuevo_phi2 <- phi_2 + sig_phi2*rnorm(M, 0, 1)
    
    # Si algun efecto aleatorio es negativo, no simular
    if(any(nuevo_phi1 <= 0) || any(nuevo_phi2 <= 0)){
      MHRatio <- rep(-Inf,M)
      
    } else {
      # Actualizar la matriz de numeros aleatorios para el filtro de particulas
      bmarrayprop=rho*bmarray+sqrt(1-rho^2)*rnorm(n_FP,0,1)
      
      # Actualizar la matriz de numeros aleatorios para el remuestreo
      uprop=pnorm(rho*qnorm(umat)+sqrt(1-rho^2)*rnorm((n-1)*M,0,1))
      
      # Calcular el logaritmo de la densidad a posteriori - efectos aleatorios
      log_post <- log_ver + log_eta(phi_1, phi_2, theta)
  
      # Calcular el logaritmo de la verosimilitud (nueva) - efectos aleatorios    
      log_ver_nuevo <- fun_log_verosimilitud(theta = theta, 
                                             ea = cbind(nuevo_phi1, 
                                                        nuevo_phi2),
                                             bmarrayprop, uprop)
      
      # Calcular el logaritmo de la densidad a posteriori (nueva) - efectos aleatorios
      log_post_nuevo <- log_ver_nuevo + log_eta(nuevo_phi1, nuevo_phi2,
                                                theta)  
      
      # Ratio de aceptacion - efectos aleatorios
      MHRatio <- log_post_nuevo - log_post
    }
    
    # Aceptar o rechazar - efectos aleatorios
    I <- log(runif(5)) < MHRatio 
    
    # Actualizar valores
    phi_1[I] <- nuevo_phi1[I]
    phi_2[I] <- nuevo_phi2[I]
    log_ver[I] <- log_ver_nuevo[I]
    tasa_ea <- tasa_ea + I
    bmarray[,,,I] <- bmarrayprop[,,,I]
    umat[,I] <- uprop[,I] 
    
    # Actualizar Bloque 1 (hiperparametros): 
    # mu_phi_1, log_tau_1,
    # mu_phi_2, log_tau_2
    
    # Hamiltonian Monte Carlo (HMC)

    # Crear una copia para la cadena de propuesta
    propuesta <- theta 
    
    # Ejecutar el algoritmo HMC
    Hamiltoniano_valores <- Hamiltoniano(theta[B1], 
                                         Sigma[B1,B1], 
                                         S_inv[B1,B1],
                                         S_inv_sqrt[B1,B1],
                                         eps,
                                         L, phi_1, phi_2)
    
    # Cadena de propuesta para los hiperparametros
    propuesta[B1] <- Hamiltoniano_valores$propueta
    
    # Restricciones
    theta_prop <- theta_ut(propuesta)
    restriccion_1 <- (theta_prop[1] > 150  & theta_prop[1] < 250)
    restriccion_2 <- (theta_prop[2] > 20  & theta_prop[2] < 40)
    restriccion_3 <- (theta_prop[3] > 300  & theta_prop[3] < 370)
    restriccion_4 <- (theta_prop[4] > 30  & theta_prop[4] < 70)
    restricciones <- restriccion_1*restriccion_2*restriccion_3*restriccion_4
    
    # Aceptar o rechazar la cadena propuesta
    if (restricciones == 0) {
      MHRatio <- -Inf
      
    } else {
      # Energia cinetica al final de la trayectoria
      log_K <- Hamiltoniano_valores$log_K
      log_K_nuevo <- Hamiltoniano_valores$log_K_nuevo
      
      # Evaluar las energias potencial de la trayectoria 
      # Logaritmo de la densidad a posteriori
      log_post <- log_eta(phi_1, phi_2, theta) + log_prior(theta)
      log_post_nuevo <- log_eta(phi_1, phi_2, propuesta) + log_prior(propuesta)
      
      # Ratio de aceptacion
      MHRatio <- min(0, log_post_nuevo - log_post + log_K - log_K_nuevo)
      }
    
    # Aceptar o rechazar - hiperparametros
    if (log(runif(1)) < MHRatio) {
      
      # Actualizar valores
      theta[B1] <- propuesta[B1]
      tasa_theta[1] <- tasa_theta[1] + 1
    }
    
    # Actualizar Bloque 2 (efectos comunes):
    # Actualizar sigma, xi
    
    # Proponer efectos comunes  
    propuesta <- theta
    propuesta[B2] <- theta[B2] + c(tcrossprod(rnorm(2), S[B2,B2]))
    
    # Restricciones
    theta_prop <- theta_ut(propuesta)
    restriccion_5 <- theta_prop[5] < 1
    restriccion_6 <- theta_prop[6] < 5
    restricciones <- restriccion_5*restriccion_6
    
    # Aceptar o rechazar la cadena propuesta
    if (restricciones == 0) {
      MHRatio <- -Inf
      
    } else {

      # Actualizar la matriz de numeros aleatorios para los efectos aleatorios 
      bmarrayprop <- rho*bmarray+sqrt(1-rho^2)*rnorm(n_FP,0,1)
      
      # Actualizar la matriz de numeros aleatorios para el filtro de particulas
      uprop <- pnorm(rho*qnorm(umat)+sqrt(1-rho^2)*rnorm((n-1)*M,0,1))    
      
      # Calcular el logaritmo de la densidad a posteriori - efectos comunes
      log_post <- sum(log_ver) + log_prior(theta)
  
      # Calcular el logaritmo de la verosimilitud (nuevo) - efectos comunes    
      log_ver_nuevo <- fun_log_verosimilitud(propuesta,
                                             cbind(phi_1, phi_2),
                                             bmarrayprop,
                                             uprop)
      
      # Calcular el logaritmo de la densidad a posteriori (nuevo) - efectos comunes
      log_post_nuevo <- sum(log_ver_nuevo) + log_prior(propuesta)
      
      # Ratio de aceptacion
      MHRatio <- min(0, log_post_nuevo - log_post)
    }

    # Aceptar o rechazar
      if (log(runif(1)) < MHRatio) {
        
        # Actualizar valores
        theta[B2] <- propuesta[B2]
        log_ver <- log_ver_nuevo
        tasa_theta[2] <- tasa_theta[2] + 1
        bmarray <- bmarrayprop
        umat <- uprop
    }
    
    # Actualizar cadena y efectos aleatorios
    cadena[i, ] <- theta
    phi1[i, ] <- phi_1 
    phi2[i, ] <- phi_2
    
  }
  
  # Tiempo final
  end.time <- Sys.time()
  
  # Tiempo del MCMC
  tiempo <- as.numeric(difftime(end.time, start.time, units = "secs"))
  
  return(list(cadena = cadena,
              tasa_theta = tasa_theta/(iteraciones+periodo_quemado),
              tasa_ea = tasa_ea/(iteraciones+periodo_quemado),
              tiempo = tiempo,
              phi1 = phi1,
              phi2 = phi2))
}
