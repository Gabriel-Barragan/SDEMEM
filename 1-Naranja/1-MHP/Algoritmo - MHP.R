MHP <- function(datos, 
                M,
                theta_0,
                x0,
                Sigma,
                iteraciones=1100,
                periodo_quemado=100,
                FP,
                D,
                N,
                L, 
                rho,
                epsilon=1e-6){
  
  # Numero de parametros
  npar <- length(theta_0)
  
  # Inicializar la cadena
  cadena <- array(0, dim = c(iteraciones+periodo_quemado, npar))
  
  # Guardar valor inicial en la cadena
  cadena[1, ] <- theta_0
  
  # Tasa de aceptacion
  tasa_theta <- 0
  
  # Descomposicion de Cholesky de la matriz de covarianza
  S <- chol(Sigma)
  
  # Parametro de actualizacion de la matriz de covarianza
  sd <- 2.4^2/npar
  
  # Logaritmo de la verosimilitud de los datos
  lv <- rep(0, iteraciones+periodo_quemado)
  
  # Generar efectos aleatorios
  # Numeros aleatorios para los efectos aleatorios para cada individuo
  n_ea <- M*L*2
  
  # Matriz de numeros aleatorios para los efectos aleatorios N[0,1]
  bmarray_ea <- array(rnorm(n_ea, 0, 1), dim = c(M*L, 2))
  
  mu_phi_1 <- theta_0[1]
  sigma_phi_1 <- exp(-0.5*theta_0[2])
  mu_phi_2 <- theta_0[3]
  sigma_phi_2 <- exp(-0.5*theta_0[4])
  
  # phi_1
  phi_1_ea <- matrix(mu_phi_1+sigma_phi_1*bmarray_ea[,1], M, L)
  # phi_2
  phi_2_ea <- matrix(mu_phi_2+sigma_phi_2*bmarray_ea[,2], M, L)
    
  # Numero de observaciones de cada individuo
  n <- dim(datos)[1]
  
  # Numeros aleatorios para el filtro de particula para cada individuo
  n_FP <- N*D*(n-1)*M
  
  # Matriz de numeros aleatorios para el filtro de particulas N[0,1]
  bmarray_FP <- array(rnorm(n_FP, 0, 1), dim=c(N, D, n-1, M)) 
  
  # Matriz de numeros aleatorios para remuestreo U[0,1]
  umat <- matrix(runif((n-1)*M), nrow=n-1, ncol=M)
  
  # Definir la funcion de verosimilitud
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
  
  # Calcular el logaritmo de la verosimilitud de cada efecto aleatorio L
  lv_ea <- array(0, dim = c(M, L))
  for (l in 1:L) {
    ea_l <- cbind(phi_1_ea[,l], phi_2_ea[,l])
    lv_ea[,l] <- fun_log_verosimilitud(cadena[1,], 
                                       ea_l,
                                       bmarray_FP,
                                       umat)
  }
  
  # Calcular el logaritmo de la verosimilitud de todos los individuos
  log_ver <- sum(apply(lv_ea, 1, logSumExp) - log(L))
  lv[1] <- log_ver
  
  # Calcular el logaritmo de la densidad a posteriori (inicial)
  log_post <- log_ver + log_prior(cadena[1,])
  
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
    
    # Actualizar de la matriz de covarianza
    if (i > periodo_quemado){
      Sigma <- (i-2)/(i-1)*Sigma + tcrossprod(cadena[i-1,] - cadena_media)/i
      cadena_media <- (i - 1)/i * cadena_media + cadena[i-1,] / i
      Sigma_act <- sd*Sigma + epsilon*diag(npar)
      S <- chol(Sigma_act)
    }
    
    # Proponer una nueva cadena
    theta_propuesta <- cadena[i-1, ] + c(tcrossprod(rnorm(npar), S))
  
    # Si algun parametro es menor que 0, no ejecutar el filtro de particulas
    # no aceptar cadena
    if (any(theta_ut(as.vector(theta_propuesta)) <= 0)){
     cadena[i, ] <- cadena[i-1, ]
     
    } else {
      
      # Restricciones
      theta_prop <- theta_ut(theta_propuesta)
      restriccion_1 <- (theta_prop[1] > 150  & theta_prop[1] < 250)
      restriccion_2 <- (theta_prop[2] > 20  & theta_prop[2] < 40)
      restriccion_3 <- (theta_prop[3] > 300  & theta_prop[3] < 370)
      restriccion_4 <- (theta_prop[4] > 30  & theta_prop[4] < 70)
      restriccion_5 <- theta_prop[5] < 1
      restriccion_6 <- theta_prop[6] < 5
      restricciones <- restriccion_1*restriccion_2*restriccion_3*restriccion_4*restriccion_5*restriccion_6
      
      # Aceptar o rechazar la cadena propuesta
      if (restricciones == 0) {
        MHRatio <- -Inf
        
      } else {
        
        # Actualizar la matriz de numeros aleatorios para los efectos aleatorios
        bmarrayprop_ea <- rho*bmarray_ea + sqrt(1-rho^2)*rnorm(n_ea, 0, 1)
        
        # Actualizar la matriz de numeros aleatorios para el filtro de particulas
        bmarrayprop_FP <- rho*bmarray_FP + sqrt(1-rho^2)*rnorm(n_FP, 0, 1)
        
        # Actualizar la matriz de numeros aleatorios para el remuestreo
        uprop <- pnorm(rho*qnorm(umat) + sqrt(1-rho^2)*rnorm((n-1)*M, 0, 1))
        
        # Generar efectos aleatorios
        # phi_1
        mu_phi_1 <- theta_propuesta[1]
        sigma_phi_1 <- exp(-0.5*theta_propuesta[2])
        
        # phi_2
        mu_phi_2 <- theta_propuesta[3]
        sigma_phi_2 <- exp(-0.5*theta_propuesta[4])
        
        phi_1_ea <- matrix(mu_phi_1 + sigma_phi_1*bmarrayprop_ea[,1], M, L)
        phi_2_ea <- matrix(mu_phi_2 + sigma_phi_2*bmarrayprop_ea[,2], M, L)
        
        # Si algun efecto aleatorio es negativo, reemplazarlo por la media
        phi_1_ea[phi_1_ea<=0] <- mu_phi_1
        phi_2_ea[phi_2_ea<=0] <- mu_phi_2
        
        # Calcular el logaritmo de la verosimilitud de cada efecto aleatorio L
        lv_ea <- array(0, dim = c(M, L))
        log_ver <- 0
        
        for (l in 1:L) {
          ea_l <- cbind(phi_1_ea[,l], phi_2_ea[,l])
          lv_ea[,l] <- fun_log_verosimilitud(theta_propuesta, ea_l,
                                             bmarrayprop_FP, umat)
        }
        
        # Calcular el logaritmo de la verosimilitud de todos los individuos
        log_ver <- sum(apply(lv_ea, 1, logSumExp) - log(L))
        lv[i] <- log_ver
        
        # Calcular el logaritmo de la densidad a posteriori (iteracion i)      
        log_post_nuevo <- log_ver + log_prior(theta_propuesta)
        
        # Ratio de aceptacion
        MHRatio <- min(0, log_post_nuevo - log_post)
        }
      
      # Aceptar o rechazar
      if (log(runif(1)) < MHRatio) {
        
        # Actualizar valores
        cadena[i, ] <- theta_propuesta
        log_post <- log_post_nuevo
        tasa_theta <- tasa_theta + 1
        bmarray_ea <- bmarrayprop_ea
        bmarray_FP <- bmarrayprop_FP
        umat <- uprop
        cat(sprintf("Aceptado: %s \n",
                    round(tasa_theta/(iteraciones+periodo_quemado)*100,2)))
        
      } else {
        cadena[i, ] <- cadena[i-1, ]
      
      }
    }
    
    # Calcular la media de la cadena (periodo de quemado)
    if (i == periodo_quemado){
      cadena_media <- colMeans(cadena[1:periodo_quemado,])
    }
  }
  
  # Tiempo final
  end.time <- Sys.time()
  
  # Tiempo del MCMC
  tiempo <- as.numeric(difftime(end.time, start.time, units = "secs"))
  
  return(list(cadena = cadena,
              tasa_theta = tasa_theta/(iteraciones+periodo_quemado),
              tiempo = tiempo,
              lv = lv,
              Sigma = Sigma))
}

