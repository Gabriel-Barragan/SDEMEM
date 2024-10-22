HMC.mixed.effects <- function(y,
                              startTheta = NULL,
                              iterations  = 10000,
                              nBI = 0 , 
                              Sigma.theta = NULL,
                              f.theta = 1,
                              startState,
                              noParticles,
                              noStep,
                              consoleUpdates=100,
                              startRe,
                              Sigma.re,
                              f.re = 1,
                              startHyper,
                              Sigma.hyper,
                              HyperPrior,
                              eps,
                              L,
                              rho=0.999) {
  noObs <- dim(y)[1]-2
  noIndividuals <- dim(y)[3]
  nb <- length(startRe)
  npar <- length(startTheta)
  lChain <- iterations
  pReChain <- array(NA_real_, dim=c(lChain - nBI, nb, noIndividuals))
  pThetaChain <- matrix(NA_real_, nrow = lChain - nBI, ncol = npar)
  pHyperChain <- array(NA_real_, dim=c(lChain - nBI, nb, 2)) # 2 corresponds to alpha and beta
  
  pReMat <- array(startRe, dim = c(nb,noIndividuals))
  pHyperMat <- startHyper
  pTheta <- startTheta
  
  #********************************************************************************
  
  bmarray <- array(rnorm(noParticles*2*noStep*noObs*noIndividuals), 
                   dim = c(noParticles,2,noStep,noObs,noIndividuals)) 
  umat <- array(runif(noObs*noIndividuals), dim=c(noObs,noIndividuals))
  
  # First call to the model. Calculate likelihood
  logLike <- logLikeCan <- rep(NA_real_, noIndividuals)
  for (m in 1:noIndividuals) {
    pRe <- pReMat[,m]
    logLike[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pRe, bmarray[,,,,m],umat[,m])
  }
  accept.prob.re <- accept.re <- rep(0,noIndividuals)
  accept.prob.theta <- 0
  accept.prob.hyper <- accept.hyper <- rep(0,nb)
  
  #********************************************************************************
  
  # Define Variance-covariance matrix (vcovProp) for proposal generation an
  scalPropRe <- f.re * 2.4^2/nb # This f is the scaling factor tuned manually
  
  scalPropTheta <- f.theta * 2.4^2/npar # This f is the scaling factor tuned manually
  covParTheta <- scalPropTheta * Sigma.theta
  
  
  #********************************************************************************
  # Build up the chain. Candidates for the parameter values (candidatepValues)
  # are assumed to stem from a multivariate normal distribution (mvrnorm) with mean
  # at the current state and covariance given by scalProp*covPar.
  #-----
  # Start time
  start.time <- Sys.time()
  for (j in 1:lChain) {
    if (j%%consoleUpdates == 0) print(c(j,accept.prob.theta/(lChain-nBI)))
    
    # Update radom effects
    for (m in 1:noIndividuals) {
      pRe <- pReMat[,m]
      pReCan <- as.vector(mvtnorm::rmvnorm(1, pRe, scalPropRe*Sigma.re[[m]]))
      accept.re[m] <- 0
      if(sum(pReCan>0)==length(pReCan)){
        
        bmarrayprop=rho*bmarray[,,,,m]+sqrt(1-rho^2)*rnorm(noParticles*2*noStep*noObs,0,1) #update bm increments
        uprop=pnorm(rho*qnorm(umat[,m])+sqrt(1-rho^2)*rnorm(noObs)) #update uniforms for resampling step
        
        # Call log likelihood
        logLikeCan[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pReCan,bmarrayprop, uprop)
        
        # Call log prior
        logPrior.Re <- logPriorRe(pRe, pHyperMat)
        logPrior.ReCan <- logPriorRe(pReCan, pHyperMat)
        
        # Calculate log posterior
        logPost <- logLike[m] + logPrior.Re
        logPostCan <- logLikeCan[m] + logPrior.ReCan
        
        # Check whether the candidates are accepted.
        log.alpha <- logPostCan - logPost
        log.u <- log(runif(1))
        if (log.u < log.alpha) {
          logLike[m] <- logLikeCan[m]
          pRe <- pReCan
          accept.re[m] <-  1
          bmarray[,,,,m] <- bmarrayprop
          umat[,m] <- uprop
        }
      }
      pReMat[,m] <- pRe
    }
    
    # Update theta
    pThetaCan <- mvtnorm::rmvnorm(1, pTheta, covParTheta)
    accept.theta <- 0
    
    if(sum(pThetaCan[-c(5,6)]>0)==length(pThetaCan[-c(5,6)])){
      for (m in 1:m) {
        pRe <- pReMat[,m]
        # Call the model and calculate the likelihood
        logLikeCan[m] <- particleFilter_b(y[,,m],pThetaCan,noParticles,startState,noStep, pRe,bmarray[,,,,m], umat[,m])
      }
      
      # Call log prior
      logPrior.Theta <- logPriorTheta(pTheta)
      logPrior.ThetaCan <- logPriorTheta(pThetaCan)
      
      # Calculate log posterior
      logPost <- sum(logLike) + logPrior.Theta
      logPostCan <- sum(logLikeCan) + logPrior.ThetaCan
      
      # Check whether the candidates are accepted.
      log.alpha <- logPostCan - logPost
      log.u <- log(runif(1))
      if (log.u < log.alpha) {
        logLike <- logLikeCan
        pTheta <- pThetaCan
        accept.theta <-  1
      }  
    }
    
    # Update hyperparameters
    for (k in 1:nb) {
      pHyper <- pHyperMat[k,]
      pRe <- pReMat[k,]
      alpha <- pHyper[1]
      beta <- pHyper[2]
      
      lambda.alpha <- lambda.beta <- HyperPrior[k,1]
      nu.alpha <- nu.beta <- HyperPrior[k,2]
      accept.hyper[k] <- 0
      
      # Auxiliary variables
      P <- c(mvtnorm::rmvnorm(1, sigma = solve(Sigma.hyper[[k]])))
      
      # Kinetic energy at the beginning of the trajectory
      logK <- c(P %*% Sigma.hyper[[k]] %*% P / 2)
      
      # Make a half step for momentum at the beginning
      pHyperCan <- pHyper
      alpha.can <- pHyperCan[1]
      beta.can <- pHyperCan[2]
      
      Pnew <- P + eps * grad.log.post(pRe, alpha.can, beta.can,
                                      lambda.alpha, lambda.beta,
                                      nu.alpha, nu.beta) / 2
      
      # Alternate full steps for position and momentum
      for (l in 1:L) {
        # Make a full step for the position
        pHyperCan <- pHyperCan + eps * c(Sigma.hyper[[k]] %*% Pnew)
        
        # Make a full step for the momentum, except at the end of trajectory
        if(l!=L){
          alpha.can <- pHyperCan[1]
          beta.can <- pHyperCan[2]
          
          Pnew <- P + eps * grad.log.post(pRe, alpha.can, beta.can,
                                          lambda.alpha, lambda.beta,
                                          nu.alpha, nu.beta)
        }
      }
      # Make a half step for momentum at the end
      alpha.can <- pHyperCan[1]
      beta.can <- pHyperCan[2]
      Pnew <- P + eps * grad.log.post(pRe, alpha.can, beta.can,
                                      lambda.alpha, lambda.beta,
                                      nu.alpha, nu.beta) / 2
      
      # Negate momentum at end of trajectory to make the proposal symmetric
      Pnew <- -Pnew
      
      # Calculate potential energy at the end of trajectory
      
      # Calculate log likelihood
      logLikeHyper <-  log.likelihood(pRe,alpha,beta)
      logLikeHyperCan <-  log.likelihood(pRe,alpha.can,beta.can)
      
      # Calculate log prior
      logPriorAlpha <- log.gamma.density(alpha, lambda.alpha, nu.alpha)
      logPriorAlphaCan <- log.gamma.density(alpha.can, lambda.alpha, nu.alpha)
      
      logPriorBeta <- log.gamma.density(beta, lambda.beta, nu.beta)
      logPriorBetaCan <- log.gamma.density(beta.can, lambda.beta, nu.beta)
      
      # Calculate log posterior
      logPostHyper <- logLikeHyper + logPriorAlpha + logPriorBeta
      logPostHyperCan <- logLikeHyperCan + logPriorAlphaCan + logPriorBetaCan
      
      # Calculate potential energy at the end of trajectory
      logKnew <- c(Pnew %*% Sigma.hyper[[k]] %*% Pnew / 2)
      
      # Check whether the candidates are accepted.
      log.alpha <-logPostHyperCan - logPostHyper + logK - logKnew 
      log.u <- log(runif(1))
      if (log.u < log.alpha) {
        pHyper <- pHyperCan
        accept.hyper[k] <-  1
      }
      
      # if (sum(pHyperCan > 0)==length(pHyperCan)){
      # }
      pHyperMat[k,] <- pHyper
      
    }
    
    if (j > nBI) {
      pReChain[j-nBI,,] <- pReMat
      pThetaChain[j-nBI,] <- pTheta
      pHyperChain[j-nBI,,] <- pHyperMat
      
      accept.prob.re <- accept.prob.re + accept.re
      accept.prob.theta <- accept.prob.theta + accept.theta
      accept.prob.hyper <- accept.prob.hyper + accept.hyper
    }
  }
  # End time
  end.time <- Sys.time()
  
  # Tiempo del MCMC
  mcmc.time <- as.numeric(difftime(end.time, start.time, units = "secs"))
  
  accept.prob.re <- accept.prob.re/(lChain-nBI)
  accept.prob.theta <- accept.prob.theta/(lChain-nBI)
  accept.prob.hyper <- accept.prob.hyper/(lChain-nBI)
  
  list(ReChain = pReChain,
       ThetaChain = pThetaChain,
       HyperChain = pHyperChain,
       accept.prob.re = accept.prob.re,
       accept.prob.theta = accept.prob.theta, 
       accept.prob.hyper = accept.prob.hyper,
       mcmc.time = mcmc.time)
}

# Run ---------------------------------------------------------------------

startTheta <- c(2,2,
                1.5,1.5,
                0.5,1,
                0.25,0.25,
                0.25,0.25)
iterations <- 5000
nBI <- floor(iterations/2)
Sigma.theta <- diag(c(0.005,0.005,
                      0.005,0.005,
                      0.005,0.005,
                      0.0001,0.0001,
                      0.0001,0.0001)
)

# Sigma.theta <- cov(res$ThetaChain)

f.theta <- 0.5
startState <- c(2,2)
noParticles <- 20
noStep <- 20
consoleUpdates <- 100
startRe <- c(1,1,1,1)
Sigma.re <- rep(list(diag(1e-5,4)),dim(y)[3])
f.re <- 1
startHyper <- rbind(c(30,30),
                    c(90,90),
                    c(110,110),
                    c(50,50))
Sigma.hyper <- rep(list(diag(5, 2)),4)

HyperPrior <- rbind(c(1,0.05),
                    c(1,0.01),
                    c(1,0.01),
                    c(1,0.02))

eps <- 0.1
L <- 10
res <- HMC.mixed.effects(y,
                         startTheta,
                         iterations,
                         nBI, 
                         Sigma.theta,
                         f.theta,
                         startState,
                         noParticles,
                         noStep,
                         consoleUpdates,
                         startRe,
                         Sigma.re,
                         f.re,
                         startHyper,
                         Sigma.hyper, 
                         HyperPrior,
                         eps,
                         L)
res$mcmc.time

res$accept.prob.re
res$accept.prob.theta
res$accept.prob.hyper

parameterNamesRe <- c(expression(b[11]), expression(b[12]),
                      expression(b[21]), expression(b[22]))

par(mfrow=c(2,2))
for (k in 1:length(startRe)) {
  boxplot(res$ReChain[,k,],main="Random effects", xlab=parameterNamesRe[k],
          ylab="", col="pink")  
}


theta
(thhat <- colMeans(res$ThetaChain)) %>% round(.,2)
(thhatSD <- colSds(res$ThetaChain)) %>% round(.,2)

parameterNames <- c(expression(beta[11]), expression(beta[12]),
                    expression(beta[12]), expression(beta[22]),
                    expression(alpha[1]), expression(alpha[2]),
                    expression(sigma[1]), expression(sigma[2]),
                    expression(sigma[obs[1]]), expression(sigma[obs[2]]))

parameterColors <- c("royalblue", "royalblue",
                     "royalblue", "royalblue",
                     "lightcoral", "lightcoral",
                     "navajowhite","navajowhite",
                     "springgreen","springgreen")

parameterACFnames <- c(expression("ACF of " *beta[11]), expression("ACF of " *beta[12]),
                       expression("ACF of " *beta[12]), expression("ACF of " *beta[22]),
                       expression("ACF of " *alpha[1]), expression("ACF of " *alpha[2]),
                       expression("ACF of " *sigma[1]), expression("ACF of " *sigma[2]),
                       expression("ACF of " *sigma[obs[1]]) ,expression("ACF of " *sigma[obs[2]])
)


parameterScales <- cbind(apply(res$ThetaChain,2,min), 
                         apply(res$ThetaChain,2,max))

iact <- c()

par(mfrow=c(3,1))
for (k in 1:10) {
  
  # Histogram of the posterior
  hist(
    res$ThetaChain[, k],
    breaks = floor(sqrt(iterations-nBI)),
    col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25),
    border = NA,
    xlab = parameterNames[k],
    ylab = "posterior estimate",
    main = "",
    xlim = parameterScales[k,],
    freq = FALSE
  )
  
  # Add lines for the kernel density estimate of the posterior
  kde <- density(res$ThetaChain[, k], kernel = "e",
                 from = parameterScales[k, 1], to = parameterScales[k, 2])
  lines(kde, lwd = 2, col = parameterColors[k])
  
  # Plot the estimate of the posterior mean
  abline(v = thhat[k], lwd = 2, lty = "dotted")
  
  # Plot the true value
  abline(v = theta[k], lwd = 2, lty = "dotted", col="red")
  
  # Plot Initial chain
  abline(v = startTheta[k], lwd = 2, lty = "dotted", col="blue")
  
  # Add lines for prior
  prior_grid <- seq(parameterScales[k, 1], parameterScales[k, 2], 0.01)
  if (k==1) {prior_values = dunif(prior_grid, 0, 4)}
  if (k==2) {prior_values = dunif(prior_grid, 0, 4)}
  if (k==3) {prior_values = dunif(prior_grid, 0, 4)}
  if (k==4) {prior_values = dunif(prior_grid, 0, 4)}
  if (k==5) {prior_values = dunif(prior_grid, -1, 2)}
  if (k==6) {prior_values = dunif(prior_grid, -1, 2)}
  if (k==7) {prior_values = dunif(prior_grid, 0, 0.75)}
  if (k==8) {prior_values = dunif(prior_grid, 0, 0.75)}
  if (k==9) {prior_values = dunif(prior_grid, 0, 0.75)}
  if (k==10) {prior_values = dunif(prior_grid, 0, 0.75)}
  lines(prior_grid, prior_values, col = "darkgrey")
  
  # Plot trace of the Markov chain
  plot(
    res$ThetaChain[, k],
    col = parameterColors[k],
    type = "l",
    xlab = "iteration",
    ylab = parameterNames[k],
    ylim = parameterScales[k,],
    bty = "n"
  )
  polygon(
    c(1:(iterations-nBI), rev(1:(iterations-nBI))),
    c(res$ThetaChain[, k], rep(-1, iterations-nBI)),
    border = NA,
    col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
  )
  abline(h = thhat[k], lwd = 2, lty = "dotted")
  # Plot the true value
  abline(h = theta[k], lwd = 2, lty = "dotted", col="red")
  # Plot Initial chain
  abline(h = startTheta[k], lwd = 2, lty = "dotted", col="blue")
  
  # Plot the autocorrelation function
  acf_res <- acf(res$ThetaChain[, k], plot = FALSE, lag.max = 100)
  plot(
    acf_res$lag,
    acf_res$acf,
    col = parameterColors[k],
    type = "l",
    xlab = "iteration",
    ylab = parameterACFnames[k],
    lwd = 2,
    ylim = c(-0.2, 1),
    bty = "n"
  )
  polygon(
    c(acf_res$lag, rev(acf_res$lag)),
    c(acf_res$acf, rep(0, length(acf_res$lag))),
    border = NA,
    col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
  )
  abline(h = 1.96 / sqrt(iterations - nBI), lty = "dotted",lwd=2)
  abline(h = -1.96 / sqrt(iterations - nBI), lty = "dotted",lwd=2)
  
  iact <- c(iact, 1 + 2 * sum(acf_res$acf))
}

iact

parameterNamesHyper <- c(expression(psi[11]), expression(psi[12]),
                         expression(psi[21]), expression(psi[22]))

par(mfrow=c(2,2))
for (k in 1:length(startRe)) {
  boxplot(res$HyperChain[,k,],main="Hyperparameters",
          xlab=parameterNamesHyper[k],
          ylab="", col=c("red","blue"))  
}

save(res, file = "rescorrHamiltonian.RData")
