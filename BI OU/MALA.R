setwd("/cloud/project/5506643/Tesis Version final 2/Version final/BI OU")
drift_fun_b <- function(beta, b, alpha, x) {
  b <- matrix(b,2,2,T)
  (beta*b)%*%(alpha-x)
}
diffusion_fun <- function(sigma){
  diag(sigma^2)
}

SDE_EM_b <- function(x0, theta, noStep, b, bmarray) {
  # bmarray[d,noStep]
  beta_11 <- theta[1]
  beta_12 <- theta[2]
  beta_21 <- theta[3]
  beta_22 <- theta[4]
  beta <- matrix(c(beta_11,beta_12,
                   beta_21, beta_22), nrow = 2, ncol = 2, T)
  
  alpha_1 <- theta[5]
  alpha_2 <- theta[6]
  alpha <- c(alpha_1,alpha_2)
  
  sigma_1 <- theta[7]
  sigma_2 <- theta[8]
  sigma <- c(sigma_1,sigma_2)
  
  x <- x0
  
  delta_t <- 1/noStep
  
  for (i in 1:noStep) {
    x <- x + drift_fun_b(beta, b, alpha, x)*delta_t + t(chol(diffusion_fun(sigma)*delta_t))%*%c(bmarray[,i])
  }
  return(x)
}

log.norm.pdf <- function(X,MU,S) {
  X <- as.vector(X)
  MU <- as.vector(MU)
  
  k <- length(MU)
  det.S <- det(S)
  inv.S <- solve(S)
  quad.form <- X-MU
  out <- -0.5*k*log(2*pi) -0.5*log(det.S) -0.5*t(quad.form)%*%inv.S%*%(quad.form)
  return(out)
}

log.uniform.density <- function(theta,a,b){
  if (theta >= a & theta <= b){
    log_p <- -log(b-a)
  } else {
    log_p <- -Inf
  }
  return(log_p)
}


logPriorTheta <- function(theta){
  beta_11 <- theta[1]
  beta_12 <- theta[2]
  beta_21 <- theta[3]
  beta_22 <- theta[4]
  alpha_1 <- theta[5]
  alpha_2 <- theta[6]
  sigma_1 <- theta[7]
  sigma_2 <- theta[8]
  sigma_obs_1 <- theta[9]
  sigma_obs_2 <- theta[10]
  lp <- sum(log.uniform.density(beta_11,0,4),
            log.uniform.density(beta_12,0,4),
            log.uniform.density(beta_21,0,4),
            log.uniform.density(beta_22,0,4),
            log.uniform.density(alpha_1,-1,2),
            log.uniform.density(alpha_2,-1,2),
            log.uniform.density(sigma_1,0,0.75),
            log.uniform.density(sigma_2,0,0.75),
            log.uniform.density(sigma_obs_1,0,0.75),
            log.uniform.density(sigma_obs_1,0,0.75)
  )
  return(lp)
}
# Resampling --------------------------------------------------------------

#Systematic resampling
sysresamp2=function(wts,N,uni)
{
  vec=rep(0,N)
  wsum=sum(wts)
  k=1
  u=uni/N
  wsumtarg=u
  wsumcurr=wts[k]/wsum
  delta=1/N
  for(i in 1:N)
  {
    while (wsumcurr<wsumtarg)
    {
      k=k+1
      wsumcurr=wsumcurr+wts[k]/wsum
    }   
    vec[i]=k 
    wsumtarg=wsumtarg+delta
  }
  return(vec)
}

# Particle filter ---------------------------------------------------------
# for noIndividuals

particleFilter_b <- function(y, 
                             theta,
                             noParticles=20,
                             initialState=0,
                             noStep=20,
                             b,
                             bmarray,
                             umat) {
  # bmarray[noParticles,d,noStep,(n-1)]
  # umat[n-1]
  y <- y[-1,]
  T <- dim(y)[1]
  
  sigma_obs_1 <- theta[9]
  sigma_obs_2 <- theta[10]
  sigma_obs <- c(sigma_obs_1,sigma_obs_2) 
  Sigma.y <- diag(sigma_obs^2)
  
  particles <- array(0, dim = c(noParticles, T,2))
  logWeights <- array(0, dim = c(noParticles, T))
  logNormalizedWeights <- array(0, dim = c(noParticles, T))
  
  logLikelihood <- 0
  
  particles[,1,] <- initialState
  
  logWeights[,1] <- mvtnorm::dmvnorm(y[1,], initialState, Sigma.y, log = T)
  
  logNormalizedWeights[,1] <- -log(noParticles)
  
  for (t in 2:T) {
    
    for (i in 1:noParticles) {
      # Propagate
      particles[i,t,] <- SDE_EM_b(particles[i,t-1,], theta, noStep, b, bmarray[i,,,t-1])
      
      # Compute weights
      logWeights[i,t] <- mvtnorm::dmvnorm(y[t,], particles[i,t,], Sigma.y, log = T)
    }
    
    maxLogWeight <- max(logWeights[,t])
    logWeights[,t] <- logWeights[,t] - maxLogWeight
    
    logSumWeights <- matrixStats::logSumExp(logWeights[,t])
    
    logNormalizedWeights[,t] <- logWeights[,t] - logSumWeights
    
    n_eff <- exp(-logSumExp(2*logNormalizedWeights[,t]))
    if (n_eff < 0.5*noParticles){
      indices <- sysresamp2(exp(logWeights[,t]), noParticles,umat[t-1])
      particles[,t,]<- particles[indices,t,]
      logNormalizedWeights[,t] <- -log(noParticles)
    }
    
    # Estimate the log-likelihood
    predictiveLikelihood <- maxLogWeight + logSumWeights - log(noParticles)
    logLikelihood <- logLikelihood + predictiveLikelihood
    
  }
  
  logLikelihood
  
}

log.gamma.density <- function(theta, lambda, nu){
  out <- (lambda-1)*log(theta) - theta*nu
  return(out)
}

logPriorRe <- function(theta, hyper){
  n <- length(theta)
  out <- 0
  for (i in 1:n) {
    par <- hyper[i,]
    out <- out + log.gamma.density(theta[i], par[1], par[2])
  }
  return(out)
}  

log.likelihood <- function(y, alpha, beta) {
  # find n from the data (ie the number of observations)
  n <- length(y)
  
  # compute the sum of the observations and the sum of the log of the observations
  sum.y <- sum(y)
  sum.log.y <- sum(log(y))
  
  out <- n*alpha*log(beta) - n*lgamma(alpha) + (alpha-1)*sum.log.y - beta*sum.y
  return(out)
} 

log.gamma.density.deriv <- function(theta, lambda, nu){
  out <- (lambda-1)/theta - nu
  return(out)
}
grad.log.post <- function(y,
                          alpha, beta,
                          lambda.alpha, lambda.beta,
                          nu.alpha, nu.beta){
  # find n from the data (ie the number of observations)
  n <- length(y)
  
  # compute the sum of the observations and the sum of the log of the observations
  sum.y <- sum(y)
  sum.log.y <- sum(log(y))
  
  grad.log.like.alpha <- n*log(beta) - n*digamma(alpha) + sum.log.y
  grad.log.like.beta <- n*alpha/beta - sum.y
  grad.log.like <- c(grad.log.like.alpha, grad.log.like.beta)
  
  grad.log.prior.alpha <- log.gamma.density.deriv(alpha, lambda.alpha, nu.alpha)
  grad.log.prior.beta <- log.gamma.density.deriv(beta, lambda.beta, nu.beta)
  grad.log.prior <- c(grad.log.prior.alpha, grad.log.prior.beta)
  
  out <- grad.log.like + grad.log.prior
  return(out)
}


MALA.mixed.effects <- function(y,
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
                              mala.eps,
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
  
  scale.mala <- mala.eps^2/2^(1/3)
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
      
      # MALA
      mu.can <- pHyper + 0.5*scale.mala * c(Sigma.hyper[[k]] %*% grad.log.post(pRe, alpha, beta,
                                                                                  lambda.alpha, lambda.beta,
                                                                                  nu.alpha, nu.beta))
      pHyperCan <- mvtnorm::rmvnorm(1, mu.can, scale.mala*Sigma.hyper[[k]]) 
      
      log.qCan <- mvtnorm::dmvnorm(pHyperCan, mu.can, scale.mala*Sigma.hyper[[k]])
      
      alpha.can <- pHyperCan[1]
      beta.can <- pHyperCan[2]
      mu <- pHyperCan + 0.5*scale.mala * c(Sigma.hyper[[k]] %*% grad.log.post(pRe, alpha.can, beta.can,
                                                                               lambda.alpha, lambda.beta,
                                                                               nu.alpha, nu.beta))
      
      log.q <- mvtnorm::dmvnorm(pHyper, mu, scale.mala*Sigma.hyper[[k]])
      
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
      
      # Check whether the candidates are accepted.
      log.alpha <-logPostHyperCan - logPostHyper + log.q - log.qCan
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
Sigma.theta <- diag(c(0.001,0.001,
                      0.001,0.001,
                      0.001,0.001,
                      0.001,0.001,
                      0.001,0.001)
)

# Sigma.theta <- cov(res$ThetaChain)

f.theta <- 0.25
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
Sigma.hyper <- rep(list(diag(15, 2)),4)

HyperPrior <- rbind(c(1,0.05),
                    c(1,0.01),
                    c(1,0.01),
                    c(1,0.02))

mala.eps <- 0.25
res <- MALA.mixed.effects(y,
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
                         mala.eps)
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

save(res, file = "rescorrMALA.RData")
