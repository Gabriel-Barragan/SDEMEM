library(ggplot2)


# OU model (unidimensional)
n <- 25
x0 <- 1
beta <- 0.5
alpha <- 1
sigma <- 1.5
sigma_obs <- 1.5
true.values <- c(beta, alpha,sigma,sigma_obs)


drift_fun <- function(beta, alpha, x) {
  beta*(alpha-x)
}

diffusion_fun <- function(sigma){
  sigma^2
}
generate.data <- function(theta, noObservations, initialStates){
  beta <- theta[1]
  alpha <- theta[2]
  sigma <- theta[3]
  sigma_obs <- theta[4]
  
  delta_t <- 1/noObservations
  time_seq <- seq(0,1, length.out=n+1)
  x <- y <- array(NA, dim=c(n+1))
  x[1] <- x0
  for (i in 1:n) {
    x[i+1] <- x[i]+drift_fun(beta, alpha, x[i])*delta_t + sqrt(diffusion_fun(sigma)*delta_t)*rnorm(1)  
    y[i+1] <- x[i+1]+sigma_obs*rnorm(1)
  }
  return(list(x=x, y=y, t=time_seq))
}

data <- generate.data(true.values,n,x0)
x <- data$x
y <- data$y
t <- data$t

# Simulation --------------------------------------------------------------

SDE_EM <- function(x0, theta, noStep) {
  x <- x0
  beta <- theta[1]
  alpha <- theta[2]
  sigma <- theta[3]
  sigma_obs <- theta[4]
  delta_t <- 1/noStep
  
  for (i in 1:noStep) {
    x <- x + drift_fun(beta, alpha, x)*delta_t + sqrt(diffusion_fun(sigma)*delta_t)*rnorm(1)
  }
  return(x)
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

particleFilter <- function(y, 
                           theta,
                           noParticles=20,
                           initialState=0,
                           noStep=20) {
  y <- y[-1]
  T <- length(y)
  sigma_obs <- theta[4]
  particles <- array(0, dim = c(noParticles, T))
  logWeights <- array(0, dim = c(noParticles, T))
  logNormalizedWeights <- array(0, dim = c(noParticles, T))

  logLikelihood <- 0

  particles[,1] <- initialState

  logWeights[,1] <- dnorm(y[1], initialState, sigma_obs, log = T)

  logNormalizedWeights[,1] <- -log(noParticles)

  for (t in 2:T) {

    for (i in 1:noParticles) {
      # Propagate
      particles[i,t] <- SDE_EM(particles[i,t-1], theta, noStep)

      # Compute weights
      logWeights[i,t] <- dnorm(y[t], particles[i,t], sigma_obs, log=T)
    }

    maxLogWeight <- max(logWeights[,t])
    logWeights[,t] <- logWeights[,t] - maxLogWeight

    logSumWeights <- matrixStats::logSumExp(logWeights[,t])

    logNormalizedWeights[,t] <- logWeights[,t] - logSumWeights

    n_eff <- exp(-logSumExp(2*logNormalizedWeights[,t]))
    if (n_eff < 0.5*noParticles){
      indices <- sysresamp2(exp(logWeights[,t]), noParticles, runif(1))
      particles[,t]<- particles[indices,t]
      logNormalizedWeights[,t] <- -log(noParticles)
    }

    # Estimate the log-likelihood
    predictiveLikelihood <- maxLogWeight + logSumWeights - log(noParticles)
    logLikelihood <- logLikelihood + predictiveLikelihood

  }

  list(logLikelihood = logLikelihood,
       particles = particles)

}

res.PF <- particleFilter(y, theta = c(0.4, 1.25, 1, 2))
x.particles <- res.PF$particles
x.filtered <- apply(x.particles, 2, mean)

#xhatSD <- apply(x.particles, 2, sd)

# xhat_upper <- x.filtered + 1.96 * xhatSD
# xhat_lower <- x.filtered - 1.96 * xhatSD
# 
# xhat_upper <- c(NA,xhat_upper)
# xhat_lower <- c(NA,xhat_lower)

x_quantiles <- apply(x.particles, 2, function(x) quantile(x, probs = c(0.05, 0.95)))
x_quantiles <- cbind(c(NA,NA),x_quantiles)

plot(t,t,bty="n", type = "n", main="Simulation", xlab="Time", ylab = "X",
     ylim = c(min(x,y,x.filtered,na.rm = T), max(x,y,x.filtered,na.rm = T)))
lines(t,x,col="blue")
lines(t,c(NA,x.filtered), col="darkgreen", lty="dotted", lwd=2)
points(t, y, col="red")

# polygon(
#   c(t, rev(t)),
#   c(xhat_upper, rev(xhat_lower)),
#   border = NA,
#   col = rgb(t(col2rgb("green")) / 256, alpha = 0.25)
# )

polygon(
  c(t, rev(t)),
  c(x_quantiles[2,], rev(x_quantiles[1,])),
  border = NA,
  col = rgb(t(col2rgb("green")) / 256, alpha = 0.25)
)


# MCMC --------------------------------------------------------------------

log.uniform.density <- function(theta,a,b){
  if (theta >= a & theta <= b){
    log_p <- -log(b-a)
  } else {
    log_p <- -Inf
  }
  return(log_p)
}

logPriorTheta <- function(theta){
  beta <- theta[1]
  alpha <- theta[2]
  sigma <- theta[3]
  sigma_obs <- theta[4]
  lp <- sum(log.uniform.density(beta,0,2),
            log.uniform.density(alpha,-2,2),
            log.uniform.density(sigma,0,3),
            log.uniform.density(sigma_obs,0,3))
  return(lp)
}

Metropolis <- function(y,
                       startValue = NULL,
                       iterations  = 10000,
                       nBI = 0 , 
                       Sigma = NULL,
                       f = 1, 
                       startState,
                       noParticles,
                       noStep,
                       consoleUpdates=100) {
  

  pValues = startValue
  lChain = iterations
  
  npar <- length(pValues)
  pChain <- matrix(NA_real_, nrow = lChain - nBI, ncol = npar)
  
  #********************************************************************************
  
  # First call to the model. Calculate likelihood and prior
  logLike <- particleFilter(y,pValues,noParticles,startState,noStep)$logLikelihood
  accept.prob <- 0
  
  #********************************************************************************
  
  # Define Variance-covariance matrix (vcovProp) for proposal generation an
  
  scalProp <- f * 2.4^2/npar # This f is the scaling factor tuned manually
  covPar <- scalProp * Sigma
  
  #********************************************************************************
  # Build up the chain. Candidates for the parameter values (candidatepValues)
  # are assumed to stem from a multivariate normal distribution (mvrnorm) with mean
  # at the current state and covariance given by scalProp*covPar.
  #-----
  
  for (j in 1:lChain) {
    if (j%%consoleUpdates == 0) print(c(j,logLike, accept.prob/(lChain-nBI)))
    candidatepValues <- as.vector(mvtnorm::rmvnorm(1, pValues, covPar))
    accept <- 0
    
    if(sum(candidatepValues[-2]>0)==length(candidatepValues[-2])){
      # Call the model and calculate the likelihood
      logLikeCan <- particleFilter(y,candidatepValues,noParticles,startState,noStep)$logLikelihood
      
      # Call log prior
      logPrior <- logPriorTheta(pValues)
      logPriorCan <- logPriorTheta(candidatepValues)
      
      # Calculate log posterior
      logPost <- logLike + logPrior
      logPostCan <- logLikeCan + logPriorCan
      
      # Check whether the candidates are accepted.
      alpha <- min(exp(logPostCan - logPost), 1)
      if (runif(1) < alpha) {
        logLike <- logLikeCan
        pValues <- candidatepValues
        accept <-  1
      }  
    }
    
    if (j > nBI) {
      pChain[j-nBI,] <- pValues
      accept.prob <- accept.prob + accept
    }
  }
  accept.prob <- accept.prob/(lChain-nBI)
  list(Draws = pChain, accept.prob = accept.prob)
}

startValue <- c(0.4, 1.25, 1, 2)
iterations <- 5000
nBI <- floor(iterations/2)
# Sigma <- diag(0.1^2, length(startValue))
(Sigma <- cov(res$Draws))
f <- 1
startState <- 2
noParticles <- 20
noStep <- 20
res <- Metropolis(y, startValue = startValue,iterations = iterations, 
                  nBI = nBI,Sigma = Sigma,f = f,startState = startState,
                  noParticles = noParticles,noStep=noStep)
res$accept.prob
true.values
(thhat <- colMeans(res$Draws)) %>% round(.,2)
(thhatSD <- colSds(res$Draws)) %>% round(.,2)
#---------------------------------------------------------------------------
# Parameter posteriors
#---------------------------------------------------------------------------
parameterNames <- c(expression(beta),
                    expression(alpha),
                    expression(sigma),
                    expression(sigma[obs]))

parameterColors <- c("azure3",
                     "lightyellow3",
                     "navajowhite",
                     "green")

parameterACFnames <- c(expression("ACF of " *beta),
                       expression("ACF of " *alpha),
                       expression("ACF of " *sigma),
                       expression("ACF of " *sigma[obs])
                       )

parameterScales <- cbind(apply(res$Draws,2,min), apply(res$Draws,2,max))

iact <- c()

par(mfrow=c(3,1))
for (k in 1:4) {
  
  # Histogram of the posterior
  hist(
    res$Draws[, k],
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
  kde <- density(res$Draws[, k], kernel = "e",
                 from = parameterScales[k, 1], to = parameterScales[k, 2])
  lines(kde, lwd = 2, col = parameterColors[k])
  
  # Plot the estimate of the posterior mean
  abline(v = thhat[k], lwd = 2, lty = "dotted")
  
  # Plot the true value
  abline(v = true.values[k], lwd = 2, lty = "dotted", col="red")
  
  # Plot Initial chain
  abline(v = startValue[k], lwd = 2, lty = "dotted", col="blue")
  
  # Add lines for prior
  prior_grid <- seq(parameterScales[k, 1], parameterScales[k, 2], 0.01)
  if (k==1) {prior_values = dunif(prior_grid, 0, 2)}
  if (k==2) {prior_values = dunif(prior_grid, -2, 2)}
  if (k==3) {prior_values = dunif(prior_grid, 0, 3)}
  if (k==4) {prior_values = dunif(prior_grid, 0, 3)}
  lines(prior_grid, prior_values, col = "darkgrey")
  
  # Plot trace of the Markov chain
  plot(
    res$Draws[, k],
    col = parameterColors[k],
    type = "l",
    xlab = "iteration",
    ylab = parameterNames[k],
    ylim = parameterScales[k,],
    bty = "n"
  )
  polygon(
    c(1:(iterations-nBI), rev(1:(iterations-nBI))),
    c(res$Draws[, k], rep(-1, iterations-nBI)),
    border = NA,
    col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
  )
  abline(h = thhat[k], lwd = 2, lty = "dotted")
  # Plot the true value
  abline(h = true.values[k], lwd = 2, lty = "dotted", col="red")
  # Plot Initial chain
  abline(h = startValue[k], lwd = 2, lty = "dotted", col="blue")
  
  # Plot the autocorrelation function
  acf_res <- acf(res$Draws[, k], plot = FALSE, lag.max = 100)
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


# Adaptive Metropolis -----------------------------------------------------

AM <- function(y,
               startValue = NULL,
               iterations  = 10000,
               nBI = 0 , 
               Sigma = NULL,
               f = 1, 
               startState,
               noParticles,
               noStep,
               consoleUpdates=100,
               eps=1e-6) {
  
  pValues = startValue
  lChain = iterations
  
  noAdapt <- 500
  n.iter <- lChain + noAdapt
  
  npar <- length(pValues)
  pChain <- matrix(NA_real_, nrow = n.iter - nBI, ncol = npar)
  
  #********************************************************************************
  
  # First call to the model. Calculate likelihood and prior
  logLike <- particleFilter(y,pValues,noParticles,startState,noStep)$logLikelihood
  accept.prob <- 0
  
  #********************************************************************************
  
  # Define Variance-covariance matrix (vcovProp) for proposal generation an
  epsDiag <- eps * diag(npar)
  scalProp <- f * 2.4^2/npar # This f is the scaling factor tuned manually
  covPar <- scalProp * Sigma
  
  #********************************************************************************
  # Build up the chain. Candidates for the parameter values (candidatepValues)
  # are assumed to stem from a multivariate normal distribution (mvrnorm) with mean
  # at the current state and covariance given by scalProp*covPar.
  #-----
  
  for (j in 1:n.iter) {
    if (j%%consoleUpdates == 0) print(c(j,logLike))
    candidatepValues <- as.vector(mvtnorm::rmvnorm(1, pValues, covPar))
    accept <- 0
    
    if(sum(candidatepValues[-2]>0)==length(candidatepValues[-2])){
      # Call the model and calculate the likelihood
      logLikeCan <- particleFilter(y,candidatepValues,noParticles,startState,noStep)$logLikelihood
      
      # Call log prior
      logPrior <- logPriorTheta(pValues)
      logPriorCan <- logPriorTheta(candidatepValues)
      
      # Calculate log posterior
      logPost <- logLike + logPrior
      logPostCan <- logLikeCan + logPriorCan
      
      # Check whether the candidates are accepted.
      alpha <- min(exp(logPostCan - logPost), 1)
      if (runif(1) < alpha) {
        logLike <- logLikeCan
        pValues <- candidatepValues
        accept <-  1
      }  
    }
    
    if (j > nBI) {
      pChain[j-nBI,] <- pValues
      accept.prob <- accept.prob + accept
    }
    if (j == (nBI + noAdapt)) {
      avePar <- apply(pChain[1:noAdapt,], 2, mean)
      covPar <- scalProp * (cov(pChain[1:noAdapt,]) + epsDiag)
    }
    if (j > (nBI + noAdapt)) {
      accept.prob <- accept.prob + accept
      t <- j - nBI
      avePar_new <- as.vector(((t-1) * avePar + pValues) / t)
      covPar_new <- ((t-2) * covPar + scalProp * ((t-1) * (avePar %o% avePar) - t * (avePar_new %o% avePar_new) + (pValues %o% pValues)) + epsDiag) / (t-1)
      avePar <- avePar_new
      covPar <- covPar_new
    }
    
  }
  accept.prob = accept.prob/(lChain-nBI)
  list(Draws = pChain[(noAdapt+1):(n.iter-nBI),], accept.prob = accept.prob)
}

startValue <- c(0.4, 1.25, 1, 2)
iterations <- 5000
nBI <- floor(iterations/2)
# Sigma <- diag(0.1^2, length(startValue))
(Sigma <- cov(res$Draws))
f <- 1
startState <- 2
noParticles <- 20
noStep <- 20
res <- AM(y, startValue = startValue,iterations = iterations, 
          nBI = nBI,Sigma = Sigma,f = f,startState = startState,
          noParticles = noParticles,noStep=noStep)
res$accept.prob
true.values
(thhat <- colMeans(res$Draws)) %>% round(.,2)
(thhatSD <- colSds(res$Draws)) %>% round(.,2)
#---------------------------------------------------------------------------
# Parameter posteriors
#---------------------------------------------------------------------------
parameterNames <- c(expression(beta),
                    expression(alpha),
                    expression(sigma),
                    expression(sigma[obs]))

parameterColors <- c("azure3",
                     "lightyellow3",
                     "navajowhite",
                     "green")

parameterACFnames <- c(expression("ACF of " *beta),
                       expression("ACF of " *alpha),
                       expression("ACF of " *sigma),
                       expression("ACF of " *sigma[obs])
)

parameterScales <- cbind(apply(res$Draws,2,min), apply(res$Draws,2,max))

iact <- c()

par(mfrow=c(3,1))
for (k in 1:4) {
  
  # Histogram of the posterior
  hist(
    res$Draws[, k],
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
  kde <- density(res$Draws[, k], kernel = "e",
                 from = parameterScales[k, 1], to = parameterScales[k, 2])
  lines(kde, lwd = 2, col = parameterColors[k])
  
  # Plot the estimate of the posterior mean
  abline(v = thhat[k], lwd = 2, lty = "dotted")
  
  # Plot the true value
  abline(v = true.values[k], lwd = 2, lty = "dotted", col="red")
  
  # Plot Initial chain
  abline(v = startValue[k], lwd = 2, lty = "dotted", col="blue")
  
  # Add lines for prior
  prior_grid <- seq(parameterScales[k, 1], parameterScales[k, 2], 0.01)
  if (k==1) {prior_values = dunif(prior_grid, 0, 2)}
  if (k==2) {prior_values = dunif(prior_grid, -2, 2)}
  if (k==3) {prior_values = dunif(prior_grid, 0, 3)}
  if (k==4) {prior_values = dunif(prior_grid, 0, 3)}
  lines(prior_grid, prior_values, col = "darkgrey")
  
  # Plot trace of the Markov chain
  plot(
    res$Draws[, k],
    col = parameterColors[k],
    type = "l",
    xlab = "iteration",
    ylab = parameterNames[k],
    ylim = parameterScales[k,],
    bty = "n"
  )
  polygon(
    c(1:(iterations-nBI), rev(1:(iterations-nBI))),
    c(res$Draws[, k], rep(-1, iterations-nBI)),
    border = NA,
    col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
  )
  abline(h = thhat[k], lwd = 2, lty = "dotted")
  # Plot the true value
  abline(h = true.values[k], lwd = 2, lty = "dotted", col="red")
  # Plot Initial chain
  abline(h = startValue[k], lwd = 2, lty = "dotted", col="blue")
  
  # Plot the autocorrelation function
  acf_res <- acf(res$Draws[, k], plot = FALSE, lag.max = 100)
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

# Mixed effects -----------------------------------------------------------
setwd("/cloud/project/5506643/Tesis Version final 2/Version final/BI OU")
drift_fun_b <- function(beta, b, alpha, x) {
  b <- matrix(b,2,2,T)
  (beta*b)%*%(alpha-x)
}
diffusion_fun <- function(sigma){
  diag(sigma^2)
}

generate.data.individuals <- function(theta=c(3, 2.5,
                                              1.8, 2,
                                              1, 1.5, 
                                              0.3, 0.5,
                                              1, 1.5),
                                      noObservations=10,
                                      trueInitialState=c(1,1),
                                      noIndividuals=50,
                                      populationEffects=c(45,100,100,25)){
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
  
  sigma_obs_1 <- theta[9]
  sigma_obs_2 <- theta[10]
  sigma_obs <- c(sigma_obs_1,sigma_obs_2) 
  
  delta_t <- 1/noObservations
  time_seq <- seq(0,1, length.out=noObservations+1)
  x <- y <- array(NA, dim=c(noObservations+1,2,noIndividuals))
  x[1,,] <- initialStates
  
  p <- length(populationEffects)
  b <- array(NA_real_,dim=c(p,noIndividuals))
  for (k in 1:p) {
    b[k,] <- rgamma(noIndividuals, populationEffects[k], populationEffects[k])
  }
  for (i in 1:noObservations) {
    for (m in 1:noIndividuals) {
      x[i+1,,m] <- x[i,,m]+drift_fun_b(beta, b[,m], alpha, x[i,,m])*delta_t + t(chol(diffusion_fun(sigma)*delta_t))%*%rnorm(2)  
      y[i+1,,m] <- x[i+1,,m]+diag(sigma_obs)%*%rnorm(2)
    }
  }
  return(list(x=x, y=y, t=time_seq))
}

# save(list=c("t","x","y","theta","populationEffects"), file = "Datos.RData")

noIndividuals <- 5
populationEffects=c(45,100,100,25)
theta=c(3, 2.5,
        1.8, 2,
        1, 1.5, 
        0.3, 0.5,
        0.1, 0.5)
trueInitialState <- c(1,1)
noObservations <- 10
data <- generate.data.individuals(theta,noObservations,trueInitialState,noIndividuals,populationEffects)
x <- data$x
y <- data$y
t <- data$t

par(mfrow=c(2,1))
plot(t,t,bty="n", type = "n", main="Simulation", xlab="Time", ylab = expression(X[1]),
     ylim = c(min(x[,1,],y[,1,], na.rm = T), 
              max(x[,1,],y[,1,], na.rm = T)))

for (m in 1:noIndividuals) {
  lines(t,x[,1,m],col="blue")
  points(t, y[,1,m], col="red")  
}

plot(t,t,bty="n", type = "n", main="Simulation", xlab="Time", ylab = expression(X[2]),
     ylim = c(min(x[,2,],y[,2,], na.rm = T), 
              max(x[,2,],y[,2,], na.rm = T)))

for (m in 1:noIndividuals) {
  lines(t,x[,2,m],col="blue")
  points(t, y[,2,m], col="red")  
}

SDE_EM_b <- function(x0, theta, noStep, b) {
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
    x <- x + drift_fun_b(beta, b, alpha, x)*delta_t + t(chol(diffusion_fun(sigma)*delta_t))%*%rnorm(2)
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
                           b) {
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
      particles[i,t,] <- SDE_EM_b(particles[i,t-1,], theta, noStep, b)
      
      # Compute weights
      logWeights[i,t] <- mvtnorm::dmvnorm(y[t,], particles[i,t,], Sigma.y, log = T)
    }
    
    maxLogWeight <- max(logWeights[,t])
    logWeights[,t] <- logWeights[,t] - maxLogWeight
    
    logSumWeights <- matrixStats::logSumExp(logWeights[,t])
    
    logNormalizedWeights[,t] <- logWeights[,t] - logSumWeights
    
    n_eff <- exp(-logSumExp(2*logNormalizedWeights[,t]))
    if (n_eff < 0.5*noParticles){
      indices <- sysresamp2(exp(logWeights[,t]), noParticles, runif(1))
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

Metropolis.mixed.effects <- function(y,
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
                       f.hype = 1, 
                       HyperPrior) {
  
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
  
  # First call to the model. Calculate likelihood
  logLike <- logLikeCan <- rep(NA_real_, noIndividuals)
  for (m in 1:noIndividuals) {
    pRe <- pReMat[,m]
    logLike[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pRe)
  }
  accept.prob.re <- accept.re <- rep(0,noIndividuals)
  accept.prob.theta <- 0
  accept.prob.hyper <- accept.hyper <- rep(0,nb)
  
  #********************************************************************************
  
  # Define Variance-covariance matrix (vcovProp) for proposal generation an
  scalPropRe <- f.re * 2.4^2/nb # This f is the scaling factor tuned manually
  
  scalPropTheta <- f.theta * 2.4^2/npar # This f is the scaling factor tuned manually
  covParTheta <- scalPropTheta * Sigma.theta
  
  scalPropHyper <- f.hyper * 2.4^2/2 # This f is the scaling factor tuned manually
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
        # Call log likelihood
        logLikeCan[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pReCan)

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
        logLikeCan[m] <- particleFilter_b(y[,,m],pThetaCan,noParticles,startState,noStep, pRe)
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
      pHyperCan <- as.vector(mvtnorm::rmvnorm(1, pHyper, scalPropHyper*Sigma.hyper[[k]]))
      accept.hyper[k] <- 0
      if (sum(pHyperCan > 0)==length(pHyperCan)){
        pRe <- pReMat[k,]
        alpha <- pHyper[1]
        beta <- pHyper[2]
        alpha.can <- pHyperCan[1]
        beta.can <- pHyperCan[2]
        lambda.alpha <- lambda.beta <- HyperPrior[k,1]
        nu.alpha <- nu.beta <- HyperPrior[k,2]

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
        log.alpha <- logPostHyperCan - logPostHyper
        log.u <- log(runif(1))
        if (log.u < log.alpha) {
          pHyper <- pHyperCan
          accept.hyper[k] <-  1
        }
      }
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
f.hyper <- 1
HyperPrior <- rbind(c(1,0.05),
                    c(1,0.01),
                    c(1,0.01),
                    c(1,0.02))
res <- Metropolis.mixed.effects(y,
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
                                f.hype, 
                                HyperPrior)
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

# save(res, file = "resHMC.RData")


# AM ----------------------------------------------------------------------

AM.mixed.effects <- function(y,
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
                                     f.hype = 1, 
                                     HyperPrior,
                             eps=1e-9) {
  
  noIndividuals <- dim(y)[3]
  nb <- length(startRe)
  npar <- length(startTheta)
  
  lChain <- iterations
  noAdapt <- 5
  n.iter <- lChain + noAdapt  
  
  pReChain <- array(NA_real_, dim=c(n.iter - nBI, nb, noIndividuals))
  pThetaChain <- matrix(NA_real_, nrow = n.iter - nBI, ncol = npar)
  pHyperChain <- array(NA_real_, dim=c(n.iter - nBI, nb, 2)) # 2 corresponds to alpha and beta
  
  pReMat <- array(startRe, dim = c(nb,noIndividuals))
  pHyperMat <- startHyper
  pTheta <- startTheta
  
  #********************************************************************************
  
  # First call to the model. Calculate likelihood
  logLike <- logLikeCan <- rep(NA_real_, noIndividuals)
  for (m in 1:noIndividuals) {
    pRe <- pReMat[,m]
    logLike[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pRe)
  }
  accept.prob.re <- accept.re <- rep(0,noIndividuals)
  accept.prob.theta <- 0
  accept.prob.hyper <- accept.hyper <- rep(0,nb)
  
  #********************************************************************************
  
  # Define Variance-covariance matrix (vcovProp) for proposal generation an
  epsDiagRe <- eps * diag(nb)
  epsDiagTheta <- eps * diag(npar)
  epsDiagHyper <- eps * diag(2)
  
  scalPropRe <- f.re * 2.4^2/nb # This f is the scaling factor tuned manually
  
  scalPropTheta <- f.theta * 2.4^2/npar # This f is the scaling factor tuned manually
  covParTheta <- scalPropTheta * Sigma.theta
  
  scalPropHyper <- f.hyper * 2.4^2/2 # This f is the scaling factor tuned manually
  
  #********************************************************************************
  # Build up the chain. Candidates for the parameter values (candidatepValues)
  # are assumed to stem from a multivariate normal distribution (mvrnorm) with mean
  # at the current state and covariance given by scalProp*covPar.
  #-----
  # Start time
  start.time <- Sys.time()
  for (j in 1:n.iter) {
    if (j%%consoleUpdates == 0) print(c(j))
    
    # Update radom effects
    for (m in 1:noIndividuals) {
      pRe <- pReMat[,m]
      if(j<=(nBI + noAdapt)){
        covParRe <- scalPropRe*Sigma.re[[m]]
      } else{
        covParRe <- covParReList[[m]]
      }
      pReCan <- as.vector(mvtnorm::rmvnorm(1, pRe, covParRe))
      accept.re[m] <- 0
      if(sum(pReCan>0)==length(pReCan)){
        # Call log likelihood
        logLikeCan[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pReCan)
        
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
        logLikeCan[m] <- particleFilter_b(y[,,m],pThetaCan,noParticles,startState,noStep, pRe)
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
      
      if(j<=(nBI + noAdapt)){
        covParHyper <- scalPropHyper*Sigma.hyper[[k]]
      } else{
        covParHyper <- covParHyperList[[k]]
      }
      
      pHyperCan <- as.vector(mvtnorm::rmvnorm(1, pHyper, covParHyper))
      accept.hyper[k] <- 0
      if (sum(pHyperCan > 0)==length(pHyperCan)){
        pRe <- pReMat[k,]
        alpha <- pHyper[1]
        beta <- pHyper[2]
        alpha.can <- pHyperCan[1]
        beta.can <- pHyperCan[2]
        lambda.alpha <- lambda.beta <- HyperPrior[k,1]
        nu.alpha <- nu.beta <- HyperPrior[k,2]
        
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
        log.alpha <- logPostHyperCan - logPostHyper
        log.u <- log(runif(1))
        if (log.u < log.alpha) {
          pHyper <- pHyperCan
          accept.hyper[k] <-  1
        }
      }
      pHyperMat[k,] <- pHyper
      
    }
    
    if (j > nBI) {
      pReChain[j-nBI,,] <- pReMat
      pThetaChain[j-nBI,] <- pTheta
      pHyperChain[j-nBI,,] <- pHyperMat
    }
    if (j == (nBI + noAdapt)){
      aveParRe <- aveParRe_new <- matrix(NA_real_, nb,noIndividuals)
      aveParHyper <- aveParHyper_new <- matrix(NA_real_, 2,nb)
      covParReList <- covParReList_new <- rep(list(diag(0,nb)),noIndividuals)
      covParHyperList <- covParHyperList_new <- rep(list(diag(0,2)),nb)
      for (m in 1:noIndividuals) {
        aveParRe[,m] <- apply(pReChain[1:noAdapt,,m],2, mean)
        covParReList[[m]] <- scalPropRe * (cov(pReChain[1:noAdapt,,m]) + epsDiagRe)
      }
      
      aveParTheta <- apply(pThetaChain[1:noAdapt,],2,mean)
      covParTheta <- scalPropTheta * (cov(pThetaChain[1:noAdapt,]) + epsDiagTheta)
      
      for (k in 1:nb) {
        aveParHyper[,k] <- apply(pHyperChain[1:noAdapt,k,],2, mean)
        covParHyperList[[k]] <- scalPropHyper * (cov(pHyperChain[1:noAdapt,k,]) + epsDiagHyper)
      }
    }
    if(j > (nBI + noAdapt)){
      accept.prob.re <- accept.prob.re + accept.re
      accept.prob.theta <- accept.prob.theta + accept.theta
      accept.prob.hyper <- accept.prob.hyper + accept.hyper
      
      t <- j - nBI
      
      for (m in 1:noIndividuals) {
        aveParRe_new[,m] <- as.vector(((t-1) * aveParRe[,m] + pReMat[,m]) / t)
        covParReList_new[[m]] <- ((t-2) * covParReList[[m]] + scalPropRe * ((t-1) * (aveParRe[,m] %o% aveParRe[,m]) - t * (aveParRe_new[,m] %o% aveParRe_new[,m]) + (pReMat[,m] %o% pReMat[,m])) + epsDiagRe) / (t-1)
        aveParRe[,m] <- aveParRe_new[,m]
        covParReList[[m]] <- covParReList_new[[m]]
      }
      
      aveParTheta_new <- as.vector(((t-1) * aveParTheta + pTheta) / t)
      covParTheta_new <- ((t-2) * covParTheta + scalPropTheta * ((t-1) * (aveParTheta %o% aveParTheta) - t * (aveParTheta_new %o% aveParTheta_new) + (c(pTheta) %o% c(pTheta))) + epsDiagTheta) / (t-1)
      aveParTheta <- aveParTheta_new
      covParTheta <- covParTheta_new
      
      for (k in 1:nb) {
        aveParHyper_new[,k] <- as.vector(((t-1) * aveParHyper[,k] + pHyperMat[k,]) / t)
        covParHyperList_new[[k]] <- ((t-2) * covParHyperList[[k]] + scalPropHyper * ((t-1) * (aveParHyper[,k] %o% aveParHyper[,k]) - t * (aveParHyper_new[,k] %o% aveParHyper_new[,k]) + (pHyperMat[k,] %o% pHyperMat[k,])) + epsDiagHyper) / (t-1)
        aveParHyper[,k] <- aveParHyper_new[,k]
        covParHyperList[[k]] <- covParHyperList_new[[k]]
      }
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
f.hyper <- 1
HyperPrior <- rbind(c(1,0.05),
                    c(1,0.01),
                    c(1,0.01),
                    c(1,0.02))
res <- AM.mixed.effects(y,
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
                                f.hype, 
                                HyperPrior)
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

# HMC ---------------------------------------------------------------------


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
                              L) {
  
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
  
  # First call to the model. Calculate likelihood
  logLike <- logLikeCan <- rep(NA_real_, noIndividuals)
  for (m in 1:noIndividuals) {
    pRe <- pReMat[,m]
    logLike[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pRe)
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
        # Call log likelihood
        logLikeCan[m] <- particleFilter_b(y[,,m],pTheta,noParticles,startState,noStep,pReCan)
        
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
        logLikeCan[m] <- particleFilter_b(y[,,m],pThetaCan,noParticles,startState,noStep, pRe)
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
