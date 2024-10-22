logGamma <- function(logY,alpha,beta){
  return(alpha*log(beta)-log(gamma(alpha)) + (alpha-1)*exp(logY) - beta*exp(logY) + logY)
}

logNorm <- function(x,m,s) {
  return(-log(s) - 0.5*log(2*pi) - 0.5*(1/s**2)*(x-m)**2)
}

logMultNorm <- function(X,M,S) {
  k <- length(X)
  return(-0.5*(k*log(2*pi) + log(det(S)) + t(X - M)%*%solve(S)%*%(X - M)))
}

model <- function (time, x, parms) {
  with(as.list(c(x, parms)), {
    dx <- 1/(phi_1*phi_2)*x*(phi_1 - x)
    return(list(dx))
  })
}

theta_ut <- function(theta) {
  if (is.null(dim(theta))) {
    theta_tmp <- theta[c(2,4,5,6)] 
    theta[c(2,4,5,6)] <- sqrt(1/exp(theta_tmp))  
  }
  else {
    theta_tmp <- theta[,c(2,4,5,6)] 
    theta[,c(2,4,5,6)] <- sqrt(1/exp(theta_tmp))  
  }
  return(theta)
}

MALA <- function(theta, precondition_Sigma, eps, phi_1, phi_2) {
  p <- length(theta)
  eps_pro <- eps**2/p**(1/3)
  log_grad <- log_gradiente(phi_1, phi_2, theta)
  mu_star <- theta +  0.5 * eps_pro * precondition_Sigma %*% log_grad
  cov <- eps_pro*precondition_Sigma
  prop <- mu_star + c(tcrossprod(rnorm(p),chol(cov)))
  # prop <- mu_star + c(crossprod(chol(cov), rnorm(p)))
  # q(theta*|theta)
  log_q <- logMultNorm(prop, mu_star, cov)
  
  # calculate q(theta|theta*)
  log_grad_nuevo <- log_gradiente(phi_1, phi_2, prop)
  mu_theta <- prop +  0.5 * eps_pro * precondition_Sigma %*% log_grad
  log_q_star <- logMultNorm(theta, mu_theta, cov)
  return(list(prop=prop, log_q=log_q, log_q_star=log_q_star))
}

Hamiltoniano <- function(theta, Sigma, S_inv, S_inv_sqrt, eps, L, phi_1, phi_2){
  d <- length(theta)
  P <- c(crossprod(S_inv_sqrt, rnorm(d))) # Variables auxiliares
  log_K <- 0.5*P %*% Sigma %*% P # Energia cinetica al inicio de la trayectoria
  
  # Al inicio, hacer un medio paso para el momento
  theta_nuevo <- theta
  P_nuevo <- P + 0.5*eps*log_gradiente(phi_1, phi_2, theta_nuevo)
  
  # Alternar los pasos completos para posicion y momento
  for (l in 1:L) {
    # Hacer un paso completo para la posicion
    theta_nuevo <- theta_nuevo + eps * c(Sigma %*% P_nuevo)
    
    # Hacer un paso completo para el momento, excepto al final de la trayectoria
    if ( l != L) P_nuevo <- P_nuevo + eps*log_gradiente(phi_1, phi_2, theta_nuevo)
  }
  # Al final, hacer un medio paso para el momento
  P_nuevo <- P_nuevo + 0.5*eps*log_gradiente(phi_1, phi_2, theta_nuevo)
  
  # Negar el momento al final de la trayectoria para hacer la propuesta simetrica
  P_nuevo <- -P_nuevo
  
  # Energia cinetica al final de la trayectoria
  log_K_nuevo <- 0.5 * P_nuevo %*% Sigma %*% P_nuevo
  
  return(list(propueta=theta_nuevo, log_K=log_K, log_K_nuevo=log_K_nuevo))
}

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


