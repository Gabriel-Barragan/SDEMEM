# SDE

n.step <- 10000
start.time <- 0
end.time <- 10
delta.t <- (end.time - start.time)/n.step 
x0 <- 0.001

drift.fun <- function(b,x){
  -b*x*log(x)
}

diffusion.fun <- function(sigma,x){
  sigma^2*x^2
}
x <- rep(0,n.step)
y <- rep(0,n.step/100)
x[1] <- x0

b <- 0.6
sigma <- 0.1
sigma_obs <- 0.02
k <- 1
for (i in 1:(n.step-1)) {
  
  if (x[i]<=0){
    x[i+1] <- x0
  } else {
    x[i+1] <- x[i] + drift.fun(b,x[i])*delta.t + sqrt(delta.t*diffusion.fun(sigma,x[i]))*rnorm(1)  
  }
  if(i%%100==0 || i==n.step-1) {
    y[k] <- x[i+1]+sqrt(sigma_obs)*rnorm(1)
    k <- k+1
  }
}
time.1 <- seq(start.time, end.time,length.out=n.step)
time.2 <- seq(100,n.step, by=100)
plot(time.1,x,type="l", ylim = c(min(x,y), max(x,y)))
points(time.2*delta.t,y,col="salmon")

theta <- c(0.45, 0.2, 0.05)

SDE <- function(b,s,x.init,no.Discretization){
  delta.t.SDE <- 1/no.Discretization
  x <- x.init
  for (i in 1:no.Discretization) {
    if(x<=0){
      x <- x.init
    } else {
      x <- x + drift.fun(b,x)*delta.t.SDE + sqrt(diffusion.fun(s, x)*delta.t.SDE)*rnorm(1)   
    }
  }
  x
}

particle.filter <- function(y, theta, no.Particles, initial.State, no.Discretization){
  obs.time <- length(y)  
  b.par <- theta[1]
  sigma.par <- theta[2]
  sigma.obs.par <- theta[3]
  
  for (t in 2:obs.time) {
    
  }
}