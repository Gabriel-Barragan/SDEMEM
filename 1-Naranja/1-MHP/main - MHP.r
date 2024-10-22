# Algoritmo Metropolis - Hastings de particulas
# Clasico

library(ggplot2)
library(dplyr)
library(deSolve)
library(coda)
library(matrixStats)

ggplot2::theme_set(theme_bw())

# Direccion de trabajo
setwd("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja")

# Obtener los datos
source("Datos_naranja.R")
naranja <- Datos_naranja(T)
M <- naranja$M
datos <- naranja$datos

# Funciones
source("log_prior.R")
source("filtro_particulas.R")
source("simulacion_trayectorias.R")
source("1-MHP/Algoritmo - MHP.R")
source("funciones.R")

# x <- c(x = 50)
# parms <- c(phi_1 = 200, phi_2 = 350)
# times <- seq(200, 260, 0.1)
# out <- ode(x, times, model, parms)
# plot(out)

# Parametros del algortimo

iteraciones <- 2000
periodo_quemado <- 0.1*iteraciones

FP <- 1
D <- 50
N <- 20
L <- 20 # Efectos aleatorios

# Parametros iniciales
# Theta_0
mu_phi_1 <- 190 
sigma_phi_1 <- 25
tau_1 <- 1/sigma_phi_1**2

mu_phi_2 <- 350 
sigma_phi_2 <- 52.5
tau_2 <- 1/sigma_phi_2**2

sig <- 0.08
tau_sigma <- 1/sig**2 

xi <- 0.5
tau_xi <- 1/xi**2

theta_0 <- c(mu_phi_1, log(tau_1),
             mu_phi_2, log(tau_2),
             log(tau_sigma), log(tau_xi))

# Estado inicial
x0 <- 30

# Periodo de quemado
burn_in <- (periodo_quemado+1):(iteraciones+periodo_quemado)

# Matriz de covarianza inicial
npar <- length(theta_0)
# Sigma <- diag(1, npar)

# EM
# load("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja/1-MHP/res_MHI_EM_1.RData")

# PDM
# load("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja/1-MHP/res_MHI_PDM_1.RData")

# PDR
# load("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja/1-MHP/res_MHI_PDR_1.RData")

Sigma <- 2.4^2/npar*var(res$cadena[burn_in,])

# Correlacion
rho <- 0.999

# Algortimo
res <- MHP(datos,
           M,
           theta_0,
           x0,
           Sigma,
           iteraciones,
           periodo_quemado,
           FP, D, N, L, rho)

# save(res, file = "1-MHP/res_MHI_EM_4.RData")

# save(res, file = "1-MHP/res_MHI_PDM_2.RData")

# save(res, file = "1-MHP/res_MHI_PDR_4.RData")

# Resultados
# Tasa de aceptacion
res$tasa_theta*100

# Tiempo de ejecucion
res$tiempo/60

# Logaritmo de verosimilitud
res$lv %>% summary(.)
res$lv %>% sd(.)
par(mfrow=c(1,1))
res$lv %>% plot(., type="l")

# Graficos
par(mfrow=c(2,3))

# Transformacion inversa de la cadena
cadena <- theta_ut(res$cadena[burn_in, ])

# Nombres de los parametros
parameterNames <- c(expression(mu[phi[1]]),
                    expression(sigma[phi[1]]),
                    expression(mu[phi[2]]),
                    expression(sigma[phi[2]]),
                    expression(sigma),
                    expression(xi))

# Colores de los parametros
parameterColors <- c("blue", "lightblue", "green", "lightgreen",
                     "purple", "darkblue")


# for (k in 1:6) {
#   # Plot 
#   plot(cadena[,k], type="l", main = parameterNames[k], xlab = "", ylab = "")
#   abline(h = mean(cadena[,k]), col='red', lwd=3)
#   abline(h = theta_ut(theta_0)[k], col='black', lwd=3)
#   
#   hist(cadena[,k], main=parameterNames[k], col = parameterColors[k], xlab = "", ylab = "")
#   abline(v=mean(cadena[,k]), col='red', lwd=3)
#   abline(v = theta_ut(theta_0)[k], col='black', lwd=3)
#   
#   boxplot(cadena[,k], main=parameterNames[k], col = parameterColors[k])
# }

# Traza
par(mfrow=c(2,3))
for (k in c(1,3,5,
            2,4,6)) {
  plot(cadena[,k], type="l", main = parameterNames[k], xlab = "", ylab = "")
  abline(h = theta_ut(theta_0)[k], col='black', lwd=3)
  abline(h = mean(cadena[,k]), col='red', lwd=3)
}

# Histograma
for (k in c(1,3,5,2,4,6)) {
  hist(cadena[,k], main=parameterNames[k], col = parameterColors[k], xlab = "", ylab = "")
  abline(v = theta_ut(theta_0)[k], col='black', lwd=3)
  abline(v=mean(cadena[,k]), col='red', lwd=3)
}

# Boxplot
for (k in c(1,3,5,2,4,6)) {
  boxplot(cadena[,k], main=parameterNames[k], col = parameterColors[k])
}

# Autocorrelacion
for (k in c(1,3,5,2,4,6)) {
  acf(cadena[,k], main = parameterNames[k])  
}

# Estadistica descriptiva
(resumen <- apply(cadena, 2, summary) %>% round(.,2))
# write.csv(resumen, "1-MHP/naranja_EM_2.csv", row.names = T)
# write.csv(resumen, "1-MHP/naranja_PDM_3.csv", row.names = T)
# write.csv(resumen, "1-MHP/naranja_PDR_2.csv", row.names = T)

apply(cadena, 2, sd) %>% round(.,2)

# Diagnostico
# https://warwick.ac.uk/fac/sci/wdsi/events/wrug/resources/mcmcse.pdf
cadena <- as.mcmc(cadena)
effectiveSize(cadena)
summary(effectiveSize(cadena)) # Effective sample size
summary(iteraciones / effectiveSize(cadena)) # Integrated autocorrelation time
summary(1 - rejectionRate(cadena)) # Acceptance rate

c(
  res$tiempo, mean(effectiveSize(cadena)),
  mean(effectiveSize(cadena)) / res$tiempo,
  1 - mean(rejectionRate(cadena))
)

### IACT
parameterNames <- c(expression(mu[phi[1]]),
                    expression(sigma[phi[1]]),
                    expression(mu[phi[2]]),
                    expression(sigma[phi[2]]),
                    expression(sigma),
                    expression(xi))
parameterACFnames <- c(expression("ACF de " * mu[phi[1]]), 
                       expression("ACF de " * sigma[phi[1]]),
                       expression("ACF de " * mu[phi[2]]),
                       expression("ACF de " * sigma[phi[2]]),
                       expression("ACF de " * sigma),
                       expression("ACF de " * xi)
                       )

parameterColors <- c("blue", "lightblue", "green", "lightgreen",
                           "purple", "darkblue")
iact <- c()
par(mfrow=c(2,3))
for (k in 1:6) {
  # Plot the autocorrelation function
  acf_res <- acf(res$cadena[, k], plot = FALSE, lag.max = 100)
  plot(
    acf_res$lag,
    acf_res$acf,
    col = parameterColors[k],
    type = "l",
    xlab = "iteracion",
    ylab = parameterACFnames[k],
    lwd = 2,
    #ylim = c(-0.2, 1),
    bty = "n"
  )
  polygon(
    c(acf_res$lag, rev(acf_res$lag)),
    c(acf_res$acf, rep(0, length(acf_res$lag))),
    border = NA,
    col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
  )
  abline(h = 1.96 / sqrt(iteraciones - periodo_quemado), lty = "dotted")
  abline(h = -1.96 / sqrt(iteraciones - periodo_quemado), lty = "dotted")
  
  iact <- c(iact, 1 + 2 * sum(acf_res$acf))
}

iact

# Actualizar matriz de covarianza
Sigma <- res$Sigma
res <- MHP(datos, M, theta_0, x0, Sigma, iteraciones, periodo_quemado,
           FP, D, N, L, rho)

# (theta_0 <- apply(res$cadena, 2, median))
# save(Sigma, file = "1-MHP/Sigma_MHI_EM.RData")
# save(res, file = "1-MHP/res_MHI_EM.RData")
