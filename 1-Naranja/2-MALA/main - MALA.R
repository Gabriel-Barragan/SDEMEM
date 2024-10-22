# Algoritmo Metropolis - Hastings de particulas
# Metrpolis - adjusted Langevin algorithm (MALA)

library(ggplot2)
library(dplyr)
library(deSolve)
library(coda)
library(matrixStats)

ggplot2::theme_set(theme_bw())

# Ubicar la direccion
setwd("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja")

# Obtener los datos
source("Datos_naranja.R")
naranja <- Datos_naranja(FALSE)
M <- naranja$M
datos <- naranja$datos

# Funciones
source("log_prior.R")
source("filtro_particulas.R")
source("simulacion_trayectorias.R")
source("2-MALA/Algoritmo - MHP_MALA.R")
source("funciones.R")
source("log_eta.R")
source("log_gradiente.R")

# x <- c(x = 50)
# parms <- c(phi_1= 200, phi_2= 350)
# times <- seq(200, 260, 0.1)
# out <- ode(x, times, model, parms)
# plot(out)

# Parametros del algortimo
iteraciones <- 2000
periodo_quemado <- 0.25*iteraciones

FP <- 2
D <- 50
N <- 20

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
Sigma <- diag(1, npar)
Sigma[6,6] <-Sigma[5,5] <- 0.1
# EM
# load("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja/2-MALA/res_MHP_MALA_EM_1.RData")

# PDM
# load("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja/2-MALA/res_MHP_MALA_PDM_1.RData")

# PDR
# load("/cloud/project/5506643/Tesis Version final 2/Version final/1-Naranja/2-MALA/res_MHP_MALA_PDR_2.RData")

# Sigma <- 2.4^2/npar*var(res$cadena[burn_in,])

# Correlacion
rho <- 0.999

# Efectos aleatorios
phi1_m <- rep(190, M)
sig_phi1 <- rep(25, M)

phi2_m <- rep(350, M)
sig_phi2 <- rep(52.5, M)

apply(res$phi1[burn_in,], 2, summary)
(phi1_m <- apply(res$phi1[burn_in,], 2, mean))
#(phi1_m <- apply(res$phi1[burn_in,], 2, median))
(sig_phi1 <- apply(res$phi1[burn_in,], 2, sd))

apply(res$phi2[burn_in,], 2, summary)
(phi2_m <- apply(res$phi2[burn_in,], 2, mean))
#(phi2_m <- apply(res$phi2[burn_in,], 2, median))
(sig_phi2 <- apply(res$phi2[burn_in,], 2, sd))

eps <- 0.91

# Algortimo
res <- MHP_MALA(datos,
                M,
                theta_0,
                x0,
                Sigma,
                iteraciones,
                periodo_quemado,
                FP,
                D,
                N,
                eps,
                phi1_m, phi2_m,
                sig_phi1, sig_phi2,
                rho)

# save(res, file = "2-MALA/res_MHP_MALA_EM_3.RData")
# save(res, file = "2-MALA/res_MHP_MALA_PDM_4.RData")
# save(res, file = "2-MALA/res_MHP_MALA_PDR_3.RData")

# Resultados
res$tasa_theta*100
res$tasa_ea*100

# Tiempo de ejecucion
res$tiempo/60

# Efectos aleatorios
# Resumen estadistico
phi_1_res <- res$phi1[burn_in,]
(resumen <- apply(phi_1_res, 2, summary) %>% round(.,2))
apply(phi_1_res, 2, mean)
apply(phi_1_res, 2, sd)
# write.csv(resumen, "2-MALA/naranja_phi_1_EM_1.csv", row.names = T)
# write.csv(resumen, "2-MALA/naranja_phi_1_PDM_1.csv", row.names = T)
# write.csv(resumen, "2-MALA/naranja_phi_1_PDR_1.csv", row.names = T)

phi_2_res <- res$phi2[burn_in,]
(resumen <- apply(phi_2_res, 2, summary) %>% round(.,2))
apply(phi_2_res, 2, mean)
apply(phi_2_res, 2, sd)
# write.csv(resumen, "2-MALA/naranja_phi_2_EM_1.csv", row.names = T)
# write.csv(resumen, "2-MALA/naranja_phi_2_PDM_1.csv", row.names = T)
# write.csv(resumen, "2-MALA/naranja_phi_2_PDR_1.csv", row.names = T)

# Graficos efectos aleatorios

par(mfrow=c(2,5))
for (k in 1:5) {
  plot(phi_1_res[,k], type="l", main = paste("m = ", k), xlab = "", ylab = "")
  abline(h = mean(phi_1_res[,k]), col='red', lwd=3)
}

for (k in 1:5) {
  hist(phi_1_res[,k], main="", col = "lightblue", xlab = "", ylab = "")
  abline(v=mean(phi_1_res[,k]), col='red', lwd=3)
}

par(mfrow=c(2,5))
for (k in 1:5) {
  plot(phi_2_res[,k], type="l", main = paste("m = ", k), xlab = "", ylab = "")
  abline(h = mean(phi_2_res[,k]), col='red', lwd=3)
}

for (k in 1:5) {
  hist(phi_2_res[,k], main="", col = "pink", xlab = "", ylab = "")
  abline(v=mean(phi_2_res[,k]), col='red', lwd=3)
}

par(mfrow=c(1,2))
phi_1_res %>% boxplot(col = "blue", main=expression(phi[1]))
phi_2_res %>% boxplot(col = "red", main=expression(phi[2]))

par(mfrow=c(2,5))
for (k in 1:5) {
  acf(phi_1_res[,k], main = expression(phi[1])) 
}

for (k in 1:5) {
  acf(phi_2_res[,k], main = expression(phi[2])) 
}


# Parametros

# Graficos 
par(mfrow=c(2,3))

cadena <- theta_ut(res$cadena[burn_in, ])
parameterNames <- c(expression(mu[phi[1]]),
                    expression(sigma[phi[1]]),
                    expression(mu[phi[2]]),
                    expression(sigma[phi[2]]),
                    expression(sigma),
                    expression(xi))
parameterColors <- c("blue", "lightblue", "green", "lightgreen",
                     "purple", "darkblue")
for (k in 1:6) {
  # Plot 
  plot(cadena[,k], type="l", main = parameterNames[k], xlab = "", ylab = "")
  abline(h = mean(cadena[,k]), col='red', lwd=3)
  abline(h=theta_ut(theta_0)[k], col="black", lwd=3)
  
  hist(cadena[,k], main=parameterNames[k], col = parameterColors[k], xlab = "", ylab = "")
  abline(v=mean(cadena[,k]), col='red', lwd=3)
  abline(v=theta_ut(theta_0)[k], col="black", lwd=3)
  
  boxplot(cadena[,k], main=parameterNames[k], col = parameterColors[k])
}

# Traza
par(mfrow=c(2,3))
for (k in c(1,3,5,
            2,4,6)) {
  # Plot 
  plot(cadena[,k], type="l", main = parameterNames[k], xlab = "", ylab = "")
  abline(h = theta_ut(theta_0)[k], col='black', lwd=3)
  abline(h = mean(cadena[,k]), col='red', lwd=3)
}

# Histograma
for (k in c(1,3,5,2,4,6)) {
  hist(cadena[,k], main=parameterNames[k], col = parameterColors[k], xlab = "", ylab = "")
  abline(v = theta_ut(theta_0)[k], col='black', lwd=3)
  abline(v = mean(cadena[,k]), col='red', lwd=3)
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
# write.csv(resumen, "2-MALA/naranja_EM_1.csv", row.names = T)
# write.csv(resumen, "2-MALA/naranja_PDM_1.csv", row.names = T)
# write.csv(resumen, "2-MALA/naranja_PDR_1.csv", row.names = T)

apply(cadena, 2, mean) %>% round(.,2)
apply(cadena, 2, median) %>% round(.,2)
apply(cadena, 2, sd) %>% round(.,2)

# Actualizar la cadena
Sigma <- 2.4^2/npar*var(res$cadena[burn_in,])

#####################
# Diagnostico

# Analisis con coda
res_coda <- as.mcmc(cadena)
plot(res_coda)
plot(res_coda[,1], main=expression(mu[phi[1]]))
plot(res_coda[,2], main=expression(sigma[phi[1]]))
plot(res_coda[,3], main=expression(mu[phi[2]]))
plot(res_coda[,4], main=expression(sigma[phi[2]]) )
plot(res_coda[,5], main=expression(sigma ))
plot(res_coda[,6], main=expression(xi))

#Autocorrelation
autocorr.plot(res_coda)
autocorr.plot(res_coda[,1], main=expression(mu[phi[1]]))
autocorr.plot(res_coda[,2], main=expression(sigma[phi[1]]) )
autocorr.plot(res_coda[,3], main=expression(mu[phi[2]]))
autocorr.plot(res_coda[,4], main=expression(sigma[phi[2]]) )
autocorr.plot(res_coda[,5], main=expression(sigma) )
autocorr.plot(res_coda[,6], main=expression(xi))
### IACT
parameterNames <- c(expression(mu[phi[1]]),
                    expression(sigma[phi[1]]),
                    expression(mu[phi[2]]),
                    expression(sigma[phi[2]]),
                    expression(sigma),
                    expression(xi))
parameterACFnames <- c(expression("ACF of " * mu[phi[1]]), 
                       expression("ACF of " * sigma[phi[1]]),
                       expression("ACF of " * mu[phi[2]]),
                       expression("ACF of " * sigma[phi[2]]),
                       expression("ACF of " * sigma),
                       expression("ACF of " * xi)
)

parameterColors <- c("blue", "lightblue", "green", "lightgreen",
                     "purple", "darkblue")
iact <- c()
for (k in 1:6) {
  # Plot the autocorrelation function
  acf_res <- acf(cadena[, k], plot = FALSE, lag.max = NULL)
  plot(
    acf_res$lag,
    acf_res$acf,
    col = parameterColors[k],
    type = "l",
    xlab = "iteration",
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
# Sigma <- 2.38**2/dim(res$cadena)[2]*var(res$cadena)
# save(Sigma, file = "2-MALA/Sigma_MHI_EM.RData")
# save(res, file = "~/Maestria/Titulación/Tesis/Gabriel Barragan/Proyecto/Version final/1-Naranja/1-MHP/res_MHI_EM.RData")

(epsilon_grid <- seq(0.1, 1, by=0.1))
n <- length(epsilon_grid)
theta_grid <- array(NA, dim=c(n, 2))
ea_grid <- array(NA, dim=c(n, 5))
for (k in 1:n) {
  res <- MHP_MALA(datos, M, theta_0, x0, Sigma, iteraciones, periodo_quemado,
                  FP, D, N, epsilon_grid[k], phi1_m, phi2_m, sig_phi1, sig_phi2, rho)
  theta_grid[k,] <- res$tasa_theta
  ea_grid[k,] <- res$tasa_ea
}

table_theta <- cbind(epsilon_grid, theta_grid)
colnames(table_theta) <- c("eps", "Pop. effects", "Fixed effects")
print(table_theta)

table_ea <- cbind(epsilon_grid, ea_grid)
colnames(table_ea) <- c("eps", "M1", "M2", "M3", "M4", "M5")
print(table_ea)

#EM: eps = [0.2, 0.4]
# write.csv(table_theta, "2-MALA/naranja_EM_tasa_theta_1.csv", row.names = T)
# write.csv(table_ea, "2-MALA/naranja_EM_tasa_ea_1.csv", row.names = T)

#PDM: eps = [0.75, 1]
# write.csv(table_theta, "2-MALA/naranja_PDM_tasa_theta_1.csv", row.names = T)
# write.csv(table_ea, "2-MALA/naranja_PDM_tasa_ea_1.csv", row.names = T)

# PDR: eps=0.5
# write.csv(table_theta, "2-MALA/naranja_PDR_tasa_theta_1.csv", row.names = T)
# write.csv(table_ea, "2-MALA/naranja_PDR_tasa_ea_1.csv", row.names = T)
