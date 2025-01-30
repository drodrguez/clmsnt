# main.R

# Cargar librerías necesarias
library(ordinal)   # Para comparar con el modelo ordinal estándar
data(wine)

# Cargar nuestras funciones implementadas
source("R/distributions.R")
source("R/model.R")
source("R/estimation.R")
source("R/methods.R")

model <- clm_stn(rating ~ temp + contact, data = wine)
fit <- fit_clm_stn(model)
print(fit)
fit$gradient

