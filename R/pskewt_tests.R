setwd("/home/drodrguez/Documentos/UC/skewCLM/skewtCLM/R")
source("skewt.R")
# 1. Verificar propiedades básicas de una CDF
test_basic_properties <- function(lambda=0, v=5) {
  # Puntos de prueba
  z <- seq(-5, 5, length.out=100)
  p <- pskewt(z, lambda=lambda, v=v)

  # Lista de pruebas
  tests <- list(
    "Monotonía" = all(diff(p) >= -1e-10),
    "Rango [0,1]" = all(p >= 0 & p <= 1),
    "Límite -Inf" = all.equal(pskewt(-Inf, lambda, v), 0),
    "Límite Inf" = all.equal(pskewt(Inf, lambda, v), 1)
  )

  return(tests)
}

# 2. Verificar que la derivada numérica coincide con la PDF
test_derivative <- function(lambda=0, v=5) {
  z <- seq(-5, 5, length.out=100)
  h <- 0.001  # paso para diferencia finita

  # Calcular derivada numérica
  p1 <- pskewt(z + h, lambda, v)
  p2 <- pskewt(z - h, lambda, v)
  numerical_pdf <- (p1 - p2)/(2*h)

  # Calcular PDF directamente
  theoretical_pdf <- dskewt(z, lambda, v)

  # Comparar
  max_diff <- max(abs(numerical_pdf - theoretical_pdf))
  return(max_diff < 1e-5)
}

# 3. Verificar caso simétrico (lambda=0)
test_symmetric_case <- function(v=5) {
  z <- seq(-5, 5, length.out=100)
  p_skewt <- pskewt(z, lambda=0, v=v)
  p_t <- pt(z, df=v)

  max_diff <- max(abs(p_skewt - p_t))
  return(max_diff < 1e-5)
}

# 4. Verificar simetría cuando lambda=0
test_symmetry <- function(v=5) {
  z <- seq(0, 5, length.out=50)
  p_pos <- pskewt(z, lambda=0, v=v)
  p_neg <- pskewt(-z, lambda=0, v=v)

  max_diff <- max(abs(p_pos + p_neg - 1))
  return(max_diff < 1e-5)
}

# Ejecutar todas las pruebas
run_all_tests <- function() {
  cat("1. Probando propiedades básicas...\n")
  print(test_basic_properties(lambda=2, v=5))

  cat("\n2. Probando relación con PDF...\n")
  print(test_derivative(lambda=2, v=5))

  cat("\n3. Probando caso simétrico...\n")
  print(test_symmetric_case(v=5))

  cat("\n4. Probando simetría cuando lambda=0...\n")
  print(test_symmetry(v=5))

  # Visualización
  z <- seq(-5, 5, length.out=100)
  df <- data.frame(
    z = rep(z, 3),
    p = c(
      pskewt(z, lambda=0, v=5),
      pskewt(z, lambda=2, v=5),
      pskewt(z, lambda=-2, v=5)
    ),
    type = rep(c("Simétrica", "Asimetría positiva", "Asimetría negativa"),
               each=length(z))
  )

  library(ggplot2)
  p <- ggplot(df, aes(x=z, y=p, color=type)) +
    geom_line() +
    theme_minimal() +
    labs(title="CDF de la distribución skew-t",
         x="z", y="F(z)",
         color="Tipo") +
    ylim(0,1)

  print(p)
}

# Ejecutar todas las pruebas
run_all_tests()
