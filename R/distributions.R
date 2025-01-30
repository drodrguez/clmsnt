# scripts for skew normal t distribution
# in this file:
# - skew normal t probability density function as dskewt
# - skew notmal t cumulative distribution function as pskewt
library("cubature")

dskewt <- function(x, lambda=0, v=5){
  # Calcula la función de densidad de probabilidad de la distribución skew-t
  # usando la representación: p = 2 * dt(x, v) * pnorm(lambda * x)
  # donde dt es la densidad t-Student y pnorm la CDF normal estándar
  #
  # Args:
  #   x: Vector numérico donde evaluar la densidad
  #   lambda: Parámetro de asimetría (default = 0)
  #   v: Grados de libertad (default = 5)
  #
  # Returns:
  #   Vector numérico con las densidades evaluadas

  # Inicializamos el vector de resultados
  p <- numeric(length(x))

  # Control de estabilidad numérica para valores extremos de x
  # Como dt usa x^2, prevenimos overflow chequeando contra sqrt(max_double)
  idx_valid <- abs(x) <= .Machine$double.xmax^0.5

  if(any(idx_valid)) {
    # Calculamos solo para valores válidos de x
    x_valid <- x[idx_valid]

    # Calculamos la densidad t-Student
    d_t <- dt(x_valid, df=v)

    # Para valores extremos de lambda*x, pnorm da efectivamente 0 o 1
    # 37 es el umbral donde pnorm alcanza la precisión de máquina
    lambda_x <- lambda * x_valid
    p_norm <- ifelse(abs(lambda_x) > 37,
                     ifelse(lambda_x < 0, 0, 1),
                     pnorm(lambda_x))

    # Calculamos la densidad final
    p[idx_valid] <- 2 * d_t * p_norm

    # Ponemos 0 donde d_t o p_norm son 0
    p[idx_valid][d_t == 0 | p_norm == 0] <- 0
  }

  return(p)
}


pskewt = function(z, lambda=0, v=5) {
  # Función que calcula la CDF de la distribución skew-t
  # Args:
  #   z: punto(s) donde evaluar la CDF
  #   lambda: parámetro de asimetría (default = 0)
  #   v: grados de libertad (default = 5)
  # Returns:
  #   Probabilidad(es) acumulada(s)

  # Inicializar vector de resultados
  p <- numeric(length(z))

  # Manejar casos infinitos y extremos de una vez
  extreme_idx <- is.infinite(z) | abs(z) > 1e10
  if(any(extreme_idx)) {
    p[extreme_idx] <- ifelse(z[extreme_idx] > 0, 1, 0)
  }

  # Procesar valores regulares de una vez
  regular_idx <- !extreme_idx
  if(any(regular_idx)) {
    z_reg <- z[regular_idx]

    # Integración numérica vectorizada
    lower_limit <- pmax(-1e10, pmin(-50, 5*z_reg))

    F_values <- mapply(
      function(zi, li) {
        F_val <- adaptIntegrate(
          f = function(x) dskewt(x, lambda=lambda, v=v),
          lowerLimit = li,
          upperLimit = zi,
          maxEval = 1000,
          absError = 1e-10,
          tol = 1e-8
        )
        return(F_val$integral)
      },
      z_reg, lower_limit
    )

    # Asignar resultados
    p[regular_idx] <- pmax(0, pmin(1, F_values))
  }

  return(p)
}
