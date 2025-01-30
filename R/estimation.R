# estimation.R

score_clm_stn <- function(theta, model) {
  # Extraer parámetros del vector theta
  n_thresholds <- model$n_levels - 1
  alpha <- theta[1:n_thresholds]              # Parámetros umbral
  beta <- theta[(n_thresholds + 1):(n_thresholds + model$n_vars)]  # Coeficientes de regresión
  lambda <- theta[length(theta)]              # Parámetro de asimetría

  # Obtener datos del modelo
  X <- model$X                                # Matriz de diseño
  y <- as.numeric(model$y)                    # Categorías como números
  n <- nrow(X)                                # Número de observaciones

  # Inicializar vectores de score
  score_alpha <- numeric(n_thresholds)
  score_beta <- numeric(model$n_vars)
  score_lambda <- 0

  # Calcular predictor lineal para cada observación
  eta <- as.vector(X %*% beta)

  # Para cada observación
  for(i in 1:n) {
    # Identificar la categoría observada
    j <- y[i]

    # Calcular las contribuciones según la categoría observada
    if(j == 1) {
      # Primera categoría (j=1)
      pi_ij <- pskewt(alpha[j] - eta[i], lambda, model$v)
      pi_ij <- max(pi_ij, .Machine$double.eps) # Evitar división por cero

      # Score para alpha (solo contribuye α₁)
      score_alpha[j] <- score_alpha[j] + dskewt(alpha[j] - eta[i], lambda, model$v)/pi_ij

      # Score para beta
      score_beta <- score_beta - X[i,] * dskewt(alpha[j] - eta[i], lambda, model$v)/pi_ij

      # Score para lambda (usando integración numérica)
      int_j <- integrate(
        function(t) 2*dt(t, df=model$v)*t*dnorm(lambda*t),
        lower=-50,
        upper=alpha[j] - eta[i],
        subdivisions=1000,
        rel.tol=1e-6,
        abs.tol=1e-6
      )$value

      score_lambda <- score_lambda + int_j/pi_ij

    } else if(j == model$n_levels) {
      # Última categoría
      pi_ij <- 1 - pskewt(alpha[j-1] - eta[i], lambda, model$v)
      pi_ij <- max(pi_ij, .Machine$double.eps)

      # Score para alpha (solo contribuye α_{j-1})
      score_alpha[j-1] <- score_alpha[j-1] -
        dskewt(alpha[j-1] - eta[i], lambda, model$v)/pi_ij

      # Score para beta
      score_beta <- score_beta + X[i,] * dskewt(alpha[j-1] - eta[i], lambda, model$v)/pi_ij

      # Score para lambda
      int_j1 <- integrate(
        function(t) 2*dt(t, df=model$v)*t*dnorm(lambda*t),
        lower=-50,
        upper=alpha[j-1] - eta[i],
        subdivisions=1000,
        rel.tol=1e-6,
        abs.tol=1e-6
      )$value

      score_lambda <- score_lambda - int_j1/pi_ij

    } else {
      # Categorías intermedias
      pi_ij <- pskewt(alpha[j] - eta[i], lambda, model$v) -
        pskewt(alpha[j-1] - eta[i], lambda, model$v)
      pi_ij <- max(pi_ij, .Machine$double.eps)

      # Score para alpha (contribuyen αj y α_{j-1})
      score_alpha[j] <- score_alpha[j] +
        dskewt(alpha[j] - eta[i], lambda, model$v)/pi_ij
      score_alpha[j-1] <- score_alpha[j-1] -
        dskewt(alpha[j-1] - eta[i], lambda, model$v)/pi_ij

      # Score para beta
      score_beta <- score_beta + X[i,] *
        (dskewt(alpha[j-1] - eta[i], lambda, model$v) -
           dskewt(alpha[j] - eta[i], lambda, model$v))/pi_ij

      # Score para lambda
      int_j <- integrate(
        function(t) 2*dt(t, df=model$v)*t*dnorm(lambda*t),
        lower=-50,
        upper=alpha[j] - eta[i],
        subdivisions=1000,
        rel.tol=1e-6,
        abs.tol=1e-6
      )$value

      int_j1 <- integrate(
        function(t) 2*dt(t, df=model$v)*t*dnorm(lambda*t),
        lower=-50,
        upper=alpha[j-1] - eta[i],
        subdivisions=1000,
        rel.tol=1e-6,
        abs.tol=1e-6
      )$value

      score_lambda <- score_lambda + (int_j - int_j1)/pi_ij
    }
  }

  # Devolver el vector de scores completo
  return(c(score_alpha, score_beta, score_lambda))
}
#' Fit CLM with skew-t-normal link using maximum likelihood
#'
#' @param model A CLM object created by clm_stn()
#' @return A fitted_clm_stn object containing the estimated parameters and model information
#' @details This function uses optim() with BFGS method to maximize the log-likelihood.
#'          The parameters are estimated in the order: thresholds (alpha), regression
#'          coefficients (beta), and asymmetry parameter (lambda).
fit_clm_stn <- function(model) {
  n_thresholds <- model$n_levels - 1

  # Valores iniciales más informativos basados en el modelo ordinal estándar
  clm_ref <- clm(rating ~ temp + contact, data=wine)
  init_alpha <- coef(clm_ref)[1:4]
  init_beta <- coef(clm_ref)[5:6]
  init_lambda <- 0.1

  init_params <- c(init_alpha, init_beta, init_lambda)

  # Valores iniciales
  n_thresholds <- model$n_levels - 1
  n_predictors <- model$n_vars

  # Escala 1.0 para umbrales, 0.1 para predictores, 0.01 para lambda
  parscale <- c(rep(1.0, n_thresholds),
                rep(0.1, n_predictors),
                0.01)

  # Optimización con más información
  opt <- optim(
    par = init_params,
    fn = loglik_clm_stn,
    gr = score_clm_stn,
    method = "BFGS",
    control = list(
      fnscale = -1,
      maxit = 2000,
      trace = 1,
      REPORT = 1,
      parscale = parscale  # Escalamiento para diferentes parámetros
    ),
    model = model
  )

  if(opt$convergence != 0) {
    warning("Optimization did not converge. Code: ", opt$convergence,
            "\nMessage: ", opt$message)
  }

  structure(
    list(
      coefficients = opt$par,
      loglik = opt$value,
      convergence = opt$convergence,
      message = opt$message,
      counts = opt$counts,
      init_values = init_params,
      gradient = score_clm_stn(opt$par, model),
      model = model,
      clm_ref = clm_ref
    ),
    class = "fitted_clm_stn"
  )
}
