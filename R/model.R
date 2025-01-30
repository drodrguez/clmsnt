# model.R

#' Create a Cumulative Link Model object with skew-t-normal link
#'
#' @param formula A formula object
#' @param data A data frame containing the variables
#' @param lambda Initial value for asymmetry parameter
#' @param v Degrees of freedom (fixed)
#' @return A CLM object

clm_stn <- function(formula, data, lambda = 0, v = 5) {
  # Extraer la variable respuesta y la matriz de diseño
  mf <- model.frame(formula, data)
  y <- model.response(mf)

  # Crear matriz de diseño sin intercepto
  X <- model.matrix(formula, data)[,-1, drop=FALSE]  # Asegurarnos de eliminar el intercepto

  # Verificar que y sea ordinal
  if (!is.ordered(y)) {
    stop("Response must be an ordered factor")
  }

  # Añadir información sobre la distribución de categorías
  category_props <- prop.table(table(y))
  init_cuts <- qnorm(cumsum(category_props)[-length(category_props)])

  structure(
    list(
      y = y,
      X = X,
      levels = levels(y),
      n_levels = length(levels(y)),
      n_obs = nrow(X),
      n_vars = ncol(X),
      lambda = lambda,
      v = v,
      init_cuts = init_cuts,  # Añadir valores iniciales informados
      formula = formula,
      call = match.call()
    ),
    class = "clm_stn"
  )
}


#' Log-likelihood function for CLM with skew-t-normal link
#'
#' @param theta Parameter vector (alpha, beta, lambda)
#' @param model A CLM object
#' @return Log-likelihood value
loglik_clm_stn <- function(theta, model) {
  # Extraer parámetros
  n_thresholds <- model$n_levels - 1
  alpha <- theta[1:n_thresholds]
  beta <- theta[(n_thresholds + 1):(n_thresholds + model$n_vars)]
  lambda <- theta[length(theta)]

  # Calcular predictor lineal
  eta <- model$X %*% beta

  # Calcular probabilidades para cada categoría
  probs <- matrix(0, nrow = model$n_obs, ncol = model$n_levels)

  for (j in 1:model$n_levels) {
    if (j == 1) {
      probs[,j] <- pskewt(alpha[j] - eta, lambda = lambda, v = model$v)
    } else if (j == model$n_levels) {
      probs[,j] <- 1 - pskewt(alpha[j-1] - eta, lambda = lambda, v = model$v)
    } else {
      probs[,j] <- pskewt(alpha[j] - eta, lambda = lambda, v = model$v) -
        pskewt(alpha[j-1] - eta, lambda = lambda, v = model$v)
    }
  }

  # Calcular log-verosimilitud
  y_mat <- model.matrix(~factor(model$y) - 1)
  loglik <- sum(y_mat * log(probs))

  return(loglik)
}
