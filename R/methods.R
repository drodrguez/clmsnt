print.fitted_clm_stn <- function(x, ...) {
  cat("CLM con enlace skew-t-normal\n\n")

  n_thresholds <- x$model$n_levels - 1
  alpha <- x$coefficients[1:n_thresholds]
  beta <- x$coefficients[(n_thresholds + 1):(n_thresholds + x$model$n_vars)]
  lambda <- x$coefficients[length(x$coefficients)]

  cat("Umbrales:\n")
  for(i in 1:n_thresholds) {
    cat(sprintf("alpha_%d: %f\n", i, alpha[i]))
  }

  cat("\nCoeficientes:\n")
  colnames <- colnames(x$model$X)
  for(i in 1:length(beta)) {
    cat(sprintf("%s: %f\n", colnames[i], beta[i]))
  }

  cat(sprintf("\nParámetro de asimetría (lambda): %f\n", lambda))

  cat(sprintf("\nLog-verosimilitud: %f\n", x$loglik))

  cat(sprintf("Convergencia: %s\n",
              ifelse(x$convergence == 0, "Exitosa", "No exitosa")))
}
