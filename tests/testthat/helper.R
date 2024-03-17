are_all_close <- function(v, w, abs_tol = 1e-6, rel_tol = 1e-6) {
  if (!(is.numeric(v) && is.numeric(w))) {
    stop("Input args must be numeric.")
  }
  if (length(v) != length(w)) {
    warning("Lengths of the input vectors do not match.")
    return(FALSE)
  }
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  are_all_within_rtol <- all(abs_diff < rel_tol * pmax(abs(v), abs(w)))
  return(are_all_within_atol && are_all_within_rtol)
}

approx_grad_via_finite_diff <- function(func, x, dx = 1e-6) {
  numerical_grad <- rep(0, length(x))
  for (i in 1:length(x)) {
    x_plus <- x 
    x_plus[i] <- x[i] + dx
    x_minus <- x
    x_minus[i] <- x[i] - dx
    numerical_grad[i] <- (func(x_plus) - func(x_minus)) / (2 * dx)
  }
  return(numerical_grad)
}

simulate_data <- function(
    n_obs, n_pred, model_name = "linear", intercept = NULL, 
    coef_true = NULL, design = NULL, seed = NULL, option = list()
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if ((model_name != "linear")  && !is.null(option$signal_to_noise)) {
    warning(paste(
      "The `signal_to_noise` option is currently unsupported for",
      "non-linear models and will be ignored."
    ))
  }
  if (is.null(coef_true)) {
    coef_true <- rnorm(n_pred, sd = 1 / sqrt(n_pred))
  } 
  if (is.null(design)) {
    design <- matrix(rnorm(n_obs * n_pred), nrow = n_obs, ncol = n_pred)
  }
  if (!is.null(intercept)) {
    if (!is.numeric(intercept)) {
      stop("The intercept argument must be numeric.")
    }
    coef_true <- c(intercept, coef_true)
    design <- cbind(rep(1, n_obs), design)
  }
  linear_pred <- as.vector(design %*% coef_true)
  if (model_name == 'linear') {
    expected_mean <- linear_pred
    signal_to_noise <- option$signal_to_noise
    if (is.null(signal_to_noise)) {
      signal_to_noise <- 0.1
    }
    noise_magnitude <- sqrt(var(expected_mean) / signal_to_noise^2)
    noise <- noise_magnitude * rnorm(n_obs)
    outcome <- expected_mean + noise 
  } else if (model_name == 'logit') {
    n_trial <- option$n_trial
    prob <- 1 / (1 + exp(-linear_pred))
    # Object type of `outcome` returned by this function is variable. One should 
    # in general be careful about introducing this type of inconsistency, but 
    # sometimes one might find it the most natural and/or reasonable thing to do.
    if (is.null(n_trial)) {
      outcome <- rbinom(n_obs, 1, prob)
    } else {
      n_success <- rbinom(n_obs, n_trial, prob)
      outcome <- list(n_success = n_success, n_trial = n_trial)
    }
  } else if (model_name == 'poisson') {
    expected_mean <- exp(linear_pred)
    outcome <- rpois(n_obs, expected_mean)
  } else stop(paste(model_name, "model is not currently supported;",
                    "use 'linear', 'logit', or 'poisson'"))
  return(list(design = design, outcome = outcome, coef_true = coef_true, model_name = model_name))
}
