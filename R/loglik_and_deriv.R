calc_loglik <- function(model, reg_coef, ...){
  UseMethod("calc_loglik")
}

calc_loglink_deriv <- function(model, reg_coef, order, ...){
  UseMethod("calc_loglink_deriv")
}

calc_grad <- function(model, reg_coef){
  design <- model$design
  loglink_deriv <- calc_loglink_deriv(model, reg_coef, order = 1)
  grad <- t(design) %*% loglink_deriv # formula for canonical links
  grad <- as.vector(grad)
  return(grad)
}

calc_hessian <- function(model, reg_coef) {
  design <- model$design
  loglink_deriv <- calc_loglink_deriv(model, reg_coef, order = 2)
  hess <- - t(design) %*% (outer(loglink_deriv, rep(1, ncol(design))) * design)
  # hess <- - t(design) %*% (loglink_deriv * design) # to do: check if this is correct (if so, use b/c more readable than above)
  return(hess)
}

calc_hessian_inverse <- function(model, reg_coef){
  if (model$name == "linear"){
    return(calc_linear_hessian_inverse(reg_coef, model$design, model$outcome, model$noise_var))
  } else{
    return(calc_nonlinear_hessian_inverse(model, reg_coef))
  }
}

calc_loglik.linear_model <- function(model, reg_coef) {
  design <- model$design
  outcome <- model$outcome
  noise_var <- model$noise_var
  predicted_val <- design %*% reg_coef
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_loglink_deriv.linear_model <- function(model, reg_coef, order){
  design <- model$design
  outcome <- model$outcome
  noise_var <- model$noise_var
  if (order == 1) {
    predicted_val <- design %*% reg_coef
    deriv <- (outcome - predicted_val) / noise_var
  } else if (order == 2) {
    deriv <- 1 / noise_var
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_linear_hessian_inverse <- function(reg_coef, design, outcome, noise_var = 1){
  R <- solve_least_sq_via_qr_cpp_eig(design, outcome)$R
  unweighted_inverse <- - invert_gram_mat_from_qr(R)
  n_obs <- nrow(design); n_pred <- ncol(design)
  noise_var <- mean((outcome - design %*% reg_coef)^2) /
    (1 - n_pred / n_obs)
  inverse <- noise_var * unweighted_inverse
  return(inverse)
}

calc_loglik.logit_model <- function(model, reg_coef) {
  design <- model$design
  outcome <- model$outcome
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }
  logit_prob <- design %*% reg_coef
  loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
    # TODO: improve numerical stability for logit_prob >> 1
  return(loglik)
}

calc_loglink_deriv.logit_model <- function(model, reg_coef, order) {
  design <- model$design
  outcome <- model$outcome
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }
  logit_prob <- as.vector(design %*% reg_coef)
  predicted_prob <- 1 / (1 + exp(-logit_prob))
  if (order == 1) {
    deriv <- n_success - n_trial * predicted_prob
  } else if (order == 2) {
    deriv <- n_trial * predicted_prob * (1 - predicted_prob)
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_nonlinear_hessian_inverse <- function(model, reg_coef) {
  design <- model$design; outcome <- model$outcome
  weight <- calc_loglink_deriv(model, reg_coef, order = 2)
  sqrt_weighted_design <- outer(sqrt(weight), rep(1, ncol(design))) * design
  # sqrt_weighted_design <- sqrt(weight) * design # to do: check if this is correct (if so, use b/c more readable than above)
  R <- qr_wrapper(sqrt_weighted_design)$R
  inverse <- - invert_gram_mat_from_qr(R)
  return(inverse)
}

calc_loglik.poisson_model <- function(model, reg_coef){
  design <- model$design
  outcome <- model$outcome
  linear_pred <- design %*% reg_coef
  predicted_counts <- exp(linear_pred)
  loglik <- sum(outcome * linear_pred - predicted_counts)
  return(loglik)
}

calc_loglink_deriv.poisson_model <- function(model, reg_coef, order){
  design <- model$design
  outcome <- model$outcome
  predicted_counts <- exp(design %*% reg_coef)
  if (order == 1){
    deriv <- outcome - predicted_counts
  } else if (order == 2){
    deriv <- predicted_counts
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}
