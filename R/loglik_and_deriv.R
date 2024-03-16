#' @export
calc_loglik <- function(reg_coef, design, outcome, model_name, ...){
  if (model_name == "linear"){
    return(calc_linear_loglik(reg_coef, design, outcome, noise_var))
  } else if (model_name == "logit"){
    return(calc_logit_loglik(reg_coef, design, outcome))
  }
}

calc_linear_loglik <- function(reg_coef, design, outcome, noise_var) {
  predicted_val <- design %*% reg_coef
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_linear_grad <- function(reg_coef, design, outcome, noise_var) {
  predicted_val <- design %*% reg_coef
  grad <- t(design) %*% (outcome - predicted_val) / noise_var
  grad <- as.vector(grad)
  return(grad)
}

calc_logit_loglik <- function(
    reg_coef, design, outcome
) {
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

calc_logit_grad <- function(reg_coef, design, outcome) {
  loglink_grad <- calc_logit_loglink_deriv(reg_coef, design, outcome, order = 1)
  grad <- t(design) %*% loglink_grad
  grad <- as.vector(grad)
  return(grad)
}

calc_logit_loglink_deriv <- function(reg_coef, design, outcome, order) {
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

calc_logit_hessian <- function(reg_coef, design, outcome) {
  weight <- calc_logit_loglink_deriv(reg_coef, design, outcome, order = 2)
  hess <- - t(design) %*% (outer(weight, rep(1, ncol(design))) * design)
  return(hess)
}

calc_logit_hessian_inverse <- function(reg_coef, design, outcome) {
  weight <- calc_logit_loglink_deriv(reg_coef, design, outcome, order = 2)
  sqrt_weighted_design <- outer(sqrt(weight), rep(1, ncol(design))) * design
  R <- qr_wrapper(sqrt_weighted_design)$R
  inverse <- - invert_gram_mat_from_qr(R)
  return(inverse)
}
