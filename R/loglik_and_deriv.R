#' @export
calc_loglik <- function(reg_coef, model){
  if (model$name == "linear"){
    return(calc_linear_loglik(reg_coef, model$design, model$outcome, model$noise_var))
  } else if (model$name == "logit"){
    return(calc_logit_loglik(reg_coef, model$design, model$outcome))
  }
}

#' @export
calc_loglink_deriv <- function(reg_coef, model, order = 1){
  if (model$name == "linear"){
    return(calc_linear_loglink_deriv(reg_coef, model$design, model$outcome, model$noise_var))
  } else if (model$name == "logit"){
    return(calc_logit_loglink_deriv(reg_coef, model$design, model$outcome, order))
  }
}

calc_grad <- function(reg_coef, model){
  design <- model$design
  loglink_deriv <- calc_loglink_deriv(reg_coef, model)
  grad <- t(design) %*% loglink_deriv # formula for canonical links
  grad <- as.vector(grad)
  return(grad)
}

calc_hessian <- function(reg_coef, model) {
  design <- model$design
  weight <- calc_loglink_deriv(reg_coef, model, order = 2)
  hess <- - t(design) %*% (outer(weight, rep(1, ncol(design))) * design)
  return(hess)
}

# calc_hessian_inverse <- function(reg_coef, model){
#   if (model$name == "linear"){
#     return(calc_linear_hessian_inverse(reg_coef, model$design, model$outcome, model$noise_var))
#   } else if (model$name == "logit"){
#     return(calc_logit_hessian_inverse(reg_coef, model$design, model$outcome))
#   }
# }

calc_linear_loglik <- function(reg_coef, design, outcome, noise_var = 1) {
  predicted_val <- design %*% reg_coef
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_linear_loglink_deriv <- function(reg_coef, design, outcome, noise_var = 1){
  predicted_val <- design %*% reg_coef
  deriv <- (outcome - predicted_val) / noise_var
  deriv <- as.vector(deriv)
  return(deriv)
}

# calc_linear_hessian_inverse <- function(reg_coef, design, outcome, noise_var = 1){
#   return(0)
# }

calc_logit_loglik <- function(reg_coef, design, outcome) {
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

calc_logit_hessian_inverse <- function(reg_coef, design, outcome) {
  weight <- calc_logit_loglink_deriv(reg_coef, design, outcome, order = 2)
  sqrt_weighted_design <- outer(sqrt(weight), rep(1, ncol(design))) * design
  R <- qr_wrapper(sqrt_weighted_design)$R
  inverse <- - invert_gram_mat_from_qr(R)
  return(inverse)
}
