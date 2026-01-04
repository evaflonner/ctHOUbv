#' Log-Likelihood for Bounded-Variation OU Model
#'
#' Computes the negative log-likelihood of the observation process
#' given the bounded-variation OU dynamics.
#'
#' @param deltaY Observed increments.
#' @param phi Estimated bounded-variation process.
#' @param h Linear drift.
#' @param dt Time increments.
#' @param sigma_s Observation noise volatility.
#'
#' @return Numeric value of the negative log-likelihood.
#'
#' @export

loglik_phi <- function(deltaY, phi, h, dt, sigma_s) {
  mu_tilde <- h * phi * dt
  -sum(dnorm(deltaY, mean = mu_tilde, sd = sigma_s * sqrt(dt), log = TRUE))
}

#' Estimate OU Parameters
#'
#' Estimates the parameters of the bounded-variation
#' OU process via likelihood maximization.
#'
#' @param deltaY Observed increments.
#' @param prob0 Mixture state 1 probabilities.
#' @param prob1 Mixture state 2 probabilities.
#' @param gamma0_init Initial value for gamma0.
#' @param gamma1_init Initial value for gamma1.
#' @param phi0_init Initial value of the process.
#' @param h Linear drift.
#' @param dt Time increment.
#' @param sigma_s Observation noise volatility.
#'
#' @return An \code{optim} object containing parameter estimates.
#'
#' @export
estimate_houbv_params <- function(
    deltaY, prob0, prob1, gamma0_init, gamma1_init, phi0_init, h, dt, sigma_s
) {
  obj <- function(par, prob0, prob1, dt, deltaY, h, sigma_s) {
    gamma0 <- par[1]
    gamma1 <- par[2]
    phi0 <- par[3]
    phi <- compute_oubv(prob0, prob1, gamma0, gamma1, phi0, dt)

    loglik_phi(deltaY, phi, h, dt, sigma_s)
  }
  
  optim(
    par = c(gamma0_init, gamma1_init, phi0_init),
    fn = obj,
    prob0 = prob0,
    prob1 = prob1,
    dt = dt,
    h = h,
    sigma_s = sigma_s,
    deltaY = deltaY,
    method = "L-BFGS-B",
    lower = c(0, 0)
  )
}
