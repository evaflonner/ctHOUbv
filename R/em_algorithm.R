#' Fit a Continuous-Time Hidden Ornstein-Uhlenbeck process of bounded variation
#'
#' Fits a continuous-time hidden OUbv model
#' using an EM-style iterative procedure.
#' The latent Markov chain is estimated in continuous time, and
#' the bounded-variation process is inferred using filtered
#' state probabilities.
#'
#' @param deltaY Numeric vector of observed increments \eqn{\Delta Y_{t_n}}.
#' @param time Numeric vector of observation times \eqn{t_n}.
#' @param h Linear drift.
#' @param trans_init Initial generator matrix for CTMC.
#' @param gamma0_init Initial value for gamma0.
#' @param gamma1_init Initial value for gamma1.
#' @param phi0_init Initial value of the process.
#' @param sigma_s Observation noise volatility.
#' @param dt Time increment.
#' @param hidden_model0 Hidden model estimated in step 1 to get prior probabilities for the first state. Should be of form hmmNorm(mean,sd).
#' @param hidden_model1 Hidden model estimated in step 1 to get prior probabilities for the second state. Should be of form hmmNorm(mean,sd).
#' @param mode Either "mixture_probs" for "MAP_probs", meaning that for each time point we consider the most likely state. 
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance for the negative log-likelihood.
#'
#' @return A list containing:
#' \describe{
#'   \item{gamma0_estim}{Estimated parameter gamma0.}
#'   \item{gamma1_estim}{Estimated parameter gamma1.}
#'   \item{phi0_estim}{Estimated initial state.}
#'   \item{phi}{Estimated hidden process.}
#'   \item{MC_res}{Estimated Markov chain.}
#' }
#' 
#' @details
#' The algorithm alternates between:
#' \enumerate{
#'   \item Estimating the continuous-time hidden Markov chain using
#'   a Gaussian emission model.
#'   \item Updating the bounded-variation process and optimizing
#'   OU parameters via likelihood maximization.
#' }
#'
#'
#' @references
#' Baum, L. E., et al. (1970). A maximization technique occurring in
#' the statistical analysis of probabilistic functions of Markov chains.
#'
#' @export
#' 
fit_cthoubv <- function(
    deltaY, time, h,
    trans_init,
    gamma0_init, gamma1_init, phi0_init,
    sigma_s, dt,
    hidden_model0,
    hidden_model1,
    mode = "mixture_probs",
    max_iter = 50, tol = 1e-4
) {
  
  if (!(mode %in% c("mixture_probs","MAP_probs"))){
    stop("Mode not recognized.")
  }
  
  gamma0_current <- gamma0_init
  gamma1_current <- gamma1_init
  phi0_current <- phi0_init
  loglik_old <- Inf
  
  for (iter in 1:max_iter) {
    
    # Step 1
    step1 <- step1_cthmm(deltaY, time, trans_init, dt, hidden_model0, hidden_model1)
    
    if(mode == "mixture_probs"){
      prob0 <- step1$vit_base[,5][,1]
      prob1 <- step1$vit_base[,5][,2]
    } else {
      prob0 <- as.integer(step1$vit_base[,4] == 1)
      prob1 <- as.integer(step1$vit_base[,4] == 2)
    }
    
    # Step 2
    
    opt <- estimate_houbv_params(
        deltaY      = deltaY,
        prob0       = prob0,
        prob1       = prob1,
        gamma0_init = gamma0_current,
        gamma1_init = gamma1_current,
        phi0_init   = phi0_current,
        h           = h,
        dt       = dt,
        sigma_s     = sigma_s
      )
    
    par_optim <- opt$par
    phi <- compute_oubv(prob0, prob1, par_optim[1], par_optim[2], par_optim[3], dt)
    
    loglik_new <- opt$value
    cat("Iter", iter, "NegLL:", loglik_new, "\n")
    
    gamma0_current <- par_optim[1]
    gamma1_current <- par_optim[2]
    phi0_current <- par_optim[3]
    
    if (abs(loglik_old - loglik_new) < tol) break
    loglik_old <- loglik_new
  }
  
  list(
    gamma0_estim = gamma0_current,
    gamma1_estim = gamma1_current,
    phi0_estim = phi0_current,
    phi = phi,
    MC_res = step1
  )
}
