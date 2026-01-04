#' Compute Bounded-Variation OU Process
#'
#' Computes the bounded-variation process \eqn{\phi_t}
#' using filtered state probabilities from a continuous-time
#' hidden Markov model.
#'
#' @param p Matrix of filtered state probabilities.
#' @param a Vector of OU drift parameters.
#' @param gamma Vector of mean-reversion parameters.
#' @param phi0 Initial value of the process.
#' @param dt Time increments.
#'
#' @return Numeric vector of estimated \eqn{\phi_t}.
#'
#' @details
#' The process evolves according to a discretized
#' OU-type dynamics driven by the hidden Markov states.
#'
#' @export
#' 
compute_oubv <- function(prob0, prob1, gamma0, gamma1, phi0, dt) {
  phi <- numeric(length(prob0))
  
  for(j in 1:length(prob0)){
    phi[j+1] <-phi[j]+prob0[j]*(gamma1*dt-gamma1*phi[j]*dt)+prob1[j]*(-gamma0*phi[j]*dt)
  }
  phi
}
