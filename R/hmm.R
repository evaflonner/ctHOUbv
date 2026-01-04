library(msm)
#' Estimate Continuous-Time Hidden Markov Model
#'
#' Estimates a two-state continuous-time hidden Markov model
#' with Gaussian emissions.
#'
#' @param deltaY Numeric vector of observed increments.
#' @param time Observation times.
#' @param trans_init Initial generator matrix.
#' @param obs_noise Observation noise.
#' @param dt Time increment.
#' @param hidden_model0 Hidden model estimated in step 1 to get prior probabilities for the first state. Should be of form hmmNorm(mean,sd).
#' @param hidden_model1 Hidden model estimated in step 1 to get prior probabilities for the second state. Should be of form hmmNorm(mean,sd).
#'
#' @return A list containing the fitted HMM, estimated generator matrix,
#' and filtered state probabilities.
#'
#' @details
#' This function relies on the \pkg{msm} package to compute
#' likelihoods for continuous-time Markov chains using matrix exponentials.
#'
#' @seealso \code{\link[msm]{msm}}
#'
#' @export
#' 
step1_cthmm <- function(deltaY, time, trans_init, dt, hidden_model0, hidden_model1) {
  
  dat <- data.frame(
    deltaY = deltaY,
    time = time
  )
  
  # --- define HMM ---
  hmm_res <- tryCatch(
    {msm(
    deltaY ~ time,
    subject = rep(1, length(time)),
    data = dat,
    qmatrix = trans_init,
    hmodel = list(hidden_model0, hidden_model1))
    },
    error = function(e){
      stop(e)
    })
  
  predict_base <- simfitted.msm(hmm_res)
  vit_base <- viterbi.msm(hmm_res)

  list(
    model = hmm_res,
    vit_base = vit_base
  )
}
