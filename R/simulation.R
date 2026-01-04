#' Simulate a Continuous-Time HMM OU Bounded-Variation Model
#'
#' Simulates data from a two-state continuous-time hidden Markov model
#' with an Ornstein--Uhlenbeck drift process of bounded variation.
#' The model is given by
#'
#' \deqn{
#'   dY_t = \phi_t \, dt + \sigma \, dW_t,
#' }
#'
#' where the latent drift \eqn{\phi_t} follows
#'
#' \deqn{
#'   d\phi_t = (a_{\zeta_t} - \gamma_{\zeta_t} \phi_t)\,dt,
#' }
#'
#' and \eqn{\zeta_t} is a two-state continuous-time Markov chain with
#' generator matrix \eqn{Q}.
#'
#' @param T Terminal time of the simulation.
#' @param N Number of time steps.
#' @param Q A \eqn{2 \times 2} generator matrix of the Markov chain.
#' @param a Numeric vector of length 2 giving the drift intercepts.
#' @param gamma Numeric vector of length 2 giving mean-reversion rates.
#' @param phi0 Initial value of the latent drift \eqn{\phi_0}.
#' @param sigma Observation noise volatility.
#'
#' @return A list containing:
#' \describe{
#'   \item{time}{Time grid of length \eqn{N+1}.}
#'   \item{dt}{Time increments.}
#'   \item{zeta}{Simulated Markov chain states.}
#'   \item{phi}{Latent drift process.}
#'   \item{dY}{Observed increments.}
#'   \item{Y}{Observed cumulative process.}
#' }
#'
#' @examples
#' Q <- matrix(c(-1, 1, 0.5, -0.5), 2, 2, byrow = TRUE)
#' sim <- simulate_ctHOUbv(
#'   T = 10, N = 1000, Q = Q,
#'   a = c(3, 0), gamma = c(3, 4),
#'   phi0 = 0.5, sigma = 0.2
#' )
#'
#' @export
simulate_ctHOUbv <- function(
    T = 10,
    N = 1000,
    Q,
    a,
    gamma,
    phi0 = 0.5,
    sigma = 0.2
) {
  # --- Time grid ---
  time  <- seq(0, T, length.out = N + 1)
  dt <- diff(time)
  
  # --- Simulate Markov chain ---
  zeta <- integer(N + 1)
  zeta[1] <- 1
  
  for (n in 2:(N + 1)) {
    i <- zeta[n - 1]
    p_switch <- 1 - exp(Q[i, i] * dt[n - 1])
    zeta[n] <- if (runif(1) < p_switch) 3 - i else i
  }
  
  # --- Simulate phi process ---
  phi <- numeric(N + 1)
  phi[1] <- phi0
  
  for (n in 2:(N + 1)) {
    i <- zeta[n - 1]
    drift <- a[i] - gamma[i] * phi[n - 1]
    phi[n] <- phi[n - 1] + drift * dt[n - 1]
  }
  
  # --- Simulate observations ---
  dY <- phi[-(N + 1)] * dt +
    sigma * sqrt(dt) * rnorm(N)
  
  Y <- cumsum(c(0, dY))
  
  list(
    time  = time,
    dt = dt,
    zeta  = zeta,
    phi   = phi,
    dY    = dY,
    Y     = Y
  )
}

#' Plot Continuous-Time Hidden OUbv Model Fit
#'
#' Visualizes observed data together with the fitted latent bounded-variation
#' OU drift component from a continuous-time hidden OU model. Optionally overlays
#' the true latent process if it is available in the simulation object.
#'
#' The observed series (increments or levels) is shown on the left axis in black.
#' The fitted latent drift is overlaid on the right axis in red (solid line).
#' If \code{sim$phi} is present, the true latent drift is also shown on the right
#' axis using a red dashed line. Vertical dashed blue lines indicate inferred
#' regime switches.
#'
#' @param sim A list containing simulated or observed data. Must include
#'   \code{time} and either \code{dY} (increments) or \code{Y} (levels).
#'   If present, \code{phi} is interpreted as the true latent drift.
#' @param fit A fitted model object containing the estimated latent drift
#'   (\code{fit$phi}) and regime sequence
#'   (\code{fit$MC_res$vit_base[,4]}).
#' @param type Character string specifying whether to plot observed
#'   \code{"increments"} or \code{"level"} data.
#'
#' @details
#' The function uses base R graphics and overlays the latent drift on a
#' secondary (right-hand) axis. Legends for the latent drift are added
#' automatically and adapt depending on whether the true latent process
#' is available.
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of
#'   producing a plot.
#'
#' @export


plot_ctHOUbv <- function(
    sim,
    fit,
    type = c("increments", "level")
) {
  type <- match.arg(type)
  
  time <- sim$time[-1]
  
  if (type == "increments") {
    y_obs <- sim$dY
    ylab  <- "Observed increments"
  } else {
    y_obs <- sim$Y[-1]
    ylab  <- "Observed level"
  }
  
  phi_hat  <- head(fit$phi, -1)
  zeta_hat <- fit$MC_res$vit_base[, 4]
  
  has_phi <- !is.null(sim$phi)
  if (has_phi) {
    phi_true <- head(sim$phi, -1)
  }
  
  plot(
    time, y_obs, type = "l", col = "black",
    xlab = "Time", ylab = ylab,
    main = "Fitted model"
  )
  
  par(new = TRUE)
  plot(
    time, phi_hat, type = "l", col = "red", lwd = 2,
    axes = FALSE, xlab = "", ylab = ""
  )
  
  if (has_phi) {
    lines(
      time, phi_true,
      col = "red", lwd = 2, lty = 2
    )
  }
  
  axis(
    side = 4, col = "red", col.axis = "red"
  )
  mtext("Latent drift", side = 4, line = 3, col = "red")
  
  # Regime switches
  switches <- which(diff(zeta_hat) != 0)
  if (length(switches) > 0) {
    abline(v = time[switches], col = "blue", lty = 2)
  }
  
  legend_labels <- "Fitted latent OUbv"
  legend_lty    <- 1
  
  if (has_phi) {
    legend_labels <- c("Fitted latent OUbv", "True latent OUbv")
    legend_lty    <- c(1, 2)
  }
  
  legend(
    "topright",
    legend = legend_labels,
    col    = "red",
    lty    = legend_lty,
    lwd    = 2,
    bty    = "n"
  )
}
