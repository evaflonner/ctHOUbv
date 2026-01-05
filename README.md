
# ctHOUbv

**Continuous-Time Ornstein Uhlenbeck process of bounded variation under
partial information**

------------------------------------------------------------------------

## Overview

`ctHOUbv` is an R package for **estimation and filtering of
continuous-time models under partial information** in which the
observation process is driven by a **Ornstein Uhlenbeck process of
bounded variation (OUbv)** whose dynamics depend on an unobserved
finite-state Markov chain.

The package is designed for applications in which:

- the latent state must remain within a **fixed interval**
  (e.g. portfolio weights, intensities, proportions),
- regime switching occurs in **continuous time**, and
- inference must be performed using **partial and noisy observations**.

The methodology combines:

1.  Continuous-time hidden Markov chain estimation,
2.  Filtering of the OUbv process,
3.  An EM-style iterative likelihood maximization procedure.

------------------------------------------------------------------------

## Model Specification

### Observation equation

Let $(Y_t)_{t \ge 0}$ be the observed process. The model assumes

$$
dY_t = \phi_t , dt + \sigma , dW_t,
$$

where:

- $W_t$ is a standard Brownian motion,
- $\sigma > 0$ is a constant volatility parameter,
- $\phi_t$ is an unobserved bounded-variation process.

------------------------------------------------------------------------

### Latent state dynamics

The latent process $(\phi_t)_{t \ge 0}$ evolves according to

$$
d\phi_t = C(\phi_t, \zeta_t),dt,
$$

where $$(\zeta_t)_{t \ge 0}$$ is a **two-state continuous-time Markov
chain** with state space $$\Sigma = {e_1, e_2}$$ and generator matrix
$Q$.

The drift function is given by

$$
C(\phi_t,\zeta_t) =
\begin{cases}
a_{e_1} - \gamma_{e_1}\phi_t,\text{ if } \zeta_t = e_1, \\
a_{e_2} - \gamma_{e_2}\phi_t, \text{ if } \zeta_t = e_2,
\end{cases}
$$

with parameters $a_{e_i}, \gamma_{e_i} \in \mathbb{R}$.

------------------------------------------------------------------------

### Bounded-variation property

If 
$$\frac{a_{e_2}}{\gamma_{e_2}} < \frac{a_{e_1}}{\gamma_{e_1}}$$
and
$$\phi_0 \in \left(\frac{a_{e_2}}{\gamma_{e_2}}, \frac{a_{e_1}}{\gamma_{e_1}}\right),$$
then the solution $\phi_t$ remains in this interval for all
$t \ge 0$.

In the default specification used for weights, $a_{e_2} = 0$, $a_{e_1} = \gamma_{e_1}$, $\phi_t \in (0,1)$.

------------------------------------------------------------------------

## Estimation Methodology

The package implements a **stepwise EM-style procedure**.

### Step 1: Continuous-time HMM estimation

A simplified Gaussian emission model is first estimated:

$$
\Delta Y_{t_n} \mid \zeta_{t_n} = i
\sim \mathcal{N}(m_i \Delta_n, \sigma^2 \Delta_n),
$$

yielding:

- an estimate of the generator matrix $Q$,
- filtered (or smoothed) state probabilities
  $\hat p_{i,t_n} = \mathbb{P}(\zeta_{t_n} = i \mid \mathcal{F}_{t_n}^Y)$.

This step uses likelihood-based methods for continuous-time Markov
chains.

------------------------------------------------------------------------

### Step 2: Filtering and parameter estimation for (\_t)

Given $\hat p_{i,t_n}$, the latent process is approximated by

$$
\hat \phi_{t_n}
\approx
\hat \phi_{t_{n-1}}
+
\sum_{i \in {0,1}}
\hat p_{i,t_n}
\left(a_i - \gamma_i \hat \phi_{t_{n-1}}\right)\Delta_n.
$$

The parameters $(\gamma_0,\gamma_1)$ are estimated by maximizing the
likelihood

$$
\prod_{n=1}^N
\mathcal{N}\left(
\Delta Y_{t_n};
\hat \phi_{t_n} \Delta_n,
\sigma^2 \Delta_n
\right).
$$

Steps 1 and 2 may be iterated until convergence of the log-likelihood.

------------------------------------------------------------------------

## Package Structure

    ctHOUbv/
    ├── R/
    │   ├── simulation.R     # Data generation and plotting functionality
    │   ├── hmm.R         # CT-HMM estimation
    │   ├── oubv.R         # OUbv filtering
    │   ├── loglik_glm.R        # Likelihood evaluation
    │   └── em_algorithm.R     # EM-style iteration
    ├── DESCRIPTION
    ├── NAMESPACE
    └── README.md

------------------------------------------------------------------------

## Installation

From GitHub:

``` r
# install.packages("devtools")
devtools::install_github("evaflonner/ctHOUbv")
```

------------------------------------------------------------------------

## Basic Usage

``` r
library(ctHOUbv)

set.seed(3012026)
Q <- matrix(c(-1, 1, 0.5, -0.5), 2, 2, byrow = TRUE)
sim <- simulate_ctHOUbv(
     T = 10, N = 1000, Q = Q,
     a = c(3, 0), gamma = c(3, 4),
     phi0 = 0.5, sigma = 0.02
   )

res <-fit_cthoubv(
  deltaY = sim$dY, time = head(sim$time,-1), 
  h = 1, trans_init = rbind( c(0, 3), c(3, 0)), 
  gamma0_init = runif(1,1,5), gamma1_init = runif(1,2,6), phi0_init = runif(1,0.3,0.7), 
  sigma_s = 0.02, dt = sim$dt[1],mu_hidden_model = 5e-4, sigma_hidden_model = 1e-2)

plot_ctHOUbv(sim,res)
```

------------------------------------------------------------------------

## Relation to Existing Methods

- Related to **continuous-time HMMs** and the Baum–Welch algorithm.

- Closely connected to **Interacting Multiple Model (IMM)** filters.

- Differs from IMM approaches in that:

  - the latent continuous state has **bounded support**,
  - no Gaussian approximation of $\phi_t$ is used.

------------------------------------------------------------------------

## Intended Use

This package is intended for:

- methodological research,
- simulation studies,
- prototyping estimation procedures for OUbv processes under partial
  information.

It is **not** optimized for high-frequency or large-scale production
use.

------------------------------------------------------------------------

## References

- Baum, L. E., et al. (1970). *A maximization technique occurring in the
  statistical analysis of probabilistic functions of Markov chains*.
- Jackson, C. (2011). *Multi-state models for panel data: the msm
  package for R*.
- Blom, H. A. P., & Bar-Shalom, Y. (1988). *The interacting multiple
  model algorithm for systems with Markovian switching coefficients*.
- Ratanov, N. (2021). *Ornstein-Uhlenbeck processes of bounded
  variation*. Methodology and Computing in Applied Probability, 23(3),
  925-946.
- Flonner, Eva, and Zehra Eksi. (2024). *Riding the Waves? Dissecting
  Momentum and Mean Reversion in Asset Prices Using Stochastic Filtering
  of the Ornstein-Uhlenbeck Process of Bounded Variation.*

------------------------------------------------------------------------

## License

MIT License.

------------------------------------------------------------------------
