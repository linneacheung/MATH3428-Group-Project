# ==============================================================================
# Course: MATH 3428 Bayesian Statistics
# Project: Bayesian Volatility Inference via MCMC (Student-t Distribution)
# Script: mcmc_sampler.R
# Role: Member 1 - Framework & Sampler Lead
# Description: Random Walk Metropolis sampler with unconstrained parameter
#              transformations, Jacobian corrections, and proposal tuning.
# ==============================================================================

#' Target Unconstrained Log-Posterior Density
#' 
#' @param theta Vector of transformed parameters: c(mu, ln(sigma), ln(nu - 2))
#' @param y Vector of observed financial log-returns
#' @return Scalar value of the unnormalized log-posterior density
log_posterior <- function(theta, y) {
  mu <- theta[1]
  sigma <- exp(theta[2])
  nu <- 2 + exp(theta[3])
  n <- length(y)
  
  # 1. Student-t Log-Likelihood (using lgamma for numerical stability)
  term1 <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * nu) - log(sigma)
  z <- (y - mu) / sigma
  log_lik <- n * term1 - ((nu + 1) / 2) * sum(log(1 + (z^2) / nu))
  
  # 2. Prior Distributions:
  # mu ~ Normal(0, 100)
  # sigma ~ Half-Cauchy(scale = 5)
  # nu ~ Gamma(shape = 2, rate = 0.1)
  log_prior_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)
  log_prior_sigma <- log(2 / (pi * 5 * (1 + (sigma / 5)^2)))
  log_prior_nu <- dgamma(nu, shape = 2, rate = 0.1, log = TRUE)
  
  # 3. Log-Jacobian adjustment for unconstrained mappings
  # d(sigma)/d(theta2) = exp(theta2), d(nu)/d(theta3) = exp(theta3)
  log_jacobian <- theta[2] + theta[3]
  
  return(log_lik + log_prior_mu + log_prior_sigma + log_prior_nu + log_jacobian)
}

#' Random Walk Metropolis Sampler for Student-t Model
#' 
#' @param y Vector of observed log-returns
#' @param n_iter Total number of MCMC iterations (default: 10000)
#' @param burn_in Number of burn-in samples to discard (default: 2000)
#' @param init_theta Initial parameter vector on transformed scale (optional)
#' @param prop_sd Vector of proposal standard deviations for random walk
#' @return A list containing posterior samples, acceptance rate, and full chain
run_student_t_mcmc <- function(y, n_iter = 10000, burn_in = 2000, 
                               init_theta = NULL, 
                               prop_sd = c(0.0006, 0.03, 0.25)) {
  
  # Set empirical initial state if not explicitly provided
  if (is.null(init_theta)) {
    init_theta <- c(mean(y), log(sd(y)), log(3))
  }
  
  samples <- matrix(NA, nrow = n_iter, ncol = 3)
  colnames(samples) <- c("mu", "sigma", "nu")
  
  curr_theta <- init_theta
  curr_log_post <- log_posterior(curr_theta, y)
  n_accepted <- 0
  
  for (i in 1:n_iter) {
    # Generate proposal from symmetric multivariate normal
    prop_theta <- rnorm(3, mean = curr_theta, sd = prop_sd)
    prop_log_post <- log_posterior(prop_theta, y)
    
    # Metropolis acceptance probability on log scale
    log_alpha <- min(0, prop_log_post - curr_log_post)
    if (log(runif(1)) < log_alpha) {
      curr_theta <- prop_theta
      curr_log_post <- prop_log_post
      n_accepted <- n_accepted + 1
    }
    
    # Transform back to original parameter space
    samples[i, 1] <- curr_theta[1]
    samples[i, 2] <- exp(curr_theta[2])
    samples[i, 3] <- 2 + exp(curr_theta[3])
  }
  
  acc_rate <- n_accepted / n_iter
  post_samples <- samples[(burn_in + 1):n_iter, ]
  
  return(list(
    post_samples = as.data.frame(post_samples),
    acceptance_rate = acc_rate,
    all_chains = samples
  ))
}

# ==============================================================================
# Validation Block (Simulated Ground Truth Recovery)
# ==============================================================================
if (sys.nframe() == 0) {
  set.seed(42)
  sim_y <- 0.001 + 0.015 * rt(1000, df = 4)
  res <- run_student_t_mcmc(sim_y, n_iter = 5000, burn_in = 1000)
  
  cat("--- MCMC Diagnostic Check ---\n")
  cat("Acceptance Rate:", round(res$acceptance_rate * 100, 2), "%\n")
  cat("Target Range: 25% - 45%\n")
}