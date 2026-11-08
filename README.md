# MATH3428 Group Project: Bayesian Volatility Modeling

Bayesian MCMC Inference for Financial Volatility (VOO vs. NVDA) using a heavy-tailed Student-$t$ framework.

---

## Team & Task Division

* **Member 1 (Framework & Sampler Lead):** CHEUNG YUK YING (21063768)
  * Custom Random Walk Metropolis MCMC sampler design with proposal tuning (`mcmc_sampler.R`).
  * Unconstrained parameter transformations ($\mu, \ln\sigma, \ln(\nu-2)$) and Jacobian corrections (`derivation_draft.md`).
  * Synthetic parameter recovery validation and baseline testing.

* **Member 2 (Market Benchmark Analysis):** [name] ([id])
  * S&P 500 ETF (VOO) log-return historical data acquisition & EDA.
  * MCMC convergence diagnostics (Trace plots, Gelman-Rubin $\hat{R}$, ACF).
  * Posterior inference and empirical QQ-plot validation for market benchmark.

* **Member 3 (High-Volatility Tech Asset & Risk Evaluation):** [name] ([id])
  * NVIDIA (NVDA) heavy-tailed volatility modeling via MCMC.
  * Comparative posterior tail analysis ($\nu_{\text{VOO}}$ vs. $\nu_{\text{NVDA}}$).
  * Bayesian Value-at-Risk (VaR) & Expected Shortfall (ES) risk estimation.

---

## Repository Structure

* `mcmc_sampler.R`: Modular Metropolis-Hastings sampling engine.
* `derivation_draft.md`: Mathematical formulation, log-likelihood, and Jacobian derivations.
* *(Upcoming)* Data scripts, convergence diagnostic reports, and final paper draft.
