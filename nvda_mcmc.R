# ============================================
# NVDA: MCMC Estimation, VaR & ES
# MATH3428 Group Project - Group 1 - member 3
# ============================================

if (!require("quantmod")) install.packages("quantmod")
if (!require("coda"))     install.packages("coda")
library(quantmod)
library(coda)

cat("Working directory:", getwd(), "\n\n")
source("mcmc_sampler.R")

# 1. Data Preparation
nvda <- getSymbols("NVDA", from = "2024-01-01", to = "2026-10-01", auto.assign = FALSE)
y_vec <- na.omit(diff(log(Ad(nvda))))
y <- as.numeric(y_vec)

# 2. Run 3 MCMC chains
cat("Running 3 MCMC chains...\n")
set.seed(123)
prop_sd_new <- c(0.0006, 0.03, 0.25) * 1.25

init_1 <- c( 0.000, log(0.005), log(5))
init_2 <- c( 0.001, log(0.020), log(15))
init_3 <- c(-0.001, log(0.010), log(40))

chain1 <- run_student_t_mcmc(y, n_iter = 20000, burn_in = 5000, init_theta = init_1, prop_sd = prop_sd_new)
chain2 <- run_student_t_mcmc(y, n_iter = 20000, burn_in = 5000, init_theta = init_2, prop_sd = prop_sd_new)
chain3 <- run_student_t_mcmc(y, n_iter = 20000, burn_in = 5000, init_theta = init_3, prop_sd = prop_sd_new)

acc <- c(chain1$acceptance_rate, chain2$acceptance_rate, chain3$acceptance_rate)
cat("Acceptance rates:", round(acc, 4), 
    "| All in range 25%-45%:", all(acc >= 0.25 & acc <= 0.45), "\n\n")

# 3. Convergence diagnostics
mcmc_list <- mcmc.list(
  mcmc(as.matrix(chain1$post_samples)),
  mcmc(as.matrix(chain2$post_samples)),
  mcmc(as.matrix(chain3$post_samples))
)

cat("--- Gelman-Rubin R-hat ---\n")
print(gelman.diag(mcmc_list, autoburnin = FALSE))
cat("\n--- Effective Sample Size (ESS) ---\n")
print(effectiveSize(mcmc_list))

# Trace plots
png("fig_nvda_trace_plot.png", width = 1000, height = 700, res = 120)
traceplot(mcmc_list, main = "Trace Plots for NVDA Student-t Parameters")
dev.off()

# Combine posterior samples
post <- do.call(rbind, lapply(list(chain1, chain2, chain3), function(ch) as.matrix(ch$post_samples)))

# 4. Posterior inference & Fitted Student-t QQ plot
nu_median <- median(post[, "nu"])
png("fig_nvda_post_nu.png", width = 800, height = 500, res = 120)
plot(density(post[, "nu"]), main = "Posterior Distribution of nu", 
     xlab = "nu", col = "blue", lwd = 2)
abline(v = nu_median, col = "red", lty = 2)
dev.off()

cat("\nPosterior median nu =", round(nu_median, 3), 
    "| 95% CI:", round(quantile(post[, "nu"], c(0.025, 0.975)), 3), "\n\n")

# Add Fitted Student-t QQ plot
mu_hat    <- mean(post[, "mu"])
sigma_hat <- mean(post[, "sigma"])
nu_hat    <- mean(post[, "nu"])
p <- ppoints(length(y))
theoretical <- mu_hat + sigma_hat * qt(p, df = nu_hat)

png("fig_nvda_qq_student_t.png", width = 900, height = 600, res = 120)
plot(theoretical, sort(y),
     xlab = "Fitted Student-t Quantiles",
     ylab = "Sorted NVDA Returns",
     main = "Fitted Student-t Q-Q Plot (NVDA)")
abline(0, 1, col = "red", lwd = 2)
dev.off()

# 5. VaR & ES calculation
alpha_95 <- 0.05
alpha_99 <- 0.01

VaR_norm_95 <- -(mean(y) + sd(y) * qnorm(alpha_95))
ES_norm_95  <- -(mean(y) - sd(y) * dnorm(qnorm(alpha_95)) / alpha_95)
VaR_norm_99 <- -(mean(y) + sd(y) * qnorm(alpha_99))
ES_norm_99  <- -(mean(y) - sd(y) * dnorm(qnorm(alpha_99)) / alpha_99)

mu_draw    <- post[, "mu"]
sigma_draw <- post[, "sigma"]
nu_draw    <- post[, "nu"]

risk_draws <- function(alpha) {
  z <- qt(alpha, df = nu_draw)
  VaR_return <- mu_draw + sigma_draw * z
  ES_return  <- mu_draw - sigma_draw * (dt(z, df = nu_draw) / alpha) * ((nu_draw + z^2) / (nu_draw - 1))
  data.frame(VaR_loss = -VaR_return, ES_loss = -ES_return)
}

summarise_risk <- function(rd) {
  c(VaR_median = median(rd$VaR_loss), VaR_q025 = unname(quantile(rd$VaR_loss, 0.025)), VaR_q975 = unname(quantile(rd$VaR_loss, 0.975)),
    ES_median  = median(rd$ES_loss),  ES_q025  = unname(quantile(rd$ES_loss, 0.025)),  ES_q975  = unname(quantile(rd$ES_loss, 0.975)))
}

risk_t_95 <- summarise_risk(risk_draws(alpha_95))
risk_t_99 <- summarise_risk(risk_draws(alpha_99))

# Kupiec test
check_violations <- function(VaR_loss, alpha, model_name) {
  viol <- sum(y < -VaR_loss); n <- length(y); rate <- viol / n
  p_value <- binom.test(viol, n, p = alpha)$p.value
  data.frame(Model = model_name, Alpha = alpha, Violations = viol, N = n, Rate = rate, Kupiec_p = p_value)
}

viol_table <- rbind(
  check_violations(VaR_norm_95, 0.05, "Normal"), 
  check_violations(VaR_norm_99, 0.01, "Normal"),
  check_violations(risk_t_95["VaR_median"], 0.05, "Student-t"), 
  check_violations(risk_t_99["VaR_median"], 0.01, "Student-t")
)

risk_table <- data.frame(
  Model = c("Normal", "Normal", "Student-t", "Student-t"),
  Alpha = c(0.05, 0.01, 0.05, 0.01),
  VaR = c(VaR_norm_95, VaR_norm_99, risk_t_95["VaR_median"], risk_t_99["VaR_median"]),
  VaR_CI_Low = c(NA, NA, risk_t_95["VaR_q025"], risk_t_99["VaR_q025"]),
  VaR_CI_High = c(NA, NA, risk_t_95["VaR_q975"], risk_t_99["VaR_q975"]),
  ES = c(ES_norm_95, ES_norm_99, risk_t_95["ES_median"], risk_t_99["ES_median"])
)

# 6. Save results
write.csv(risk_table, "fig_nvda_risk_comparison.csv", row.names = FALSE)
write.csv(viol_table, "fig_nvda_var_violations.csv", row.names = FALSE)

cat("--- Risk comparison table ---\n")
print(risk_table, digits = 4)

cat("\n--- Kupiec violations table ---\n")
print(viol_table, digits = 4)

cat("\nNVDA_MCMC done! \n")