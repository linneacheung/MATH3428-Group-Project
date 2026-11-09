# ==============================================================
# MATH3428 Group Project - Group 1 (21172181)
# VOO EDA
# Data range: 2024-01-01 to 2026-10-01
# ==============================================================

install.packages("quantmod")
install.packages("moments")
install.packages("coda")

library("coda")
library(quantmod)
library(moments)

getSymbols("VOO", from = "2024-01-01", to = "2026-10-01")
price <- Ad(VOO)
y <- as.numeric(diff(log(price)))
y <- y[!is.na(y)]


stats <- c(Mean = mean(y),
           SD = sd(y),
           Skewness = skewness(y),
           Kurtosis = kurtosis(y),
           exc_Kurtosis = kurtosis(y) -3,
           Number = length(y))
round(stats, 5)
#histogram
hist(y, breaks = 50, freq = FALSE,
     main = "VOO Daily Log Returns", xlab = "Log return")
curve(dnorm(x, mean(y), sd(y)), add = TRUE, lwd = 2)

#QQ plot of EDA
qqnorm(y, main = "Normal Q-Q Plot of VOO Returns")
qqline(y)

#use mcmc_sampler.R
setwd("C:/Users/vicky/OneDrive/文件/voo_group_project")

getwd() 
source("mcmc_sampler.R")
length(y)
mean(y)
#get 3 mcmc chain
set.seed(123)
c1 <- run_student_t_mcmc(y, n_iter = 20000, burn_in = 5000,
                         init_theta = c( 0.000, log(0.005), log(5)))
c2 <- run_student_t_mcmc(y, n_iter = 20000, burn_in = 5000,
                         init_theta = c( 0.001, log(0.020), log(15)))
c3 <- run_student_t_mcmc(y, n_iter = 20000, burn_in = 5000,
                         init_theta = c(-0.001, log(0.010), log(40)))
                         
#3 mcmc chains on VOO
str(c1)
str(c2)
str(c3)
# acceptance rates of all 3 chains
c(c1$acceptance_rate, c2$acceptance_rate, c3$acceptance_rate)

# build the chains
chains <- mcmc.list(mcmc(as.matrix(c1$post_samples)),
                    mcmc(as.matrix(c2$post_samples)),
                    mcmc(as.matrix(c3$post_samples)))

# trace plots
traceplot(chains)

# Gelman-Rubin R-hat
gelman.diag(chains)

#ESS
effectiveSize(chains)

#ACF
autocorr.plot(chains)
autocorr.diag(chains)

#posterior inference code
post <- do.call(rbind, lapply(list(c1, c2, c3),
                              function(ch) as.matrix(ch$post_samples)))

colMeans(post)
apply(post, 2, median)
apply(post, 2, sd)
apply(post, 2, quantile, probs = c(0.025, 0.975))

#student-t qq-plot

mu_hat    <- mean(post[, "mu"])
sigma_hat <- mean(post[, "sigma"])
nu_hat    <- mean(post[, "nu"])

p <- ppoints(length(y))
theoretical <- mu_hat + sigma_hat * qt(p, df = nu_hat)

plot(theoretical, sort(y),
     xlab = "Fitted Student-t quantiles",
     ylab = "Sorted VOO returns",
     main = "Fitted Student-t QQ plot (VOO)")
abline(0, 1, col = "red")
