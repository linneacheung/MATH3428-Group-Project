# ============================================
# NVDA EDA
# MATH3428 Group Project - Group 1 - XIE YI QI
# Data range: 2024-01-01 to 2026-09-30
# ============================================

# install
if (!require("quantmod")) install.packages("quantmod")
if (!require("moments"))  install.packages("moments")
if (!require("ggplot2"))  install.packages("ggplot2")

library(quantmod)
library(moments)
library(ggplot2)

# NVDA
nvda <- getSymbols("NVDA", src = "yahoo",
                   from = "2024-01-01",
                   to   = "2026-10-01",
                   auto.assign = FALSE)

price <- Ad(nvda)          # Adjusted Close
price <- na.omit(price)

# daily log return y_t = ln(P_t / P_{t-1})
y <- na.omit(diff(log(price)))
y_vec <- as.numeric(y)

cat("number of observations:", length(y_vec), "\n")

# descriptive statistics
stats_table <- data.frame(
  Mean            = mean(y_vec),
  SD              = sd(y_vec),
  Skewness        = skewness(y_vec),
  Kurtosis        = kurtosis(y_vec),          # raw kurtosis, normal = 3
  Excess_Kurtosis = kurtosis(y_vec) - 3,
  Min             = min(y_vec),
  Max             = max(y_vec),
  N               = length(y_vec)
)
print(stats_table)

# histogram with fitted normal curve
p1 <- ggplot(data.frame(y = y_vec), aes(x = y)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = "grey70", colour = "white") +
  stat_function(fun = dnorm,
                args = list(mean = mean(y_vec), sd = sd(y_vec)),
                colour = "red", linewidth = 1) +
  labs(title = "NVDA Daily Log Returns vs Normal Fit",
       x = "Log return", y = "Density") +
  theme_minimal()

print(p1)

# Q-Q plot
png("nvda_qqplot.png", width = 800, height = 600)
qqnorm(y_vec, pch = 16, cex = 0.6,
       main = "NVDA Log Returns: Normal Q-Q Plot")
qqline(y_vec, col = "red", lwd = 2)
dev.off()

# save histogram
ggsave("nvda_hist.png", p1, width = 8, height = 5)


# Quick summary
cat("\n--- Summary ---\n")
cat("Kurtosis (normal = 3):", kurtosis(y_vec), "\n")
cat("Skewness:", skewness(y_vec), "\n")
cat("Daily SD:", sd(y_vec), "\n")
cat("Annualized SD (approx):", sd(y_vec) * sqrt(252), "\n")