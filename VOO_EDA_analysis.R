# ==============================================================
# MATH3428 Group Project - Group 1 (21172181)
# VOO EDA
# Data range: 2024-01-01 to 2026-10-01
# ==============================================================

install.packages("quantmod")
install.packages("moments")

library(quantmod)
getSymbols("VOO", from = "2024-01-01", to = "2026-10-01")
price <- Ad(VOO)
y <- as.numeric(diff(log(price)))
y <- y[!is.na(y)]


library(moments)

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

#QQ plot
qqnorm(y, main = "Normal Q-Q Plot of VOO Returns")
qqline(y)
