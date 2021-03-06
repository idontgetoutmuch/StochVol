rm(list=ls())
set.seed(123)
library("stochvol")
data("exrates")
ret <- logret(exrates$USD, demean = TRUE)

png(file="diagrams/eurusd.png")
par(mfrow = c(2, 1))
plot(exrates$date, exrates$USD, type = "l",main = "Price of 1 EUR in USD", xlab="Year", ylab="Rate")
plot(exrates$date[-1], ret, type = "l", main = "De-meaned log returns", xlab="Year", ylab="Log Return")
dev.off()
