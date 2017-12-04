##
## Script for reading in & plotting implied vols from C++ code
##

setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

mar_dat <- read.csv("./data/output/q2_imp_vols_mar.csv", stringsAsFactors = F)
apr_dat <- read.csv("./data/output/q2_imp_vols_apr.csv", stringsAsFactors = F)
may_dat <- read.csv("./data/output/q2_imp_vols_may.csv", stringsAsFactors = F)

### Plot march smile
pdf("./plots/q2/aapl_mar_smile.pdf", height = 4.75, width = 4.75)
plot(mar_dat$imp_vol~mar_dat$moneyness, ylim = c(0.3,0.6),
     main = "AAPL Volatility Smile",
     ylab = "Implied Volatility", xlab = "Moneyness = S/K")
mtext(paste0(paste0("S = $502.12, r = 1.75%, ", mar_dat$date[1]), " Maturity"))
dev.off()

### Plot april smile
pdf("./plots/q2/aapl_apr_smile.pdf", height = 4.75, width = 4.75)
plot(apr_dat$imp_vol~apr_dat$moneyness, ylim = c(0.3,0.6),
     main = "AAPL Volatility Smile",
     ylab = "Implied Volatility", xlab = "Moneyness = S/K")
mtext(paste0(paste0("S = $502.12, r = 1.75%, ", apr_dat$date[1]), " Maturity"))
dev.off()

### Plot may smile
pdf("./plots/q2/aapl_may_smile.pdf", height = 4.75, width = 4.75)
plot(may_dat$imp_vol~may_dat$moneyness, ylim = c(0.3,0.6),
       main = "AAPL Volatility Smile",
       ylab = "Implied Volatility", xlab = "Moneyness = S/K")
mtext(paste0(paste0("S = $502.12, r = 1.75%, ", may_dat$date[1]), " Maturity"))
dev.off()

### Plot all three smiles in one figure
pdf("./plots/q2/aapl_smile.pdf", height = 6, width = 4.75)
par(mfrow = c(3,1), oma = c(0,0,2,0))
plot(mar_dat$imp_vol~mar_dat$moneyness, ylim = c(0.3,0.6),
     main = paste0(mar_dat$date[1], " Maturity"),
     ylab = "Implied Volatility", xlab = "Moneyness = S/K")

plot(apr_dat$imp_vol~apr_dat$moneyness, ylim = c(0.3,0.6),
     main = paste0(apr_dat$date[1], " Maturity"),
     ylab = "Implied Volatility", xlab = "Moneyness = S/K")

plot(may_dat$imp_vol~may_dat$moneyness, ylim = c(0.3,0.6),
     main = paste0(may_dat$date[1], " Maturity"),
     ylab = "Implied Volatility", xlab = "Moneyness = S/K")

mtext("AAPL Volatility Smile: S = $502.12, r = 1.75%" , outer = T, cex = 1.25)
dev.off()