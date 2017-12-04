##
## Script for plotting the implied volatilities output by the corresponding C++ program
## Most of this trash is written for aesthetic purposes
##

setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment1/")
implied_vols <- read.csv("./data/implied_vols.csv", stringsAsFactors = F)
times <- read.csv("./data/times.csv", stringsAsFactors = F)

## get volatility data only for the largest N
maxN_callvols <- implied_vols$call[implied_vols$N == max(implied_vols$N)]
maxN_putvols <- implied_vols$put[implied_vols$N == max(implied_vols$N)]

## begin plotting smile
pdf("./plots/smile.pdf", height = 4.75, width = 9.5)
par(mfrow = c(1,2), oma = c(0, 0, 2, 0)) # 
## OTM OPTIONS: plot volatility smile for the largest N computed
otm_callvols <- maxN_callvols
otm_callvols[which(unique(implied_vols$moneyness) <= 1)] <- NaN
otm_putvols <- maxN_putvols
otm_putvols[which(unique(implied_vols$moneyness) >= 1)] <- NaN
# endpoints for y axis, consider the highest/lowest volatilities & round them up/down
lo <- floor(min(otm_callvols, otm_putvols, na.rm = T)*100)/100
hi <- ceiling(max(otm_callvols, otm_putvols, na.rm = T)*100)/100 

plot(otm_callvols ~ unique(implied_vols$strike), 
     lwd = 2, type = 'b', ylim = c(lo, hi),
     xlab = "Strike", ylab = "Implied Volatility", main = "Volatilities: OTM Options",
     xaxt = "n", yaxt = "n") # xaxt = 'n', yaxt = 'n', suppress drawing axis values
axis(2, at = seq(lo, hi, 0.05), las = 1) # draw y axis ticks from lo to hi every 0.05
axis(1, at = unique(implied_vols$strike), las = 1) # draw y axis ticks from lo to hi every 0.05
lines(otm_putvols ~ unique(implied_vols$strike), lwd = 2, pch = 24, type = 'b', yaxt = "n")
legend("topright", legend = c("Call","Put"), pch = c(21,24), text.font = 2)

## plot ALL strikes regardless of moneyness
# endpoints for y axis, consider the highest and lowest volatilities
lo <- floor(min(maxN_callvols, maxN_putvols, na.rm = T)*100)/100
hi <- ceiling(max(maxN_callvols, maxN_putvols, na.rm = T)*100)/100

## ALL OPTIONS: smile for largest N
plot(maxN_callvols ~ unique(implied_vols$strike), 
     lwd = 2, type = 'b', ylim = c(lo, hi),
     xlab = "Strike", ylab = "Implied Volatility", main = "Volatilities: All Options",
     xaxt = "n", yaxt = "n") # xaxt = 'n', yaxt = 'n', suppress drawing axis values
axis(2, at = seq(lo, hi, 0.05), las = 1) # draw y axis ticks from lo to hi every 0.05
axis(1, at = unique(implied_vols$strike), las = 1) # draw y axis ticks from lo to hi every 0.05
lines(maxN_putvols ~ unique(implied_vols$strike), lwd = 2, pch = 24, type = 'b', yaxt = "n")
legend("topright", legend = c("Call","Put"), pch = c(21,24), text.font = 2)
mtext(paste0("Volatility Smile for N-steps = ", max(implied_vols$N)), outer = T, cex = 1.25)
dev.off() 

## calls
## plot convergence of implied volatility for increasing N
pdf("./plots/call_convergence.pdf", height = 6.5, width = 5.5)
# get unique strike values of observed call options
call_strikes <- unique(implied_vols$strike[!is.nan(implied_vols$call)])
par(mfrow = c(3,2), oma = c(0, 0, 2, 0)) # create a frame with 3x2 boxes

for (i in 1:length(call_strikes)) { # loop over all strike prices
  vol_subset <- implied_vols[implied_vols$strike == call_strikes[i],]
  lo <- floor(min(vol_subset$call)*100)/100 # round y lower axis endpoints
  hi <- ceiling(max(vol_subset$call)*100)/100 # round y upper axis endpoints
  
  plot(vol_subset$call ~ log10(unique(implied_vols$N)), 
       type = 'b', ylim = c(lo, hi),
       xlab = 'log10(N-steps)', ylab = 'Implied Volatility', 
       main = paste0("Strike = ",sprintf("%.2f",call_strikes[i])),
       yaxt = "n", lwd = 1.5) # yaxt = 'n' suppress drawing the y axis values
  axis(2, at = seq(lo, hi, by = 0.01), las = 1) # draw y axis ticks from lo to hi every 0.01

}
mtext("Implied Volatility Convergence: Call Options", outer = T, cex = 1)
dev.off()

## puts
## plot convergence of implied volatility for increasing N
pdf("./plots/put_convergence.pdf", height = 6.5, width = 5.5)
# get unique strike values of observed put options 
put_strikes <- unique(implied_vols$strike[!is.nan(implied_vols$put)])
par(mfrow = c(3,2), oma = c(0, 0, 2, 0)) # create a frame with 3x2 boxes

for (i in 1:length(put_strikes)) { # loop over all strike prices
  vol_subset <- implied_vols[implied_vols$strike == put_strikes[i],]
  lo <- floor(min(vol_subset$put)*100)/100 # round y axis lower endpoint
  hi <- ceiling(max(vol_subset$put)*100)/100 # round y axis upper endpoint
  
  plot(vol_subset$put ~ log10(unique(implied_vols$N)), 
       type = 'b', ylim = c(lo, hi),
       xlab = 'log10(N-steps)', ylab = 'Implied Volatility', 
       main = paste0("Strike = ",sprintf("%.2f",put_strikes[i])),
       yaxt = "n", lwd = 1.5) # yaxt = 'n' suppress drawing the y axis values
  axis(2, at = seq(lo, hi, by = 0.01), las = 1) # draw y axis ticks from lo to hi every 0.01
}
mtext("Implied Volatility Convergence: Put Options", outer = T, cex = 1)
dev.off()

## plot computation time as a fxn of N
pdf("./plots/times.pdf", height = 4.75, width = 9.25)
par(mfrow = c(1,2), oma = c(0, 0, 2, 0)) # create a frame with 1x2 boxes

plot(times$seconds ~ log10(times$N), type = 'b',
     ylab = "Seconds", xlab = "log10(N-steps)",
     main = "Seconds by log10 Steps", las = 1, lwd = 2)
plot(log10(times$seconds) ~ log10(times$N), type = 'b',
     ylab = "log10(Seconds)", xlab = "log10(N-steps)",
     main = "log10 Seconds by log10 Steps", las = 1, lwd = 2)
mtext("Implied Volatility Computation Time", outer = T, cex = 1.25)
dev.off()