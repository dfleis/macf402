##
## Script for plotting fitted curves/surfaces
##
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

#==== LOAD DATA ====#
mar_data <- read.csv("./data/output/q2_imp_vols_mar.csv", stringsAsFactors = F)
fit_mar_data <- read.csv("./data/output/q2_fitted_1d_mar_data.csv", stringsAsFactors = F)
apr_data <- read.csv("./data/output/q2_imp_vols_apr.csv", stringsAsFactors = F)
fit_apr_data <- read.csv("./data/output/q2_fitted_1d_apr_data.csv", stringsAsFactors = F)
may_data <- read.csv("./data/output/q2_imp_vols_may.csv", stringsAsFactors = F)
fit_may_data <- read.csv("./data/output/q2_fitted_1d_may_data.csv", stringsAsFactors = F)

fit_2d_data <- read.csv("./data/output/q2_fitted_2d_data.csv", stringsAsFactors = F)

#==== PLOT 1-D FIT ====#
pdf("./plots/q2/fitted_1d.pdf", height = 6, width = 4.75)
par(mfrow = c(3,1), oma = c(0, 0, 2, 0)) # set up 1x3 box for 1-D smoothing plot

#==== MARCH DATA ====#
plot(mar_data$imp_vol ~ mar_data$moneyness,
     xlab = "Moneyness = S/K", ylab = "Implied Volatility",
     main = paste0("Maturity = ", mar_data$date[1]))
lines(fit_mar_data$fitted_vol ~ fit_mar_data$moneyness, type = 'l', lwd = 2)

#==== APRIL DATA ====#
plot(apr_data$imp_vol ~ apr_data$moneyness,
     xlab = "Moneyness = S/K", ylab = "Implied Volatility",
     main = paste0("Maturity = ", apr_data$date[1]))
lines(fit_apr_data$fitted_vol ~ fit_apr_data$moneyness, type = 'l', lwd = 2)

#==== MAY DATA ====#
plot(may_data$imp_vol ~ may_data$moneyness,
     xlab = "Moneyness = S/K", ylab = "Implied Volatility",
     main = paste0("Maturity = ", may_data$date[1]))
lines(fit_may_data$fitted_vol ~ fit_may_data$moneyness, type = 'l', lwd = 2)

# global title
mtext("Fitted Smile: Options on AAPL", outer = T, cex = 1.25)
dev.off()


#==== PLOT 2-D FIT ====#
# append all 1-D observed data
obs_2d_data <- rbind(mar_data, apr_data, may_data)

moneyness <- unique(fit_2d_data$moneyness)
maturity <- unique(fit_2d_data$maturity)
implied_vol <- matrix(fit_2d_data$fitted_vol, ncol = length(moneyness), nrow = length(maturity), byrow = T)


closest_fit <- data.frame(matrix(nrow = nrow(obs_2d_data), ncol = 3))
colnames(closest_fit) <- c('moneyness', 'maturity', 'fitted_vol')
for (i in 1:nrow(closest_fit)) {
  moneyness_err <- abs(fit_2d_data$moneyness - obs_2d_data$moneyness[i])
  maturity_err <- abs(fit_2d_data$maturity - obs_2d_data$maturity[i])
  idx <- which( (moneyness_err == min(moneyness_err)) & (maturity_err == min(maturity_err)))
  closest_fit[i,] <- as.numeric(fit_2d_data[idx,])
}

pdf("./plots/q2/fitted_2d_new_grid.pdf", height = 6, width = 7)
par(mfrow = c(1,1)) # create 1x1 frame to place plot
# create 2-D plot
p <- persp(moneyness, maturity, implied_vol, phi = 15, theta = -30,
           col = "gray80", border = "gray30", expand = 0.5, ticktype = "detailed",
           xlab = "Moneyness = S/K", ylab = "Maturity", zlab = "Implied Volatility",
           cex.lab = 0.75, cex.axis = 0.75, zlim = c(min(obs_2d_data$imp_vol), max(obs_2d_data$imp_vol)),
           main = "Fitted Surface with Grid Min/Max Corresponding to Data Range\nh1 = 0.05, h2 = 0.05")
obs <- trans3d(obs_2d_data$moneyness, obs_2d_data$maturity, obs_2d_data$imp_vol, pmat = p)
pred <- trans3d(closest_fit$moneyness, closest_fit$maturity, closest_fit$fitted_vol, pmat = p)
points(obs, col = 'black', pch = 1)
segments(x0 = obs$x, y0 = obs$y, x1 = pred$x, y1 = pred$y)
dev.off()



























