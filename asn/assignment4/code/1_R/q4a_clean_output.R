setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude_n7 <- read.csv("./data/q4_heston_crude_est_-7.csv")
dat_crude_0 <- read.csv("./data/q4_heston_crude_est_0.csv")
dat_crude_7 <- read.csv("./data/q4_heston_crude_est_7.csv")

## Useful functions
removeCol <- function(df, colnames) {
  df[,!(names(df) %in% colnames)]
}
prependCol <- function(df, colname) {
  df[c(colname, names(df)[!(names(df) %in% colname)])]
}
moveCol <- function(df, colname, idx) {
  ## should check if idx is beyond min/max of df cols
  df <- prependCol(df, colname)
  for (i in 2:idx) {
    df <- prependCol(df, colnames(df)[i])
  }
  df
}

## Add ID variable
dat_crude_n7$rho <- -0.7
dat_crude_0$rho <- 0
dat_crude_7$rho <- 0.7

## Combine data frames
dat <- rbind(dat_crude_n7, dat_crude_0, dat_crude_7)

## Compute CI Width
dat$diff <- dat$CI_hi - dat$CI_lo
dat$eff <- 1/dat$var

## Clean
dat <- moveCol(dat, "diff", 4)
dat <- removeCol(dat, c("X"))
dat <- prependCol(dat, "rho")
dat <- prependCol(dat, "N")
names(dat) <- c("N Sim.", "rho", "CI High", "Mean", "CI Low", "CI Width", "Var.", "Time", "Eff.")
dat <- format(dat, nsmall = 4) # control digits

## Write csv
write.csv(dat, './data/q4a_clean_heston_crude_est.csv', quote = F, row.names = F)

