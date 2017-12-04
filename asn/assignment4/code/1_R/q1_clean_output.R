setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude <- read.csv("./data/q1_crude_est.csv")
dat_anti <- read.csv("./data/q1_anti_est.csv")

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
dat_crude$Estimator <- "Crude"
dat_anti$Estimator <- "Antithetic"

## Combine data frames
dat <- rbind(dat_crude, dat_anti)

## Compute CI Width & eff
dat$diff <- dat$CI_hi - dat$CI_lo
dat$eff <- 1/dat$var
dat$time <- floor(dat$time * 10^4)/10^4

## Clean
dat <- moveCol(dat, "diff", 4)
dat <- removeCol(dat, c("X"))
dat <- prependCol(dat, "Estimator")
dat <- prependCol(dat, "N")
names(dat) <- c("N Sim.", "Est.", "CI High", "Mean", "CI Low", "CI Width", "Var.", "Time", "Eff.")
dat <- format(dat, nsmall = 4) # control digits

## Write csv
write.csv(dat, './data/q1_clean_estimators.csv', quote = F, row.names = F)

