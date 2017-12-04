setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment4/")

## Load data
dat_crude <- read.csv("./data/q3_crude_est.csv")
dat_anti <- read.csv("./data/q3_anti_est.csv")
dat_cont <- read.csv("./data/q3_cont_est.csv")

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
dat_anti$Estimator <- "Anti."
dat_cont$Estimator <- "Control"

## Combine data frames
dat <- rbind(dat_crude, dat_anti, dat_cont)

## Compute CI Width
dat$diff <- dat$CI_hi - dat$CI_lo
dat$eff <- 1/dat$var
dat$eff <- floor(dat$eff * 10^5)/10^5 # round
dat$time <- floor(dat$time * 10^4)/10^4 # round

## Clean
dat <- moveCol(dat, "diff", 4)
dat <- removeCol(dat, c("X"))
dat <- prependCol(dat, "Estimator")
dat <- prependCol(dat, "N")
names(dat) <- c("N Sim.", "Est.", "CI High", "Mean", "CI Low", "CI Width", "Var.", "Geom. Price", "Time", "Eff.")
dat <- format(dat, nsmall = 4) # control digits

## Write csv
write.csv(dat, './data/q3_clean_estimators.csv', quote = F, row.names = F)

