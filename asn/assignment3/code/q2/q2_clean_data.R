##
## Script for cleaning the AAPL stock data 
##
setwd("~/drive/concordia/2015-2016/1_fall2015/macf402/assignment3/")

mar2012 <- read.csv("./data/raw_mar2012.csv", stringsAsFactors = F)
apr2012 <- read.csv("./data/raw_apr2012.csv", stringsAsFactors = F)
may2012 <- read.csv("./data/raw_may2012.csv", stringsAsFactors = F)

mar2012$date <- as.Date("17-Mar-2012", "%d-%b-%Y")
apr2012$date <- as.Date("17-Apr-2012", "%d-%b-%Y")
may2012$date <- as.Date("17-May-2012", "%d-%b-%Y")

dat <- rbind(mar2012, apr2012, may2012)
dat <- dat[c("date","Last","Vol","Strike", "type")]

dat$Last <- as.numeric(dat$Last)
dat$Vol <- as.numeric(dat$Vol)
dat$Strike <- as.numeric(dat$Strike)

write.csv(dat, "./data/cleaned/cleaned_opt_data.csv", row.names = F)
write.csv(dat[dat$date == "2012-03-17",], "./data/cleaned/cleaned_mar2012.csv", row.names = F)
write.csv(dat[dat$date == "2012-04-17",], "./data/cleaned/cleaned_apr2012.csv", row.names = F)
write.csv(dat[dat$date == "2012-05-17",], "./data/cleaned/cleaned_may2012.csv", row.names = F)
