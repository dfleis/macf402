library(ggplot2)
setwd("~/Desktop/macf401_a5_test/")
dat_list <- list()
my_files <- list.files("./data/out")
idx <- 1
for (file in my_files) {
  dat_list[[idx]] <- read.csv(paste0("./data/out/",file), stringsAsFactors = F)
  idx <- idx + 1
}

ymin <- min(sapply(dat_list, FUN = function(x) min(x$diff)))
ymax <- max(sapply(dat_list, FUN = function(x) max(x$diff)))

ymax_last <- max(dat_list[[length(dat_list)]]$diff)
ymin_last <- min(dat_list[[length(dat_list)]]$diff)

dmax <- (ymax - ymax_last)/length(dat_list)
dmin <- (ymin - ymin_last)/length(dat_list)

my_plots <- list()
for (i in 1:length(dat_list)) {
  
  ymax_dx <- ymax - i * dmax
  ymin_dx <- ymin - i * dmin
  
  print(c(ymax_dx, ymin_dx, ymax_dx - ymin_dx))
  
  dat <- dat_list[[i]]
  N <- dat$N[1]
  plot_title <- paste("American Option Prices: N = ", N, "\nBinomial Model Price vs. Observed Price", sep = "")
  dat$Type <- "Call"
  dat$Type[dat$type == "P"] <- "Put"
  
  my_plots[[i]] <- ggplot(dat, aes(x = seq(1:length(dat$price)), y = diff, shape = Type)) + 
    labs(title = plot_title) +
    ylab("Binomial - Observed Price Difference") +
    scale_y_continuous(limits = c(ymin_dx, ymax_dx)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(1,2)) +  
    geom_hline(yintercept = 0) +
    facet_grid(. ~ Type, scales = "free") + 
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank(),
          legend.key = element_blank())
  
  plot_filename <- paste("./plots/am_price_diff_N", sprintf("%04d",N), ".png", sep = "")
  png(filename = plot_filename, width = 800, height = 500)
  print(my_plots[[i]])
  dev.off()
}
setwd("~/Desktop/macf401_a5_test/plots/")
system("convert *.png -delay 1 -loop 0 0_am_price_diff.gif")
setwd("~/Desktop/macf401_a5_test/")

x <- cbind(c(-1,-1),c(1,1))
