setwd("~/drive/concordia/2015-2016/2_winter2016/macf401/assignments/a5/")
dat_vols <- read.csv("./data/implied_vols_googl_NMAX.csv", stringsAsFactors = F)
dat_vols <- dat_vols[dat_vols$N == max(dat_vols$N),]

p <- 0.5 # global variable
q <- 0.5 # global variable

u_fxn <- function(R, sigma, dt) {
  return ( exp((R - 0.5 * sigma^2) * dt + sigma * sqrt(dt)) )
}
d_fxn <- function(R, sigma, dt) {
  return ( exp((R - 0.5 * sigma^2) * dt - sigma * sqrt(dt)) )
}
asset <- function(S0, n, h, R, sigma, dt) {
  u <- u_fxn(R, sigma, dt)
  d <- d_fxn(R, sigma, dt)
  S0 * u^h * d^(n - h)
}
payoff <- function(S, K, type) {
  if (!(type %in% c('c','p'))) {
    warning("Warning in price_option: Invalid type.")
    return (0)
  } else if (type == 'c') {
    return (max(S - K, 0))
  } else
    return (max(K - S, 0))
}
price_option_am_r <- function(S0, K, n, h, N, R, sigma, dt, type) {
  if (h > n) {
    warning("Warning in price_option: h > n.")
  }
  if (!(type %in% c('c','p'))) {
    warning("Error in price_option: Invalid type.")
    return (0)
  }
  S <- asset(S0, n, h, R, sigma, dt)
  G_n <- payoff(S, K, type)
  if (n == N) {
    return ( G_n )
  } else {
    Vd <- price_option_am_r(S0, K, n + 1, h, N, R, sigma, dt, type)
    Vu <- price_option_am_r(S0, K, n + 1, h + 1, N, R, sigma, dt, type)
    return ( max( G_n, (1/(1 + R)^dt * (p * Vu + q * Vd)) ) )
  }
}
price_option_am_i <- function(S0, K, n, h, N, R, sigma, dt, type) {
  if (h > n) {
    warning("Warning in price_option: h > n.")
    return (0)
  }
  if (!(type %in% c('c','p'))) {
    warning("Error in price_option: Invalid type.")
    return (0)
  }
  
  prices <- double(N + 1)
  for (l in 0:N) { # final step will have N + 1 nodes
    # final step at node (N, i), Nth step, l heads
    S <- asset(S0, N, l, R, sigma, dt) 
    prices[l + 1] <- payoff(S, K, type)
  }
  for (k in (N - 1):n) { # move backwards through the tree
    for (l in 1:(k + 1)) {
      S <- asset(S0, k, l - 1, R, sigma, dt) # price asset at the k-th step given l - 1 heads
      prices[l] <- max( payoff(S, K, type), 1/(1 + R)^dt * (p * prices[l + 1] + q * prices[l]) )
    }
  }
  return (prices[1])
}


N <- 1000
R <- 0.005
prices <- NA
for (i in 1:nrow(dat_vols)) {
  print(i)
  type <- as.character(tolower(dat_vols$type[i]))
  dt <- dat_vols$tau[i]/N
  prices[i] <- price_option_am_i(S0 = dat_vols$spots[i], K = dat_vols$strike[i], 
                                 n = 0, h = 0, N = N, R = R, sigma = dat_vols$vol[i], dt = dt, type = type)
}

dat_vols$prices <- prices
dat_vols$diff <- dat_vols$prices - dat_vols$ask
dat_vols$Type <- "Call"
dat_vols$Type[dat_vols$type == "P"] <- "Put"

ggplot(dat_vols, aes(x = seq(1:length(dat_vols$prices)), y = diff, shape = Type)) + 
  labs(title = "American Option Prices\nBinomial Model Price vs. Ask Price") + 
  ylab("Binomial - Ask Price Difference") +
  geom_point(size = 5) + 
  geom_hline(yintercept = 0) +
  facet_grid(. ~ Type, scales = "free") + 
  theme_bw(15) + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) 


dat_vols







