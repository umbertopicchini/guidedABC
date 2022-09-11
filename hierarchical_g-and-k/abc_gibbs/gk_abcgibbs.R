# This is an adaptation of the code https://github.com/GClarte/ABCG which supports the paper Clarté, G., Robert, C. P., Ryder, R. J., & Stoehr, J. (2021). Componentwise approximate Bayesian computation via Gibbs-like steps. Biometrika, 108(3), 591-607.

# We stripped away samplers that are not of ABC-Gibbs type and loaded data made of 20 g-and-k datasets each having length 1000.
# These are loaded below via X <- read.table('yobs.dat')

library(gk)
library(ggplot2)


source("functions_gk_simple.R")

hyper <- 5.7072  # this is the alpha parameter ("hyperparameter")
B <- 0.1915
g <- 0.6221
k <- 0.4377
As <- c(7.0242,
        5.3015,
        5.2623,
        7.0356,
        6.5410,
        6.3116,
        5.6005,
        5.7102,
        6.1643,
        6.6289,
        5.5098,
        5.7827,
        5.7161,
        6.5460,
        5.1648,
        6.2747,
        2.8875,
        5.4935,
        6.8803,
        6.1476)

X <- read.table('yobs.dat') # these data have been generated within Matlab using the true parameters above

#for (i in 1:ncol(X)) {
#  X[, i] <- rgk(nrow(X), As[i], B, g, k)
#}




#Statstar <- numeric()
#for (i in 1:ncol(X)) {
#  Statstar <- c(Statstar, quantile(X[, i], (0:8) / 8))
#}

# définition du modèle pour SMCmaison
#rdist <- function(par, Statstar) {
#  stat <- numeric()
#  for (i in 1:ncol(X)) {
#    stat <- c(stat, quantile(rgk(nrow(X), par[1 + i], B, g, k), (0:8) / 8))
#  }
#  return(sum(abs(stat - Statstar)))
#}

#rprior <- function() {
#  P <- numeric(ncol(X)+1)
#  P[1] <- runif(1, -10, 10)
#  P[2:ncol(X)+1] <- rnorm(ncol(X), P[1], 1)
#  return(P)
#}

#dprior <- function(P, ...) {
#  dunif(P[1], -10, 10) * prod(dnorm(P[2:ncol(X)+1], P[1], 1))
#}

#model <- rdist

#monprior <- list(density = dprior, simu = rprior)


outgib <- gibbs(c(100, 50, 50, 50, 50), X, 10000)
write.table(outgib[[1]], file = paste("ABCgibbs_out.dat", sep = ''), sep = "\t", row.names = F, col.names = F)



Dat <- data.frame(
  value = c(
    outgib[[1]][, 1],
    outgib[[1]][, 2],
    outgib[[1]][, 3],
    outgib[[1]][, 4],
    outgib[[1]][, 5]
  ),
  Method = rep(c(rep("ABC Gibbs", 1000)), 5),
  Parameter = c(
    rep("Hyperparameter", 3000),
    rep("mu1", 3000), rep("mu2", 3000), rep("mu3", 3000), rep("mu4", 3000)
  )
)

Dat$Parameter <- factor(Dat$Parameter,
                        levels = c("Hyperparameter", "mu1", "mu2", "mu3", "mu4"),
                        labels = c(
                          "Hyperparameter" = expression(alpha),
                          "mu1" = expression(mu[1]),
                          "mu2" = expression(mu[2]),
                          "mu3" = expression(mu[3]),
                          "mu4" = expression(mu[4])
                        )
)

theta_star <- data.frame(
  value = c(hyper, As[1:4]),
  Parameter = levels(Dat$Parameter)
)


f <- ggplot(data = Dat) + stat_density(aes(x = value, linetype = Method),
                                       geom = "line", position = "identity"
) +
  geom_vline(data = theta_star, aes(xintercept = value)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + coord_cartesian(
    xlim = c(8, 12), ylim = c(0, 10), expand = TRUE,
    default = FALSE, clip = "on"
  )
# geom_density(data=Dat,aes(x=value,color=Method))

f <- f + facet_wrap(Parameter ~ ., labeller = label_parsed, ncol = 3) + scale_colour_viridis_d()
f
