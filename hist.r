library(latex2exp)
library(magick)
set.seed(353535)

# === BART SIMPSON DENSITY ===
# http://www.math.yorku.ca/~hkj/Teaching/4230/Coverage/density.R
dbart <- function(x) {
    d1 <- dnorm(x, 0, 1)
    d2 <- dnorm(x, 0 / 2 - 1, 1 / 10)
    d3 <- dnorm(x, 1 / 2 - 1, 1 / 10)
    d4 <- dnorm(x, 2 / 2 - 1, 1 / 10)
    d5 <- dnorm(x, 3 / 2 - 1, 1 / 10)
    d6 <- dnorm(x, 4 / 2 - 1, 1 / 10)
    f  <- 0.5 * d1 + (1 / 10) * (d2 + d3 + d4 + d5 + d6)
    return(f)
}


rbart <- function(n) {
    ind <- sample(1:6, n, replace = TRUE, prob = c(5, rep(1,5)) / 10)
    mean <- (ind - 2) / 2 - 1
    mean[ind == 1] <- 0
    sd <- rep(0.1, n) + rep(0.9, n) * (ind == 1)
    x <- rnorm(n, mean, sd)
    return(x)
}

