library(tidyverse)
library(ggthemes)

# bart simpson density
# http://www.math.yorku.ca/~hkj/Teaching/4230/Coverage/density.R
dbart <- function(x) {
    d1 <- dnorm(x, 0, 1)
    d2 <- dnorm(x, 0 / 2 - 1, 1 / 10)
    d3 <- dnorm(x, 1 / 2 - 1, 1 / 10)
    d4 <- dnorm(x, 2 / 2 - 1, 1 / 10)
    d5 <- dnorm(x, 3 / 2 - 1, 1 / 10)
    d6 <- dnorm(x, 4 / 2 - 1, 1 / 10)
    return(d1 + (1 / 10) * (d2 + d3 + d4 + d5 + d6))
}

# randomly draw from the bart simpson density
rbart <- function(n) {
    ind          <- sample(1:6, n, replace = TRUE, prob = c(5, rep(1, 5)) / 10)
    mu           <- (ind - 2) / 2 - 1
    mu[ind == 1] <- 0
    sg           <- rep(0.1, n) + rep(0.9, n) * (ind == 1)
    return(rnorm(n, mu, sg))
}

# leave-one-out cross validation scores for various bandwidtbws,
# for a given sample of xs, using the gaussian kernel
gaussian_loocv <- function(xs, bws) {
    n      <- length(xs)  # number of xs sampled
    k      <- length(bws) # number of bandwidtbws trialed
    scores <- numeric(length = n)
    for (h in seq_len(k)) {
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                den1 <- dnorm((xs[i] - xs[j]) / bws[h], sd = sqrt(2))
                den2 <- dnorm((xs[i] - xs[j]) / bws[h])
                scores[h] <- scores[h] + den1 - 2 * den2
            }
        }
        scores[h] <- scores[h] / (bws[h] * n^2) + 2 * dnorm(0) / (n * bws[h])
    }
    return(scores)
}

# determine which bandwidth produced the minimum loocv score
optimal_h <- function(scores, bws) {
    return(bws[which(scores == min(scores))])
}

# given a sample, and bandwidths to trial, produce a kernel density estimate
produce_kde <- function(xs, bw) {
    return(
        sapply(xs, function(xi) sum(dnorm((xi - xs) / bw) / (length(xs) * bw)))
    )
}

# TODO: Can this be functionalized?
sample_sz <- 1000
sample_xs <- sort(rbart(sample_sz))
test_bws  <- seq(from = 0.01, to = 0.20, length = 1000)
loocv_scr <- gaussian_loocv(sample_xs, test_bws)

optimal_h <- max(which(loocv_scr == min(loocv_scr))) # max is just a tie-breaker
kde_fits  <- sapply(test_bws, function(h) produce_kde(sample_xs, h))
true_den  <- dbart(
    seq(from = min(sample_xs), to = max(sample_xs), length = sample_sz)
)
