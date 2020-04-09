library(latex2exp)
# library(magick)
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


# === DENSITY ESIMATION ===
# Gaussian kernel
k <- function(x) return(dnorm(x))


# Estimate the density at point x
fhat <- function(x, xs, h) {
    n <- length(xs)
    return(sum(k(((x - xs) / h))) / (n * h))
}


# Calculate the leave-one-out cross-validation score for various bandwidths
loocv <- function(xs, hs) {
    n <- length(xs)
    scores <- vector()
    for (h in hs) {
        score <- 0
        for (xi in xs) {
            for (xj in xs) {
                score <- score + dnorm((xi-xj)/h, 0, sqrt(2)) - 2 * k((xi-xj)/h)
            }
        }
        score <- score / (h * n^2)
        score <- score + 2 * k(0) / (n * h)
        scores <- c(scores, score)
    }
    return(scores)
}


# Determine the bandwidth which corresponds to minimum loocv score
optimal_h <- function(loocv_scores, hs) {
    min_cv <- min(loocv_scores)
    indexes <- which(loocv_scores %in% min_cv)
    index <- min(indexes)
    return(hs[index])
}


# Produce estimates of the density at all values of x
kde <- function(xs, h) {
    vals <- vector()
    for (x in xs) {
        vals <- c(vals, fhat(x, xs, h))
    }
    return(vals)
}


# === PLOT OPTIMAL DENSITY ESTIMATES ===

gamma_data <- rgamma(250, 3, 1)
gamma_data <- sort(gamma_data)
gamma_x_axis <- seq(from = min(gamma_data), to = max(gamma_data), length = 250)
gamma_test_bandwidths <- seq(from = 0.10, to = 0.90, length = 250)
gamma_loocv_scores <- loocv(gamma_data, gamma_test_bandwidths)
gamma_optimal_bandwidth <- optimal_h(gamma_loocv_scores, gamma_test_bandwidths)
gamma_density_estimate <- kde(gamma_data, gamma_optimal_bandwidth)

bart_data <- rbart(250)
bart_data <- sort(bart_data)
bart_x_axis <- seq(from = -3, to = 3, length = 250)
bart_test_bandwidths <- seq(from = 0.01, to = 0.20, length = 250)
bart_loocv_scores <- loocv(bart_data, bart_test_bandwidths)
bart_optimal_bandwidth <- optimal_h(bart_loocv_scores, bart_test_bandwidths)
bart_density_estimate <- kde(bart_data, bart_optimal_bandwidth)

png("./png/gamma_loocv.png")
plot(gamma_test_bandwidths, gamma_loocv_scores, pch = "•", cex = 1, 
     col = "orange", xlab = TeX("$h$"), ylab = TeX("LOOCV Score"),
     main = TeX("LOOCV Score as a Function of $h$"))
abline(v = gamma_optimal_bandwidth, lwd = 3, col = "red")
dev.off()

png("./png/gamma_kde.png")
plot(gamma_data, gamma_density_estimate, type = "l", col = "darkgreen", lwd = 3,
     xlab = "x", ylab = TeX("$\\widehat{f}(x)"),
     main = paste("Density Estimate: h =", signif(gamma_optimal_bandwidth)),
     ylim = c(0, 0.3))
lines(gamma_x_axis, dgamma(gamma_x_axis, 3, 1), col = "blue", lty = 3, lwd = 2)
dev.off()

png("./png/bart_loocv.png")
plot(bart_test_bandwidths, bart_loocv_scores, pch = "•", cex = 1, 
     col = "orange", xlab = TeX("$h$"), ylab = TeX("LOOCV Score"),
     main = TeX("LOOCV Score as a Function of $h$"))
abline(v = bart_optimal_bandwidth, lwd = 3, col = "red")
dev.off()

png("./png/bart_kde.png")
plot(bart_data, bart_density_estimate, type = "l", col = "darkgreen", lwd = 3,
     xlab = "x", ylab = TeX("$\\widehat{f}(x)"),
     main = paste("Density Estimate: h =", signif(bart_optimal_bandwidth)),
     ylim = c(0, 0.7))
lines(bart_x_axis, dbart(bart_x_axis), col = "blue", lty = 3, lwd = 2)
dev.off()


# === MAKE GIF ===
frame_num <- 0
gamma_gif <- function(delay) {
    for (h in gamma_test_bandwidths) {
        if (h == gamma_optimal_bandwidth) {
            for (i in 1:100) {
                if (frame_num < 10) {
                    png(paste("./png/gamma/kde_000", frame_num, ".png", 
                              sep = ""))
                } else if (frame_num < 100) {
                    png(paste("./png/gamma/kde_00", frame_num, ".png", 
                              sep = ""))
                } else {
                    png(paste("./png/gamma/kde_0", frame_num, ".png", 
                              sep = ""))
                }
                plot(gamma_data, kde(gamma_data, h), type = "l", 
                     col = "darkgreen", lwd = 3, xlab = "x", 
                     ylab = TeX("$\\widehat{f}(x)"),
                     main = paste("Density Estimate: h =", signif(h)),
                     ylim = c(0, 0.3))
                lines(gamma_x_axis, dgamma(gamma_x_axis, 3, 1), col = "blue", 
                      lty = 3, lwd = 2)
                dev.off()
                frame_num <- frame_num + 1
            }
        } else {
            if (frame_num < 10) {
                png(paste("./png/gamma/kde_000", frame_num, ".png", sep = ""))
            } else if (frame_num < 100) {
                png(paste("./png/gamma/kde_00", frame_num, ".png", sep = ""))
            } else {
                png(paste("./png/gamma/kde_0", frame_num, ".png", sep = ""))
            }
            plot(gamma_data, kde(gamma_data, h), type = "l", col = "red", 
                 lwd = 3, xlab = "x", ylab = TeX("$\\widehat{f}(x)"),
                 main = paste("Density Estimate: h =", signif(h)),
                 ylim = c(0, 0.3))
            lines(gamma_x_axis, dgamma(gamma_x_axis, 3, 1), col = "blue", 
                  lty = 3, lwd = 2)
            dev.off()
            frame_num <- frame_num + 1
        }
    }
    system(paste("convert -delay ", delay, 
                 " ./png/gamma/*.png -loop 0 ./gamma.gif", sep = ""))
}


bart_gif <- function(delay) {
    for (h in bart_test_bandwidths) {
        if (h == bart_optimal_bandwidth) {
            for (i in 1:100) {
                if (frame_num < 10) {
                    png(paste("./png/bart/kde_000", frame_num, ".png", 
                              sep = ""))
                } else if (frame_num < 100) {
                    png(paste("./png/bart/kde_00", frame_num, ".png", sep = ""))
                } else {
                    png(paste("./png/bart/kde_0", frame_num, ".png", sep = ""))
                }
                plot(bart_data, kde(bart_data, h), type = "l", 
                     col = "darkgreen", lwd = 3, xlab = "x", 
                     ylab = TeX("$\\widehat{f}(x)"),
                     main = paste("Density Estimate: h =", signif(h)),
                     ylim = c(0, 0.7))
                lines(bart_x_axis, dbart(bart_x_axis), col = "blue", lty = 3, 
                      lwd = 2)
                dev.off()
                frame_num <- frame_num + 1
            }
        } else {
            if (frame_num < 10) {
                png(paste("./png/bart/kde_000", frame_num, ".png", sep = ""))
            } else if (frame_num < 100) {
                png(paste("./png/bart/kde_00", frame_num, ".png", sep = ""))
            } else {
                png(paste("./png/bart/kde_0", frame_num, ".png", sep = ""))
            }
            plot(bart_data, kde(bart_data, h), type = "l", col = "red", 
                 lwd = 3, xlab = "x", ylab = TeX("$\\widehat{f}(x)"),
                 main = paste("Density Estimate: h =", signif(h)),
                 ylim = c(0, 0.7))
            lines(bart_x_axis, dbart(bart_x_axis), col = "blue", lty = 3, 
                  lwd = 2)
            dev.off()
            frame_num <- frame_num + 1
        }
    }
    system(paste("convert -delay ", delay, 
                 " ./png/bart/*.png -loop 0 ./bart.gif", sep = ""))
}


gamma_gif(4)
bart_gif(4)
