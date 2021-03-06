% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_ss}
\alias{get_ss}
\title{Get sufficient statistics}
\usage{
get_ss(dat)
}
\arguments{
\item{dat}{Data set with items as columns and examinees as rows. Missing
responses should be coded as NA.}
}
\value{
Vector of sufficient statistics
}
\description{
Find sufficient statistics from a data set
}
\examples{

set.seed(2524)

diff_1 <- rnorm(10)
diff_2 <- rnorm(10)

N <- 500

th <- MASS::mvrnorm(N, mu = c(0, -1),
                    Sigma = matrix(c(1, .5 * 2, .5 * 2, 4), nrow = 2))

probs_1 <- 1 / (1 + exp(-outer(th[, 1], diff_1, "-")))
probs_2 <- 1 / (1 + exp(-outer(th[, 2], diff_2, "-")))

probs <- cbind(probs_1, probs_2)

dat <- apply(probs, 2, function(p) as.numeric(p > runif(N)))

get_ss(dat)

}
