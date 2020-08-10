#' Scale Alignment Wrapper Function for 'TAM' Objects
#'
#'
#' Apply scale alignment methods to models previously fit with tam.mml in the
#' 'TAM' package.
#'
#' @param mod Fitted model of class tam.mml. Importantly, mod$irtmodel must be
#' either "1PL", "PCM", or "PCM2"
#' @param method Either "DDA1", "DDA2", "LRA", or "best", see details
#' @param refdim Which is the reference dimension (unchanged during alignment)
#'
#' @return Aligned tam.mml object with the following added list items:
#' \item{method}{Alignment method: "DDA1", "DDA2", or "LRA"}
#' \item{rhat}{Vector of estimated scaling parameters r, see details}
#' \item{shat}{Vector of estimated shift parameters s, see details}
#' \item{cor_before}{Kendall's rank-order correlation between sufficient
#' statistics and Thurstone thresholds before alignment}
#' \item{cor_after}{Kendall's rank-order correlation between sufficient
#' statistics and Thurstone thresholds after alignment}
#'
#' @details Scales can be said to be aligned if the item sufficient statistics
#' imply the same item parameter estimates, regardless of dimension. Scale
#' alignment is currently defined only for Rasch family models with
#' between-items multidimensionality (i.e., each scored item belongs to exactly
#' one dimension).
#'
#' MODEL PARAMETERIZATIONS
#'
#' The partial credit model is a general Rasch family model for polytomous item
#' responses. Within 'TAM', the partial credit model can be parameterized in two
#' ways. If a 'TAM' model is fit with the option irtmodel = "PCM", then the
#' following model is specified for an item with \eqn{m + 1} response categories:
#'
#' \deqn{\log(\frac{P(x | \theta_d)}{P(x-1 | \theta_d)}) = \alpha_d \theta_d - \xi_{i(d)x}}{log(P(x | \theta_d)/(P(x - 1 | \theta_d))) = \alpha_d \theta_d - \xi_i(d)x}
#' for response category \eqn{x = 1,...,m}, and
#'
#' \deqn{P(x = 0 | \theta_d) = \frac{1}{\sum_{j=0}^{m}\exp \sum_{k=0}^j (\alpha_d \theta_d -
#' \xi_{i(d)k})}}{P(x = 0 | \theta_d) = 1 / (\sum_{j=0}^{m}\exp \sum_{k=0}^j (\alpha_d \theta_d -
#'  \xi_i(d)k))}
#' for response category \eqn{x=0}. \eqn{\alpha_d} is a dimension
#' steepness parameter, typically fixed to 1, \eqn{\theta_d} is a latent
#' variable on dimension \eqn{d}, and \eqn{\xi_i(d)x} is a step parameter for item
#' step \eqn{x} on item \eqn{i} belonging to dimension \eqn{d}.
#'
#' If instead a TAM model is fit with the option irtmodel = "PCM2", the model is
#' specified as
#'
#' \deqn{\log(\frac{P(x | \theta_d)}{P(x-1 | \theta_d)}) = \alpha_d \theta_d - \delta_{i(d)} +
#' \tau_{i(d)x}.}{log(P(x | \theta_d)/(P(x - 1 | \theta_d))) = \alpha_d \theta_d - \delta_i(d) +
#' \tau_i(d)x.}
#'
#' MODEL TRANSFORMATIONS
#'
#' Under Rasch family models, the latent trait metric can be linearly
#' transformed. For each dimension \eqn{d} the parameters on the transformed
#' metric (denoted by the \eqn{\sim}{'} symbol) are found through the transformation
#' parameters \eqn{r_d} and \eqn{s_d} as described by the following equations:
#'
#' \deqn{\tilde{\theta}_d = r_{d} \theta_{d} + s_{d}}{\theta'_d = r_{d}
#' \theta_{d} + s_{d}}
#'
#' \deqn{\tilde{\alpha}_d = \alpha_{d} / r_{d}}{\alpha'_d = \alpha_d / r_d}
#'
#' \deqn{\tilde{\xi}_{i(d)x} = \xi_{i(d)x} + \alpha_d s_d / r_d}{\xi'_i(d)x = \xi_i(d)x + \alpha_d * s_d / r_d}
#'
#' \deqn{\tilde{\delta}_{i(d)} = \delta_{i(d)} + \alpha_d s_d / r_d}{\delta'_i(d) = \delta_i(d) + \alpha_d * s_d / r_d}
#'
#' \deqn{\tilde{\tau}_{i(d)x} = \tau_{i(d)x}}{\tau'_i(d)x = \tau_{i(d)x}}
#'
#'
#' SUFFICIENT STATISTICS
#'
#' Under Rasch family models, the item sufficient statistics are the number of
#' examinees that score in response category \eqn{x} or higher,
#' \eqn{x = 1,...,m}. For the purpose of scale alignment, we consider sufficient
#' statistics to be the proportion of examinees that score in response category
#' \eqn{x} or higher. This definition allows for scale alignment in the
#' presence of missing data.
#'
#' THURSTONE THRESHOLDS
#'
#' Scales are aligned if the same sufficient statistics imply the
#' same item parameters, regardless of dimension. The success of scale alignment
#' is difficult to assess because the item sufficient statistics typically
#' differ across items and dimensions. Under the Rasch model for binary item
#' responses, the success of scale alignment can be assessed by looking at the
#' rank-order correlation (e.g., Kendall's tau) between item sufficient
#' statistics and item parameter estimates.
#'
#' However, under the partial credit model, item sufficient statistics need not
#' be monotonically related to estimated item parameters. Under this model,
#' we can assess the quality of scale alignment by taking the rank-order
#' correlation between item sufficient statistics and Thurstone thresholds.
#' Thurstone thresholds are defined as the \eqn{\theta} value at which the
#' probability of responding in category \eqn{x} or higher equals .5. Thurstone
#' thresholds, in most cases, will be monotonically related to item sufficient
#' statistics (within dimensions). Note that the item difficulty estimates
#' under the Rasch model for binary items are also Thurstone thresholds.
#'
#' ALIGNMENT METHODS
#'
#' Two types of scale alignment methods have been developed.
#'
#' The first class of methods, historically called delta-dimensional alignment
#' (DDA), requires fitting both a multidimenisonal model and a model in which
#' all items belong to a single dimension. With these two sets of parameter
#' estimates, the transformation parameters \eqn{r_d} and \eqn{s_d} are then
#' found so that, for each dimension, the means and standard deviation of
#' parameters from the transformed multidimensional models equal the means and
#' standard deviations of parameters from the unidimensional model. Under the
#' ordinary Rasch model, the estimated item difficulties can be used for
#' transformation (which is done if either method "DDA1" or "DDA2" is selected).
#' Under the partial credit model, either the \eqn{\delta} parameters or the
#' Thurstone thresholds from the two models may be used within the DDA (note
#' that DDA using item \eqn{\xi} parameters tends to be unsuccessful). Method
#' "DDA1" uses the item \eqn{\delta} parameters, and method "DDA2" uses the
#' Thurstone thresholds. If all items are binary, "DDA1" and "DDA2" are identical.
#'
#' The second class of methods, called logistic regression alignment (LRA),
#' requires fitting a logistic regression between item sufficient statistics
#' and Thurstone thresholds for each dimension. The fitted logistic regression
#' coefficients can then be used to estimate \eqn{r_d} and \eqn{s_d} so that the
#' same logistic regression curve expresses the relationship between sufficient
#' statistics and Thurstone thresholds for all dimensions.
#'
#' For either the DDA or LRA method, a reference dimension (by default, the
#' first dimension) is specified such that \eqn{r_d = 1} and \eqn{s_d = 0} for
#' the reference dimension.
#'
#' @references
#' Feuerstahler, L. M., & Wilson, M. (2019). Scale alignment in between-item multidimensional Rasch models. Journal of Educational Measurement, 56(2), 280--301. <doi: 10.1111/jedm.12209>
#'
#' Feuerstahler, L. M., & Wilson, M. (under review). Scale alignment in the between-items multidimensional partial credit model.
#'
#' @examples
#'
#' ## Example 1: binary item response data
#'
#' ## generate data for a 2-dimensional model with 10 items on each dimension
#'
#' if(require(TAM)){
#'
#' set.seed(2524)
#'
#' diff_1 <- rnorm(10)
#' diff_2 <- rnorm(10)
#'
#' N <- 500
#'
#' th <- MASS::mvrnorm(N, mu = c(0, -1),
#'                     Sigma = matrix(c(1, .5 * 2, .5 * 2, 4), nrow = 2))
#'
#' probs_1 <- 1 / (1 + exp(-outer(th[, 1], diff_1, "-")))
#' probs_2 <- 1 / (1 + exp(-outer(th[, 2], diff_2, "-")))
#'
#' probs <- cbind(probs_1, probs_2)
#'
#' dat <- apply(probs, 2, function(p) as.numeric(p > runif(N)))
#'
#' Q <- cbind(c(rep(1, 10), rep(0, 10)),
#'             c(rep(0, 10), rep(1, 10)))
#'
#' # fit the model
#'
#' mod <- TAM::tam.mml(resp = dat, irtmodel = "1PL", Q = Q)
#'
#' # align the model
#'
#' mod_aligned <- align(mod)
#'
#' ## check alignment success
#'
#' mod_aligned$cor_before
#' mod_aligned$cor_after
#'
#' ## view "best" alignment method
#' mod_aligned$method
#'
#' ## view alignment parameters
#'
#' mod_aligned$rhat
#' mod_aligned$shat
#'
#' }
#'
#' \donttest{
#'
#' ## Example 2: Partial Credit Model
#'
#' # generate 3-category data for a 2-dimensional model with 5 items on each dimension
#'
#' set.seed(8491)
#'
#' N <- 500
#'
#' th <- MASS::mvrnorm(N, mu = c(0, 0),
#'                     Sigma = matrix(c(1, .5 * 2, .5 * 2, 4), nrow = 2))
#'
#' xi_1 <- rnorm(5)
#' xi_1 <- cbind(xi_1, xi_1 + rnorm(5, mean = 1, sd = .5))
#'
#' xi_2 <- rnorm(5)
#' xi_2 <- cbind(xi_2, xi_2 + rnorm(5, mean = 1, sd = .5))
#'
#' dat1 <- catR::genPattern(th[, 1], it = xi_1, model = "PCM")
#' dat2 <- catR::genPattern(th[, 2], it = xi_2, model = "PCM")
#'
#' dat <- cbind(dat1, dat2)
#'
#' Q <- cbind(c(rep(1, 5), rep(0, 5)),
#'            c(rep(0, 5), rep(1, 5)))
#'
#' ## fit the model using both parameterizations
#'
#' mod1 <- TAM::tam.mml(resp = dat, irtmodel = "PCM", Q = Q)
#' mod2 <- TAM::tam.mml(resp = dat, irtmodel = "PCM2", Q = Q)
#'
#' ## align the models
#' mod1_aligned <- align(mod1)
#' mod2_aligned <- align(mod2)
#'
#' ## check alignment success
#'
#' mod1_aligned$cor_before
#' mod1_aligned$cor_after
#'
#' mod2_aligned$cor_before
#' mod2_aligned$cor_after
#'
#' ## view "best" alignment method
#' mod1_aligned$method
#' mod2_aligned$method
#'
#' ## view alignment parameters
#'
#' mod1_aligned$rhat
#' mod1_aligned$shat
#' mod2_aligned$rhat
#' mod2_aligned$shat
#'
#' }
#'
#' @export


align <- function(mod, method = "best", refdim = 1){

  # check for between-items multidimensionality
  if (!all(rowSums(mod$Q) == 1) | mod$ndim == 1)
    stop("Model must be between-items multidimensional")

  # make indicators

  # which dimension corresponds to each item
  dim_ind_i <- as.numeric(mod$Q %*% 1:ncol(mod$Q))

  # which item corresponds to each parameter
  item_ind <- as.numeric(apply(mod$A, 3, function(p)
    which(rowSums(p, na.rm = TRUE) != 0)))

  # find alphas for the original model
  alpha <- apply(mod$B[, 2, ], 2, max)


  if (method %in% c("DDA1", "DDA2", "best")){
    mod_uni <- TAM::tam.mml(resp = mod$resp,
                            irtmodel = mod$irtmodel, ndim = 1,
                            verbose = FALSE)

    thresh_u <- get_thresh(pars = mod_uni$xsi$xsi,
                           itemtype = mod_uni$irtmodel,
                           item_ind = item_ind,
                           alpha = 1)
  }

  ncats <- apply(mod$resp, 2, function(x) length(unique(x[!is.na(x)])))
  nobs <- apply(mod$resp, 2, function(x) sum(!is.na(x)))

  # find thresholds for the original model
  thresh_m <- get_thresh(pars = mod$xsi$xsi, itemtype = mod$irtmodel,
                         item_ind = item_ind, alpha = alpha[dim_ind_i])

  # get sufficient statistics
  ss <- get_ss(dat = mod$resp)

  cor_before <- -stats::cor(thresh_m, ss, method = "kendall")

  if (method %in% c("DDA1", "best")){

    ## DDA on item deltas
    dda1_out <- dda1(multi_pars = mod$xsi$xsi, uni_pars = mod_uni$xsi$xsi,
                     itemtype = mod$irtmodel, item_ind = item_ind,
                     dim_ind_i = dim_ind_i, refdim = refdim, alpha = alpha)
    dda1_out$cor_after <- -stats::cor(dda1_out$thresh, ss, method = "kendall")

  }

  if (method %in% c("DDA2", "best")){
    dda2_out <- dda2(multi_pars = mod$xsi$xsi, uni_pars = mod_uni$xsi$xsi,
                     itemtype = mod$irtmodel, item_ind = item_ind,
                     dim_ind_i = dim_ind_i, refdim = refdim, alpha = alpha)
    dda2_out$cor_after <- -stats::cor(dda2_out$thresh, ss, method = "kendall")
  }

  if (method %in% c("LRA", "best")){

    lra_out <- lra(multi_pars = mod$xsi$xsi, itemtype = mod$irtmodel, ss = ss,
                   nobs = nobs, ncats = ncats, thresh_m = thresh_m,
                   item_ind = item_ind, dim_ind_i = dim_ind_i, refdim = refdim,
                   alpha = alpha)
    lra_out$cor_after <- -stats::cor(lra_out$thresh, ss, method = "kendall")
  }

  if (method == "DDA1") res <- dda1_out
  if (method == "DDA2") res <- dda2_out
  if (method == "LRA") res <- lra_out
  if (method == "best"){
    best <- which.max(c(dda1_out$cor_after, dda2_out$cor_after,
                        lra_out$cor_after))
    res <- switch(best, dda1_out, dda2_out, lra_out)
    method <- switch(best, "DDA1", "DDA2", "LRA")
  }




  ## refit the model
  xsi.fixed <- mod$xsi.fixed.estimated
  xsi.fixed[, 2] <- res$new_pars
  B.fixed <- mod$B.fixed.estimated
  for (d in 1:mod$ndim){
    B.fixed[B.fixed[, 3] == d, 4] <- B.fixed[B.fixed[, 3] == d, 4] *
      res$alphatilde[d]
  }

  variance.matrix <- diag(res$rhat) %*% mod$variance %*% diag(res$rhat)
  variance.fixed <- as.matrix(expand.grid(1:mod$ndim, 1:mod$ndim, 0))
  variance.fixed <- variance.fixed[variance.fixed[, 2] >= variance.fixed[, 1], ]
  for (i in 1:nrow(variance.fixed)){
    variance.fixed[i, 3] <- variance.matrix[variance.fixed[i, 1],
                                            variance.fixed[i, 2]]
  }


  beta.fixed <- mod$beta.fixed
  beta.fixed[, 3] <- res$rhat * mod$beta.fixed[, 3] + res$shat


  mod_aligned <- TAM::tam.mml.2pl(resp = mod$resp,
                                  group = mod$group,
                                  irtmodel = "GPCM",
                                  ndim = mod$ndim, pid = mod$pid,
                                  xsi.fixed = xsi.fixed,
                                  xsi.inits = xsi.fixed,
                                  beta.fixed = beta.fixed,
                                  beta.inits = beta.fixed,
                                  variance.fixed = variance.fixed,
                                  variance.inits = variance.matrix,
                                  A = mod$A,
                                  B.fixed = B.fixed,
                                  Q = mod$Q,
                                  pweights = mod$pweights,
                                  control = mod$control,
                                  verbose = FALSE)

  mod_aligned$method <- method
  mod_aligned$rhat <- res$rhat
  mod_aligned$shat <- res$shat
  mod_aligned$cor_before <- cor_before
  mod_aligned$cor_after <- res$cor_after

  mod_aligned
}
