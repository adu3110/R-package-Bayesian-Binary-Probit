
#' Samples From Gaussian Distribution
#'
#' Generate samples from Gaussian distribution
#' Uses central limit theorem to generate gaussian samples from uniform samples
#' @param num_samples: Number of samples to be generated
#' @param mean_norm : Mean of Gaussian distribution
#' @param sd_norm : Standard Deviation of Gaussian distribution
#' @param num_uniform_samples: Number of Uniform samples to be used for generating Gaussian sample by averaging
#' @return samples: Gaussian samples
#' @export
gaussiansamplesbyCLT <- function(num_samples, mean_norm = 0, sd_norm = 1,
                                 num_uniform_samples = 50){

  if(num_samples < 1){
    stop("Number of samples must be at least one")
  }

  sampling_result <- .C('RandNorm',
                        num_samples = as.integer(num_samples),
                        num_unif_samples = as.integer(num_uniform_samples),
                        mean_norm = as.double(mean_norm),
                        var_norm = as.double(sd_norm^2),
                        samples = as.double(matrix(data=0, nrow = num_samples, ncol = 1)))

  samples <- c(sampling_result$samples)

  return(samples)

}


#' Samples From One Sided Gaussian Distribution
#'
#' Generate samples from one sided Gaussian distribution
#' Uses central limit theorem to generate gaussian samples from uniform samples
#' @param num_samples: Number of samples to be generated
#' @param mean_norm : Mean of Gaussian distribution
#' @param sd_norm : Standard Deviation of Gaussian distribution
#' @param cut_off : Cut off value for truncation
#' @param side : right or left based on which sideddistribution to be generated
#' @param num_uniform_samples: Number of Uniform samples to be used for generating Gaussian sample by averaging
#' @return samples: one sided Gaussian samples
#' @export

samplegaussianonesided <- function(num_samples, mean_norm = 0, sd_norm = 1,
                                   cut_off = 0, side = c("right", "left"),
                                   num_uniform_samples = 50){

  if(num_samples < 1){
    stop("Number of samples must be at least one")
  }

  side <- match.arg(side)

  greater <- as.integer(1)
  if(side == "left"){
    greater <- as.integer(0)
  }

  sampling_result <- .C('RandTruncNorm',
                        num_samples = as.integer(num_samples),
                        num_unif_samples = as.integer(num_uniform_samples),
                        mean_norm = as.double(mean_norm),
                        var_norm = as.double(sd_norm^2),
                        cutoff = as.double(cut_off),
                        greater = greater,
                        samples = as.double(matrix(data=0, nrow = num_samples, ncol = 1)),
                        attempts_failed = as.integer(0))

  if(sampling_result$attempts_failed == 1){
    stop("Truncated Normal sampling failed after 100000 attempts,
         probably the cutoff is too extreme")
  }

  samples <- c(sampling_result$samples)

  return(samples)

  }


#' Matrix Inversion
#'
#' Invert a matrix by Gauss Elimination method
#' @param mat: a square matrix
#' @return matinv: inverse of mat
#' @export
gausseliminationinverse <- function(mat){

  if(!is.matrix(mat)){
    stop("Input is not a matrix")
  }

  nrow_mat <- nrow(mat)
  ncol_mat <- nrow(mat)

  if(nrow_mat != ncol_mat){
    stop("Input is not a square matrix")
  }

  inversion_result <- .C('GaussEliminationInverse',
                         A = as.double(mat),
                         num_rows = as.integer(nrow_mat),
                         Ainv = as.double(matrix(data = 0, nrow = nrow_mat, ncol = nrow_mat)),
                         inverse_failed = as.integer(0))

  if(inversion_result$inverse_failed == 1){
    stop("Inversion failed due to leading 0")
  }

  matinv <- matrix(data = inversion_result$Ainv, nrow = nrow_mat, ncol = nrow_mat)

  return(matinv)

}


#' Cholesky Decomposition
#'
#' Lower Triangular Cholesky decomposition of a matrix
#' @param mat: a square matrix
#' @return L: cholesky lower traingular decomposition of mat
#' @export

choleskylower <- function(mat){

  if(!is.matrix(mat)){
    stop("Input is not a matrix")
  }

  nrow_mat <- nrow(mat)
  ncol_mat <- nrow(mat)

  if(nrow_mat != ncol_mat){
    stop("Input is not a square matrix")
  }

  cholesky_result <- .C('CholeskyDecomposition',
                        A = as.double(mat),
                        L = as.double(matrix(data = 0, nrow = nrow_mat, ncol = nrow_mat)),
                        num_rows = as.integer(nrow_mat),
                        positive_deifinte_failed = as.integer(0))

  if(cholesky_result$positive_deifinte_failed == 1){
    stop("Matrix is not positive definite")
  }

  L <- matrix(data = cholesky_result$L, nrow = nrow_mat, ncol = nrow_mat)

  return(L)

}


#' MCMC chain for Bayesian Binary Probit Regression
#'
#' Run MCMC chain and get estimates for Logistic Regression Coefficients
#' @param y : vector with either 0 or 1
#' @param X : matrix with independent variables as columns, and first column is a vector of 1s
#' @param covar_prior: covariance matrix of 0 mean gaussian priors for coefficients
#' @param num_mcmc: Number of MCMC samples
#' @param burn_in: Number of burn in samples
#' @param thinning: Thinning interval
#' @return coeffs: data frame containing coeffcient estimates and their standard deviation
#' @export

bayesprobitmcmc <- function(y, X, covar_prior, num_mcmc=10000, burn_in=1000, thinning = 1){


  if(!(is.matrix(y) || is.vector(y))){
    stop("y must be a vector or a matrix")
  }

  if(!is.matrix(X)){
    stop("X must be a matrix")
  }

  if(any(is.na(y)) || any(is.infinite(y))){
    stop("na or infinite values in y")
  }

  if(any(is.na(X)) || any(is.infinite(X))){
    stop("na or infinite values in X")
  }


  y_vals <- unique(y)

  if(length(y_vals) > 2){
    stop("Only Binary values are supported for respose")
  }

  if(!(all(y_vals %in% c(0, 1)))){
    stop("Allowed values for y are 0's and 1's")
  }

  if(!(all(is.numeric(X)))){
    stop("Only numeric values are supported for X")
  }

  if(length(y) != nrow(X)){
    stop("length of y should be equal to number of rows in X")
  }

  if(!is.matrix(covar_prior)){
    stop("covar_matrix must be a matrix")
  }

  if(nrow(covar_prior) != ncol(covar_prior)){
    stop("covar_matrix must be a square matrix")
  }

  if(nrow(covar_prior) != ncol(X)){
    stop("Number of rows covar_matrix must be equal to number of independent variables + 1 (Intercept)")
  }

  if(!(all(is.numeric(covar_prior)))){
    stop("Only numeric values are supported for covar_prior")
  }

  if(burn_in >= num_mcmc){
    stop("Burn in should be less than Number of MCMC samples")
  }

  if(thinning >= (num_mcmc - burn_in)){
    stop("Thinning should be smaller than the number of
         MCMC samples after removing burn in")
  }

  if(thinning < 1){
    stop("Thinning should be at least 1")
  }

  num_observations <- nrow(X)
  num_variables <- ncol(X)

  coeff_samples <- matrix(data=0, nrow = num_variables, ncol = num_mcmc)

  bayes_probit_samples <- .C('BayesProbitC',
                             y = as.integer(y),
                             X = as.double(X),
                             covar_prior = as.double(covar_prior),
                             num_variables = as.integer(num_variables),
                             num_observations = as.integer(num_observations),
                             num_mcmc = as.integer(num_mcmc),
                             coeff_samples = as.double(coeff_samples),
                             positive_deifinte_failed = as.integer(0),
                             sample_attemps_failed = as.integer(0),
                             inverse_failed = as.integer(0))


  if(bayes_probit_samples$inverse_failed == 1){
    stop("Matrix Inversion by Gauss Elimination failed due to leading 0")
  }

  if(bayes_probit_samples$positive_deifinte_failed == 1){
    stop("Cannot Do Cholesky Decomposition for Non-Positive Definite Matrix")
  }

  if(bayes_probit_samples$sample_attemps_failed == 1){
    stop("Truncated Normal sampling failed after 100000 attempts,
         probably the cutoff is too extreme")
  }

  coeff_samples <- matrix(data = bayes_probit_samples$coeff_samples,
                          nrow = num_variables, ncol = num_mcmc)

  colvals_thining <- seq((burn_in+1), num_mcmc, by = thinning)

  coeffs_mean <- apply(coeff_samples[, colvals_thining], 1, mean)
  coeffs_sd <- apply(coeff_samples[, colvals_thining], 1, sd)

  coeffs <- data.frame(estimates = coeffs_mean,
                       standard_deviation = coeffs_sd)

  rownames(coeffs) <- colnames(X)

  return(coeffs)

  }



#' Binary Probit Regression
#'
#' Binary Probit Regression function
#' @param formula : formula for the model
#' @param data : data frame containing input data
#' @param covar_prior: covariance matrix of 0 mean gaussian priors for coefficients
#' @param num_mcmc: Number of MCMC samples
#' @param burn_in: Number of burn in samples
#' @param thinning: Thinning interval
#' @return coeffs: data frame containing coeffcient estimates and their standard deviation
#' @export

bayesianbinaryprobit <- function(formula, data, covar_prior,
                                 num_mcmc=10000, burn_in=1000, thinning = 1){

  coeffs <- NULL
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)

  mf <- eval(expr = mf, envir = parent.frame())

  X_prime <- as.matrix(mf[,2:ncol(mf)])

  X <- cbind(matrix(data = 1, nrow = nrow(mf), ncol = 1), X_prime)

  colnames(X)[1] <- "Intercept"

  y_original <- as.matrix(mf[, 1])
  y_factors <- as.factor(y_original)

  y_levels <- levels(y_factors)

  if(length(y_levels) != 2){
    stop("Only two unique values allowed for y")
  }

  y <- matrix(data = 0, nrow = nrow(X), ncol = 1)

  y[y_original == y_levels[2]] <- 1

  coeffs <- bayesprobitmcmc(y, X, covar_prior, num_mcmc, burn_in, thinning)

  rownames(coeffs) <- paste(rownames(coeffs), y_levels[2], sep = '.')

  return(coeffs)

}


