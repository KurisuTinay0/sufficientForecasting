#' Conformal inference of the sufficient forecasting
#'
#' @param y Response, T by 1 matrix
#' @param X Predictors, p by T matrix
#' @param newX New predictors, a vector contains p entries (or \code{NULL})
#' @param type \code{LM} or \code{LLM} (default = \code{LM})
#' @param K The number of common factors (default = obtained by \code{\link{getK}})
#' @param L The number of predictive indices, L is required to be no greater than K (default = 1)
#' @param alpha Mis-coverage rate
#' @param discretization Hyperparameter in SIR (default = \code{TRUE})
#' @param nslices Hyperparameter in SIR (default = 10)
#'
#' @return A list with components
#' \describe{
#'   \item{yhat}{Out-of-sample forecast for \code{newX}; or in-sample forecast for the last observed data point if \code{newX} is \code{NULL}}
#'   \item{ci_lower}{Lower bound of conformal interval}
#'   \item{ci_upper}{Upper bound of conformal interval}
#'   \item{psi_transpose_ft}{Matrix of sufficient predictive indices ψ'f_t (T × L)}
#'   \item{Phi}{Matrix of sufficient dimension reduction directions (K × L)}
#'   \item{factors}{Matrix of estimated factor scores (T × K)}
#' }
#' @export
#' @references
#' Yu, X., Yao, J. and Xue, L. (2022), Nonparametric estimation and conformal inference
#' of the sufficient forecasting with a diverging number of factors,
#' \emph{Journal of Business & Economic Statistics} 40(1), 342–354.
#' @examples
#' utils::data(dataExample, package = "sufficientForecasting")
#' result <- SF.CI(dataExample$y, dataExample$X, type = "LM", alpha = 0.05)
#' # Access sufficient predictive indices
#' psi_ft <- result$psi_transpose_ft
#' # Access dimension reduction directions
#' Phi_matrix <- result$Phi
#' # Access estimated factors
#' factor_scores <- result$factors

SF.CI <- function(
    y, X, newX = NULL, type = "LM", K = "default", L = 1, alpha = 0.1,
    discretization = TRUE, nslices = 10) {
  
  # Load required packages
  if (!requireNamespace("dr", quietly = TRUE)) {
    stop("Package 'dr' is required for sir.cov function. Please install it.")
  }
  
  # Default K determination
  if (K == "default") {
    K <- getK(y, X, 12)
  }
  
  # Input validation and dimension checking
  pp <- nrow(X)
  TT <- ncol(X)
  
  # Warning and error messages
  if (!is.matrix(X) | !is.matrix(y)) {
    stop("X and y must be matrices")
  }
  
  if (dim(y)[1] != dim(X)[2]) {
    stop("X must be a P by T matrix and y must be a T by 1 matrix")
  }
  
  if (!is.null(newX)) {
    if (!is.vector(newX) | length(newX) != pp) {
      stop("new predictors must be a vector containing p entries")
    }
  }
  
  if (L > K | L < 1 | L %% 1 != 0) {
    stop("invalid L: L must be an integer and must be not smaller than 1 and not greater than K")
  }
  
  if (!is.numeric(K)) {
    stop("invalid K: try K = 'default'")
  }
  
  if (K < 1 | K %% 1 != 0) {
    stop("invalid K: K must be an integer and not smaller than 1")
  }
  
  maxi <- max(L, 2)
  if (nslices < maxi | nslices %% 1 != 0) {
    stop("invalid nslices: nslices must be an integer and >= max{L,2} is required")
  }
  
  if (alpha >= 1 | alpha <= 0) {
    stop("invalid alpha value")
  }
  
  if (type != "LM" & type != "LLM") {
    stop("type must be one of 'LM' and 'LLM'")
  }
  
  if (alpha >= 1 | alpha <= 0) {
    stop("alpha must be (0,1)")
  }
  
  # Helper function for linear model fitting
  SF_lm_fit <- function(yy, Predictor) {
    T0 <- length(yy)
    Predictor <- as.matrix(Predictor)
    xtemp <- Predictor
    ytemp <- yy
    beta <- solve(t(xtemp) %*% xtemp) %*% (t(xtemp) %*% ytemp)
    eps <- yy - Predictor %*% beta
    return(eps)
  }
  
  # Helper function for local linear regression fitting
  SF_LLR_fit <- function(yy, Predictor) {
    T0 <- length(yy)
    Predictor <- data.frame(Predictor)
    xtemp <- Predictor
    ytemp <- yy
    LLR.fit <- loess(ytemp ~ ., data = xtemp, span = 0.9, degree = 1,
                     control = loess.control(surface = "direct"))
    eps <- LLR.fit$residuals
    return(eps)
  }
  
  # Main function to calculate residuals and sufficient predictive indices
  hateps <- function(yy, XX) {
    tt <- dim(XX)[2]
    
    # PCA for factors and loadings
    PCA <- eigen(t(XX) %*% XX)
    hFF <- as.matrix(PCA$vectors[, 1:K] * sqrt(tt))   # T × K factor scores
    hBB <- XX %*% hFF / tt                           # p × K factor loadings
    
    # Sliced inverse regression for sufficient directions
    hFF.cov <- dr::sir.cov(hFF, yy, discretization, nslices)
    Phi.h <- eigen(hFF.cov)$vectors[, 1:L]           # K × L sufficient directions
    Predictor <- hFF %*% Phi.h                       # T × L predictive indices
    
    # Calculate residuals
    eps_SF_LML <- SF_lm_fit(yy, Predictor)
    eps_SF_LLRL <- SF_LLR_fit(yy, Predictor)
    
    epsmat <- cbind(eps_SF_LML, eps_SF_LLRL)
    colnames(epsmat) <- c("SF_LML", "SF_LLRL")
    
    return(list(epsmat = epsmat, 
                psi_transpose_ft = Predictor,  # ψ'f_t matrix
                hFF = hFF,                     # factor scores
                Phi.h = Phi.h))                # sufficient directions
  }
  
  # Function to calculate p-values for conformal inference
  pyhat <- function(epsvec) {
    Tlast <- length(epsvec)
    pvalue <- mean(abs(epsvec) >= abs(epsvec)[Tlast])
    return(pvalue)
  }
  
  # Main conformal inference function
  SF_conformal <- function(yy, XX) {
    Tlast <- length(yy)
    
    # Get residuals and sufficient predictive indices
    eps_result <- hateps(yy, XX)
    eps.hat <- eps_result$epsmat
    psi_transpose_ft <- eps_result$psi_transpose_ft
    
    # Get point forecast
    yhat <- SF.SIR(y = y, X = X, newX = newX, type = type,
                   K = K, L = L, discretization = discretization,
                   nslices = nslices)
    
    # Conformal inference for confidence intervals
    y.grid <- seq(min(yy) - 2 * sd(yy), max(yy) + 2 * sd(yy), length = 200)
    p.vec <- matrix(NA, length(y.grid), 2)
    colnames(p.vec) <- c("SF_LML", "SF_LLRL")
    
    for (i in 1:length(y.grid)) {
      yy[Tlast] <- y.grid[i]
      eps.hat.mat <- hateps(yy, XX)$epsmat
      p.vec[i, ] <- apply(eps.hat.mat, 2, pyhat)
    }
    
    ci_SF_LML <- y.grid[p.vec[, "SF_LML"] > alpha]
    ci_SF_LLRL <- y.grid[p.vec[, "SF_LLRL"] > alpha]
    
    ci_lower <- c(min(ci_SF_LML), min(ci_SF_LLRL))
    ci_upper <- c(max(ci_SF_LML), max(ci_SF_LLRL))
    resultmat <- rbind(yhat, ci_lower, ci_upper)
    
    return(list(result = round(resultmat, 4),
                psi_transpose_ft = psi_transpose_ft,
                Phi.h = eps_result$Phi.h,
                hFF = eps_result$hFF))
  }
  
  # In-sample forecasting
  if (is.null(newX)) {
    out <- SF_conformal(y, X)
    if (type == "LM") {
      return(list(yhat = out$result[, 1], 
                  ci_lower = out$result[2, 1], 
                  ci_upper = out$result[3, 1],
                  psi_transpose_ft = out$psi_transpose_ft,
                  Phi = out$Phi.h,
                  factors = out$hFF))
    } else if (type == "LLM") {
      return(list(yhat = out$result[, 2], 
                  ci_lower = out$result[2, 2], 
                  ci_upper = out$result[3, 2],
                  psi_transpose_ft = out$psi_transpose_ft,
                  Phi = out$Phi.h,
                  factors = out$hFF))
    }
  }
  
  # Out-of-sample forecasting
  if (!is.null(newX)) {
    cX <- cbind(X, newX)
    cy <- c(y, mean(y))
    out <- SF_conformal(cy, cX)
    if (type == "LM") {
      return(list(yhat = out$result[, 1], 
                  ci_lower = out$result[2, 1], 
                  ci_upper = out$result[3, 1],
                  psi_transpose_ft = out$psi_transpose_ft,
                  Phi = out$Phi.h,
                  factors = out$hFF))
    } else if (type == "LLM") {
      return(list(yhat = out$result[, 2], 
                  ci_lower = out$result[2, 2], 
                  ci_upper = out$result[3, 2],
                  psi_transpose_ft = out$psi_transpose_ft,
                  Phi = out$Phi.h,
                  factors = out$hFF))
    }
  }
}

# Helper function to determine optimal number of factors K
getK <- function(y, X, maxK = 12) {
  # Eigenvalue ratio method to determine number of factors
  TT <- ncol(X)
  PCA <- eigen(t(X) %*% X)
  eigenvalues <- PCA$values[1:maxK]
  
  ratios <- eigenvalues[1:(maxK-1)] / eigenvalues[2:maxK]
  K <- which.max(ratios)
  
  return(max(1, min(K, maxK)))
}

# SF.SIR function for sufficient forecasting
SF.SIR <- function(y, X, newX = NULL, type = "LM", K = 3, L = 1, 
                   discretization = TRUE, nslices = 10) {
  
  # Implementation of sufficient forecasting using SIR
  TT <- length(y)
  pp <- nrow(X)
  
  # PCA factor extraction
  PCA <- eigen(t(X) %*% X)
  hFF <- as.matrix(PCA$vectors[, 1:K] * sqrt(TT))
  
  # Sliced inverse regression
  hFF.cov <- dr::sir.cov(hFF, y, discretization, nslices)
  Phi.h <- eigen(hFF.cov)$vectors[, 1:L]
  
  # Predictive indices
  Predictor <- hFF %*% Phi.h
  
  if (type == "LM") {
    # Linear model forecasting
    model <- lm(y[-1] ~ Predictor[-TT, ])
    if (is.null(newX)) {
      # In-sample forecast for last observation
      yhat <- predict(model, newdata = data.frame(Predictor[TT, , drop = FALSE]))
    } else {
      # Out-of-sample forecast
      new_factors <- t(newX) %*% (X %*% hFF / TT) / (TT * pp)
      yhat <- predict(model, newdata = data.frame(t(new_factors %*% Phi.h)))
    }
  } else if (type == "LLM") {
    # Local linear model forecasting
    model <- loess(y[-1] ~ Predictor[-TT, 1], 
                   span = 0.9, degree = 1,
                   control = loess.control(surface = "direct"))
    if (is.null(newX)) {
      yhat <- predict(model, newdata = data.frame(Predictor[TT, 1, drop = FALSE]))
    } else {
      new_factors <- t(newX) %*% (X %*% hFF / TT) / (TT * pp)
      yhat <- predict(model, newdata = data.frame(new_factors %*% Phi.h[, 1]))
    }
  }
  
  return(yhat)
}
