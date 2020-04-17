#' Solve optimization algorithm for one value of lambda
#'
#' @param X 
#' @param lambda 
#' @param pvec 
#' @param k_max 
#' @param eps 
#' @param Ustart 
#'
#' @return
#' @export
#'
#' @examples
solve_optim <- function(X, lambda, pvec, k_max = 1000, eps = 1e-06, Ustart = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(Ustart)) {
    # Perform svd of X to initialize U and V
    X_svd <- svd(X)
    rank_total <- sum(X_svd$d > eps)
    U <- X_svd$u[, 1:rank_total]
    V <- X_svd$v[, 1:rank_total] %*% diag(X_svd$d[1:rank_total])
  } else {
    U = Ustart
  }
  # Call to C function
  out = solve_optimC(X, U, lambda, pvec, k_max, eps)
  if (out$k + 1 >= k_max) {
    warning(paste("Algorithm didn't converge in ", k_max, " iterations!", sep = ""))
    return(list(U = out$U, V = out$V, k = out$k, error = out$error, f = out$f))
  }else{
    return(list(U = out$U, V = out$V, k = out$k, error = out$error[1:(out$k+1)], f = out$f[1:(out$k+1)]))
  }
}

apply_slide_givenS <- function(X, pvec, S, Ustart = NULL, eps = 1e-06, k_max = 1000, standardized = F) {
  r <- ncol(S)
  if (standardized == F) {
    out <- standardizeX(X, pvec, center = T)
    X <- out$X
  }
  
  # Initialize U
  if (is.null(Ustart)) {
    outs <- svd(X, nu = r)
    U <- outs$u
  } else if ((ncol(Ustart) != r) | (nrow(Ustart) != nrow(X))) {
    stop("Supplied dimensions of Ustart don't match the dimensions in X and S")
  }else{
    U <- Ustart
  }
  
  out = slide_givenS_C(X, pvec, S, U, eps, k_max)
  if (out$k + 1 >= k_max) {
    warning(paste("Algorithm didn't converge in ", k_max, " iterations!", sep = ""))
    return(list(U = out$U, V = out$V, k = out$k, error = out$error))
  }else{
    return(list(U = out$U, V = out$V, k = out$k, error = out$error[1:(out$k+1)]))
  }
}

