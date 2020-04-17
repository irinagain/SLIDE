#' Fits SLIDE model for X using the pre-specified structure S.
#'
#' @param X A n x p concatenated data matrix of views X_1,...,X_d.
#' @param pvec A vector of values p_1,....,p_d corresponding to the number of measurements within each  view.
#' @param S A binary matrix with nonzero columns of size d x r.
#' @param Ustart An n X r optional starting value for U, if not supplied the first r left singular vectors of X are used.
#' @param eps A convergence tolerance criterion, the default value is 1e-6.
#' @param k_max A maximal number of allowable iterations, the default values is 1000.
#' @param standardized A logical indicator of whether X is centered and standardized. The default value is FALSE and the standardization is performed within the function.
#'
#' @return A list with the elements
#'   \item{U}{A n x r score matrix for the SLIDE model.}
#'   \item{V}{A p x r loadings matrix for the SLIDE model with sparsity pattern according to S.}
#'   \item{error}{Tolerance value at convergence.}
#' @export
#' @examples
#' n = 100
#' p1 = 25
#' p2 = 25
#' data = generateModel1(n = n, pvec = c(p1, p2))
#'
#' # Specify binary structure
#' S = matrix(c(1,1,0,1,0,1),nrow = 2, ncol = 3)
#'
#' # Unstandardized
#' fit_slide = slide_givenS(data$X, pvec = c(p1,p2), S = S)
#'
#' # Standardized
#' out = standardizeX(data$X, pvec = c(p1,p2))
#' fit_slide = slide_givenS(out$X, pvec = c(p1,p2), S = S, standardized = TRUE)
slide_givenS <- function(X, pvec, S, Ustart = NULL, eps = 1e-06, k_max = 1000, standardized = F) {
    d <- length(pvec)
    pcum <- c(0, cumsum(pvec))
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

    # Initialize V
    V <- matrix(0, sum(pvec), r)
    for (i in 1:d) {
        # Index the measurements corresponding to ith dataset
        index <- c((pcum[i] + 1):pcum[i + 1])
        # Identify columns in V that are present for the dataset i
        nonzero <- c(1:r)[S[i, ] == 1]
        if (length(nonzero) > 0) {
            V[index, nonzero] <- crossprod(X[, index], U[, nonzero])
        }
    }

    error <- 1000
    k <- 1
    while ((k < k_max) & (error[k] > eps)) {
        k <- k + 1
        UVold <- tcrossprod(U, V)
        # Refit U
        XV_svd <- svd(X %*% V)
        U <- tcrossprod(XV_svd$u, XV_svd$v)

        # Refit V
        for (i in 1:d) {
            # Index the measurements corresponding to ith dataset
            index <- c((pcum[i] + 1):pcum[i + 1])
            # Identify columns in V that are present for the dataset i
            nonzero <- c(1:r)[S[i, ] == 1]
            V[index, nonzero] <- crossprod(X[, index], U[, nonzero])
        }
        # Calculate the difference due to refitting
        error[k] <- sum((tcrossprod(U, V) - UVold)^2)
    }
    if (k == k_max) {
        warning(paste("Algorithm didn't converge in ", k_max, " iterations!", sep = ""))
    }
    return(list(U = U, V = V, error = error))
}
