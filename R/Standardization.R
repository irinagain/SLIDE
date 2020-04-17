
#' Standardization of multi-view data.
#'
#' Performs column-centering and standardization so that the Frobenius norm within each view is equal to one.
#'
#' @param X A n x p concatenated data matrix of views X_1,...,X_d.
#' @param pvec A vector of values p_1,....,p_d corresponding to the number of measurements within each  view.
#' @param center A logical indicator of whether the columns of X should be centered. The default value is TRUE.
#'
#' @return A list with the elements
#'   \item{X}{A n x p concatenated data matrix that has been column-centered and standardized so that Frobenius norm within each view is equal to one.}
#'   \item{svec}{A vector of largest singular values for each view after centering and standardization.}
#'   \item{norms}{A vector of Frobenius norms for each view before scaling.}
#'   \item{Xmean}{A vector of column means of X before centering. A zero vector is returned if \code{center = F}.}
#' @export
#' @examples
#' n = 100
#' p1 = 40
#' p2 = 60
#' X1 = matrix(rnorm(n*p1), n, p1)
#' X2 = matrix(rnorm(n*p2), n, p2)
#' X = cbind(X1, X2)
#' out = standardizeX(X, pvec = c(p1,p2))
standardizeX <- function(X, pvec, center = T) {
    d <- length(pvec)
    norms <- rep(0, d)
    svec <- rep(0, d)
    pcum <- c(0, cumsum(pvec))
    # Center each column
    if (center) {
        Xmean <- colMeans(X)
        X <- X - matrix(Xmean, nrow(X), ncol(X), byrow = T)
    } else {
        Xmean <- rep(0, ncol(X))
    }
    for (i in 1:d) {
        index <- c((pcum[i] + 1):pcum[i + 1])
        norms[i] <- sum(X[, index]^2)
        # Scale the dataset
        X[, index] <- X[, index]/sqrt(norms[i])
        # Calculate largest singular value
        svec[i] <- max(svd(X[, index], nu = 1)$d)
    }
    return(list(X = X, svec = svec, norms = norms, Xmean = Xmean))
}

