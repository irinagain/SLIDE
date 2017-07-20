updateVoptim1 <- function(XU, lambda, pvec) {
    d <- length(pvec)
    pcum <- c(0, cumsum(pvec))
    V <- matrix(0, nrow(XU), ncol(XU))
    for (i in 1:d) {
        # Index the measurements corresponding to ith dataset
        index <- c((pcum[i] + 1):pcum[i + 1])
        # Calculate norm of each column in XU
        norms <- sqrt(colSums(XU[index, ]^2))
        # Perform soft-thresholding for each column - only those that have norms > lambda
        indexr = norms > lambda
        if (sum(indexr) > 1) {
            V[index, indexr] = XU[index, indexr] %*% diag(1 - lambda/norms[indexr])
        } else if (sum(indexr) == 1) {
            V[index, indexr] = XU[index, indexr] * (1 - lambda/norms[indexr])
        }
    }
    return(V)
}

# Value of the objective function for the penalized optimization problem at the
# current iterate
evaluatef_optim1 <- function(XU, V, lambda, pvec) {
    f <- sum(V^2)/2 - sum(diag(crossprod(V, XU)))
    d <- length(pvec)
    pcum <- c(0, cumsum(pvec))
    for (i in 1:d) {
        index <- c((pcum[i] + 1):pcum[i + 1])
        f <- f + lambda * sum(sqrt(colSums(V[index, ]^2)))
    }
    return(f)
}

solve_optim1 <- function(X, lambda, pvec, k_max = 1000, eps = 1e-06,
     Ustart = NULL) {
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
        V = crossprod(X, U)
    }
    # Current function value
    f <- evaluatef_optim1(crossprod(X, U), V, lambda, pvec)
    k <- 1
    error <- 100
    while ((k < k_max) & (error[k] > eps)) {
        k <- k + 1
        # Update V
        V <- updateVoptim1(crossprod(X, U), lambda, pvec)
        # Update U
        nonzero = c(1:ncol(V))[colSums(abs(V)) > 0]
        if (length(nonzero) == 0) {
            # V is exactly zero, terminate
            f[k] = 0
            error[k] <- f[k - 1] - f[k]
            return(list(U = U, V = V, k = k, error = error, f = f))
        } else {
            XV_svd <- svd(X %*% V)
            U <- tcrossprod(XV_svd$u, XV_svd$v)
        }

        # Current function value
        f[k] <- evaluatef_optim1(crossprod(X, U), V, lambda, pvec)

        # Current difference in function values
        error[k] <- f[k - 1] - f[k]
    }
    if (k == k_max) {
        warning(paste("Algorithm didn't converge in ", k_max, " iterations!", sep = ""))
    }
    return(list(U = U, V = as.matrix(V), k = k, error = error, f = f))
}

# Iterative algorithm for penalized matrix factorization problem for the range of lambda values
# X - n x p concatenated data matrix, standardized lambdavec - vector of
# nonnegative tuning parameters pvec - values p_1,....,p_d corresponding to the
# number of measurements within each data type k_max - maximal number of
# iterations allowed eps - convergence tolerance as measured by the difference in
# objective function values reduced - if true, the rank of U is allowed to
# decrease with iterations (by taking only nonzero columns of V) rank_total -
# starting rank
#' Solves penalized matrix factorization problem for the range of lambda values
#'
#' @param X A n x p concatenated data matrix of views X_1,...,X_d.
#' @param pvec A vector of values p_1,....,p_d corresponding to the number of measurements within each data view.
#' @param lambda_seq An optional sequence of tuning parameters for the penalized matrix decomposition problem. By default, the algorithm generates its own sequence based on supplied values of \code{n_lambda}, \code{lambda_min} and \code{lambda_max}.
#' @param n_lambda A length of tuning parameter sequence. The default value is 50. It is only used when \code{lambda_seq = NULL}.
#' @param lambda_max A maximal value for tuning parameter. The default value is 1. If X is already standardized, it is recommended to set \code{lambda_max} to the largest singular value within the view.
#' @param lambda_min A minimal tuning parameter to be considered, the default value is 0.1
#' @param k_max A maximal number of allowed iterations, the default value is 1000.
#' @param eps A convergence tolerance criterion as measured by the differene in objective functions at successive iterations, the default value is 1e-06.
#' @export
#' @return A list with the elements
#' \item{lambda}{A sequence of tuning parameters used.}
#' \item{param}{A list with estimates of \code{U} and \code{V} obtained from solving the penalized matrix factorization problem with corresponding values of tuning parameter.}
#' @examples
#' n = 50
#' p1 = 20
#' p2 = 40
#' X1 = matrix(rnorm(n*p1), n, p1)
#' X2 = matrix(rnorm(n*p2), n, p2)
#' X = cbind(X1, X2)
#' out = standardizeX(X, pvec = c(p1,p2))
#' out_solve = solve_optim1_seq(X = out$X, pvec = c(p1,p2), lambda_max = max(out$svec), n_lambda = 30)
solve_optim1_seq <- function(X, pvec, lambda_seq = NULL,
     n_lambda = 50, lambda_max = 1, lambda_min = 0.01, k_max = 1000, eps = 1e-06) {
    # Check for user supplied lambda_sequence
    if (!is.null(lambda_seq)) {
        # Truncate the sequnce by lambda_max
        lambda_seq <- lambda_seq[lambda_seq <= lambda_max]
        # Order the sequence from largest to smallest
        lambda_seq <- sort(lambda_seq, decreasing = T)
        # Update n_lambda value
        n_lambda <- length(lambda_seq)
    } else {
        # generate the sequence of tuning parameters
        lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = n_lambda))
    }

    # Generate the same Ustart for each, so no need to repeat SVD Perform svd of X to
    # initialize U and V
    X_svd <- svd(X)
    rank_total <- sum(X_svd$d > eps)
    Ustart <- X_svd$u[, 1:rank_total]

    # Solve for each lambda from smallest to largest (but return from largest to
    # smallest)
    param = list()
    for (l in n_lambda:1) {
        # Use neighboring U as a new starting point
        if (l < n_lambda) {
            U_start = param[[l + 1]]$U
        }
        # Solve group lasso problem
        out <- solve_optim1(X, lambda = lambda_seq[l], pvec = pvec, k_max = k_max,
            eps = eps, Ustart = Ustart)
        param[[l]] <- list(U = out$U, V = out$V, fmin = min(out$f))
    }
    return(list(param = param, lambda = lambda_seq))
}

create_structure_from_solve <- function(out_solve, pvec, ratio_max = 2/3) {
    n_lambda <- length(out_solve$lambda)
    n <- nrow(out_solve$param[[1]]$U)
    # Keeps track of which structures were considered, i.e. if s_used = 1, then the
    # structure is used, if s_used = 0, then the structure is not used
    s_used <- rep(0, n_lambda)
    d <- length(pvec)
    pcum <- c(0, cumsum(pvec))
    Slist <- list()
    iter <- 0
    # Sold keeps track of the previous structure so that the comparison can be made
    Sold <- NULL
    rall <- ncol(out_solve$param[[1]]$V)
    for (t in 1:n_lambda) {
        # Calculate the number of nonzero columns in V
        r <- sum(colSums(abs(out_solve$param[[t]]$V)) > 0)
        if (r > 0) {
            if (r > ratio_max * min(n, sum(pvec))) {
                # Total rank is larger than 2/3 of possible rank
                return(list(Slist = Slist, s_used = s_used))
            }
            # Figure out current structure
            S <- matrix(0, d, rall)
            for (j in 1:d) {
                index <- c((pcum[j] + 1):pcum[j + 1])
                S[j, ] = (colSums(abs(out_solve$param[[t]]$V[index, ])) > 0)
                # Check whether the maximal rank is exceeded
                if (sum(S[j, ]) > ratio_max * min(n, pvec[j])) {
                  return(list(Slist = Slist, s_used = s_used))
                }
            }
            # Reduce to only have non-zero columns
            S = S[, colSums(S) > 0, drop = FALSE]
            # Compare with what is already there to see whether it is new or not
            if (is.null(Sold)) {
                # Have not seen this structure before
                iter = iter + 1
                s_used[t] = 1
                Slist[[iter]] <- S
                Sold <- S
            } else if (ncol(Sold) != ncol(S)) {
                # Have not seen this structure before
                iter = iter + 1
                s_used[t] = 1
                Slist[[iter]] <- S
                Sold <- S
            } else if (max(abs(Sold - S)) != 0) {
                # Have not seen this structure before
                iter = iter + 1
                s_used[t] = 1
                Slist[[iter]] <- S
                Sold <- S
            }
        }
    }  #end for lambda_seq
    return(list(Slist = Slist, s_used = s_used))
}

#' Creates a list of candidate structures for the SLIDE model
#'
#' Creates a list of candidate structures for the SLIDE model based on the solution path of penalized matrix factorization problem.
#'
#' The function solves the penalized matrix factorization problem for each value of from the sequence of tuning parameters lambda_1,..., lambda_m. The block-sparsity pattern of the resulting V(lambda_1), ..., V(lambda_m) is used to generate the sequence of binary structures S(lambda_1),..., S(lambda_m). This sequence is further trimmed to remove any repetitions, and to only keep structures with the number of nonzero columns being less than \code{ratio_max * min(n, sum(pvec))}.
#'
#'
#' @param X A n x p concatenated data matrix of views X_1,...,X_d.
#' @param pvec A vector of values p_1,....,p_d corresponding to the number of measurements within each  view.
#' @param lambda_seq An optional sequence of tuning parameters for the penalized matrix decomposition problem. By default, the algorithm generates its own sequence based on supplied values of \code{n_lambda}, \code{lambda_min} and \code{lambda_max}.
#' @param n_lambda A length of tuning parameter sequence. The default value is 50. It is only used when \code{lambda_seq = NULL}.
#' @param lambda_min A minimal value for tuning parameter. The default value is 0.01.
#' @param lambda_max A maximal value for tuning parameter. The default value is 1. If X is already standardized, it is recommended to set \code{lambda_max} to the largest singular value within the view. If \code{standardized = F}, this choice is made automatically by the algorithm.
#' @param ratio_max A maximal allowable rank of the binary structure as a ratio of maximal rank of X, the default value is 2/3.
#' @param standardized A logical indicator of whether X is centered and standardized. The default value is FALSE and the standardization is performed within the function.
#' @param eps A convergence tolerance criterion as measured by the differene in objective functions at successive iterations, the default value is 1e-06.
#' @param k_max A maximal number of allowed iterations, the default value is 1000.
#'
#' @return A list with the elements
#'   \item{lambda_seq}{A sequence of tuning parameters used to generate the structure sequence.}
#'   \item{Slist}{A list of distinct binary structures which forms a candidate set for SLIDE model.}
#'   \item{id_used}{A binary vector indicating the correspondence between \code{lambda_seq} and \code{Slist}. The zero values indicate the tuning parameters that either lead to the same structures, or resulted in structures with the rank larger than allowed by \code{ratio_max}.}
#' @export
#' @examples
#' n = 50
#' p1 = 20
#' p2 = 40
#' X1 = matrix(rnorm(n*p1), n, p1)
#' X2 = matrix(rnorm(n*p2), n, p2)
#' X = cbind(X1, X2)
#' pvec = c(p1,p2)
#'
#' # Unstandardized
#' slist = create_structure_list(X, pvec)
#'
#' # Already standardized
#' out = standardizeX(X, pvec = c(p1,p2))
#' slist = create_structure_list(out$X, pvec, lambda_max = max(out$svec), standardized = TRUE)
create_structure_list <- function(X, pvec, lambda_seq = NULL, n_lambda = 50, lambda_min = 0.01, lambda_max = 1, ratio_max = 2/3, standardized = F, eps = 1e-6, k_max = 1000) {
  if (standardized == F) {
    out <- standardizeX(X, pvec, center = T)
    lambda_max <- min(lambda_max, max(out$svec))
    X <- out$X
  }

  # Solve penalized matrix optimizaiton problem for the sequence of lambda values
  if (is.null(lambda_seq)){
    out <- solve_optim1_seq(X = X, pvec = pvec, n_lambda = n_lambda, lambda_max = lambda_max,
                          lambda_min = lambda_min, eps = eps, k_max = k_max)
  }else {
    lambda_seq <- sort(lambda_seq[(lambda_seq <= lambda_max)&(lambda_seq >= lambda_min)], decreasing = T)
    if (length(lambda_seq) < 1){
      stop(paste("Supplied sequence of tuning parameters is outside the range of lamdba_min", lambda_min,"and lambda_max", lambda_max, sep=" "))
    }
    out <- solve_optim1_seq(X = X, pvec = pvec, lambda_seq = lambda_seq, eps = eps, k_max = k_max)
  }

  # Form the list of candidate structures
  out_struct <- create_structure_from_solve(out, pvec, ratio_max = ratio_max)

  return(list(lambda_seq = out$lambda, Slist = out_struct$Slist, id_used = out_struct$s_used))
}


