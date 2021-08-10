
#' Selects binary structure for SLIDE model using bi-cross-validation
#'
#'Selects binary structure for SLIDE model using bi-cross-validation based on the supplied set of candidate structures
#'
#' @param X A n x p concatenated data matrix of views X_1,...,X_d.
#' @param pvec A vector of values p_1,....,p_d corresponding to the number of measurements within each data view.
#' @param structure_list A list of distinct binary structures which forms a candidate set for SLIDE model.
#' @param n_fold A number of folds for rows, the default value is 3.
#' @param p_fold A number of folds for columns, the default value is 3.
#' @param k_max A maximal number of allowed iterations for SLIDE model fitting, the default value is 1000.
#' @param eps A convergence tolerance criterion, the default value is 1e-6.
#' @param center A logical indicator of whether the data is centered, the default value is TRUE.
#'
#' @return A list with the elements
#'   \item{structure_list}{A list of candidate structures, same as the input \code{structure_list}.}
#'   \item{structure_min}{A binary matrix corresponding to the structure with minimal error across the folds.}
#'   \item{error_mean}{Total error values across the folds for each structure from the list.}
#'   \item{error}{Matrix of error values for each fold and each structure.}
#' @export
#'
#' @examples
#' n = 100
#' p1 = 25
#' p2 = 25
#' data = generateModel1(n = n, pvec = c(p1, p2))
#'
#' # Standardize the data
#' out = standardizeX(data$X, pvec = c(p1,p2))
#'
#' # Create a list of candidate structures
#' out_struct <- create_structure_list(out$X, data$pvec, lambda_max = max(out$svec),
#'  standardized = TRUE)
#'
#' # Select a structure from the list using bcv procedure
#' out_bcv <- slide_BCV(out$X, data$pvec, structure_list = out_struct$Slist)
#' out_bcv$structure_min
slide_BCV <- function(X, pvec, structure_list, n_fold = 3, p_fold = 3, k_max = 1000,
    eps = 1e-06, center = T) {
    d <- length(pvec)
    n <- nrow(X)
    p <- ncol(X)
    pcum <- c(0, cumsum(pvec))

    # Get length of structure list
    n_structure = length(structure_list)

    # Get the folds For samples
    fold_id_n <- sample(rep(seq_len(n_fold), length.out = n))
    # For measurements
    fold_id_p <- c()
    for (i in 1:d) {
        fold_id_p <- c(fold_id_p, sample(rep(seq_len(p_fold), length.out = pvec[i])))
    }
    
    # Since bcv submatrices have lower dimensions than original ones, need to ensure that the tried ranks (number of columns in structure_list) do not exceed the dimensions of subfolds

    # Save prediction error
    error <- matrix(0, n_fold * p_fold, n_structure)
    iter <- 0
    for (k in 1:n_fold) {
        for (j in 1:p_fold) {
            iter <- iter + 1
            # Update the pvec
            pvec_foldj = pvec
            for (i in 1:d) {
                index <- c((pcum[i] + 1):pcum[i + 1])
                pvec_foldj[i] = sum(fold_id_p[index] != j)
            }
            pcum_fold = c(0, cumsum(pvec_foldj))
            # Standardize each dataset
            out_s <- standardizeX(X = X[fold_id_n != k, fold_id_p != j], pvec = pvec_foldj,
                center = center)
            Xfold <- out_s$X

            # Consider all structure from the list, fit the model, evaluate BCV error
            for (l in 1:n_structure) {
                # Find U and V based on the given structure
                out <- slide_givenS(X = Xfold, pvec = pvec_foldj, S = structure_list[[l]],
                  k_max = k_max, eps = eps)
                Vadj <- out$V
                for (i in 1:d) {
                  index <- c((pcum_fold[i] + 1):pcum_fold[i + 1])
                  Vadj[index, ] <- out$V[index, ] * sqrt(out_s$norms[i])
                }
                # Perform SVD on the output
                out_svd <- svd(tcrossprod(out$U, Vadj))
                if (max(out_svd$d) < eps) {
                  error[iter, l] <- d
                } else {
                  # because of est_givenranks, only have nonzero columns already
                  if (center == F) {
                    # No centering was done
                    tildeV = crossprod(X[fold_id_n != k, fold_id_p == j], out$U)
                    tildeU = X[fold_id_n == k, fold_id_p != j] %*% Vadj %*% solve(crossprod(Vadj))
                    # Calculate the error
                    pcum_notfold <- pcum - pcum_fold
                    for (i in 1:d) {
                      index <- c((pcum_notfold[i] + 1):pcum_notfold[i + 1])
                      error[iter, l] = error[iter, l] + sum((X[fold_id_n == k, which(fold_id_p == j)[index]] - tcrossprod(tildeU, tildeV[index, ]))^2)/sum(X[fold_id_n ==k, which(fold_id_p == j)[index]]^2)
                    }
                  } else {
                    n_newf = sum(fold_id_n != k)  # number of samples in the submatrix for model fitting
                    # Centering was done, stored in out_s$Xmean
                    tildeV = cbind(colMeans(X[fold_id_n != k, fold_id_p == j]) * sqrt(n_newf),
                      crossprod(X[fold_id_n != k, fold_id_p == j], out$U))
                    Vnew = cbind(sqrt(n_newf) * out_s$Xmean, Vadj)
                    # Adjust for potential singularity here
                    Vsvd = svd(Vnew)
                    Vproj = Vsvd$u[, Vsvd$d > eps] %*% diag(1/Vsvd$d[Vsvd$d > eps]) %*%
                      t(Vsvd$v[, Vsvd$d > eps])
                    tildeU = X[fold_id_n == k, fold_id_p != j] %*% Vproj
                    # Calculate the error
                    pcum_notfold <- pcum - pcum_fold
                    for (i in 1:d) {
                      index <- c((pcum_notfold[i] + 1):pcum_notfold[i + 1])
                      error[iter, l] = error[iter, l] + sum((X[fold_id_n == k, which(fold_id_p ==
                        j)[index]] - tcrossprod(tildeU, tildeV[index, ]))^2)/sum(scale(X[fold_id_n ==
                        k, which(fold_id_p == j)[index]], scale = F)^2)
                    }
                  }
                }
            }  # end for lamda_seq
        }  # end for p folds
    }  # end for n folds
    # Calculate average prediction error for each tuning parameter
    error_sum <- colSums(error)
    # Find largest lambda corresponding to minimal average error, here assume that
    # lambda_seq is sorted from largest to smallest
    id_min <- min(c(1:n_structure)[error_sum == min(error_sum)])
    structure_min <- structure_list[[id_min]]
    return(list(structure_list = structure_list, structure_min = structure_min, error_sum = error_sum,
        error = error))
}
