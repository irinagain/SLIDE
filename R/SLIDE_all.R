#' Fits SLIDE model to the multi-view data X
#'
#' Fits SLIDE model to the multi-vew data X
#'
#' The function automatically column-centers and standardizes X so that the Frobenius norm within each view is equal to one (using function \code{standardizeX}). The function then performs 3 steps of SLIDE workflow:
#' \enumerate{
#' \item Create a candidate list of distinct binary structures S_1,...,S_m (using function \code{create_structure_list}).
#' \item Select one structure S from the list using the bi-cross-validation (using function \code{slideBCV}).
#' \item Fit SLIDE model with selected structure S (using function \code{slide_givenS}).
#' }
#'
#' @param X A n x p concatenated data matrix of views X_1,...,X_d.
#' @param pvec A vector of values p_1,....,p_d corresponding to the number of measurements within each  view.
#' @param n_lambda A length of tuning parameter sequence used for generation of candidate structures. The default value is 50.
#' @param lambda_min A minimal value for tuning parameter. The default value is 0.01.
#' @param n_fold A number of folds for rows in BCV procedure, the default value is 3.
#' @param p_fold A number of folds for columns in BCV procedure, the default value is 3.
#' @param center A logical indicator of whether the data should be centered, the default value is TRUE.
#' @param k_max A maximal number of allowed iterations for SLIDE model fitting, the default value is 1000.
#' @param eps A convergence tolerance criterion, the default value is 1e-6.
#'
#' @return A list with the elements
#'   \item{out_s}{Output from \code{standardizedX} applied to X.}
#'   \item{structure_list}{A list of considered candidate structures.}
#'   \item{S}{A selected binary structure for the SLIDE model.}
#'   \item{model}{Slide model fitted on selected \code{S}, output from \code{slide_givenS}.}
#' @export
#'
#' @examples
#' n = 100
#'
#' # Two matched datasets
#' p1 = 25
#' p2 = 25
#' data = generateModel1(n = n, pvec = c(p1, p2))
#' out_slide = slide(X = data$X, pvec = c(p1,p2))
#' out_slide$S
#'
#' # Three matched datasets
#' p3 = 40
#' data = generateModel2(n = n, pvec = c(p1, p2, p3))
#' out_slide = slide(X = data$X, pvec = c(p1, p2, p3))
#' out_slide$S
slide <- function(X, pvec, n_lambda = 50, lambda_min = 0.01, n_fold = 3, p_fold = 3, center = T, k_max = 5000, eps = 1e-06) {
    # Center and scale the original dataset
    out_s <- standardizeX(X, pvec, center = center)

    # Form the list of candidate structures
    out_struct <- create_structure_list(X = out_s$X, pvec = pvec, lambda_max = max(out_s$svec), standardized = T, lambda_min = lambda_min, n_lambda = n_lambda)

    # Select the structure from the list using bcv procedure
    out_bcv <- slide_BCV(X = out_s$X, pvec = pvec, structure_list = out_struct$Slist, n_fold = n_fold,
        p_fold = p_fold, k_max = k_max, eps = eps, center = center)

    # Fit slide model on a selected structure S
    out_slide <- slide_givenS(X = out_s$X, pvec = pvec, S = out_bcv$structure_min, k_max = k_max,
        eps = eps, standardized = T)

    return(list(out_s = out_s, structure_list = out_struct$Slist, S = out_bcv$structure_min, model = out_slide))
}
