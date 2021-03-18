#' Percent of variance explained by each component of fitted SLIDE model
#'
#' @param slide Output from \code{slide} function
#' @inheritParams slide_givenS
#'
#' @return A list with the elements
#'   \item{var_total}{A vector of length r (number of components in fitted SLIDE model), each element is a percentage variance explained in the total concatenated dataset by that component}
#'   \item{var_vew}{A d x r matrix, each element is a percentage variance explained in the respective view by that component.}
#' @export
#'
#' @examples
#' n = 25
#' # Two matched datasets
#' p1 = p2 = 10
#' data = generateModel1(n = n, pvec = c(p1, p2))
#' out_slide = slide(X = data$X, pvec = c(p1, p2))
#' var_explained(out_slide, pvec = c(p1, p2))
#' 
var_explained <- function(slide, pvec){
  # Total number of views
  d <- length(pvec)
  pcum <- c(0, cumsum(pvec))
  # how many total components did the model select?
  n_components <- ncol(slide$S)
  # how much variance is explained by each component: (a) on the whole dataset; (b) on each view individually
  var_total = rep(NA, n_components)
  var_view = matrix(NA, d, n_components)
  for (j in 1:n_components){
    # Standardized total dataset always has frobenius norm d
    var_total[j] = sum((slide$out_s$X - tcrossprod(slide$model$U[, j, drop = F], slide$model$V[, j, drop = F]))^2)/d
    for (l in 1:d){
      index <- c((pcum[l] + 1):pcum[l + 1])
      # Each standardized individual dataset has frobenius norm 1
      var_view[l, j] = sum((slide$out_s$X[, index] - tcrossprod(slide$model$U[, j, drop = F], slide$model$V[index, j, drop = F]))^2)
    }
  }
  
  return(list(var_total = round((1-var_total)*100, 2), var_view = round((1-var_view) * 100, 2)))
}

#' Calculate scores U on a new dataset based on supplied V from SLIDE model
#'
#' @param V A n x p loadings matrix from fitted SLIDE model, part of the output of \code{slide} and \code{slide_givenS}
#' @inheritParams slide_givenS
#' @param Xtest - A ntest x p concatenated matrix of d views for which the scores U are desired
#'
#' @return A ntest x r matrix of orthogonal scores for Xtest, where r is the number of components (columns) in V
#' @export
#'
#' @examples
#' n = 25; p1 = p2 = 25
#' data = generateModel1(n = n, pvec = c(p1, p2))
#' # Specify binary structure
#' S = matrix(c(1,1,0,1,0,1),nrow = 2, ncol = 3)
#' fit_slide = slide_givenS(data$X, pvec = c(p1,p2), S = S)
#' 
#' # Generate new testing data
#' ntest = 30
#' data_test = generateModel1(n = ntest, pvec = c(p1, p2))
#' # Calculate scores on test data
#' U_test = scores_test(fit_slide$V, pvec = c(p1, p2), Xtest = data_test$X)
#' 
scores_test <- function(V, pvec, Xtest){
  # Standardize Xtest 
  test_out <- standardizeX(Xtest, pvec)
  
  # Calculate scores U based on supplied V (modification of Algorithm 3)
  XV_svd <- svd(test_out$X %*% V)
  U <- tcrossprod(XV_svd$u, XV_svd$v)
  
  return(U)
}