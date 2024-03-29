% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SLIDE_all.R
\name{slide}
\alias{slide}
\title{Fits SLIDE model to the multi-view data X}
\usage{
slide(
  X,
  pvec,
  n_lambda = 50,
  lambda_min = 0.01,
  n_fold = 3,
  p_fold = 3,
  center = T,
  k_max = 5000,
  eps = 1e-06,
  ratio_max = NULL
)
}
\arguments{
\item{X}{A n x p concatenated data matrix of views X_1,...,X_d.}

\item{pvec}{A vector of values p_1,....,p_d corresponding to the number of measurements within each  view.}

\item{n_lambda}{A length of tuning parameter sequence used for generation of candidate structures. The default value is 50.}

\item{lambda_min}{A minimal value for tuning parameter. The default value is 0.01.}

\item{n_fold}{A number of folds for rows in BCV procedure, the default value is 3.}

\item{p_fold}{A number of folds for columns in BCV procedure, the default value is 3.}

\item{center}{A logical indicator of whether the data should be centered, the default value is TRUE.}

\item{k_max}{A maximal number of allowed iterations for SLIDE model fitting, the default value is 1000.}

\item{eps}{A convergence tolerance criterion, the default value is 1e-6.}

\item{ratio_max}{A maximal allowable rank of the binary structure as a ratio of maximal rank of X, the default value is set based on bi-cross-validation folds; with 3 folds corresponding to \code{ratio_max = 2/3} and 2 folds corresponding to \code{ratio_max = 1/2}.}
}
\value{
A list with the elements
  \item{out_s}{Output from \code{standardizedX} applied to X.}
  \item{structure_list}{A list of considered candidate structures.}
  \item{S}{A selected binary structure for the SLIDE model.}
  \item{model}{Slide model fitted on selected \code{S}, output from \code{slide_givenS}.}
}
\description{
Fits SLIDE model to the multi-vew data X
}
\details{
The function automatically column-centers and standardizes X so that the Frobenius norm within each view is equal to one (using function \code{standardizeX}). The function then performs 3 steps of SLIDE workflow:
\enumerate{
\item Create a candidate list of distinct binary structures S_1,...,S_m (using function \code{create_structure_list}).
\item Select one structure S from the list using the bi-cross-validation (using function \code{slideBCV}).
\item Fit SLIDE model with selected structure S (using function \code{slide_givenS}).
}
}
\examples{
n = 25

# Two matched datasets
p1 = 10
p2 = 10
data = generateModel1(n = n, pvec = c(p1, p2))
out_slide = slide(X = data$X, pvec = c(p1,p2))
out_slide$S

# Three matched datasets
p3 = 15
data = generateModel2(n = n, pvec = c(p1, p2, p3))
out_slide = slide(X = data$X, pvec = c(p1, p2, p3))
out_slide$S
}
