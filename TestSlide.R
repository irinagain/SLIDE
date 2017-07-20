# Test

# generate the model
data <- generateModel1(n = 100, pvec = c(25,25), cvec = c(1,1), snr = 1, orthogonalV = T)

X = data$X
pvec = data$pvec
center = T
eps = 1e-6
k_max = 1000
lambda_min = 0.05
nl = 100
n_fold = 3
p_fold = 3

#  # Scale the original dataset
out_s <- standardizeX(X, pvec, center = center)
X <- out_s$X
svec <- out_s$svec

# Solve penalized matrix optimizaiton problem for the sequence of lambda values
out <- solve_optim1_seq(X = X, pvec = pvec, n_lambda = nl, lambda_max = max(svec), lambda_min = lambda_min, eps = eps, k_max = k_max)

# Form the list of candidate structures
out_struct <- create_structure_list(out, pvec)

# Select the structure from the list using bcv procedure
outbcv <- slide_BCV(X, pvec = data$pvec, structure_list = out_struct$Slist, n_fold = n_fold, p_fold=p_fold, k_max = k_max, eps = eps, center = center)
# Fit slide model on a selected structure S
out_slide <- slide_givenS(X = X, pvec = pvec, S = outbcv$structure_min, k_max = k_max, eps = eps)
