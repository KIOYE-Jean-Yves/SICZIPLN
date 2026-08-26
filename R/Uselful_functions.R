#' @importFrom MASS ginv
#' @importFrom parallel makeCluster stopCluster detectCores parLapply clusterExport clusterEvalQ
#' @importFrom stats model.matrix optim
#' @importFrom PLNmodels prepare_data ZIPLN
#' Compute the product x by logarithm x
#'
#' This function computes the product  x by logarithm x.
#'
#' @param x A scalar, vector, or matrix. If `x` is a vector or matrix,
#'   the product is computed element-wise.
#'
#' @returns
#' A scalar if `x` is a scalar, a vector if `x` is a vector,
#' and a matrix if `x` is a matrix. The returned object has the
#' same dimensions as `x`.
#'
#' @encoding UTF-8
#'
#' @examples
#' data <- c(10, 20, -99)
#' xlogx(data)
#'
#' @export
xlogx <- function(x){
  ifelse(x < .Machine$double.eps, 0, x*log(x))
}
#########################################################
#' Compute a numerically stable logarithm
#'
#' Computes the natural logarithm after truncating the input values
#' to the interval [`eps`, `1 - eps`]. This prevents numerical issues
#' caused by values equal or very close to 0 or 1.
#'
#' @param tens A numeric scalar, vector, matrix, or array.
#' @param eps A small positive constant used to bound the values before
#'   applying the logarithm. Default is `1e-16`.
#'
#' @returns
#' An object with the same dimensions as `tens`, containing the
#' logarithm of the truncated values.
#'
#' @details
#' Each element of `tens` is first truncated to the interval
#' [`eps`, `1 - eps`] using:
#' \deqn{\max(\min(x, 1-\varepsilon), \varepsilon)}
#' and the natural logarithm is then applied.
#'
#' This function is useful when the probabilities are close to 0 or 1 may
#' lead to numerical instability log(0)= -Inf
#'
#' @examples
#' trunc_log(c(0, 1e-20, 0.2, 0.8, 1))
#'
#' @export
trunc_log <- function(tens, eps = 1e-16) {
  # Limiter les valeurs entre eps et 1 - eps
  integer <- pmin(pmax(tens, eps), 1 - eps)

  # Appliquer le logarithme
  return(log(integer))
}
#################################################################
#' Approximate the logarithm of a factorial
#'
#' Computes an approximation of \eqn{\log(n!)} using Ramanujan's formula.
#' This approximation is highly accurate, especially for moderate and
#' large values of `n`.
#'
#' For numerical convenience, values equal to 0 are replaced by 1,
#' since \eqn{0! = 1! = 1} and therefore
#' \eqn{\log(0!) = \log(1!) = 0}.
#'
#' @param n A non-negative numeric vector or matrix.
#'
#' @returns
#' An object with the same dimensions as `n` containing an approximation
#' of \eqn{\log(n!)} for each element.
#'
#' @details
#' The approximation is based on Ramanujan's formula:
#' \deqn{
#' \log(n!) \approx
#' n\log(n) - n +
#' \frac{1}{6}\log\left(8n^3 + 4n^2 + n + \frac{1}{30}\right)
#' + \frac{1}{2}\log(\pi).
#' }
#'
#' @examples
#' logfactorial(5)
#' logfactorial(0:10)
#'
#' @export
## focntion de calcul du log factorial
logfactorial <- function(n) { # Ramanujan's formula
  ## Handle 0 separately since 0! = 1 and log(0!) = 0
  n[n == 0] <- 1 ## 0! = 1!
  n*log(n) - n + log(8*n^3 + 4*n^2 + n + 1/30)/6 + log(pi)/2
}
################################################################################
#' Truncated logit transformation
#'
#' Applies a logit transformation to the input values and then applies
#' the \code{trunc_log()} function to the result.
#'
#' The logit transformation is defined as:
#'
#' \deqn{\log\left(\frac{p}{1-p}\right)}
#'
#' where \eqn{p} is a probability between 0 and 1.
#'
#' @param params A numeric scalar, vector or matrix of probabilities strictly between 0 and 1.
#'
#' @return A  numeric scalar, vector or matrix of the same dimension as \code{params}, corresponding
#'   to the transformed values after applying the logit and then
#'   \code{trunc_log()}.
#'
#' @details
#' This function computes the logit transformation of each element of
#' \code{params} using \code{trunc_log()}.
#'
#' @examples
#' logit(c(0.2, 0.5, 0.8))
#'
#' @export
logit <- function(params) {
  return(trunc_log(params / (1 - params)))
}
################################################################""
#' Dirac-like transformation
#'
#' Transforms a numeric scalar vector into a two-level representation:
#' values equal to 0 are kept as 0, and all non-zero values are set to -100.
#'
#' @param params A numeric scalar, vector or matrix.
#'
#' @return A  numeric scalar, vector or matrix of the same dimesion as \code{params}, where:
#' \itemize{
#'   \item 0 values remain 0
#'   \item all non-zero values are replaced by -100
#' }
#'
#' @details
#' This function performs a simple thresholding operation. It can be interpreted
#' as a discrete indicator-like transformation.
#'
#' @examples
#' dirac(c(0, 1, 2, 0, -3))
#' # returns: 0 -100 -100 0 -100
#'
#' @export
dirac_old <- function(params) {
  params[params == 0] <- 0
  params[params != 0] <- -100
  return(params)
}
dirac <- function(params) {
  params[] <- ifelse(params == 0, 1, 0)
  params
}
####################################################################"
#' Logarithm with Zero Handling
#'
#' Computes the natural logarithm of a numeric vector, matrix, or array.
#' Values equal to \code{-Inf}, typically resulting from \code{log(0)},
#' are replaced with \code{0}.
#'
#' @param mat A numeric vector, matrix, or array.
#'
#' @return An object with the same dimensions as \code{mat}, containing
#' the natural logarithm of each element, with \code{-Inf} values replaced
#' by \code{0}.
#'
#' @details
#' This function provides a convenient way to avoid negative infinite values
#' when zero entries are present in the input data.
#'
#' @examples
#' x <- c(1, exp(1), 0)
#' log_with_zero(x)
#'
#' @export
log_with_zero <- function(mat){
  res_log <- log(mat)
  res_log[res_log == -Inf] <- 0
  return(res_log)
}
##############################################################################
#' Matrix Product with Infinite Value Handling
#'
#' Computes the matrix product of two numeric matrices using the
#' matrix multiplication operator \code{\%*\%}. Any \code{-Inf}
#' values in the resulting matrix are replaced with \code{0}.
#'
#' @param mat1 A numeric matrix.
#' @param mat2 A numeric matrix compatible with \code{mat1} for matrix multiplication.
#'
#' @return A numeric matrix corresponding to the matrix product of
#' \code{mat1} and \code{mat2}, where all \code{-Inf} values have been
#' replaced by \code{0}.
#'
#' @examples
#' mat1 <- matrix(c(1, 2, 3, 4), nrow = 2)
#' mat2 <- matrix(c(5, 6, 7, 8), nrow = 2)
#' product_dirac(mat1, mat2)
#'
#' @export
product_dirac <- function(mat1, mat2){
  prod_mat <- mat1 %*% mat2
  prod_mat[prod_mat == -Inf] <- 0
  return(prod_mat)
}
#############################################

#############################################
#' Evidence Lower Bound for a Covariate-Dependent Zero-Inflated
#' Poisson Log-Normal Model
#'
#' Computes the variational evidence lower bound (ELBO) for a multivariate
#' zero-inflated Poisson log-normal model where the zero-inflation probability
#' depends on covariates through a logistic regression model.
#'
#' @param X A numeric matrix of covariates for the latent Gaussian component.
#' @param Y A matrix of count observations.
#' @param O A matrix of offsets on the log scale.
#' @param B A matrix of regression coefficients for the latent Gaussian model.
#' @param Sigma Covariance matrix of the latent Gaussian distribution.
#' @param Pi Initial or estimated zero-inflation probabilities.
#' @param R Matrix of posterior probabilities of structural zeros.
#' @param M Matrix of variational means of the latent Gaussian variables.
#' @param S Matrix of variational standard deviations of the latent Gaussian variables.
#' @param X_zero Design matrix for the zero-inflation component.
#' @param B_zero Matrix of regression coefficients for the latent zero-inflation model.
#'
#' @return A numeric value corresponding to the variational evidence lower
#'   bound (ELBO).
#'
#' @details
#' The zero-inflation probabilities are modeled through a logistic regression:
#' \deqn{\pi_{ij} = \frac{\exp(X^{(0)}_{ij} B^{(0)})}
#'                     {1 + \exp(X^{(0)}_{ij} B^{(0)})}.}
#'
#' The ELBO combines:
#' \itemize{
#'   \item The expected Poisson log-likelihood.
#'   \item The expected log-likelihood of the zero-inflation model.
#'   \item The entropy of the latent Bernoulli indicators.
#'   \item The variational Gaussian contribution.
#' }
#'
#' @examples
#' # Example with simulated data
#' # elbo <- vloglik_col(
#' #   X, Y, O, B, Sigma, Pi, R, M, S,
#' #   X_zero, B_zero
#' # )
#'
#' @export
vloglik_col_old<-function(X,Y,O,B,Sigma,Pi,R,M,S,X_zero,B_zero){
  P<-R
  p<-ncol(Y)
  Omega<-solve(Sigma)
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  Q<-1-P
  term1<-sum(diag(t(Q) %*% (Y*(O+M) -A - logfactorial(Y))+ product_dirac(t(P),dirac(Y))))

  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  u_0<-X_zero%*%B_zero#logit(Pi)
  term2<-sum(diag(t(P)%*%u_0-t(I_n_p)%*%log_with_zero(1+exp(u_0))))
  term3<--sum(diag(t(P)%*%log_with_zero(P)+t(Q)%*%log_with_zero(Q)))
  term4<-(1/2)*sum(diag(t(I_n_p)%*%log_with_zero(S2)))+
    (nrow(Y)/2)*log(det(Omega))-
    (1/2)*sum(diag( Omega %*% (diag( c(t(rep(1,nrow(Y)))%*%S2 ) ) +t(M-X%*%B)%*%(M-X%*%B))))+
    (nrow(Y)*p)/2
  elbo<-term1+term2+term3+term4
  return(elbo)
}

vloglik_col <- function(X, Y, O, B, Sigma, Pi, R, M, S, X_zero, B_zero) {
  n <- nrow(Y)
  p <- ncol(Y)

  P  <- R
  Q  <- 1 - P
  S2 <- S^2
  OM <- O + M
  A  <- exp(OM + 0.5 * S2)

  ## Cholesky unique : inverse + log-det de Sigma en une seule factorisation
  cholSigma   <- chol(Sigma)
  Omega       <- chol2inv(cholSigma)
  logdetOmega <- -2 * sum(log(diag(cholSigma)))

  ## term1 : trace(t(Q)%*%M1) = sum(Q*M1)
  M1 <- Y * OM - A - logfactorial(Y)
  term1 <- sum(Q * M1) + sum(diag(product_dirac(t(P), dirac(Y))))
  # si product_dirac(A, B) n'est qu'un produit matriciel A %*% B, alors
  # trace(t(P)%*%dirac(Y)) = sum(P * dirac(Y)), remplacer la ligne ci-dessus par:
  # term1 <- sum(Q * M1) + sum(P * dirac(Y))

  ## term2 : u_0 déjà donné directement
  u_0 <- X_zero %*% B_zero
  term2 <- sum(P * u_0) - sum(log_with_zero(1 + exp(u_0)))

  ## term3
  term3 <- -sum(P * log_with_zero(P)) - sum(Q * log_with_zero(Q))

  ## term4
  v    <- colSums(S2)
  Res  <- M - X %*% B
  ResO <- Res %*% Omega                       # O(n*p^2), au lieu de O(p^3)
  trace_term <- sum(diag(Omega) * v) + sum(Res * ResO)

  term4 <- 0.5 * sum(log_with_zero(S2)) +
    (n / 2) * logdetOmega -
    0.5 * trace_term +
    (n * p) / 2

  term1 + term2 + term3 + term4
}
#######################################################################
#' Evidence Lower Bound for a Multivariate Zero-Inflated Poisson Log-Normal Model
#'
#' Computes the variational evidence lower bound (ELBO) for a multivariate
#' zero-inflated Poisson log-normal model with a constant zero-inflation
#' probability.
#'
#' @param X A numeric matrix of covariates for the latent Gaussian component.
#' @param Y A matrix of count observations.
#' @param O A matrix of offsets on the log scale.
#' @param B A matrix of regression coefficients for the latent Gaussian model.
#' @param Sigma Covariance matrix of the latent Gaussian distribution.
#' @param Pi Scalar zero-inflation probability.
#' @param R Matrix of posterior probabilities of structural zeros.
#' @param M Matrix of variational means of the latent Gaussian variables.
#' @param S Matrix of variational standard deviations of the latent Gaussian variables.
#'
#' @return A numeric value corresponding to the variational evidence lower
#'   bound (ELBO).
#'
#' @details
#' The ELBO is composed of four terms:
#' \itemize{
#'   \item Expected complete-data log-likelihood of the Poisson component.
#'   \item Expected log-likelihood of the zero-inflation component.
#'   \item Entropy of the latent Bernoulli variables.
#'   \item Variational contribution of the latent Gaussian variables.
#' }
#'
#' @examples
#' # Example with simulated data
#' # elbo <- vloglik(X, Y, O, B, Sigma, Pi, R, M, S)
#'
#' @export
vloglik_old<-function(X,Y,O,B,Sigma,Pi,R,M,S){
  p<-ncol(Y)
  P<-R
  Omega<-solve(Sigma)
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  Q<-1-P
  term1<-sum(diag(t(Q) %*% (Y*(O+M) -A - logfactorial(Y))+ product_dirac(t(P),dirac(Y))))

  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  u_0<-I_n_p*log(Pi/(1-Pi))#logit(Pi)
  term2<-sum(diag(t(P)%*%u_0-t(I_n_p)%*%log_with_zero(1+exp(u_0))))
  term3<--sum(diag(t(P)%*%log_with_zero(P)+t(Q)%*%log_with_zero(Q)))
  # term4<-(1/2)*sum(diag(t(I_n_p)%*%log_with_zero(S2)))+
  term4<-(1/2)*sum(diag(colSums(log_with_zero(S2))))+
    ( nrow(Y)/2)*log(det(Omega))-
    (1/2)*sum(diag( Omega %*% (diag( c(t(rep(1, nrow(Y)))%*%S2 ) ) +t(M-X%*%B)%*%(M-X%*%B))))+
    ( nrow(Y)*p)/2
  elbo<-term1+term2+term3+term4
  return(elbo)
}
vloglik <- function(X, Y, O, B, Sigma, Pi, R_ind, M, S) {
  n <- nrow(Y)
  p <- ncol(Y)

  P  <- R_ind
  Q  <- 1 - P
  S2 <- S^2
  OM <- O + M
  A  <- exp(OM + 0.5 * S2)

  ## Cholesky unique : Omega = Sigma^-1 ET log(det(Omega)), plus rapide/stable
  cholSigma   <- chol(Sigma)
  Omega       <- chol2inv(cholSigma)
  logdetOmega <- -2 * sum(log(diag(cholSigma)))

  ## term1 : trace(t(Q)%*%M1) = sum(Q*M1)
  M1 <- Y * OM - A - logfactorial(Y)
  term1 <- sum(Q * M1) + sum(diag(product_dirac(t(P), dirac(Y))))
  # si product_dirac(A,B) fait juste A %*% B, remplacer la ligne au-dessus par :
  # term1 <- sum(Q * M1) + sum(P * dirac(Y))

  ## term2 : trace(t(P)%*%u0)=sum(P*u0) ; trace(t(I)%*%X)=sum(X)
  u_0 <- matrix(log(Pi / (1 - Pi)), nrow = n, ncol = p, byrow = TRUE)
  term2 <- sum(P * u_0) - sum(log_with_zero(1 + exp(u_0)))

  ## term3 : trace(t(P)%*%f(P)) = sum(P*f(P))
  term3 <- -sum(P * log_with_zero(P)) - sum(Q * log_with_zero(Q))

  ## term4
  v    <- colSums(S2)          # remplace t(rep(1,n)) %*% S2
  Res  <- M - X %*% B          # calculé une seule fois
  ResO <- Res %*% Omega        # O(n*p^2)
  trace_term <- sum(diag(Omega) * v) + sum(Res * ResO)  # évite le O(p^3)

  term4 <- 0.5 * sum(log_with_zero(S2)) +
    (n / 2) * logdetOmega -
    0.5 * trace_term +
    (n * p) / 2

  term1 + term2 + term3 + term4
}
########################################################""
#' Objective Function for the Variational E-Step
#'
#' Computes the negative penalized evidence lower bound (ELBO) used during
#' the variational E-step of a multivariate zero-inflated Poisson log-normal
#' model. The function is intended to be minimized with respect to the
#' variational parameters.
#'
#' @param X Numeric matrix of covariates for the latent Gaussian component
#'   (\eqn{n \times d}).
#' @param Y Matrix of observed counts (\eqn{n \times p}).
#' @param O Matrix of offsets on the log scale (\eqn{n \times p}).
#' @param params Vector containing the variational parameters stacked in the
#'   following order: means (\code{M}), standard deviations (\code{S}), and
#'   posterior probabilities of structural zeros (\code{R}).
#' @param B Matrix of regression coefficients for the latent Gaussian model.
#' @param Omega Precision matrix of the latent Gaussian distribution
#'   (\eqn{\Omega = \Sigma^{-1}}).
#' @param X_zero Design matrix for the zero-inflation component.
#' @param B_zero Matrix of regression coefficients for the zero-inflation model.
#' @param lambda Non-negative regularization parameter controlling the SIC
#'   penalty.
#' @param epsilon Small positive constant used in the SIC to
#'   avoid numerical instability.
#'
#' @return A numeric value corresponding to the negative penalized ELBO.
#'   The value is negated so that standard minimization algorithms can be used.
#'
#' @details
#' The ELBO consists of:
#' \itemize{
#'   \item the expected log-likelihood of the Poisson component,
#'   \item the expected log-likelihood of the zero-inflation component,
#'   \item the entropy of the latent Bernoulli indicators,
#'   \item the variational Gaussian contribution.
#' }
#'
#' A smooth SIC-type penalty is added on the regression coefficients:
#' \deqn{
#' \frac{\lambda}{2}
#' \sum_j \frac{B_j^2}{B_j^2 + \epsilon^2},
#' }
#' promoting sparse solutions while maintaining differentiability.
#'
#' The function returns the negative penalized ELBO because most optimization
#' routines in R perform minimization rather than maximization.
#'
#' @examples
#' # Negative ELBO at current variational parameters
#' # objective_E_step(
#' #   X = X,
#' #   Y = Y,
#' #   O = O,
#' #   params = params,
#' #   B = B,
#' #   Omega = Omega,
#' #   X_zero = X_zero,
#' #   B_zero = B_zero,
#' #   lambda = 1,
#' #   epsilon = 1e-4
#' # )
#'
#' @export
# R,M,S=params
objective_E_step_a_supprimer<-function(X,Y,O=O,params,B,Omega,X_zero,B_zero,lambda, epsilon){
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)
  X_zero<-as.matrix(X_zero)
  # Recuperation des paramètres variationnels dan params
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p*d)+p)#nb_params<-log(n)*((((p+1)*p)/2)+(p*d))/2
  # params<-c(M,S,R)
  M <- matrix(params[1:(n*p)],n,p)  # (n,p)
  S <-matrix(params[(n*p+1):((n*p)+(n*p))],n,p)
  R<-matrix(params[((n*p)+(n*p)+1):length(params)],n,p)
  P<-R
  Sigma<-((Omega))
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  Q<-1-P

  Pi=(1/(1+exp(-X_zero%*%B_zero)))
  term1<-sum(diag(t(Q) %*% (Y*(O+M) -A - logfactorial(Y))+ product_dirac(t(P),dirac(Y))))
  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  u_0<-X_zero%*%B_zero#logit(Pi)
  term2<-sum(diag(t(P)%*%u_0-t(I_n_p)%*%log_with_zero(1+exp(u_0))))
  term3<--sum(diag(t(P)%*%log_with_zero(P)+t(Q)%*%log_with_zero(Q)))
  term4<-(1/2)*sum(diag(t(I_n_p)%*%log_with_zero(S2)))+
    (nrow(Y)/2)*log(det(Omega))-
    (1/2)*sum(diag( Omega %*% (diag( c(t(rep(1,nrow(Y)))%*%S2 ) ) +t(M-X%*%B)%*%(M-X%*%B))))+
    (nrow(Y)*p)/2
  elbo<-term1+term2+term3+term4
  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  # SIC_penalty[1,]<-0
  elbo<-elbo-(sum(SIC_penalty)-nb_params)
  return(-elbo)
}
objective_E_step <- function(X, Y, O, params, R, B, Omega, X_zero, B_zero, lambda, epsilon) {
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  M <- matrix(params[1:(n * p)], n, p)
  S <- matrix(params[(n * p + 1):(2 * n * p)], n, p)
  P <- R
  S2 <- S^2
  A <- exp(O + M + 0.5 * S2)
  Q <- 1 - P

  term1 <- sum(Q * (Y * (O + M) - A - logfactorial(Y))) + sum(P * dirac(Y))

  u_0 <- X_zero %*% B_zero
  term2 <- sum(P * u_0) - sum(log_with_zero(1 + exp(u_0)))
  term3 <- -sum(P * log_with_zero(P) + Q * log_with_zero(Q))

  Sbar <- colSums(S2)
  resid <- M - X %*% B
  term4 <- 0.5 * sum(log_with_zero(S2)) +
    (n / 2) * log(det(Omega)) -
    0.5 * (sum(diag(Omega) * Sbar) + sum((resid %*% Omega) * resid)) +
    (n * p) / 2
  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p)+d*p+d)
  elbo <- term1 + term2 + term3 + term4-sum(SIC_penalty)-nb_params
  # nb_params <- (lambda / 2) * ((((p + 1) * p) / 2) + (p * d) + p)
  # SIC_penalty <- (lambda / 2) * (B^2 / (B^2 + epsilon^2))
  # -(elbo - (sum(SIC_penalty) - nb_params))
  -(elbo)
}
################################################################
#' Gradient of the Variational E-Step Objective Function
#'
#' Computes the gradient of the variational E-step objective (negative ELBO)
#' for a multivariate zero-inflated Poisson log-normal model.
#' This gradient is used for optimization of the variational parameters
#' \code{M}, \code{S}, and \code{R}.
#'
#' @param X Numeric matrix of covariates for the latent Gaussian component
#'   (\eqn{n \times d}).
#' @param Y Matrix of observed counts (\eqn{n \times p}).
#' @param O Matrix of offsets on the log scale (\eqn{n \times p}).
#' @param params Numeric vector containing stacked variational parameters in
#'   the order: means (\code{M}), standard deviations (\code{S}),
#'   and posterior probabilities of structural zeros (\code{R}).
#' @param B Matrix of regression coefficients for the latent Gaussian model.
#' @param Omega Precision matrix of the latent Gaussian distribution
#'   (\eqn{\Omega = \Sigma^{-1}}).
#' @param X_zero Design matrix for the zero-inflation component.
#' @param B_zero Matrix of regression coefficients for the zero-inflation model.
#' @param lambda Regularization parameter used in the penalized objective
#'   (not directly used in the gradient expression but part of the full model).
#' @param epsilon Small positive constant used in the SIC penalty approximation.
#'
#' @return A numeric vector containing the gradient of the objective function
#'   with respect to the stacked variational parameters
#'   (\code{M}, \code{S}, \code{R}).
#'
#' @details
#' The gradient is derived from the variational lower bound of a
#' zero-inflated Poisson log-normal model. It consists of three blocks:
#' \itemize{
#'   \item gradient with respect to \code{M} (variational means),
#'   \item gradient with respect to \code{S} (variational standard deviations),
#'   \item gradient with respect to \code{R} (posterior zero-inflation probabilities).
#' }
#'
#' These quantities are stacked into a single vector to be used by
#' gradient-based optimization algorithms (e.g., \code{optim}).
#'
#' @note
#' The function assumes that \code{params} is correctly ordered and has
#' length equal to \eqn{2np + np}, where \eqn{n} is the number of samples
#' and \eqn{p} the number of responses.
#'
#' @examples
#' # Gradient evaluation at current parameters
#' # g <- grad_E_step(
#' #   X, Y, O, params, B, Omega,
#' #   X_zero, B_zero, lambda = 1, epsilon = 1e-4
#' # )
#'
#' @export
grad_E_step_a_supprimer <-function(X,Y,O=O,params,B,Omega,X_zero,B_zero,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<- ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)

  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O) # offsets (n,p)
  X_zero<-as.matrix(X_zero)
  # Recuperation des paramètres variationnels dan params

  # params<-c(M,S,R)
  M <- matrix(params[1:(nrow(Y)*p)],nrow(Y),p)  # (n,p)
  S <-matrix(params[(nrow(Y)*p+1):((nrow(Y)*p)+(nrow(Y)*p))],nrow(Y),p)
  R<-matrix(params[((nrow(Y)*p)+(nrow(Y)*p)+1):length(params)],nrow(Y),p)
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  P<-R#(1/(1+exp(A+X_zero%*%B_zero)))*dirac(Y)#

  # B<-matrix(params[1:(d*p)],d,p)
  # matrix(params_init_E_step[(n*p+1):length(params_init_E_step)],n,p)
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  B<-B
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  Omega <-Omega
  Sigma<-((Omega))
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)
  # S2 <- S^2
  Pi=(1/(1+exp(-X_zero%*%B_zero)))

  # Omega<-solve(Sigma)
  # nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty
  I_n_p<-matrix(1,nrow =nrow(Y),ncol = ncol(Y))
  u_0<-X_zero%*%B_zero#logit(Pi)
  grad_M <- ((I_n_p-P)*(Y-A)-(M-X%*%B)%*%Omega)

  grad_S <- ((1/S)-(I_n_p-P)*S*A-S*rep(1,nrow(Y))%*% t(diag(Omega)))
  grad_R <- (P*(A+u_0-logit(P))-log_with_zero(1-P))

  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  gradient <- c(as.vector(grad_M), as.vector(grad_S),as.vector(grad_R))
  return(gradient)
}

grad_E_step <-function(X,Y,O=O,params,R,B,Omega,X_zero,B_zero,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<- ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)

  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O) # offsets (n,p)
  X_zero<-as.matrix(X_zero)
  # Recuperation des paramètres variationnels dan params

  # params<-c(M,S,R)
  M <- matrix(params[1:(nrow(Y)*p)],nrow(Y),p)  # (n,p)
  S <-matrix(params[(nrow(Y)*p+1):((nrow(Y)*p)+(nrow(Y)*p))],nrow(Y),p)
  # R<-matrix(params[((nrow(Y)*p)+(nrow(Y)*p)+1):length(params)],nrow(Y),p)
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  P<-R#(1/(1+exp(A+X_zero%*%B_zero)))*dirac(Y)#

  # B<-matrix(params[1:(d*p)],d,p)
  # matrix(params_init_E_step[(n*p+1):length(params_init_E_step)],n,p)
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  B<-B
  # Omega<-matrix(params[((d*p)+1):((d*p)+(p*p))],p,p)
  Omega <-Omega
  # Sigma<-((Omega))
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  # M <- matrix(prams_vec_E_step[1:(n * p)], n, p)
  # S <- matrix(prams_vec_E_step[(n * p + 1):(n * p + n * p)], n, p)
  # S2 <- S^2
  Pi=(1/(1+exp(-X_zero%*%B_zero)))

  # Omega<-solve(Sigma)
  # nSigma <- t(M) %*% (M * w) + diag(w %*% S2)
  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty
  I_n_p<-matrix(1,nrow =nrow(Y),ncol = ncol(Y))
  u_0<-X_zero%*%B_zero#logit(Pi)
  grad_M <- ((I_n_p-P)*(Y-A)-(M-X%*%B)%*%Omega)

  grad_S <-  ((1/S)-(I_n_p-P)*S*A-S*rep(1,nrow(Y))%*% t(diag(Omega)))

  # grad_R <- (P*(A+u_0-logit(P))-log_with_zero(1-P))

  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  gradient <-  c(as.vector(grad_M), as.vector(grad_S))#c(as.vector(grad_M), as.vector(grad_S),as.vector(grad_R))
  return(-gradient)
}
#######################################################################
#' Objective Function for the Variational M-Step
#'
#' Computes the negative penalized evidence lower bound (ELBO) used in the
#' variational M-step of a multivariate zero-inflated Poisson log-normal model.
#' The function is optimized with respect to the model parameters
#' \code{B} and \code{B_zero}, while variational parameters are held fixed.
#'
#' @param X Numeric matrix of covariates for the latent Gaussian component
#'   (\eqn{n \times d}).
#' @param Y Matrix of observed count data (\eqn{n \times p}).
#' @param O Matrix of offsets on the log scale (\eqn{n \times p}).
#' @param paramsM Numeric vector containing model parameters stacked as:
#'   regression coefficients for the latent Gaussian component (\code{B})
#'   followed by coefficients for the zero-inflation model (\code{B_zero}).
#' @param Omega Precision matrix of the latent Gaussian distribution
#'   (\eqn{\Omega = \Sigma^{-1}}).
#' @param M Matrix of variational means (\eqn{n \times p}).
#' @param S Matrix of variational standard deviations (\eqn{n \times p}).
#' @param R Matrix of posterior probabilities of structural zeros (\eqn{n \times p}).
#' @param X_zero Design matrix for the zero-inflation component.
#' @param lambda Regularization parameter controlling the SIC-type penalty.
#' @param epsilon Small positive constant used for numerical stabilization
#'   in the penalty term.
#'
#' @return A numeric value corresponding to the negative penalized ELBO.
#'   The function is designed for minimization in numerical optimization
#'   procedures.
#'
#' @details
#' The objective function is composed of:
#' \itemize{
#'   \item Expected Poisson log-likelihood under the variational distribution,
#'   \item Zero-inflation logistic likelihood,
#'   \item Entropy of latent Bernoulli variables,
#'   \item Gaussian variational contribution,
#'   \item SIC-type sparsity penalty on regression coefficients.
#' }
#'
#' The penalty term encourages sparsity in \code{B} while remaining smooth
#' for gradient-based optimization.
#'
#' @note
#' The function returns the negative ELBO to be compatible with standard
#' minimization routines such as \code{optim}.
#'
#' @examples
#' # Evaluate objective at current parameters
#' # obj <- objective_M_step(
#' #   X, Y, O,
#' #   paramsM,
#' #   Omega,
#' #   M, S, R,
#' #   X_zero,
#' #   lambda = 1,
#' #   epsilon = 1e-4
#' # )
#'
#' @export
objective_M_step_a_supprimer<-function(X,Y,O=O,paramsM,Omega,M,S,R,X_zero,lambda, epsilon){
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]
  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)
  X_zero<-as.matrix(X_zero)
  # Paramètres variationnels
  M<-M
  S<-S
  R<-R
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p*d)+p)#log(n)*(((p+1)*p)/2)+(p*d)
  # Recuperation des paramètres du modèle dans params
  # params<-c(M,S,R)
  B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)
  # Omega<-matrix(paramsM[(d*p+1):((d*p)+(p*p))],p,p)
  B_zero<-matrix(paramsM[((d*p)+1):length(paramsM)],d,p)
  P<-R
  Omega<-Omega
  Sigma<-((Omega))
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  Pi=(1/(1+exp(-X_zero%*%B_zero)))
  Q<-1-P
  term1<-sum(diag(t(Q) %*% (Y*(O+M) -A - logfactorial(Y))+ product_dirac(t(P),dirac(Y))))

  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  u_0<-X_zero%*%B_zero#I_n_p*log_with_zero(Pi/(1-Pi))#logit(Pi)
  term2<-sum(diag(t(P)%*%u_0-t(I_n_p)%*%log_with_zero(1+exp(u_0))))
  term3<--sum(diag(t(P)%*%log_with_zero(P)+t(Q)%*%log_with_zero(Q)))
  term4<-(1/2)*sum(diag(t(I_n_p)%*%log_with_zero(S2)))+
    (nrow(Y)/2)*log(det(Omega))-
    (1/2)*sum(diag( Omega %*% (diag( c(t(rep(1,nrow(Y)))%*%S2 ) ) +t(M-X%*%B)%*%(M-X%*%B))))+
    (nrow(Y)*p)/2
  elbo<-term1+term2+term3+term4
  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  elbo<-elbo-(sum(SIC_penalty)-nb_params)
  return(-elbo)
}

objective_M_step <- function(X, Y, O, paramsM, Omega,B, M, S, R, X_zero, lambda, epsilon) {
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  B_zero <- matrix(paramsM, d, p)
  # B_zero<- matrix(paramsM[1:(d*p)],d,p)

  B<-B#matrix(paramsM[((d*p)+1):length(paramsM)],d,p)
  P <- R
  Q <- 1 - P
  S2 <- S^2
  A <- exp(O + M + 0.5 * S2)

  term1 <- sum(Q * (Y * (O + M) - A - logfactorial(Y))) + sum(P * dirac(Y))

  u_0 <- X_zero %*% B_zero
  term2 <- sum(P * u_0) - sum(log_with_zero(1 + exp(u_0)))
  term3 <- -sum(P * log_with_zero(P) + Q * log_with_zero(Q))

  Sbar <- colSums(S2)
  resid <- M - X %*% B
  term4 <- 0.5 * sum(log_with_zero(S2)) +
    (n / 2) * log(det(Omega)) -
    0.5 * (sum(diag(Omega) * Sbar) + sum((resid %*% Omega) * resid)) +
    (n * p) / 2

  elbo <- term1 + term2 + term3 + term4
  # nb_params <- (lambda / 2) * ((((p + 1) * p) / 2) + (p * d) + p)
  # SIC_penalty <- (lambda / 2) * (B^2 / (B^2 + epsilon^2))
  # SIC_penalty[1, ] <- 0
  # -(elbo - (sum(SIC_penalty) - nb_params))
  -(elbo)
}
objective_M_step_with_grad_B <- function(X, Y, O, paramsM, Omega, M, S, R, X_zero, lambda, epsilon) {
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  B_zero <-matrix(paramsM[((d*p)+1):length(paramsM)],d,p)# matrix(paramsM, d, p)
  # B_zero<- matrix(paramsM[1:(d*p)],d,p)

  B<-matrix(paramsM[1:(d*p)],d,p) #B
  P <- R
  Q <- 1 - P
  S2 <- S^2
  A <- exp(O + M + 0.5 * S2)

  term1 <- sum(Q * (Y * (O + M) - A - logfactorial(Y))) + sum(P * dirac(Y))

  u_0 <- X_zero %*% B_zero
  term2 <- sum(P * u_0) - sum(log_with_zero(1 + exp(u_0)))
  term3 <- -sum(P * log_with_zero(P) + Q * log_with_zero(Q))

  Sbar <- colSums(S2)
  resid <- M - X %*% B
  term4 <- 0.5 * sum(log_with_zero(S2)) +
    (n / 2) * log(det(Omega)) -
    0.5 * (sum(diag(Omega) * Sbar) + sum((resid %*% Omega) * resid)) +
    (n * p) / 2
  SIC_penalty<-(lambda/2) *(B^2 / (B^2 + epsilon^2))
  SIC_penalty[1,]<-0
  nb_params<-(lambda/2)*((((p+1)*p)/2)+(p)+d*p+d)
  elbo <- term1 + term2 + term3 + term4-sum(SIC_penalty)-nb_params
  # nb_params <- (lambda / 2) * ((((p + 1) * p) / 2) + (p * d) + p)
  # SIC_penalty <- (lambda / 2) * (B^2 / (B^2 + epsilon^2))
  # SIC_penalty[1, ] <- 0
  # -(elbo - (sum(SIC_penalty) - nb_params))
  (-elbo)
}
###############################################################"
#' Gradient of the Variational M-Step Objective Function
#'
#' Computes the gradient of the negative penalized ELBO with respect to the
#' model parameters in the M-step of a multivariate zero-inflated Poisson
#' log-normal model.
#'
#' The gradient is evaluated with respect to the regression coefficients
#' \code{B} (latent Gaussian component) and \code{B_zero} (zero-inflation
#' component), while variational parameters are held fixed.
#'
#' @param X Numeric matrix of covariates for the latent Gaussian component
#'   (\eqn{n \times d}).
#' @param Y Matrix of observed count data (\eqn{n \times p}).
#' @param O Matrix of offsets on the log scale (\eqn{n \times p}).
#' @param paramsM Numeric vector containing stacked model parameters:
#'   regression coefficients for the latent Gaussian model (\code{B}) followed
#'   by coefficients for the zero-inflation model (\code{B_zero}).
#' @param Omega Precision matrix of the latent Gaussian distribution
#'   (\eqn{\Omega = \Sigma^{-1}}).
#' @param M Matrix of variational means (\eqn{n \times p}).
#' @param S Matrix of variational standard deviations (\eqn{n \times p}).
#' @param R Matrix of posterior probabilities of structural zeros (\eqn{n \times p}).
#' @param X_zero Design matrix for the zero-inflation component.
#' @param lambda Regularization parameter controlling the SIC-type penalty.
#' @param epsilon Small positive constant used for numerical stabilization
#'   in the penalty term.
#'
#' @return A numeric vector containing the concatenated gradients with respect
#'   to \code{B} and \code{B_zero}.
#'
#' @details
#' The gradient consists of two blocks:
#' \itemize{
#'   \item Gradient with respect to \code{B}, including the Gaussian likelihood
#'   contribution and the SIC-type sparsity penalty,
#'   \item Gradient with respect to \code{B_zero}, corresponding to the
#'   logistic zero-inflation model.
#' }
#'
#' The SIC penalty is a smooth sparsity-inducing penalty defined as:
#' \deqn{
#' \frac{\lambda}{2} \frac{2 B \epsilon^2}{(B^2 + \epsilon^2)^2}.
#' }
#'
#' The function returns a vectorized gradient suitable for optimization
#' routines such as \code{optim}.
#'
#' @note
#' The gradient assumes that variational parameters (\code{M}, \code{S}, \code{R})
#' are fixed during the M-step optimization.
#'
#' @examples
#' # Compute gradient at current M-step parameters
#' # g <- grad_M_step(
#' #   X, Y, O,
#' #   paramsM,
#' #   Omega,
#' #   M, S, R,
#' #   X_zero,
#' #   lambda = 1,
#' #   epsilon = 1e-4
#' # )
#'
#' @export
grad_M_step_a_supprimer <- function(X,Y,O=O,paramsM,Omega,M,S,R,X_zero,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]

  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)

  # Paramètres variationnels
  M<-M
  S<-S
  R<-R

  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  # Recuperation des paramètres du modèle dans params
  # params<-c(M,S,R)
  B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)
  # Omega<-matrix(paramsM[(d*p+1):((d*p)+(p*p))],p,p)
  Omega<-Omega
  B_zero<-matrix(paramsM[((d*p)+1):length(paramsM)],d,p)
  Pi=(1/(1+exp(-X_zero%*%B_zero)))
  P<-R#(1/(1+exp(A+X_zero%*%B_zero)))*dirac(Y)#
  Sigma<-((Omega))
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  u_0<-X_zero%*%B_zero

  # SIC_penalty<-lambda * sum(beta^2 / (beta^2 + a^2))
  # objective <- sum(w * (A - Y * Z - 0.5 * log(S2))) + 0.5 * sum(diag(Omega %*% nSigma))+SIC_penalty
  SIC_deriv<-(lambda/2) * (2*(B*epsilon^2) / (B^2 + epsilon^2)^2)
  # print(SIC_deriv)
  SIC_deriv[1,]<-0
  grad_B <-((t(X)%*%X%*%B%*%Omega)-t(X)%*%M%*%Omega +SIC_deriv) #(t(X) %*% (w * (A - Y)) - (lambda/2) * (B^3) / (B^2 + epsilon^2)^2)
  # grad_B[1,]<-0 diag( c(t(rep(1,n))%*%S2 ) )
  S2_bar<-diag(c(t(rep(1,nrow(Y)))%*%S2))#matrix((c(t(rep(1,n))%*%S2)),p,p,byrow = TRUE)
  grad_Omega <-((nrow(Y)/2)*(Omega)-(1/2)*(t(M-X%*%B)%*%(M-X%*%B)+S2_bar))
  # grad_Omega <-diag(0,p)#(-nSigma+solve(nSigma))/2
  # list(objective = objective, gradient = c(as.vector(grad_M), as.vector(grad_S)))
  # grad_R <- (P*(A-logit(P))-log_with_zero(1-P))
  grad_B_zero<-(t(X_zero)%*%P-t(X_zero)%*%(exp(u_0)/(1+exp(u_0))))
  gradient_M <- c(as.vector(grad_B),as.vector(grad_B_zero))
  return(gradient_M)
}

grad_M_step <- function(X,Y,O=O,paramsM,Omega,B,M,S,R,X_zero,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]

  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)

  # Paramètres variationnels
  M<-M
  S<-S
  R<-R

  # Les paramètres qui seront mise à jour après l'étape E de façon analytique

  B<-B
  Omega<-Omega

  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  # Recuperation des paramètres du modèle dans params
  # params<-c(M,S,R)
  # B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)

  B_zero<-matrix(paramsM,d,p)  # matrix(paramsM[((d*p)+1):length(paramsM)],d,p)
  Pi=(1/(1+exp(-X_zero%*%B_zero)))
  P<-R
  # Sigma<-((Omega))
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  u_0<-X_zero%*%B_zero


  # Grad B
  # grad_B <-((t(X)%*%X%*%B%*%Omega)-t(X)%*%M%*%Omega +SIC_deriv)
  S2_bar<-diag(c(t(rep(1,nrow(Y)))%*%S2))

  # Grad Omega
  # grad_Omega <-((nrow(Y)/2)*(Omega)-(1/2)*(t(M-X%*%B)%*%(M-X%*%B)+S2_bar))

  # Grad B_zero
  grad_B_zero<-(t(X_zero)%*%P-t(X_zero)%*%(exp(u_0)/(1+exp(u_0))))
  gradient_M <-  as.vector(grad_B_zero)

  #c(as.vector(grad_B),as.vector(grad_B_zero))
  return(-gradient_M)
}
# B n'est pas mis à jour analytiquement
grad_M_step_with_grad_B  <- function(X,Y,O=O,paramsM,Omega,M,S,R,X_zero,lambda, epsilon) {
  n <- nrow(Y)
  p <- ncol(Y)
  d<-ncol(X)
  Y <- as.matrix(Y) # réponses (n,p)
  X <- as.matrix(X) # covariables (n,d)
  # X<-X[,-2]

  O <- as.matrix(O)#matrix(0,n,p) # offsets (n,p)

  # Paramètres variationnels
  M<-M
  S<-S
  R<-R

  # Les paramètres qui seront mise à jour après l'étape E de façon analytique

  # B<-B
  Omega<-Omega

  I_n_p<-matrix(1,nrow = nrow(Y),ncol = ncol(Y))
  # Recuperation des paramètres du modèle dans params
  # paramsM<-c(B,B_zero)
   B <- matrix(paramsM[1:(d*p)],d,p)  # (n,p)

  B_zero<- matrix(paramsM[((d*p)+1):length(paramsM)],d,p) # matrix(paramsM,d,p)  #
  Pi=(1/(1+exp(-X_zero%*%B_zero)))
  P<-R
  # Sigma<-solve((Omega))
  S2<-S^2
  A<-exp(O+M+0.5*S2)
  u_0<-X_zero%*%B_zero


  # Grad B
  #  lambda<-2000
  # epsilon<-0.00001
  SIC_deriv<-(lambda/2) * (2*(B*epsilon^2) / (B^2 + epsilon^2)^2)
  SIC_deriv[1,]<-0
  # déjà gradient de -elbo ?
  grad_B <-((t(X)%*%X%*%B%*%Omega)-t(X)%*%M%*%Omega +SIC_deriv)
  S2_bar<-diag(c(t(rep(1,nrow(Y)))%*%S2))

  # Grad Omega
  # grad_Omega <-((nrow(Y)/2)*(Omega)-(1/2)*(t(M-X%*%B)%*%(M-X%*%B)+S2_bar))

  # Grad B_zero
  grad_B_zero<-(t(X_zero)%*%P-t(X_zero)%*%(exp(u_0)/(1+exp(u_0))))
  gradient_M <-  c( c(grad_B), -c(grad_B_zero))

  #c(as.vector(grad_B),as.vector(grad_B_zero))
  return(gradient_M)
}

####################################################################
# Initialisation ZIPLN
#' Helper function for ZIPLN initialization.
#'
#' @description
#' Fast LM-based starting point for ZIPLN: one multivariate `lm.fit` for the PLN
#' component and empirical zero rates / binomial GLMs for the ZI component.
#' Replaces the previous per-species `pscl::zeroinfl` loop.
#'
#' @param Y Response count matrix (n × p)
#' @param X Design matrix for the PLN component (n × d)
#' @param X0 Design matrix for the ZI component (n × d0, empty `matrix(NA,0,0)` when unused)
#' @param O Offset matrix in log-scale (n × p)
#' @param w Weight vector of length n (defaults to uniform weights)
#' @return Named list: `B` (d × p), `M` (n × p), `S2` (n × p), `R` (n × p), `B0` (d0 × p)
#'
#' @importFrom stats lm.fit glm.fit binomial PLNmodels
#' @export
compute_ZIPLN_starting_point <- function(Y, X, X0, O, w = NULL) {
  if (is.null(w)) w <- rep(1.0, nrow(Y))
  n <- nrow(Y); p <- ncol(Y); d0 <- ncol(X0)
  if (ncol(X) == 0) stop("PLN component requires at least one covariate (or an intercept) to fit ZIPLN.")

  ## PLN component: fast multivariate LM (identical to compute_PLN_starting_point "LM")
  sp <- PLNmodels::compute_PLN_starting_point(Y, X, O, w)

  ## ZI component: empirical per-species zero rates
  zero_ind <- (Y == 0) * 1.0
  R <- matrix(colMeans(zero_ind), n, p, byrow = TRUE)

  ## B0: p binomial GLMs on zero indicator vs X0 (only for "covar" ziparam where d0 > 0)
  B0 <- if (!is.null(d0) && d0 > 0) {
    binom_fam <- binomial()
    vapply(seq_len(p), function(j)
      suppressWarnings(glm.fit(X0, zero_ind[, j], family = binom_fam))$coefficients,
      numeric(d0))
  } else {
    matrix(0.0, 0L, p)
  }

  list(B = sp$B, M = sp$M, S = sp$S, R = R, B0 = B0)
}
####################################################################
##########################################################################
#' Optimization of the SIC-ZIPLN Model via Variational EM
#'
#' Performs variational EM optimization for a multivariate zero-inflated
#' Poisson log-normal (ZIPLN) model with SIC-type sparsity regularization.
#' The procedure alternates between a variational E-step and a penalized
#' M-step with a telescoping annealing scheme on \code{epsilon}.
#'
#' @param X Numeric matrix of covariates for the latent Gaussian component
#'   (\eqn{n \times d}).
#' @param Y Matrix of observed count data (\eqn{n \times p}).
#' @param X_zero Design matrix for the zero-inflation component.
#' @param offset Logical; if \code{TRUE}, offsets are included in the model.
#' @param lambda Fixed regularization parameter controlling SIC penalization.
#'   Default is \code{log(n) * ncol(X)}.
#' @param optim_method Optimization method passed to \code{optim}
#'   (e.g., \code{"BFGS"}).
#' @param max_it Maximum number of EM iterations per epsilon level.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{B} regression coefficients for the latent Gaussian model,
#'   \item \code{model_par} list of model parameters (\code{B}, \code{Omega},
#'     \code{Sigma}, \code{Pi}, \code{B_zero}),
#'   \item \code{var_par} variational parameters (\code{M}, \code{S}, \code{R}),
#'   \item \code{SICprediction} predicted latent intensities,
#'   \item \code{data} original data used in the model,
#'   \item \code{v_loglik} log-likelihood and ELBO diagnostics for different models,
#'   \item \code{BIC_SICZIPLN} Bayesian Information Criterion for SIC-ZIPLN,
#'   \item \code{BIC_SICZIPLN_approximer} BIC with sparsity approximation,
#'   \item \code{BIC_ZIPLN} BIC of the non-penalized ZIPLN model,
#'   \item \code{elbo} final ELBO value,
#'   \item \code{res_ZIPLN} fitted ZIPLN initialization model,
#'   \item \code{ICL} Integrated Completed Likelihood approximation,
#'   \item \code{ICL_approximer} alternative ICL approximation
#' }
#'
#' @details
#' The algorithm is based on a variational EM framework:
#'
#' \itemize{
#'   \item E-step: updates variational parameters (\code{M}, \code{S}, \code{R})
#'     by maximizing the ELBO using gradient-based optimization.
#'   \item M-step: updates model parameters (\code{B}, \code{B_zero}) under
#'     SIC-type sparsity penalization.
#' }
#'
#' A telescoping annealing scheme is used on \code{epsilon} to stabilize
#' optimization and gradually enforce sparsity.
#'
#' Model selection is performed using BIC and ICL criteria.
#'
#' @note
#' The function relies on an initial fit obtained from a ZIPLN model
#' (e.g., \code{ZIPLN::ZIPLN}). Convergence is assessed via changes in
#' regression coefficients.
#'
#' @examples
#' # Example (pseudo-code)
#' # fit <- SICZIPLN_optim(X, Y, X_zero)
#'
#' @export
SICZIPLN_optim<-function(X,Y,X_zero,offset=FALSE,lambda_fixed=log(n)*(ncol(X)),optim_method="BFGS",max_it=200){
  prep_pln<-prepare_data(Y,X)
  Y<-prep_pln$Abundance
  X<-prep_pln[,-c(1,length(prep_pln))]
  d<-ncol(X)
  p<-ncol(Y)
  n<-nrow(X)
  lambda<-lambda_fixed
  # if(lambda_fixed==TRUE){lambda<-log(nrow(X))*(ncol(X))}
  # if(lambda_fixed!=TRUE){lambda<-log(nrow(X))*(ncol(X))}
  if(offset==TRUE){
    O<-matrix(rep(log(prep_pln$Offset),p),ncol=p,byrow = FALSE)
  }
  if(offset==FALSE){
    O<-matrix(0,nrow=n,ncol=p)
  }
  # X<-as.matrix(X)
  # X_zero<-as.matrix(X_zero)
  X<-model.matrix(~.,as.data.frame(X))
  X_zero<-model.matrix(~.,as.data.frame(X_zero))

  # Initialisation ZIPLN
   res_PLN<-ZIPLN_init <- ZIPLN(Y~X[,-1]+offset(O[,1])|X_zero[,-1], zi = "col")
  #
  #
  B<-res_PLN$model_par$B
  Omega<-res_PLN$model_par$Omega
  Sigma<-res_PLN$model_par$Sigma
  B_zero<-res_PLN$model_par$B0
  Pi<-res_PLN$model_par$Pi
  S<-res_PLN$var_par$S
  M<-res_PLN$var_par$M
  R<-res_PLN$var_par$R
  X<-as.matrix(X)
  X_zero<-as.matrix(X_zero)
  Y<-as.matrix(Y)
  d<-ncol(X)
  n<-nrow(Y)

  #Initialisation avec starting point de ZIPLN
  # res_init_ZI<-compute_ZIPLN_starting_point(Y=Y, X=X, X0=X_zero, O=O, w = NULL)
  # B<-res_init_ZI$B
  #
  # # Omega<-res_PLN$model_par$Omega
  # # Sigma<-res_PLN$model_par$Sigma
  # B_zero<-res_init_ZI$B0
  # Pi<-res_init_ZI$R
  # S<-(res_init_ZI$S)
  # S2<-(res_init_ZI$S)^2
  # M<-res_init_ZI$M
  # R<-res_init_ZI$R
  # Omega<-n*solve((t(M-X%*%B)%*%(M-X%*%B))+diag(c(t(rep(1,nrow(Y)))%*%S2)))
  # Sigma<-solve(Omega)
  # X<-as.matrix(X)
  # X_zero<-as.matrix(X_zero)
  # Y<-as.matrix(Y)
  # d<-ncol(X)
  # n<-nrow(Y)

  params_init_E_step<-c(as.vector(M),as.vector(S),as.vector(R))
  # v <- res_PLN$model_par$B
  # #c(0.44844978, -0.08463290, -0.31860861, 0.50380695, 0.96334493, -0.18933897,
  #        0.50788630, -0.02754789, 3.14563421)
  # v<-rep(3,9)
  params_init_M_step<-c(as.vector(B),as.vector(B_zero))
  # param_initpln$B
  # B<-matrix(rnorm(d*p),d,p)
  # B_zero<-matrix(rnorm(d*p),d,p)
  # M<-matrix(1,n,p)
  # S<-matrix(0.1,n,p)
  # R<-matrix(0.1,n,p)
  # Pi<-matrix(0.2,n,p)
  #  Omega<-diag(1,p)
  # #  # En mettant d'autres pâramètres
  #  params_init_E_step<-c(as.vector(M),as.vector(S),as.vector(R))
  #  params_init_M_step<-c(as.vector(B),as.vector(B_zero))
  length_epsilon<-100
  # stock_vloglik_espsilon<-list(rep(NA,length_epsilon))

  E<-c()
  e1<-10
  E[1]<-e1
  for(t in 2:length_epsilon){
    E[t]=e1*(0.87)^(t-1)
  }
  iter<-0
  for(epsilon_val in E){
    iter<-iter+1
    v_loglik_tmp<-vloglik_col(X,Y,O,B,Sigma,Pi,R,M,S,X_zero,B_zero)
    cat("SICZIPLN: Epsilon value for telescoping :\t",iter,epsilon_val,"\t",v_loglik_tmp,"\n")
    for(i in 1:max_it){
      result_E_step <- optim(
        fn = objective_E_step,
        gr = grad_E_step,
        Y=Y,
        X=(X),
        O=O,
        par = params_init_E_step,
        B=B,
        Omega=Omega,
        X_zero=X_zero,
        B_zero=B_zero,
        lambda=lambda,
        epsilon=epsilon_val,
        method = optim_method,
        control = list(
          maxit = 300
          # reltol = 1e-5,
          # trace=FALSE,
          #fnscale= 1
          # gradtol = 1e-8
        )
      )
      # cat("E step",result_E_step$value,"\t", "indice eps",iter,"\n")
      # M <- (matrix(result_E_step$par[1:(n*p)],n,p))  # (n,p)
      # S <-  matrix(result_E_step$par[(n*p+1):length(result_E_step$par)],n,p)
      M <- matrix(result_E_step$par[1:(n*p)],n,p)  # (n,p)
      S <-matrix(result_E_step$par[(n*p+1):((n*p)+(n*p))],n,p)
      R<-matrix(result_E_step$par[((n*p)+(n*p)+1):length(result_E_step$par)],n,p)

      params_init_E_step<-c(as.vector(M),as.vector(S),as.vector(R))
      # omega_test=n*solve((t(M-X%*%B)%*%(M-X%*%B))+diag(c(t(rep(1,n))%*%S2)))
      result_M_step <- optim(
        fn = objective_M_step,
        gr = grad_M_step,
        Y=Y,
        X=(X),
        O=O,
        par = params_init_M_step,
        Omega=Omega,
        M=M,
        S=S,
        R=R,
        X_zero=X_zero,
        lambda=lambda,
        epsilon=epsilon_val,
        method = optim_method,
        control = list(
          maxit = 300
          #fnscale= 1
          # reltol = 1e-5,
          # trace=FALSE
          # gradtol = 1e-8
        )
      )
      # result_M_step$value
      # cat("M step",result_M_step$value,"\t", "indice eps",iter,"\n")
      B_old<-B
      B <- matrix(result_M_step$par[1:(d*p)],d,p)  # (n,p)
      # Omega<-matrix(result_M_step$par[(d*p+1):((d*p)+(p*p))],p,p)
      S2<-S^2
      Omega<-n*solve((t(M-X%*%B)%*%(M-X%*%B))+diag(c(t(rep(1,nrow(Y)))%*%S2)))
      Sigma<-solve(Omega)
      B_zero<-matrix(result_M_step$par[((d*p)+1):length(result_M_step$par)],d,p)
      # B <- matrix(result_M_step$par[1:(d*p)],d,p)  # (n,p)
      # Omega<- res_PLN$model_par$Omega#matrix(result_M_step$par[(d*p+1):length(result_M_step$par)],p,p)

      params_init_M_step<-c(as.vector(B),as.vector(B_zero))
      #     convergence<-abs(result_E_step$value-result_M_step$value)
      convergence<-sum(abs(B-B_old))
      if(convergence<=1e-8){
        cat("SICZIPLN:epsilon telescoping : ",epsilon_val,"Convergence avant maxit\n")
        break}
    }

  }
  # Post traitement
  B[abs(B)<=1e-5]<-0
  # Calcul de la vraissemblance en Omettant l'intercept car il n'est pas pénalisé

  # Calcul avec SICZIPLN

  v_loglik_SICZIPLN_no_intercept<-vloglik_col(X[,-1],Y,O,B[-1,],Sigma,Pi,R,M,S,X_zero,B_zero)

  v_loglik_SICZIPLN_with_intercept<-vloglik_col(X,Y,O,B,Sigma,Pi,R,M,S,X_zero,B_zero)
  #calcul avec ZIPLN
  v_loglik_ZIPLN_no_intercept<-vloglik_col(X[,-1],Y,O,res_PLN$model_par$B[-1,],res_PLN$model_par$Sigma,res_PLN$model_par$Pi,res_PLN$var_par$R,res_PLN$var_par$M,res_PLN$var_par$S,X_zero,B_zero)
  # Calcul du BIC en omettant l'intercept dans le nombre de paramètres
  BIC_SICZIPLN_no_intercept<-v_loglik_SICZIPLN_with_intercept-(res_PLN$nb_param)*0.5*log(n)
  BIC_SICZIPLN_approximer<-v_loglik_SICZIPLN_with_intercept-(res_PLN$nb_param-length(which(c(B[-1,])==0)))*0.5*log(n)
  # -length(c(B)[c(B)==0])
  # BIC ZIPLN new
  BIC_SICZIPLN<-v_loglik_ZIPLN_no_intercept-(res_PLN$nb_param)*0.5*log(n)

  # BIC_SICZIPLN<-list(BIC_SICZIPLN=BIC_SICZIPLN_no_intercept,BIC_SICZIPLN_approximer=BIC_SICZIPLN_approximer)

  BIC_ZIPLN<-res_PLN$BIC
  # list(BIC_ZIPLN_no_intercept=BIC_ZIPLN_no_intercept,BIC_ZIPLN_with_intercept=res_PLN$BIC)
  model_par<-list(B=B,Omega=Omega,Sigma=Sigma,Pi=Pi,B_zero=B_zero)
  var_par<-list(M=M,S=S,R=R,S2=S^2)

  v_loglik<-list(v_loglik_SICZIPLN_with_intercept=v_loglik_SICZIPLN_with_intercept,v_loglik_SICZIPLN_no_intercept=v_loglik_SICZIPLN_no_intercept,v_loglik_ZIPLN_no_intercept=v_loglik_ZIPLN_no_intercept,v_loglik_ZIPLN_with_intercept=res_PLN$loglik)
  SICprediction<-(1-R)*exp(M+0.5*S^2)
  # ICL
  q<-p
  entropy_ZI  <-  -sum(xlogx(1-R)) - sum(xlogx(R))
  entropy_PLN <- .5 * (n * p * log(2*pi*exp(1)) + sum(log(S^2)))
  entropy     <-  entropy_ZI + entropy_PLN
  ICL_approximer        <- BIC_SICZIPLN_approximer - entropy
  ICL<- BIC_SICZIPLN_no_intercept - entropy
  data<-list(X=X,Y=Y,O=O,X_zero=X_zero)
  return(list(B=B,model_par=model_par,var_par=var_par,SICprediction=SICprediction,data=data,v_loglik=v_loglik,BIC_SICZIPLN=BIC_SICZIPLN,BIC_SICZIPLN_approximer=BIC_SICZIPLN_approximer,BIC_ZIPLN=BIC_ZIPLN,elbo=result_M_step$value,res_ZIPLN=res_PLN,ICL=ICL,ICL_approximer=ICL_approximer))
}

#' @export
SICZIPLN_NLOPTR <- function(X, Y, X_zero, offset = FALSE,
                            lambda_fixed = log(n) * (ncol(X)),
                            optim_method = "NLOPT_LD_LBFGS",
                            max_it = max_it,initialisation_method){

  # ---- Petit mappage pour rester compatible avec les anciens noms optim() ----
  # Si l'utilisateur passe encore "CG", "BFGS", "L-BFGS-B", "Nelder-Mead", etc.,
  # on les convertit vers un algorithme nloptr equivalent. Comme des gradients
  # analytiques (grad_E_step, grad_M_step) sont fournis, un algorithme basé
  # gradient est privilégié par défaut.
  .map_to_nloptr_algo <- function(method){
    if (grepl("^NLOPT_", method)) return(method)   # déjà un nom nloptr valide
    switch(method,
           "CG"          = "NLOPT_LD_LBFGS",
           "BFGS"        = "NLOPT_LD_LBFGS",
           "L-BFGS-B"    = "NLOPT_LD_LBFGS",
           "Nelder-Mead" = "NLOPT_LN_NELDERMEAD",
           "NLOPT_LD_LBFGS")  # valeur par défaut
  }
  nloptr_algo <- .map_to_nloptr_algo(optim_method)

  prep_pln <- prepare_data(Y, X)
  Y <- prep_pln$Abundance
  X <- prep_pln[, -c(1, length(prep_pln))]
  p <- ncol(Y)
  n <- nrow(X)
  lambda <- lambda_fixed

  if (offset) {
    O <- matrix(rep(log(prep_pln$Offset), p), ncol = p, byrow = FALSE)
  } else {
    O <- matrix(0, nrow = n, ncol = p)
  }

  X      <- model.matrix(~., as.data.frame(X))
  X_zero <- model.matrix(~., as.data.frame(X_zero))

  # Choix de l'initialisation par ZIPLN ou LM (initialisation utilisé dans ZIPLN)

  # ---- Initialisation ZIPLN (inchangé) ----
  #
  res_PLN <- ZIPLN(
    Y ~ X[, -1] + offset(O[, 1]) | X_zero[, -1],
    zi = "col",
    control = ZIPLN_param(
      config_optim = list(
        algorithm ="CCSAQ",    # remplace CCSAQ par défaut
        maxeval   = 2000,       # réduit la limite max d'itérations (défaut 10000)
        ftol_rel  = 1e-6,       # tolérance moins stricte = arrêt plus rapide
        xtol_rel  = 1e-4
      ),
      trace = 0                 # supprime les logs verbeux de PLNmodels (gain marginal mais utile)
    )
  )

  res_init_ZI <- compute_ZIPLN_starting_point(Y = Y, X = X, X0 = X_zero, O = O, w = NULL)
  if(initialisation_method=="LM"){
    cat("***** Initialisation with LM *****", "\n")
  B      <-res_init_ZI$B##res_PLN$model_par$B#res_PLN$model_par$B+rnorm(1)#res_init_ZI$B+rnorm(1)#
  B_zero <-res_init_ZI$B0#  res_PLN$model_par$B0#
  Pi     <- res_init_ZI$R #res_PLN$model_par$Pi#
  S      <- res_init_ZI$S #res_PLN$var_par$S#
  S2     <- S^2
  M      <- res_init_ZI$M #res_PLN$var_par$M#
  R      <- res_init_ZI$R #res_PLN$var_par$R#
  } else  if(initialisation_method=="ZIPLN"){
    cat("***** Initialisation with ZIPLN *****", "\n")
    # res_ZISIC$res_ZIPLN$var_par$S
    # res_ZISIC$res_ZIPLN$var_par$R
    B      <-res_PLN$model_par$B#res_init_ZI$B#res_PLN$model_par$B#res_PLN$model_par$B+rnorm(1)#res_init_ZI$B+rnorm(1)#
    B_zero <-  res_PLN$model_par$B0#-res_init_ZI$B0#
    Pi     <- res_PLN$model_par$Pi#res_init_ZI$R #
    S      <- res_PLN$var_par$S#res_init_ZI$S #
    S2     <- S^2
    M      <- res_PLN$var_par$M#res_init_ZI$M
    R      <- res_PLN$var_par$R#res_init_ZI$R #

  } else{
    cat("***** Choose initialisation method between LM or ZIPLN *****", "\n")
  }
  # ---- Précalculs hors boucle (gain majeur) ----
  X      <- as.matrix(X)
  X_zero <- as.matrix(X_zero)
  Y      <- as.matrix(Y)
  d      <- ncol(X)
  n      <- nrow(Y)
  tX        <- t(X)
  tX_zero   <- t(X_zero)
  XtX       <- crossprod(X)                 # t(X) %*% X, une seule fois
  XtX_inv_Xt <-solve(XtX) %*% tX
  # XtX_chol  <- chol(XtX)
  # XtX_inv_Xt <- chol2inv(XtX_chol) %*% tX    # remplace solve(t(X)%*%X)%*%t(X) à chaque itération


  Omega <- omega_update(M, B, S2, X)
  Sigma <- solve(Omega)#chol2inv(chol(Omega))

  params_init_E_step <- c(as.vector(M), as.vector(S))
  params_init_M_step <- c(B,B_zero)

  # iter <- 0
  # epsilon_val <- 0.00005
  # lambda <- 0
  length_epsilon<-100
  E<-c()
  e1<-10
  E[1]<-e1
  for(t in 2:length_epsilon){
    E[t]=e1*(0.87)^(t-1)
  }
  iter<-0
  ELBO_epsilon <- rep(NA, length_epsilon)
  for(epsilon_val in E){
    iter<-iter+1
  v_loglik_tmp <- vloglik_col(X, Y, O, B, Sigma, Pi, R, M, S, X_zero, B_zero)
  cat("SICZIPLN: Epsilon value for telescoping :\t", iter, epsilon_val, "\t", v_loglik_tmp, "\n")
  ELBO_epsilon[iter]<-v_loglik_tmp
  ELBO <- c(v_loglik_tmp, rep(NA, max_it))

  for (i in 1:max_it) {

    # ---- Wrappers locaux : évite les soucis de .checkfunargs() de nloptr ----
    # On capture Y, X, O, R, B, Omega, X_zero, B_zero, lambda, epsilon_val
    # directement depuis l'environnement de la boucle plutôt que de les
    # repasser en arguments nommés à nloptr() (ce que nloptr valide plus
    # strictement que optim()).
    obj_E_wrap <- function(par) {
      objective_E_step(par, Y = Y, X = X, O = O, R = R, B = B, Omega = Omega,
                       X_zero = X_zero, B_zero = B_zero,
                       lambda = lambda, epsilon = epsilon_val)
    }
    grad_E_wrap <- function(par) {
      grad_E_step(par, Y = Y, X = X, O = O, R = R, B = B, Omega = Omega,
                  X_zero = X_zero, B_zero = B_zero,
                  lambda = lambda, epsilon = epsilon_val)
    }

    result_E_step <- nloptr(
      x0 = params_init_E_step,
      eval_f = obj_E_wrap,
      eval_grad_f = grad_E_wrap,
      opts = list(
        algorithm = nloptr_algo,
        maxeval   = 2000,
        maxtime  = 200,
        xtol_rel  = 1e-5,
        xtol_abs  =0,
        ftol_rel  = 1e-6,
        ftol_abs  = 1e-6

      )
    )
    # cat("Status E-step:", result_E_step$status, result_E_step$message, "\n")
    # cat("objective E_step : ",result_E_step$objective,"\n")
    M <- matrix(result_E_step$solution[1:(n * p)], n, p)
    # cat("sum of M : ",(sum(M)),"\n")
    S <- matrix(result_E_step$solution[(n * p + 1):(2 * n * p)], n, p)
    # cat("sum of S: ",(sum(S)),"\n")
    params_init_E_step <- c(as.vector(M), as.vector(S))

    S2 <- S^2

    A <- exp(O + M + 0.5 * S2)
    R <- logistic(A + X_zero %*% B_zero) * dirac(Y)

    Omega <- omega_update(M, B, S2, X)
    # Sigma <- solve(Omega)#chol2inv(chol(Omega))

    # obj_M_wrap <- function(par) {
    #   objective_M_step(par, Y = Y, X = X, O = O,
    #                    Omega = Omega, B = B, M = M, S = S, R = R,
    #                    X_zero = X_zero,
    #                    lambda = lambda, epsilon = epsilon_val)
    # }
    obj_M_wrap_with_grad_B<- function(par) {
      objective_M_step_with_grad_B(par, Y = Y, X = X, O = O,
                       Omega = Omega, M = M, S = S, R = R,
                       X_zero = X_zero,
                       lambda = lambda, epsilon = epsilon_val)
    }
    # grad_M_wrap_with_grad_B <- function(par) {
    #   objective_M_step_with_grad_B(par, Y = Y, X = X, O = O,
    #               Omega = Omega, M = M, S = S, R = R,
    #               X_zero = X_zero,
    #               lambda = lambda, epsilon = epsilon_val)
    # }
# Avec mise à jour B non analytique
    grad_M_wrap_with_grad_B <- function(par) {
      grad_M_step_with_grad_B(par, Y = Y, X = X, O = O,
                   Omega = Omega, M = M, S = S, R = R,
                   X_zero = X_zero,
                   lambda = lambda, epsilon = epsilon_val)
    }
    result_M_step <- nloptr(
      x0 = params_init_M_step,
      eval_f = obj_M_wrap_with_grad_B,
      eval_grad_f = grad_M_wrap_with_grad_B,
      opts = list(
        algorithm = nloptr_algo,
        maxeval   = 2000,
        maxtime  = 200,
        xtol_rel  = 1e-5,
        xtol_abs  =0,
        ftol_rel  =1e-6,
        ftol_abs  =1e-6
      )
    )
    # cat("objective M_step : ",result_M_step$objective,"\n")
    # ---- Mises à jour analytiques (accélérées) ----
    B_old <- B
    B  <- matrix(result_M_step$solution[1:(d*p)],d,p) #XtX_inv_Xt %*% M                  # plus de solve() recalculé ici
    # cat("sum of B : \n",sum((B)),"\n")
    Omega <- omega_update(M, B, S2, X)
    Sigma <- chol2inv(chol(Omega))

    B_zero <- matrix(result_M_step$solution[((d*p)+1):length(result_M_step$solution)],d,p)#matrix(result_M_step$solution, d, p)
    Pi<-(1/(1+exp(-X_zero%*%B_zero)))
    # cat("sum of B_zero : ",(sum(B_zero)),"\n")
    params_init_M_step <- c(B,B_zero)

    ELBO[i + 1] <- vloglik_col(X, Y, O, B, Sigma, Pi, R, M, S, X_zero, B_zero)
    # Test de la convergence
    convergence <- sum(abs(B - B_old))
    diff_ELBO   <- abs(ELBO[i] - ELBO[i + 1])
    # cat("diff_ELBO : ", diff_ELBO, "\n")
    # cat("convergence : ", convergence, "\n")
    # if (convergence <= 1e-5 || diff_ELBO <= 1e-5) {
      if (diff_ELBO <= 1e-5) {

      cat("SICZIPLN:epsilon telescoping : ", epsilon_val, "Convergence avant maxit\n")
      break
    }
  }
}
  # ---- Post-traitement ----
    B[abs(B)<=1e-5]<-0
  if (!is.null(colnames(Y))) {
    colnames(Omega)<-colnames(Y)
    rownames(Omega)<-colnames(Y)
    colnames(Sigma)<-colnames(Y)
    rownames(Sigma)<-colnames(Y)
    colnames(Pi)<-colnames(Y)
    colnames(R)<-colnames(Y)
    colnames(M)<-colnames(Y)
    colnames(S)<-colnames(Y)
    colnames(B)<-colnames(Y)
    colnames(B_zero)<-colnames(Y)
  }
  if(!is.null(rownames(Y))) {
    rownames(M)<-rownames(Y)
    rownames(S)<-rownames(Y)
  }
  # X
  if (!is.null(colnames(X))) {
    rownames(B)<-colnames(X)
  }
  if (!is.null(colnames(X_zero))) {
    rownames(B_zero)<-colnames(X_zero)
  }
  v_loglik_SICZIPLN_no_intercept   <- vloglik_col(X[, -1], Y, O, B[-1, ], Sigma, Pi, R, M, S, X_zero, B_zero)
  v_loglik_SICZIPLN_with_intercept <- vloglik_col(X, Y, O, B, Sigma, Pi, R, M, S, X_zero, B_zero)
  v_loglik_ZIPLN_no_intercept <- vloglik_col(X[, -1], Y, O, res_PLN$model_par$B[-1, ],
                                             res_PLN$model_par$Sigma, res_PLN$model_par$Pi,
                                             res_PLN$var_par$R, res_PLN$var_par$M, res_PLN$var_par$S,
                                             X_zero, B_zero)

  BIC_SICZIPLN_no_intercept <- v_loglik_SICZIPLN_with_intercept - (res_PLN$nb_param) * 0.5 * log(n)
  BIC_SICZIPLN_approximer   <- v_loglik_SICZIPLN_with_intercept -
    (res_PLN$nb_param - length(which(c(B[-1, ]) == 0))) * 0.5 * log(n)
  BIC_SICZIPLN <- v_loglik_SICZIPLN_with_intercept - (res_PLN$nb_param) * 0.5 * log(n)
  BIC_ZIPLN <- res_PLN$BIC

  model_par <- list(B = B, Omega = Omega, Sigma = Sigma, Pi = Pi, B_zero = B_zero)
  var_par   <- list(M = M, S = S, R = R, S2 = S^2)
  v_loglik  <- list(
    v_loglik_SICZIPLN_with_intercept = v_loglik_SICZIPLN_with_intercept,
    v_loglik_SICZIPLN_no_intercept   = v_loglik_SICZIPLN_no_intercept,
    v_loglik_ZIPLN_no_intercept      = v_loglik_ZIPLN_no_intercept,
    v_loglik_ZIPLN_with_intercept    = res_PLN$loglik
  )
  SICprediction <- (1 - R) * exp(M + 0.5 * S^2)

  entropy_ZI  <- -sum(xlogx(1 - R)) - sum(xlogx(R))
  entropy_PLN <- .5 * (n * p * log(2 * pi * exp(1)) + sum(log(S^2)))
  entropy     <- entropy_ZI + entropy_PLN
  ICL_approximer <- BIC_SICZIPLN_approximer - entropy
  ICL <- BIC_SICZIPLN_no_intercept - entropy

  data <- list(X = X, Y = Y, O = O, X_zero = X_zero)

  return(list(ELBO_epsilon=ELBO_epsilon,ELBO = na.omit(ELBO), res_init_ZI = res_init_ZI,
              B = B, model_par = model_par, var_par = var_par,
              SICprediction = SICprediction, data = data, v_loglik = v_loglik,
              BIC_SICZIPLN_with_Intercept = BIC_SICZIPLN, BIC_SICZIPLN_approximer = BIC_SICZIPLN_approximer,
              BIC_ZIPLN = BIC_ZIPLN, BIC_SICZIPLN_no_intercept = BIC_SICZIPLN_no_intercept,
              elbo = result_M_step$objective, res_ZIPLN = res_PLN,
              ICL = ICL, ICL_approximer = ICL_approximer
  ))
  }

#####################################
#' SIC-Regularized ZIPLN Model with Optional Grid Search
#'
#' Fits a Sparse Information Criterion (SIC) penalized
#' Zero-Inflated Poisson Log-Normal (ZIPLN) model using a variational EM
#' algorithm. The procedure estimates model parameters and variational
#' parameters jointly, and optionally performs a grid search over the
#' regularization parameter \code{lambda}.
#'
#' @param X Numeric matrix of covariates for the latent Gaussian component
#'   (\eqn{n \times d}).
#' @param Y Matrix of observed count data (\eqn{n \times p}).
#' @param X_zero Design matrix for the zero-inflation component.
#' @param offset Logical; if \code{TRUE}, offsets are included in the model.
#' @param lambda_fixed Fixed value of the SIC regularization parameter.
#'   Default is \code{log(n) * ncol(X)}.
#' @param optim_method Optimization method passed to \code{optim}
#'   (e.g., \code{"BFGS"}).
#' @param max_it Maximum number of EM iterations used in the optimization
#'   procedure.
#' @param length_lambda Number of grid points used for lambda search when
#'   \code{grid_search = TRUE}.
#' @param grid_search Logical; if \code{TRUE}, performs a parallel grid search
#'   over a range of \code{lambda} values and returns all fitted models.
#'
#' @return If \code{grid_search = FALSE}, returns a fitted SIC-ZIPLN model
#'   as a list (output of \code{SICZIPLN_optim}).
#'
#'   If \code{grid_search = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{solution}: list of fitted models for each lambda value,
#'   \item \code{seq_lambda}: sequence of tested lambda values,
#'   \item \code{BIC_all_lambda}: BIC values for each lambda,
#'   \item \code{best_lambda_indice}: index of selected lambda based on BIC,
#'   \item \code{best_lambda_indice_BIC_approx}: index based on approximate BIC.
#' }
#'
#' @details
#' The function operates in two modes:
#'
#' \itemize{
#'   \item \strong{Single fit mode:} runs a single variational EM optimization
#'   using a fixed \code{lambda}.
#'   \item \strong{Grid search mode:} explores multiple values of \code{lambda}
#'   using parallel computing (via \code{parallel::parLapply}).
#' }
#'
#' The grid includes:
#' \itemize{
#'   \item heuristic values based on \eqn{\log(n)},
#'   \item logarithmically spaced values up to a maximum threshold.
#' }
#'
#' Model selection is performed using BIC criteria computed from fitted models.
#'
#' Parallel computation is used to speed up the grid search by distributing
#' independent model fits across multiple CPU cores.
#'
#' @note
#' This function depends on:
#' \itemize{
#'   \item \code{SICZIPLN_optim} (main optimizer),
#'   \item \code{prepare_data},
#'   \item the \code{parallel} package.
#' }
#'
#' @examples
#' # Fit single model
#' # fit <- SICZIPLN(X, Y, X_zero)
#'
#' # Grid search version
#' # fit_grid <- SICZIPLN(X, Y, X_zero, grid_search = TRUE)
#'
#' @export
SICZIPLN<-function(X,Y,X_zero,offset=FALSE,lambda_fixed=(log(nrow(X))*ncol(X)),optim_method="BFGS",max_it=200,length_lambda=100,grid_search=FALSE,initialisation_method="LM"){
  prep_pln<-prepare_data(Y,X)
  p<-ncol(prep_pln$Abundance)
  n<-nrow(prep_pln$Abundance)
  if(offset==TRUE){
    O<-matrix(rep(log(prep_pln$Offset),p),ncol=p,byrow = FALSE)
  }
  if(offset==FALSE){
    O<-matrix(0,nrow=n,ncol=p)
  }
  # solution<-SICZIPLN_optim(X,Y,X_zero,offset=offset,lambda_fixed=log(nrow(X))*(ncol(X)),optim_method=optim_method,max_it=200)
    solution<-SICZIPLN_NLOPTR(X,Y,X_zero,offset=offset,lambda_fixed=lambda_fixed,optim_method=optim_method,max_it=max_it,initialisation_method=initialisation_method)
  BIC_ZIPLN<-solution$res_ZIPLN$BIC
  if(grid_search==FALSE){
    solution<-solution
    return(solution)
  }
  if(grid_search==TRUE){
    # Creation de la grille
    length_lambda=length_lambda
    lambda_min=0#(log(nrow(X))*(ncol(X)*ncol(Y))+ncol(Y))/10
    lambda_max=3*(log(nrow(X))*(ncol(X)))
    # Charger le package parallel


    # Calculer la valeur log(n) une seule fois pour l'utiliser plus tard
    log_n <- log(nrow(X))

    # Définir les valeurs spécifiques de lambda
    lambda_min <- 0
    lambda_max <- 5 * log_n * ncol(X)*ncol(Y)
    lambda_sic1 <- log_n
    lambda_sic2 <- log_n * ncol(X)
    lambda_sic3 <- log_n * ncol(X) / 2
    lambda_sic4 <- log_n * ncol(X) * ncol(Y)

    # Créer la séquence de valeurs lambda
    seq_lambda <- c(0, lambda_sic1, lambda_sic2, lambda_sic3, lambda_sic4,
                    10^(seq(0, log(lambda_max), length.out = length_lambda)))

    # Créer un cluster avec un nombre de cœurs spécifié
    n_cores <- detectCores() - 2  # Utiliser tous les cœurs sauf un pour ne pas saturer le système
    cl <- makeCluster(n_cores)
    # clusterExport(cl, list("SICPLN_optim", "X", "Y"))
    # Exporter toutes les fonctions et objets de l'environnement global vers le cluster
    clusterExport(cl, varlist = ls(envir = .GlobalEnv))  #
    # clusterEvalQ(cl, library(PLNmodels))
    # Initialiser une liste pour stocker les résultats et un vecteur pour le BIC
    res_fs <- list()
    BIC_lambda <- numeric(length(seq_lambda))

    # Paralléliser l'exécution de la boucle avec parLapply
    res_fs <- parLapply(cl, seq_lambda, function(lambda) {
      # Exécuter l'optimisation pour chaque valeur de lambda
      # result <- SICZIPLN_optim(X,Y,X_zero,offset=FALSE,lambda_fixed=lambda,optim_method=optim_method,max_it=200)
      result <- SICZIPLN_NLOPTR(X,Y,X_zero,offset=FALSE,lambda_fixed=lambda,optim_method=optim_method,max_it=200)
      # Retourner le résultat et la valeur du BIC
      list(result = result, BIC= result[7],BIC_approximer = result[8])  # Assumons que le BIC est à la 10e position
    })

    # Extraire les BIC et les résultats
    BIC_lambda <- sapply(res_fs, function(x) x$BIC)
    BIC_lambda_approximer <- sapply(res_fs, function(x) x$BIC_approximer)
    res_fs <- lapply(res_fs, function(x) x$result)

    # Nommer les éléments dans la liste des résultats
    names(res_fs) <- paste0("SICZIPLN_lambda_", round(1:length(seq_lambda), 0))

    # Trouver le meilleur lambda basé sur le BIC (plus bas BIC)
    best_lambda_indice <- which(BIC_ZIPLN <= max(unlist(BIC_lambda)))+1
    best_lambda_indice_BIC_approx <- which(BIC_ZIPLN <= max(unlist(BIC_lambda_approximer)))+1

    # Créer la solution finale
    solution <- list(solution = res_fs, seq_lambda = seq_lambda,
                     BIC_all_lambda = unlist(BIC_lambda), best_lambda_indice = best_lambda_indice,best_lambda_indice_BIC_approx=best_lambda_indice_BIC_approx)

    # Arrêter le cluster après l'exécution
    stopCluster(cl)

  }
  return(solution)
}

#' @export
#'
omega_update <- function(M, B, S2,X) {
  ones_n    <- rep(1, n)
  M_resid <- M - X %*% B
  mat <- crossprod(M_resid) + diag(c(crossprod(ones_n, S2)))
  n * chol2inv(chol(mat))                 # plus rapide + plus stable que solve()
}
