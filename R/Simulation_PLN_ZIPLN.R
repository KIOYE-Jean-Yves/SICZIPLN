# library(tidyverse)
# library(parallel)
# Fonction de generation de données zéro inflaté
# Paramètres d'entrée :
#   - n : taille des données
#   - mu=B : matrice qui encode les coefficients représentant l'influence de chaque variable explicative sur les comptages

# Sigma : matrice qui encode les rélations de dependances entre les comptages
# Pi: matrice qui encode la probabilité d'inflation de zero (elle est calulée dans la fonction Simulation_ZIPLN par une régression logistic)
# Retour :
#   Y: matrice de comptage
#
rZIPLN <- function(n     = 200,
                   mu    = rep(0, ncol(Sigma)),
                   Sigma = diag(1, 5, 5),
                   Pi    = matrix(1, n, ncol(Sigma)),
                   depths = rep(10, n))  {
  p <- ncol(Sigma)
  if (any(is.vector(mu), ncol(mu) == 1)) {
    mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE)
  }
  if (length(depths) != n) {
    depths <- rep(depths[1], n)
  }
  # adjust depths
  exp_depths <-rowSums(exp(rep(1, n) %o% diag(Sigma)/2 + mu)) ## sample-wise expected depths
  offsets <- 0#log(depths %o% rep(1, p)) - log(exp_depths) #log(1+matrix(rep(depths,p),n,p))
  Z <- mu + MASS::mvrnorm(n, rep(0, ncol(Sigma)), as.matrix(Sigma)) + offsets
  W <- matrix(rbinom(n * p, 1,  prob = Pi), n, p)
  Y <- matrix(rpois(n * p, as.vector(exp(Z))), n, p) * (1 - W)
  dimnames(Y) <- list(paste0("S", 1:n), paste0("Y", 1:p))
  return(list(Y=Y,offsets=exp(offsets)))
}
# logit    <- function(x) log(x/(1-x))

# Fonction pou calculer la probabilité d'inflation de zero
# Pi: matrice qui encode la probabilité d'inflation de zero (elle est calulée dans la fonction Simulation_ZIPLN par une régression logistic)
# Paramètres d'entrée :
# matrice x
# Retour:
# logistic(x)
logistic <- function(x) 1 / (1 + exp(-x))

# Fonction de génération de la matrice B
#
# Description :
#   Cette fonction construit une matrice B dont les dimensions sont déterminées par le nombre de colonnes de X (ncol_covar) et le nombre de colonnes de la matrice de comptages (ncol_especes). Elle initialise d'abord l'ensemble de la matrice avec la valeur 1, puis remplace aléatoirement certaines entrées par 0 afin d'atteindre selon la proportion de zéros spécifiée (proprortion_zero)
# Paramètres d'entrée :
#   - ncol_covar (int)   : Nombre de colonnes souhaité pour X (variables explicatives).
#   - ncol_especes (int) : Nombre de colonnes de la matrice de comptages (espèces).
#   - prop_zero (float)  : Proportion de zéros souhaitée dans la matrice B.
#
# Retour :
#   - Matrice B générée selon la proportion de zéros spécifiée.

Generation_matrice_B_control_entree_zero<-function(ncol_covar=6,ncol_especes=4,proprortion_zero=0.25){
  B_star  <- rep(1,ncol_especes*ncol_covar)
  taille_B<-length(B_star)
  zero_entry<-sample(1:taille_B,proprortion_zero*taille_B)
  B_star[zero_entry]<-0
  B_star  <- matrix(B_star,ncol=ncol_especes,nrow=ncol_covar,byrow = TRUE)
  return(B_star)
}

# Fonction simulation ZIPLN
#
# Description :
#   Cette fonction génère des données de comptage simulées suivant le modèle ZIPLN (Zero-Inflated Poisson Log-Normal). Elle genère aussi les matrices des variables X (associée aux comptages) et la matrice de variables X_zero associer à la proba d'inflation de zéro.
#
# Paramètres d'entrée :
#   - n (int) : Taille de l'échantillon.
#   - ncol_X (int) : Nombre de variables explicatives associées aux comptages.
#   - ncol_inflation (int) : Nombre de variables explicatives influençant l'inflation de zéros.
#   - B (float) : Matrice B (B_star) qui encode les coefficients représentant l'influence de chaque variable explicative sur les comptages.
#   - B0_star (float) : Matrice B0_star qui encode l'intensité des variables explicatives responsables de la probabilité d'inflation de zéros.
# Sigma_star: martice qui encode les rélations de dependances entre les comptages.
#
# Retour :
# Y: Matrice de comptage.
# X: matrice de variables explicatives associées aux comptages.
# X_zero: matrice de variables explicatives associées à l'inflation de zero.
# B: martice qui encode les coefficients représentant l'influence de chaque variable explicative sur les comptages.
# Sigma: martice qui encode les rélations de dependances entre les comptages.
# B0: Matrice qui encode l'intensité des variables explicatives responsables de la probabilité d'inflation de zéros.
# Pi Matrice qui encode la probabilité d'inflation de zéros.
#' @export
Simulation_ZIPLN<-function(n=100,ncol_X=5,ncol_X_zero=5,p=4,B_star,Sigma_star,B0_star,depths = rep(10000, n)){
  # generation des covariables
  # nombre de covariables environnementale
  n_=50*n
  d=ncol_X
  d_zero<-d
  data_evironement <- runif(n_*d,0,1.5)#rnorm(n, 2, 1)
  X<-matrix(data_evironement,nrow=n_,ncol=d)
  X <-cbind(rep(1,n_),X)
  # generation de covariable X_zero
  # nombre de covariables X_zero
  # d=4
  # data_evironement_X_zero <- runif(n*d_zero,0,1.5)#rnorm(n, 2, 1)
  # X_zero<-matrix(data_evironement_X_zero,nrow=n,ncol=d_zero)
  X_zero<-X

  # Paramettre de moyenne couche latente
  mu_star <- X %*% B_star

  # Pi pour composante de zero
  Pi_star <- logistic(X_zero %*% B0_star)#matrix(runif(n*p,0,1), n, p)#
  offsets=0

  # Z <- mu_star + mvrnorm(n, rep(0, ncol(Sigma_star)), as.matrix(Sigma_star)) + offsets
  # W <- matrix(rbinom(n * p, 1,  prob = Pi_star), n, p)
  # Y <- matrix(rpois(n * p, as.vector(exp(Z))), n, p) * (1 - W)
  # dimnames(Y) <- list(paste0("S", 1:n), paste0("Y", 1:p))
  # Y
  #
  Count_Offset<-rZIPLN(n_, mu = mu_star, Sigma = Sigma_star, Pi = Pi_star,depths = depths)
  Y <- Count_Offset$Y
  Offsets<-Count_Offset$offsets
  X <- X[rowSums(Y) > 0,]
  Pi_star<-Pi_star[rowSums(Y) > 0,]
  X_zero<-X_zero[rowSums(Y) > 0,]
  Y <- Y[rowSums(Y) > 0, ]
  X <- X[1:n,]
  Pi_star<-Pi_star[1:n,]
  X_zero<-X_zero[1:n,]
  Y <- Y[1:n, ]
  dimnames(Pi_star) <- list(paste0("obs_", 1:nrow(Y)), paste0("Pi_star_", 1:p))
  dimnames(X) <- list(paste0("obs_", 1:nrow(Y)), paste0("X_", 1:(d+1)))
  dimnames(X_zero) <- list(paste0("obs_", 1:nrow(Y)), paste0("X_zero_", 1:(d_zero+1)))
  dimnames(Y) <- list(paste0("obs_", 1:nrow(Y)), paste0("Y_", 1:p))
  return(data=list(X=X,Y=Y,X_zero=X_zero,Sigma_simulatation=Sigma_star,B_simulatation=B_star,Pi_simulatation=Pi_star,B0_simulatation=B0_star,depths = depths,Offsets=Offsets))
}

# Fonction de calcul du taux de Vrai positif & du taux de vrai negatif
# Paramètres d'entrée :
#   - vrai_beta : Vrai matrice
#   - beta_estime : matrice estimé
# Retour :
#   - TNR : Taux de vrai négatif
#   - TPR : Taux de vrai positif
sparsity_recognition<-function(vrai_beta, beta_estime){
  TN=0
  TP=0
  for(i in 1:nrow(vrai_beta)){
    for(j in 1:ncol(vrai_beta)){
      if(vrai_beta[i,j]!=0 & beta_estime[i,j]!=0){TP=TP+1}

      if(vrai_beta[i,j]==0 & beta_estime[i,j]==0){
        # print(c(i,j))
        TN=TN+1
      }

      # print(TN)
    }
  }

  beta0=c(vrai_beta==0)
  beta_non_zero=c(vrai_beta!=0)

  # Taux de non zero estimée estimée zero
  TPR=TP/(length(beta_non_zero[beta_non_zero==TRUE]))
  # Taux de zero estimée zéro
  TNR= TN/(length(beta0[beta0==TRUE]))
  # precision=TP/(TP+abs(length(beta_non_zero[beta_non_zero==T])-TP))
  # recall=TP/(TP+1)
  sol=list(TNR=TNR,TPR=TPR)
  return(sol)
}

# Generation de la matrice Sigma
# Sigma est generé suivant une loi uniforme [-1,1]: matrix(runif(p*p,-1,1),p,p)
# Paramètres d'entrée :
#   - p : nombre d'espèce
#   - structure_dependance ("full", ou "diag") : Structure de dependance souhaité. Si structure_dependance="full" un matrice de covariance pleine est géneré, si structure_dependance="diag" un matrice de covariance diagonale est généré.
# Retour :
#   Sigma: matrice de covariance
generate_Sigma<-function(p,structure_dependance="full"){
  if(structure_dependance=="full"){
    mcov=matrix(runif(p*p,-1,1),p,p)#1
    Sigma_star=(mcov%*%t(mcov))
  }
  if(structure_dependance=="diag"){
    mcov=matrix(runif(p*p,-1,1),p,p)#1
    Sigma_star=diag(diag((mcov%*%t(mcov))))
  }
  return(Sigma_star)
}
# Application de la méthode de selection de variables
# source("/home/tjkioye/Documents/These/2024_2025/SICZIPLN/Optim_ZI_PLN_col_dep.R")
# source("~/Documents/These/2024_2025/SICZIPLN/SICPLN/Good_SICPLN_optim2.R")
