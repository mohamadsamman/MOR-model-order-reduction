install.packages("dtwclust", dependencies=TRUE)
install.packages("sensitivity", dependencies=TRUE)
install.packages("MASS", dependencies=TRUE)
install.packages("stats", dependencies=TRUE)
install.packages("graphics", dependencies=TRUE)

library(dtwclust)
library(sensitivity)
library(MASS)
library(stats)
library(graphics)


#set seed 
set.seed(123456) 

##################################################################################

borehole <- function(x)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # x = c(rw, riw, r, Tu, Hu, Tum, Hum, Tlm, Hlm, Tl, Hl, L, Kw)
  #
  ##########################################################################
  
  xx <- matrix(x, ncol=13)
  
  rw  <- xx[,1]
  riw <- xx[,2]
  r   <- xx[,3]
  Tu  <- xx[,4]
  Hu  <- xx[,5]
  Tum <- xx[,6]
  Hum <- xx[,7]
  Tlm <- xx[,8]
  Hlm <- xx[,9]
  Tl  <- xx[,10]
  Hl  <- xx[,11]
  L   <- xx[,12]
  Kw  <- xx[,13]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  frac11 <- (Tum - Tlm) * (Hum - Hlm)
  frac2 <- frac1 / frac11
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

##################################################################################

EchantBorehole <- function(N){
  
  # Description de la fonction :
  # Cette fonction genere un échantillon (pour le modele borhole) de taille N.
  #
  # Ces 13 variables sont supposées statistiquement independantes
  #
  # Entrées de la fonction :
  # - taille = taille de l'échantillon
  #
  # Description des réponses : 
  # - X = matrice de taille N x 13
  
  X = matrix(NA, N, 13)
  
  X[,1] <- rnorm(N, 0.1, 0.015)
  X[,2] <- rnorm(N, 0.05, 0.01)
  X[,3] <- rlnorm(N, 7.71, 1.0056)
  X[,4] <- runif(N, 63100, 116000)
  X[,5] <- runif(N, 1000, 1100)
  X[,6] <- runif(N, 6310, 11600)
  X[,7] <- runif(N, 900, 1000)
  X[,8] <- runif(N, 631, 1160)
  X[,9] <- runif(N, 800, 900)
  X[,10] <- runif(N, 63.1, 116)
  X[,11] <- runif(N, 700, 800)
  X[,12] <- runif(N, 1120, 1680)
  X[,13] <- runif(N, 3000, 12000)
  
  colnames(X) <- c("rw","riw","r","Tu","Hu","Tum","Hum","Tlm","Hlm","Tl","Hl","L","Kw") # noms de variables aleatoires
  
  return(X)
  
}#end function

# Valeurs minimales et maximales pour les variables
valeurs_min <- c(0.05, 0.02, 100, 63100, 1000, 6310, 900, 631, 800, 63.1, 700, 1120, 3000)
valeurs_max <- c(0.15, 0.08, 50000, 116000, 1100, 11600, 1000, 1160, 900, 116, 800, 1680, 12000)

# Moyenne de la lognormale pour r (μ = 7.71, σ = 1.0056)
mu_r <- 7.71
sigma_r <- 1.0056
moyenne_r <- exp(mu_r + (sigma_r^2)/2)

# Valeurs moyennes pour les autres variables
valeurs_moy <- c(0.10, 0.05, moyenne_r, 93150, 1050, 9235, 950, 895.5, 850, 89.55, 750, 1400, 7500)

# Calcul du débit pour ces valeurs
debit_min <- borehole(valeurs_min)
debit_max <- borehole(valeurs_max)
debit_moy <- borehole(valeurs_moy)

print(c(debit_min, debit_max, debit_moy))

# Générer un échantillon Monte Carlo de taille 1000
N <- 1000
X_mc <- EchantBorehole(N)
y_mc <- apply(X_mc, 1, borehole)

# Moyenne et variance de la sortie
moyenne_mc <- mean(y_mc)
variance_mc <- var(y_mc)

# Histogramme de la sortie
hist(y_mc, main="Histogramme des sorties", xlab="Débit (m3/an)", col="lightblue", border="black")

# Sortie évaluée avec la moyenne des entrées
debit_moy_entrees <- borehole(valeurs_moy)

# Comparaison entre la moyenne Monte Carlo et la moyenne des entrées
print(c(moyenne_mc,variance_mc ,debit_moy_entrees))

# Quantile d'ordre 95%
quantile_95 <- quantile(y_mc, 0.95)
print(quantile_95)

# Intervalle de confiance à 95%
n_repetitions <- 100
ic_95 <- replicate(n_repetitions, {
  y_mc_repet <- apply(EchantBorehole(N), 1, borehole)
  quantile(y_mc_repet, c(0.025, 0.975))
})

moyenne_ic <- colMeans(ic_95)
moyenne_ic

# Affichage du quantile et de l'intervalle de confiance
print(quantile_95)
print(moyenne_ic)

# Définir la taille d'échantillon initiale
taille_echantillon <- 2 * 10^6
proba <- numeric()

# Itérer en augmentant la taille de l'échantillon
for (n in seq(2*10^6, 10^7, by=2*10^6)) {
  X_mc <- EchantBorehole(n)
  y_mc <- apply(X_mc, 1, borehole)
  proba <- c(proba, mean(y_mc > 250))
}

# Affichage de la convergence
plot(seq(2*10^6, 10^7, by=2*10^6), proba, type="l", xlab="Taille de l'échantillon", ylab="Probabilité (y > 250)")

# Appliquer la méthode de Morris avec un nombre d'évaluations inférieur à 100
N_morris <- 100
X_morris <- EchantBorehole(N_morris)

# Calcul des indices de Morris
morris_results <- morris(model=borehole, factors=colnames(X_morris), X1=X_morris, r=10, design="oat")

# Résultats du criblage
print(morris_results)

# Sélectionner les variables influentes (en supposant qu'on ait fait ce tri)
X_selected <- X_morris[, c(1, 2, 3)]  # Exemple avec quelques variables

# Régression linéaire
lm_model <- lm(y_mc ~ X_selected)
summary(lm_model)

# Calcul des indices de Sobol
sobol_results <- sobol(model=borehole, X1=X_mc, X2=X_mc, order=2, nboot=100)

# Affichage des résultats
print(sobol_results)


