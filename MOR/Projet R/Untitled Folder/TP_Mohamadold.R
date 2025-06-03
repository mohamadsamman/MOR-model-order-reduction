library(dtwclust)
library(sensitivity)
library(MASS)
library(stats)
library(graphics)

#set seed 
set.seed(10571) 

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

################################################################################"
#non testé

# library(ggplot2)

# # Paramètres de la simulation
# taille_initiale <- 2 * 10^6
# pas_incrementation <- 2 * 10^6
# seuil_debit <- 250
# erreur_relative_cible <- 0.1  # 10 %

# taille <- taille_initiale
# erreur_relative <- 1  # Erreur initialisée à 100%
# probabilites <- c()
# taille_echantillons <- c()

# while (erreur_relative > erreur_relative_cible) {
  
#   # Générer un échantillon de taille 'taille'
#   X_mc <- EchantBorehole(taille)
#   y_mc <- apply(X_mc, 1, borehole)
  
#   # Calculer la probabilité de dépasser le seuil
#   prob_sup_250 <- mean(y_mc > seuil_debit)
  
#   # Stocker les résultats
#   probabilites <- c(probabilites, prob_sup_250)
#   taille_echantillons <- c(taille_echantillons, taille)
  
#   # Calculer l'erreur relative si taille > 1 (pour éviter division par zéro)
#   if (length(probabilites) > 1) {
#     erreur_relative <- abs(probabilites[length(probabilites)] - probabilites[length(probabilites) - 1]) / 
#       probabilites[length(probabilites)]
#   }
  
#   # Incrémenter la taille de l'échantillon
#   taille <- taille + pas_incrementation
# }

# # Taille d'échantillon nécessaire pour atteindre l'erreur relative cible
# taille_finale <- taille_echantillons[length(taille_echantillons)]
# cat("La taille minimale d'échantillon nécessaire est :", taille_finale, "\n")

# # Graphe de convergence
# df <- data.frame(Taille = taille_echantillons, Probabilite = probabilites)

# ggplot(df, aes(x = Taille, y = Probabilite)) +
#   geom_line(color = "blue") +
#   geom_point(color = "red") +
#   xlab("Taille de l'échantillon") +
#   ylab("Probabilité (Débit > 250 m³/an)") +
#   ggtitle("Convergence de la probabilité par Monte Carlo") +
#   theme_minimal()

#3a

# # Noms des variables d'entrée
# names_vars <- c("rw", "riw", "r", "Tu", "Hu", "Tum", "Hum", "Tlm", "Hlm", "Tl", "Hl", "L", "Kw")

# # Configuration de la méthode de Morris
# r <- 7   # Nombre de trajectoires (adapté pour moins de 100 évaluations)
# levels <- 4  # Nombre de niveaux

# # Bornes des variables d'entrée
# bornes_min <- c(0.05, 0.02, 100, 63100, 1000, 6310, 900, 631, 800, 63.1, 700, 1120, 3000)
# bornes_max <- c(0.15, 0.08, 50000, 116000, 1100, 11600, 1000, 1160, 900, 116, 800, 1680, 12000)

# # Fonction de wrapper pour adapter la fonction borehole à morris()
# borehole_wrapper <- function(X) {
#   apply(X, 1, borehole)
# }

# # Plan d'échantillonnage de Morris
# morris_result <- morris(
#   model = borehole_wrapper,
#   factors = names_vars,
#   r = r,
#   design = list(type = "oat", levels = levels, grid.jump = 1),
#   binf = bornes_min,
#   bsup = bornes_max
# )

# # Affichage des résultats
# print(morris_result)

# # Graphique des résultats
# plot(morris_result, xlab = expression(mu*), ylab = expression(sigma),
#      main = "Criblage des variables d'entrée - Méthode de Morris")

#3b
# # Étape 1 : Variables influentes d'après le criblage de Morris
# variables_influentes <- c("Tu", "Hu", "L")  # Par exemple, à adapter selon les résultats du criblage

# # Moyenne des variables non influentes
# moyennes_non_influentes <- c(
#   rw = 0.10, riw = 0.05, r = exp(7.71 + (1.0056^2)/2),
#   Tum = 9235, Hum = 950, Tlm = 895.5, Hlm = 850,
#   Tl = 89.55, Hl = 750, Kw = 7500
# )

# # Taille de l'échantillon Monte Carlo
# N <- 100

# # Génération des variables influentes
# X_influentes <- data.frame(
#   Tu = runif(N, 63100, 116000),
#   Hu = runif(N, 1000, 1100),
#   L = runif(N, 1120, 1680)
# )

# # Ajout des variables non influentes fixées à leur moyenne
# for (var in names(moyennes_non_influentes)) {
#   X_influentes[[var]] <- moyennes_non_influentes[[var]]
# }

# # Calcul du débit d'eau (sortie)
# debit <- apply(X_influentes, 1, borehole)

# # Étape 3 : Scatterplots
# pairs(data.frame(X_influentes[, variables_influentes], Debit = debit),
#       main = "Scatterplots entre les variables influentes et la sortie (Débit)")

# # Étape 4 : Modèle de régression linéaire
# modele_lm <- lm(debit ~ ., data = X_influentes[, variables_influentes])
# summary(modele_lm)

# # Indices de sensibilité SRC²
# SRC2 <- summary(modele_lm)$coefficients[-1, 1]^2
# SRC2 <- SRC2 / sum(SRC2)  # Normalisation des indices de sensibilité
# print(SRC2)


 # Taille de l'échantillon pour chaque bloc de variables
 n <- 10000

 # Noms des variables d'entrée
 names_vars <- c("rw", "riw", "r", "Tu", "Hu", "Tum", "Hum", "Tlm", "Hlm", "Tl", "Hl", "L", "Kw")

 # Génération des échantillons X1 et X2
 X1 <- data.frame(
   rw = runif(n, 0.05, 0.15),
   riw = runif(n, 0.02, 0.08),
   r = rlnorm(n, 7.71, 1.0056),
   Tu = runif(n, 63100, 116000),
   Hu = runif(n, 1000, 1100),
   Tum = runif(n, 6310, 11600),
   Hum = runif(n, 900, 1000),
   Tlm = runif(n, 631, 1160),
   Hlm = runif(n, 800, 900),
   Tl = runif(n, 63.1, 116),
   Hl = runif(n, 700, 800),
   L = runif(n, 1120, 1680),
   Kw = runif(n, 3000, 12000)
 )

 X2 <- data.frame(
   rw = runif(n, 0.05, 0.15),
   riw = runif(n, 0.02, 0.08),
   r = rlnorm(n, 7.71, 1.0056),
   Tu = runif(n, 63100, 116000),
   Hu = runif(n, 1000, 1100),
   Tum = runif(n, 6310, 11600),
   Hum = runif(n, 900, 1000),
   Tlm = runif(n, 631, 1160),
   Hlm = runif(n, 800, 900),
   Tl = runif(n, 63.1, 116),
   Hl = runif(n, 700, 800),
   L = runif(n, 1120, 1680),
   Kw = runif(n, 3000, 12000)
 )

 # Wrapper pour adapter la fonction borehole à sobol()
 borehole_wrapper <- function(X) {
   apply(X, 1, borehole)
 }

 # Calcul des indices de Sobol
 sobol_result <- sobol(
   model = borehole_wrapper,
   X1 = X1,
   X2 = X2,
   order = 2,
   nboot = 100
 )

 # Affichage des résultats
 print(sobol_result)

 # Graphique des indices de Sobol
 plot(sobol_result, main = "Indices de Sobol - Analyse de sensibilité", xlab = "Variables d'entrée")
