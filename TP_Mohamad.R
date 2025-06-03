library(dtwclust)
library(sensitivity)
library(MASS)
library(stats)
library(graphics)

#set seed 
set.seed(105) 

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

print(c(debit_min, debit_moy, debit_max))

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
n_repetitions <- 1000
ic_95 <- replicate(n_repetitions, {
  y_mc_repet <- apply(EchantBorehole(N), 1, borehole)
  quantile(y_mc_repet, c(0.025, 0.975))
})

moyenne_ic <- rowMeans(ic_95)
moyenne_ic

library(ggplot2)
#########################################
#Monte-Carlo
#########################################
# Paramètres de la simulation
taille_initiale <- 2 * 10^6
pas_incrementation <- 2 * 10^6
seuil_debit <- 250
erreur_relative_cible <- 0.1  # 10 %

taille <- taille_initiale
erreur_relative <- 1  # Erreur initialisée à 100%
probabilites <- c()
taille_echantillons <- c()

# Initialiser un compteur pour les itérations
iteration <- 0
max_iterations <- 20

while (erreur_relative > erreur_relative_cible && iteration < max_iterations) {
  print(iteration)
  # Générer un échantillon de taille 'taille'
  X_mc <- EchantBorehole(taille)
  y_mc <- apply(X_mc, 1, borehole)
  
  # Calculer la probabilité de dépasser le seuil
  prob_sup_250 <- mean(y_mc > seuil_debit)
  
  # Stocker les résultats
  probabilites <- c(probabilites, prob_sup_250)
  taille_echantillons <- c(taille_echantillons, taille)
  
  # Calculer l'erreur relative si taille > 1 (pour éviter division par zéro)
  if (length(probabilites) > 1) {
    erreur_relative <- abs(probabilites[length(probabilites)] - probabilites[length(probabilites) - 1]) / 
      probabilites[length(probabilites)]
  }
  
  # Incrémenter la taille de l'échantillon
  taille <- taille + pas_incrementation
  
  # Incrémenter le compteur d'itérations
  iteration <- iteration + 1
}

# Taille d'échantillon nécessaire pour atteindre l'erreur relative cible
taille_finale <- taille_echantillons[length(taille_echantillons)]
cat("La taille minimale d'échantillon nécessaire est :", taille_finale, "\n")

# Graphe de convergence
df <- data.frame(Taille = taille_echantillons, Probabilite = probabilites)

ggplot(df, aes(x = Taille, y = Probabilite)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  xlab("Taille de l'échantillon") +
  ylab("Probabilité (Débit > 250 m³/an)") +
  ggtitle("Convergence de la probabilité par Monte Carlo") +
  theme_minimal()

#########################################
#morris
#########################################
#3a

# Noms des variables d'entrée
names_vars <- c("rw", "riw", "r", "Tu", "Hu", "Tum", "Hum", "Tlm", "Hlm", "Tl", "Hl", "L", "Kw")

# Configuration de la méthode de Morris
r <- 7   # Nombre de trajectoires (adapté pour moins de 100 évaluations)
levels <- 4  # Nombre de niveaux

# Bornes des variables d'entrée
bornes_min <- c(0.05, 0.02, 100, 63100, 1000, 6310, 900, 631, 800, 63.1, 700, 1120, 3000)
bornes_max <- c(0.15, 0.08, 50000, 116000, 1100, 11600, 1000, 1160, 900, 116, 800, 1680, 12000)

# Fonction de wrapper pour adapter la fonction borehole à morris()
borehole_wrapper <- function(X) {
  apply(X, 1, borehole)
}

# Plan d'échantillonnage de Morris
morris_result <- morris(
  model = borehole_wrapper,
  factors = names_vars,
  r = r,
  design = list(type = "oat", levels = levels, grid.jump = 1),
  binf = bornes_min,
  bsup = bornes_max
)

# Affichage des résultats
print(morris_result)

# Graphique des résultats
plot(morris_result, xlim=c(0, 120), main = "Criblage des variables d'entrée - Méthode de Morris")


#3b
# Variables influentes d'après le criblage de Morris
variables_influentes <- c("Hl", "Hu", "L", "Kw", "rw")

# Moyenne des variables non influentes fixées
moyennes_non_influentes <- c(
  riw = 0.05, r = exp(7.71 + (1.0056^2)/2),  # moyenne de la lognormale
  Tum = 9235, Hum = 950, Tlm = 895.5, Hlm = 850, Tl = 89.55
)

# Taille de l'échantillon Monte Carlo
N <- 100

# Génération des variables influentes
X_influentes <- data.frame(
  Hl = runif(N, 700, 800),
  Hu = runif(N, 1000, 1100),
  L = runif(N, 1120, 1680),
  Kw = runif(N, 3000, 12000),
  rw = runif(N, 0.05, 0.15)
)

# Fixation des variables non influentes à leur moyenne
for (var in names(moyennes_non_influentes)) {
  X_influentes[[var]] <- moyennes_non_influentes[[var]]
}

# Calcul du débit d'eau (sortie)
debit <- apply(X_influentes, 1, borehole)

# Étape 3 : Scatterplots entre les variables influentes et la sortie (Débit)
pairs(data.frame(X_influentes[, variables_influentes], Debit = debit),
      main = "Scatterplots entre les variables influentes et la sortie (Débit)",
      col = "dodgerblue", pch = 16)

# Étape 4 : Régression linéaire
modele_lm <- lm(debit ~ ., data = X_influentes[, variables_influentes])
summary(modele_lm)

# Calcul des indices de sensibilité SRC²
SRC2 <- summary(modele_lm)$coefficients[-1, 1]^2
SRC2 <- SRC2 / sum(SRC2)  # Normalisation des indices SRC²
names(SRC2) <- variables_influentes
print("Indices de sensibilité SRC² normalisés :")
print(SRC2)

#########################################
#Sobol
#########################################
# Paramètres de simulation
taille_initiale <- 5000
pas_incrementation <- 5000
erreur_relative_cible <- 0.05  # Erreur cible à 5 %
taille <- taille_initiale
erreur_relative <- 1
indices_sobol <- list()
taille_echantillons <- c()
iteration <- 0
max_iterations <- 10

# Fonction borehole wrapper
borehole_wrapper <- function(X) {
  apply(X, 1, borehole)
}

while (erreur_relative > erreur_relative_cible && iteration < max_iterations) {
  print(paste("Itération:", iteration + 1, " - Taille d'échantillon:", taille))
  
  # Génération des échantillons
  X1 <- EchantBorehole(taille)
  X2 <- EchantBorehole(taille)
  
  # Analyse de Sobol
  sobol_result <- sobol(model = borehole_wrapper, X1 = X1, X2 = X2, order = 1, nboot = 100)
  # Stockage des résultats
  indices_sobol[[length(indices_sobol) + 1]] <- sobol_result$S$original
  taille_echantillons <- c(taille_echantillons, taille)
  
  # Calcul de l'erreur relative sur le premier indice de Sobol (ex. "rw")
  if (length(indices_sobol) > 1) {
    erreur_relative <- abs(indices_sobol[[length(indices_sobol)]][1] - indices_sobol[[length(indices_sobol) - 1]][1]) / 
      abs(indices_sobol[[length(indices_sobol)]][1])
  }
  
  # Augmenter la taille de l'échantillon pour la prochaine itération
  taille <- taille + pas_incrementation
  iteration <- iteration + 1
}

# Résultats finaux
taille_finale <- taille_echantillons[length(taille_echantillons)]
cat("Taille d'échantillon finale pour atteindre l'erreur cible :", taille_finale, "\n")

# Visualisation des indices de Sobol pour la dernière itération
df_sobol <- data.frame(
  Variable = c("rw", "riw", "r", "Tu", "Hu", "Tum", "Hum", "Tlm", "Hlm", "Tl", "Hl", "L", "Kw"),
  Indice_Sobol = indices_sobol[[length(indices_sobol)]]
)

ggplot(df_sobol, aes(x = Variable, y = Indice_Sobol)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  xlab("Variables d'entrée") +
  ylab("Indice de Sobol") +
  ggtitle("Indices de Sobol - Analyse de sensibilité (convergence)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

