# Taille de la population
N <- 20000

# Générer les variables nécessaires au calcul des x selon les distributions spécifiées
z1 <- rbinom(N, size = 1, prob = 0.3)                # z1i ∼ Bernoulli(0.3)
z2 <- runif(N, min = 0, max = 2)                     # z2i ∼ Uniforme(0, 2)
z3 <- rexp(N, rate = 1)                              # z3i ∼ Exponentielle(1)
epsilon <- rnorm(N, mean = 0, sd = 1)                # εi ∼ Normale(0, 1)

# Calculer les variables x selon les relations données
x1 <- z1                                             
x2 <- z2 + 0.3 * x1
x3 <- z3 + 0.2 * (x1 + x2)

x_beta <- 2 + x1 + x2 + x3

# Définir la corrélation entre les variables auxiliaires et la variable d'intérêt
# 0.3 / 0.5 / 0.8
rho_target <- 0.5
sigma <- sqrt(var(x_beta) * (1 / rho_target^2 - 1))

# Calculer la variable réponse y selon le modèle
y <- x_beta + sigma * epsilon         # yi = 2 + x1i + x2i + x3i + σ*εi
# Proche de la corrélation cible
cor(y, x_beta)

pop <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2,
  x3 = x3
)

vrai_tot <- round(sum(pop$y),1)
vrai_tot

vrai_moy <- round(mean(pop$y),3)
vrai_moy

# Paramètres de la simulation

# Taille de l'échantillon probabiliste
n_prob <- 500
# Taille de l'échantillon non probabiliste
n_non_prob <- 500

# Paramètres pour la définition des vrai proba d'inclusion (inconnu) pour l'échantillon non proba
theta1 <- 0.1
theta2 <- 0.2
# 0.1 ou 0.3
theta3 <- 0.3

# Paramètre du modèle de Frank
a <- 10
# Nombres de groupes homogènes de réponse
nb_GHR <- 10
