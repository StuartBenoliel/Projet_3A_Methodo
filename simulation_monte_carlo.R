library(sampling)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

rm(list=ls())
set.seed(1)

source(file = "fonction.R")

# Taille de la population
N <- 20000

# Générer les variables auxiliaires selon les distributions spécifiées
z1 <- rbinom(N, size = 1, prob = 0.5)                # z1i ∼ Bernoulli(0.5)
z2 <- runif(N, min = 0, max = 2)                     # z2i ∼ Uniform(0, 2)
z3 <- rexp(N, rate = 1)                              # z3i ∼ Exponential(1)
z4 <- rchisq(N, df = 4)                              # z4i ∼ χ^2(4)

# Calculer les variables x selon les relations données
x1 <- z1                                             # x1i = z1i
x2 <- z2 + 0.3 * x1                                  # x2i = z2i + 0.3 * x1i
x3 <- z3 + 0.2 * (x1 + x2)                           # x3i = z3i + 0.2 * (x1i + x2i)
x4 <- z4 + 0.1 * (x1 + x2 + x3)                      # x4i = z4i + 0.1 * (x1i + x2i + x3i)

# Calcul de Var(eta)
x_beta <- 2 + x1 + x2 + x3 + x4

# Définir la corrélation cible
rho_target <- 0.5

# Calcul de sigma
sigma <- sqrt(var(x_beta) * (1 / rho_target^2 - 1))
sigma

# Générer les erreurs aléatoires normales εi ∼ Normal(0, 1)
epsilon <- rnorm(N, mean = 0, sd = 1)

# Calculer la variable réponse y selon le modèle
y <- x_beta + sigma * epsilon         # yi = 2 + x1i + x2i + x3i + x4i + σεi
cor(y, x_beta)

# Mettre les données dans un data.frame pour plus de clarté
pop <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  x4 = x4
)

vrai_total <- round(sum(pop$y),1)
vrai_total

vrai_moyenne <- round(mean(pop$y),3)
vrai_moyenne

df <- data.frame(
  total_y_c = numeric(0), total_y_inc = numeric(0), total_y_naif = numeric(0),
  biais_relatif_total_c = numeric(0), biais_relatif_total_inc = numeric(0),
  biais_relatif_total_naif = numeric(0),
  moyenne_y_c = numeric(0), moyenne_y_inc = numeric(0), moyenne_y_naif = numeric(0), 
  biais_relatif_moyenne_c = numeric(0), biais_relatif_moyenne_inc = numeric(0), 
  biais_relatif_moyenne_naif = numeric(0)
)

n_simu <- 10000
n_prob <- 500
n_non_prob <- 500
theta1 <- 0.1
theta2 <- 0.2
theta3 <- 0.1
theta4 <- 0.2
a <- 10
nb_GHR <- 3

for (i in 1:n_simu){
  vec <- fonction_simulation(n_prob = n_prob, n_non_prob = n_non_prob,
                                theta1 = theta1, theta2 = theta2, 
                                theta3 = theta3, a = a, nb_GHR = nb_GHR)
  df <- rbind(df, data.frame(
    total_y_c = vec[1], total_y_inc = vec[2], total_y_naif = vec[3],
    biais_relatif_total_c = vec[4], biais_relatif_total_inc = vec[5],
    biais_relatif_total_naif = vec[6],
    moyenne_y_c = vec[7], moyenne_y_inc = vec[8], moyenne_y_naif = vec[9], 
    biais_relatif_moyenne_c = vec[10], biais_relatif_moyenne_inc = vec[11], 
    biais_relatif_moyenne_naif = vec[12]
  ))
}

df_long <- df %>% 
  select(biais_relatif_total_c, biais_relatif_total_inc, 
         biais_relatif_total_naif) %>% 
  pivot_longer(cols = c(biais_relatif_total_c, biais_relatif_total_inc, 
                        biais_relatif_total_naif),
               names_to = "variable", values_to = "value")

ggplot(df_long, aes(x="", y= value, fill = variable)) +
  geom_violin(adjust = 1L, scale = "area", width = 0.8) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  theme_bw() + scale_fill_hue(direction = 1) +
  facet_wrap(~ variable) +
  scale_fill_manual(values = c("biais_relatif_total_c" = "lightgreen", 
                               "biais_relatif_total_inc" = "#286AC7", 
                               "biais_relatif_total_naif" = "red"),
                    guide = "none") +
  labs(x = "", y = "")

df_long <- df %>% 
  select(biais_relatif_moyenne_c, biais_relatif_moyenne_inc, 
         biais_relatif_moyenne_naif) %>% 
  pivot_longer(cols = c(biais_relatif_moyenne_c, biais_relatif_moyenne_inc, 
                        biais_relatif_moyenne_naif),
               names_to = "variable", values_to = "value")

ggplot(df_long, aes(x="", y= value, fill = variable)) +
  geom_violin(adjust = 1L, scale = "area", width = 0.8) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  theme_bw() + scale_fill_hue(direction = 1) +
  facet_wrap(~ variable) +
  scale_fill_manual(values = c("biais_relatif_moyenne_c" = "lightgreen", 
                               "biais_relatif_moyenne_inc" = "#286AC7",
                               "biais_relatif_moyenne_naif" = "red"),
                    guide = "none") +
  labs(x = "", y = "")
