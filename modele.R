library(sampling)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(pracma)
# install.packages(pkgs = c("devtools"))
# devtools::install_github("StatCan/nppR")
library(nppR)
library(ggrepel)

rm(list=ls())
set.seed(1)

source(file = "fonction.R")
source(file = "parametre.R")

pop <- genere_pop(rho = rho, N=N)

vrai_tot <- sum(pop$y)
vrai_tot

vrai_moy <- mean(pop$y)
vrai_moy

# Pour l'échantillon non probabiliste

# Tirage de Poisson avec pik selon un modèle logistique
# pi = 1 / (1 + exp(-(theta0 + theta1*x1 + theta2*x2 + theta3*x3)))
ech_non_prob <- tirage_non_proba()
summary(ech_non_prob$Prob)
max(ech_non_prob$Prob) / min(ech_non_prob$Prob)

# Pour l'échantillon probabiliste

# STSRS par rapport à x1
type_tirage <- 'Stratifié'
ech_prob <- tirage_proba(type = type_tirage)
summary(ech_prob$Prob)
max(ech_prob$Prob) / min(ech_prob$Prob)

tot_np <- ech_non_prob %>% 
  summarise(tot_1 = sum(x1),
            tot_2 = sum(x2),
            tot_3 = sum(x3)) %>%
  as.numeric()
tot_np

alpha_complet <- newtonsys(U_complet, c(0, 0, 0), 
                           data = ech_prob, tot_np = tot_np)
alpha_incomplet <- newtonsys(U_incomplet, c(0, 0), 
                             data = ech_prob, tot_np = tot_np)

ech_prob <- ech_prob %>% 
  mutate(prob_particip_c =
           1 / (1 + exp(-as.matrix(ech_prob[, c("x1", "x2", "x3")]) %*% 
                          alpha_complet$zero)),
         prob_particip_inc =
           1 / (1 + exp(-as.matrix(ech_prob[, c("x1", "x2")]) %*% 
                          alpha_incomplet$zero)),
         rang_c = rank(prob_particip_c),
         rang_inc = rank(prob_particip_inc))

ech_non_prob <- ech_non_prob %>% 
  mutate(prob_particip_c =
           1 / (1 + exp(-as.matrix(ech_non_prob[, c("x1", "x2", "x3")]) %*% 
                          alpha_complet$zero)),
         prob_particip_inc =
           1 / (1 + exp(-as.matrix(ech_non_prob[, c("x1", "x2")]) %*% 
                          alpha_incomplet$zero)),
         rang_c = rank(prob_particip_c),
         rang_inc = rank(prob_particip_inc))

# Méthode nppCART

cart_c <- nppCART(np.data = ech_non_prob, 
                  p.data = ech_prob %>% mutate(poids_sondage = 1/Prob), 
                  sampling.weight = 'poids_sondage', predictors = c("x1", "x2", "x3"))

cart_inc <- nppCART(np.data = ech_non_prob, 
                    p.data = ech_prob %>% mutate(poids_sondage = 1/Prob), 
                    sampling.weight = 'poids_sondage', predictors = c("x1", "x2"))

### Extract the nppCART-estimated propensities
ech_non_prob <- ech_non_prob %>% 
  mutate(propensity_c = cart_c$get_npdata_with_propensity()$propensity,
         propensity_inc = cart_inc$get_npdata_with_propensity()$propensity)

# On suppose 4 GHR

# Méthode de Frank : f(r_k) = log(1 + a*r_k / n_non_prob)
log_a <- log(1 + a)

ech_non_prob <- ech_non_prob %>% 
  mutate(
    frank_c = log(1 + a * rang_c / nrow(ech_non_prob)),
    frank_inc = log(1 + a * rang_inc / nrow(ech_non_prob))) %>% 
  mutate(
    GHR_c = cut(frank_c, 
                breaks = seq(0, log_a, length.out = nb_GHR + 1), 
                labels = 1:nb_GHR, 
                include.lowest = TRUE) %>% as.integer(),
    
    GHR_inc = cut(frank_inc, 
                  breaks = seq(0, log_a, length.out = nb_GHR + 1), 
                  labels = 1:nb_GHR, 
                  include.lowest = TRUE) %>% as.integer()
  )

table(ech_non_prob$GHR_c)
table(ech_non_prob$GHR_inc)

# Calcul des min, max et des points milieux dans ech_non_prob
group_limits_c <- ech_non_prob %>%
  group_by(GHR_c) %>%
  summarise(
    min_prob = min(prob_particip_c),
    max_prob = max(prob_particip_c)
  ) %>%
  mutate(
    pt = lead(min_prob),
    midpoint_next = (max_prob + lead(min_prob)) / 2
  )

group_limits_inc <- ech_non_prob %>%
  group_by(GHR_inc) %>%
  summarise(
    min_prob = min(prob_particip_inc),
    max_prob = max(prob_particip_inc)
  ) %>%
  mutate(
    pt = lead(min_prob),
    midpoint_next = (max_prob + lead(min_prob)) / 2
  )

ech_prob <- ech_prob %>%
  mutate(
    GHR_c = case_when(
      prob_particip_c < group_limits_c$midpoint_next[1] ~ 1,
      prob_particip_c < group_limits_c$midpoint_next[2] ~ 2,
      prob_particip_c < group_limits_c$midpoint_next[3] ~ 3,
      prob_particip_c >= group_limits_c$midpoint_next[3] ~ 4
    ),
    GHR_inc = case_when(
      prob_particip_inc < group_limits_inc$midpoint_next[1] ~ 1,
      prob_particip_inc < group_limits_inc$midpoint_next[2] ~ 2,
      prob_particip_inc < group_limits_inc$midpoint_next[3] ~ 3,
      prob_particip_inc >= group_limits_inc$midpoint_next[3] ~ 4
    ),
    produit = y / Prob
  )

table(ech_prob$GHR_c)
table(ech_prob$GHR_inc)

Ng_c <- ech_prob %>% group_by(GHR_c) %>% summarise(Ng = sum(1/Prob))
Ng_inc <- ech_prob %>% group_by(GHR_inc) %>% summarise(Ng = sum(1/Prob))

ech_non_prob <- ech_non_prob %>% 
  mutate(
    poids_frank_c = case_when(
      GHR_c == 1 ~ Ng_c$Ng[1] / nrow(ech_non_prob[ech_non_prob$GHR_c == 1, ]),
      GHR_c == 2 ~ Ng_c$Ng[2] / nrow(ech_non_prob[ech_non_prob$GHR_c == 2, ]),
      GHR_c == 3 ~ Ng_c$Ng[3] / nrow(ech_non_prob[ech_non_prob$GHR_c == 3, ]),
      GHR_c == 4 ~ Ng_c$Ng[4] / nrow(ech_non_prob[ech_non_prob$GHR_c == 4, ])
    ),
    poids_frank_inc = case_when(
      GHR_inc == 1 ~ Ng_inc$Ng[1] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 1, ]),
      GHR_inc == 2 ~ Ng_inc$Ng[2] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 2, ]),
      GHR_inc == 3 ~ Ng_inc$Ng[3] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 3, ]),
      GHR_inc == 4 ~ Ng_inc$Ng[4] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 4, ])
    )
  ) %>% 
  mutate(
    produit_logit_c = y / prob_particip_c,
    produit_logit_inc = y / prob_particip_inc,
    produit_frank_c = poids_frank_c * y,
    produit_frank_inc = poids_frank_inc * y,
    produit_cart_c = y / propensity_c,
    produit_cart_inc = y / propensity_inc
  )

table_poids <- ech_non_prob %>% 
  mutate(poids_particip_c = 1/prob_particip_c,
         poids_particip_inc = 1/prob_particip_inc,
         poids_cart_c = 1/propensity_c,
         poids_cart_inc = 1/propensity_inc) %>% 
  select(poids_particip_c, poids_particip_inc,
         poids_frank_c, poids_frank_inc, poids_cart_c, poids_cart_inc)


ech_prob$source <- "Probabiliste"
ech_non_prob$source <- "Non probabiliste"

# Fusionner les deux jeux de données
data <- rbind(ech_prob %>% 
                select(rang_c, rang_inc, prob_particip_c,
                       prob_particip_inc, GHR_c, GHR_inc, source),
              ech_non_prob %>% 
                select(rang_c, rang_inc, prob_particip_c,
                       prob_particip_inc, GHR_c, GHR_inc, source))

data <- data %>%
  arrange(prob_particip_c)

# Récupérer les indices pour les valeurs à annoter
indices <- c(1, round(nrow(data)/4), round(2*nrow(data)/4), 
             round(3*nrow(data)/4), nrow(data))

# Ajouter les labels arrondis à 2 chiffres significatifs
selected_values <- data[indices, ] %>%
  mutate(label = format(prob_particip_c, digits = 2, scientific = TRUE))

# Graphique avec annotations des valeurs sélectionnées
plot <- ggplot(data, aes(x = rang_c, y = prob_particip_c,  
                         color = as.factor(GHR_c), shape = source)) + 
  geom_point(alpha = 0.6, size = 2.5) +
  geom_text_repel(data = selected_values, aes(label = label), 
                  size = 4, show.legend = FALSE, color = "black",
                  nudge_y = 0.05,  # Décalage vertical pour éviter les chevauchements
                  box.padding = 4,  # Ajoute de l'espace autour du texte
                  max.overlaps = 10) +
  scale_y_log10() +  # Échelle logarithmique
  labs(title = "", 
       x = "Rang", 
       y = "Probabilité de participation estimée", 
       color = "GHR :", 
       shape = "Échantillon :") + 
  theme_minimal() + 
  scale_color_manual(values = c("1" = "#2674DD", 
                                "2" = "#42BB65", 
                                "3" = "#FFC300",  
                                "4" = "#E91422")) + 
  scale_shape_manual(values = c("Probabiliste" = 1, 
                                "Non probabiliste" = 4)) + 
  theme(legend.position = "top")

plot

ggsave(paste0("png/distrib_proba_selon_rang_c.png"), 
       plot = plot, width = 8, height = 6, dpi = 300, bg= "white")

data <- data %>%
  arrange(prob_particip_inc)

# Ajouter les labels arrondis à 2 chiffres significatifs
selected_values <- data[indices, ] %>%
  mutate(label = format(prob_particip_inc, digits = 2, scientific = TRUE))

plot <- ggplot(data, aes(x = rang_inc, y = prob_particip_inc, 
                         color = as.factor(GHR_inc), shape = source)) +
  geom_point(alpha = 0.6, size= 2.5) +
  geom_text_repel(data = selected_values, aes(label = label), 
                  size = 4, show.legend = FALSE, color = "black",
                  nudge_y = 0.05,  # Décalage vertical pour éviter les chevauchements
                  box.padding = 3.5,  # Ajoute de l'espace autour du texte
                  max.overlaps = 10) +
  scale_y_log10() +
  labs(title = "",
       x = "Rang",
       y = "Probabilité de participation estimée",
       color = "GHR :",
       shape = "Échantillon :") +
  theme_minimal() +
  scale_color_manual(values = c("1" = "#2674DD",
                                "2" = "#42BB65",
                                "3" = "#FFC300",  
                                "4" = "#E91422")) +
  scale_shape_manual(values = c("Probabiliste" = 1, 
                                "Non probabiliste" = 4)) + 
  theme(legend.position = "top")

plot
ggsave(paste0("png/distrib_proba_selon_rang_inc.png"), 
       plot = plot, width = 8, height = 6, dpi = 300, bg= "white")

table_deciles <- table_poids %>%
  reframe(across(everything(), list(decile = ~quantile(., probs = seq(0, 1, 0.1))))) %>%
  mutate(index = row_number() - 1) %>%  # Création d'index en premier
  mutate(across(where(is.numeric), ~round(., 1))) %>%  # Arrondi des colonnes numériques
  relocate(index)  # Remet index en première colonne
table_deciles 


# Calcul des estimateurs

tot_prob <- round(sum(ech_prob$produit), 1)
tot_naif <- round(sum(ech_non_prob$y)* N/n_non_prob, 1)

tot_logit_c <- round(sum(ech_non_prob$produit_logit_c),1)
tot_logit_inc <- round(sum(ech_non_prob$produit_logit_inc),1)
tot_frank_c <- round(sum(ech_non_prob$produit_frank_c),1)
tot_frank_inc <- round(sum(ech_non_prob$produit_frank_inc),1)
tot_cart_c <- round(sum(ech_non_prob$produit_cart_c),1)
tot_cart_inc <- round(sum(ech_non_prob$produit_cart_inc),1)

vrai_tot
tot_prob
tot_naif
tot_logit_c
tot_logit_inc
tot_frank_c
tot_frank_inc
tot_cart_c 
tot_cart_inc 

biais_r_tot_prob <- round(100*(tot_prob - vrai_tot)/vrai_tot,3)
biais_r_tot_naif <- round(100*(tot_naif - vrai_tot)/vrai_tot,3)

biais_r_tot_logit_c <- round(100*(tot_logit_c - vrai_tot)/vrai_tot,3)
biais_r_tot_logit_inc <- round(100*(tot_logit_inc - vrai_tot)/vrai_tot,3)
biais_r_tot_frank_c <- round(100*(tot_frank_c - vrai_tot)/vrai_tot,3)
biais_r_tot_frank_inc <- round(100*(tot_frank_inc - vrai_tot)/vrai_tot,3)
biais_r_tot_cart_c <- round(100*(tot_cart_c - vrai_tot)/vrai_tot,3)
biais_r_tot_cart_inc <- round(100*(tot_cart_inc - vrai_tot)/vrai_tot,3)

biais_r_tot_prob
biais_r_tot_naif

biais_r_tot_logit_c
biais_r_tot_logit_inc
biais_r_tot_frank_c
biais_r_tot_frank_inc
biais_r_tot_cart_c
biais_r_tot_cart_inc

moy_prob <- round(sum(ech_prob$produit) / sum(1/ech_prob$Prob),3)
moy_naif <- mean(ech_non_prob$y)

moy_logit_c <- 
  round(sum(ech_non_prob$produit_logit_c) / sum(ech_non_prob$prob_particip_c),3)
moy_logit_inc <- 
  round(sum(ech_non_prob$produit_logit_inc) / sum(ech_non_prob$prob_particip_inc),3)
moy_frank_c <- round(sum(ech_non_prob$produit_frank_c) / sum(ech_non_prob$poids_frank_c),3)
moy_frank_inc <- round(sum(ech_non_prob$produit_frank_inc) / sum(ech_non_prob$poids_frank_inc),3)
moy_cart_c <- round(sum(ech_non_prob$produit_cart_c) / sum(1/ech_non_prob$propensity_c),3)
moy_cart_inc <- round(sum(ech_non_prob$produit_cart_inc) /  sum(1/ech_non_prob$propensity_inc),3)

vrai_moy
moy_prob
moy_naif

moy_logit_c
moy_logit_inc
moy_frank_c
moy_frank_inc
moy_cart_c
moy_cart_inc

biais_r_moy_prob <- round(100*(moy_prob - vrai_moy)/vrai_moy,3)
biais_r_moy_naif <- round(100*(moy_naif - vrai_moy)/vrai_moy,3)

biais_r_moy_logit_c <- round(100*(moy_logit_c - vrai_moy)/vrai_moy,3)
biais_r_moy_logit_inc <- round(100*(moy_logit_inc - vrai_moy)/vrai_moy,3)
biais_r_moy_frank_c <- round(100*(moy_frank_c - vrai_moy)/vrai_moy,3)
biais_r_moy_frank_inc <- round(100*(moy_frank_inc - vrai_moy)/vrai_moy,3)
biais_r_moy_cart_c <- round(100*(moy_cart_c - vrai_moy)/vrai_moy,3)
biais_r_moy_cart_inc <- round(100*(moy_cart_inc - vrai_moy)/vrai_moy,3)

biais_r_moy_prob
biais_r_moy_naif

biais_r_moy_logit_c
biais_r_moy_logit_inc
biais_r_moy_frank_c
biais_r_moy_frank_inc
biais_r_moy_cart_c
biais_r_moy_cart_inc


