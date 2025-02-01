library(sampling)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(nppR)

rm(list=ls())
set.seed(1)

source(file = "parametre.R")
source(file = "fonction.R")

if (rho_target == 0.3){
  dossier_cor <- "cor_faible"
} else if (rho_target == 0.5){
  if (theta3 == 0.3){
    dossier_cor <- "cor_moyenne/effet_theta_faible/"
  } else if (theta3 == 0.6){
    dossier_cor <- "cor_moyenne/effet_theta_fort/"
  } 
} else if (rho_target == 0.8){
  dossier_cor <- "cor_elevee"
}

# dataframe contenant les résultats des simulations
df <- data.frame(
  biais_r_tot_prob = numeric(0), biais_r_tot_naif = numeric(0),
  biais_r_tot_frank_c = numeric(0), biais_r_tot_frank_inc = numeric(0),
  biais_r_tot_cart_c = numeric(0), biais_r_tot_cart_inc = numeric(0),
  biais_r_moy_frank_prob = numeric(0), biais_r_moy_frank_naif = numeric(0),
  biais_r_moy_cart_c = numeric(0), biais_r_moy_cart_inc = numeric(0)
)

# nombre de simulations
n_simu <- 100

for (i in 1:n_simu){
  vec <- fonction_simulation(n_prob = n_prob, n_non_prob = n_non_prob,
                             theta1 = theta1, theta2 = theta2, 
                             theta3 = theta3, a = a, nb_GHR = nb_GHR)
  
  df <- rbind(df, data.frame(
    biais_r_tot_prob = vec[1], biais_r_tot_naif = vec[2],
    biais_r_tot_frank_c = vec[3], biais_r_tot_frank_inc = vec[4],
    biais_r_tot_cart_c = vec[5], biais_r_tot_cart_inc = vec[6],
    biais_r_moy_prob = vec[7], biais_r_moy_naif = vec[8],
    biais_r_moy_frank_c = vec[9], biais_r_moy_frank_inc = vec[10],
    biais_r_moy_cart_c = vec[11], biais_r_moy_cart_inc = vec[12]
  ))
}

df_long <- df %>% 
  select(biais_r_tot_frank_c, biais_r_tot_frank_inc,
         biais_r_tot_cart_c, biais_r_tot_cart_inc,
         biais_r_tot_naif, biais_r_tot_prob) %>% 
  pivot_longer(cols = c(biais_r_tot_frank_c, 
                        biais_r_tot_frank_inc, 
                        biais_r_tot_cart_c, 
                        biais_r_tot_cart_inc, 
                        biais_r_tot_naif, 
                        biais_r_tot_prob),
               names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable, 
                           levels = c("biais_r_tot_prob", 
                                      "biais_r_tot_frank_c", 
                                      "biais_r_tot_cart_c",
                                      "biais_r_tot_frank_inc", 
                                      "biais_r_tot_cart_inc", 
                                      "biais_r_tot_naif")))

plot <- ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_violin(adjust = 1L, scale = "area", width = 0.8) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = 
                      c("biais_r_tot_prob" = "lightgreen",
                        "biais_r_tot_frank_c" = "#286AC7",
                        "biais_r_tot_cart_c" = "purple",
                        "biais_r_tot_frank_inc" = "lightblue",
                        "biais_r_tot_cart_inc" = "lightpink",
                        "biais_r_tot_naif" = "red"),
                    guide = "none") +  # Supprime la légende
  scale_x_discrete(labels = c("biais_r_tot_prob" = "Probabiliste",
                              "biais_r_tot_frank_c" = "Frank \n (complet)",
                              "biais_r_tot_cart_c" = "CART \n (complet)",
                              "biais_r_tot_frank_inc" = "Frank \n (incomplet)",
                              "biais_r_tot_cart_inc" = "CART \n (incomplet)",
                              "biais_r_tot_naif" = "Naif")) +
  labs(x = "", y = "Biais relatif des estimateurs du total de y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot
ggsave(paste0("png/",dossier_cor,"/biais_total.png"), 
       plot = plot, width = 8, height = 6, dpi = 300)


df_long <- df %>% 
  select(biais_r_moy_frank_c, biais_r_moy_frank_inc, 
         biais_r_moy_cart_c, biais_r_moy_cart_inc, 
         biais_r_moy_naif, biais_r_moy_prob) %>% 
  pivot_longer(cols = c(biais_r_moy_prob,
                        biais_r_moy_frank_c,
                        biais_r_moy_frank_inc,
                        biais_r_moy_cart_c,
                        biais_r_moy_cart_inc, 
                        biais_r_moy_naif),
               names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable, 
                           levels = c("biais_r_moy_prob", 
                                      "biais_r_moy_frank_c", 
                                      "biais_r_moy_cart_c", 
                                      "biais_r_moy_frank_inc",
                                      "biais_r_moy_cart_inc", 
                                      "biais_r_moy_naif")))

plot <- ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_violin(adjust = 1L, scale = "area", width = 0.8) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = 
                      c("biais_r_moy_prob" = "lightgreen",
                        "biais_r_moy_frank_c" = "#286AC7",
                        "biais_r_moy_cart_c" = "purple",
                        "biais_r_moy_frank_inc" = "lightblue",
                        "biais_r_moy_cart_inc" = "lightpink",
                        "biais_r_moy_naif" = "red"),
                    guide = "none") +  # Supprime la légende
  scale_x_discrete(labels = c("biais_r_moy_prob" = "Probabiliste",
                              "biais_r_moy_frank_c" = "Frank \n (complet)",
                              "biais_r_moy_cart_c" = "CART \n (complet)",
                              "biais_r_moy_frank_inc" = "Frank \n (incomplet)",
                              "biais_r_moy_cart_inc" = "CART \n (incomplet)",
                              "biais_r_moy_naif" = "Naif")) +
  labs(x = "", y = "Biais relatif des estimateurs de la moyenne de y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot
ggsave(paste0("png/",dossier_cor,"/biais_moyenne.png"), 
       plot = plot, width = 8, height = 6, dpi = 300)

df %>%
  summarise(across(everything(), 
                   list(mean = ~ mean(.), sd = ~ sd(.)))) %>%
  pivot_longer(cols = everything(), 
               names_to = c("variable", ".value"), 
               names_pattern = "^(.*)_(mean|sd)$")
