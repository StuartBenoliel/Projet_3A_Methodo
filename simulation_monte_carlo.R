library(sampling)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

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
  biais_relatif_total_prob = numeric(0), biais_relatif_total_naif = numeric(0),
  biais_relatif_total_c = numeric(0), biais_relatif_total_inc = numeric(0),
  biais_relatif_moyenne_prob = numeric(0), 
  biais_relatif_moyenne_naif = numeric(0),
  biais_relatif_moyenne_c = numeric(0), biais_relatif_moyenne_inc = numeric(0)
)

# nombre de simulations
n_simu <- 100

for (i in 1:n_simu){
  vec <- fonction_simulation(n_prob = n_prob, n_non_prob = n_non_prob,
                             theta1 = theta1, theta2 = theta2, 
                             theta3 = theta3, a = a, nb_GHR = nb_GHR)
  df <- rbind(df, data.frame(
    biais_relatif_total_prob = vec[1], biais_relatif_total_naif = vec[2],
    biais_relatif_total_c = vec[3], biais_relatif_total_inc = vec[4],
    biais_relatif_moyenne_prob = vec[5], biais_relatif_moyenne_naif = vec[6],
    biais_relatif_moyenne_c = vec[7], biais_relatif_moyenne_inc = vec[8]
  ))
}

df_long <- df %>% 
  select(biais_relatif_total_c, biais_relatif_total_inc, 
         biais_relatif_total_naif, biais_relatif_total_prob) %>% 
  pivot_longer(cols = c(biais_relatif_total_c, 
                        biais_relatif_total_inc, 
                        biais_relatif_total_naif, 
                        biais_relatif_total_prob),
               names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable, 
                           levels = c("biais_relatif_total_prob", 
                                      "biais_relatif_total_c", 
                                      "biais_relatif_total_inc", 
                                      "biais_relatif_total_naif")))

plot <- ggplot(df_long, aes(x=variable, y= value, fill = variable)) +
  geom_violin(adjust = 1L, scale = "area", width = 0.8) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = 
                      c("biais_relatif_total_prob" = "lightgreen", 
                        "biais_relatif_total_c" = "#286AC7",
                        "biais_relatif_total_inc" = "yellow",
                        "biais_relatif_total_naif" = "red"),
                    guide = "none") +
  labs(x = "", y = "")
plot
ggsave(paste0("png/",dossier_cor,"/biais_total.png"), 
       plot = plot, width = 8, height = 6, dpi = 300)


df_long <- df %>% 
  select(biais_relatif_moyenne_c, biais_relatif_moyenne_inc, 
         biais_relatif_moyenne_naif, biais_relatif_moyenne_prob) %>% 
  pivot_longer(cols = c(biais_relatif_moyenne_prob,
                        biais_relatif_moyenne_c,
                        biais_relatif_moyenne_inc, 
                        biais_relatif_moyenne_naif),
               names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable, 
                           levels = c("biais_relatif_moyenne_prob", 
                                      "biais_relatif_moyenne_c", 
                                      "biais_relatif_moyenne_inc", 
                                      "biais_relatif_moyenne_naif")))

plot <- ggplot(df_long, aes(x=variable, y= value, fill = variable)) +
  geom_violin(adjust = 1L, scale = "area", width = 0.8) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = 
                      c("biais_relatif_moyenne_prob" = "lightgreen",
                        "biais_relatif_moyenne_c" = "#286AC7",
                        "biais_relatif_moyenne_inc" = "yellow",
                        "biais_relatif_moyenne_naif" = "red"),
                    guide = "none") +
  labs(x = "", y = "")
plot
ggsave(paste0("png/",dossier_cor,"/biais_moyenne.png"), 
       plot = plot, width = 8, height = 6, dpi = 300)


df %>%
  summarise(across(everything(), 
                   list(mean = ~ mean(.), sd = ~ sd(.)))) %>%
  pivot_longer(cols = everything(), 
               names_to = c("variable", ".value"), 
               names_pattern = "^(.*)_(mean|sd)$")
