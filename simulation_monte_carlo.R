library(sampling)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

rm(list=ls())
set.seed(1)

source(file = "parametre.R")
source(file = "fonction.R")

df <- data.frame(
  total_y_c = numeric(0), total_y_inc = numeric(0), total_y_naif = numeric(0),
  biais_relatif_total_c = numeric(0), biais_relatif_total_inc = numeric(0),
  biais_relatif_total_naif = numeric(0),
  moyenne_y_c = numeric(0), moyenne_y_inc = numeric(0), moyenne_y_naif = numeric(0), 
  biais_relatif_moyenne_c = numeric(0), biais_relatif_moyenne_inc = numeric(0), 
  biais_relatif_moyenne_naif = numeric(0)
)

n_simu <- 10000

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
