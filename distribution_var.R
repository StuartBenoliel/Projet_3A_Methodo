library(ggplot2)
library(dplyr)

rm(list = ls())
set.seed(1)
source(file = "parametre.R") 

# Définir plusieurs valeurs de rho
rho_values <- c(0.3, 0.5, 0.8)

# Générer les y pour chaque rho
data <- lapply(rho_values, function(rho) {
  sigma <- sqrt(var(x_beta) * (1 / rho^2 - 1))
  y <- x_beta + sigma * epsilon
  data.frame(y = y, rho = as.factor(rho))
}) %>% bind_rows()

# Visualisation avec ggplot2
plot <- ggplot(data, aes(x = y, fill = rho, color = rho)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "",
       x = "Variable d'intérêt y",
       y = "Densité",
       fill = expression(rho),
       color = expression(rho)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

plot
ggsave(paste0("png/distrib_y_selon_correlation.png"), 
       plot = plot, width = 8, height = 6, dpi = 300, bg= "white")

#########

# 0.6 ou 0.3
theta3_values <- c(0.3, 0.6)

sum_pik <- function(theta0, theta3) {
  eta <- theta0 + theta1*x1 + theta2*x2 + theta3*x3
  piA <- 1 / (1 + exp(-eta))
  sum(piA)
}

theta0_values <- sapply(theta3_values, function(theta3) {
  uniroot(function(theta) sum_pik(theta, theta3) - n_non_prob, 
          interval = c(-50, 50))$root
})

pop_long <- do.call(rbind, lapply(1:length(theta3_values), function(i) {
  theta0 <- theta0_values[i]
  theta3 <- theta3_values[i]
  
  pop$Prob <- 1 / (1 + exp(-(theta0 + theta1*x1 + theta2*x2 + theta3*x3)))
  pop$theta3 <- factor(theta3, levels = theta3_values)  # Ajouter une colonne pour distinguer les valeurs de theta3
  
  return(pop)
}))

library(gridExtra)
# Plage de données (limites de l'axe x)
x_limits <- c(0, 0.06)

# Calculer la largeur des bins en fonction de la plage des données
bin_width <- diff(x_limits) / 20  # Diviser l'intervalle par le nombre de bins souhaité

# Graphique pour theta3 = 0.3
plot1 <- ggplot(pop_long[pop_long$theta3 == 0.3, ], aes(x = Prob)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))),
                 binwidth = bin_width, fill = "#69b3a2", color = "white", alpha = 0.8) + 
  coord_cartesian(xlim = x_limits) + 
  labs(title = expression(theta[3] == 0.3), 
       x = "", 
       y = "Fréquence")+ 
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank())

# Graphique pour theta3 = 0.6
plot2 <- ggplot(pop_long[pop_long$theta3 == 0.6, ], aes(x = Prob)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))),
                 binwidth = bin_width, fill = "#69b3a2", color = "white", alpha = 0.8) + 
  coord_cartesian(xlim = x_limits) + 
  labs(title = expression(theta[3] == 0.6), 
       x = "Probabilité de participation à l'échantillon non probabiliste pour la population U",
       y = "Fréquence") + 
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank())

# Utiliser grid.arrange pour organiser les deux graphiques l'un en dessous de l'autre
plot <- grid.arrange(plot1, plot2, nrow = 2)
ggsave(paste0("png/distrib_proba_np.png"), 
       plot = plot, width = 8, height = 8, dpi = 300, bg= "white")
