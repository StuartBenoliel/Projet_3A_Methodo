library(ggplot2)
library(dplyr)

rm(list = ls())
set.seed(1)
source(file = "parametre.R") 
source(file = "fonction.R")
pop <- genere_pop(rho = rho, N=N)

# Définir plusieurs valeurs de rho
rho_values <- c(0.3, 0.5, 0.8)

# Générer les y pour chaque rho
data <- lapply(rho_values, function(rho) {
  z1 <- rbinom(N, size = 1, prob = 0.3)       
  z2 <- runif(N, min = 0, max = 2)                   
  z3 <- rexp(N, rate = 1)                             
  epsilon <- rnorm(N, mean = 0, sd = 1)                
  x1 <- z1                                             
  x2 <- z2 + 0.3 * x1
  x3 <- z3 + 0.2 * (x1 + x2)
  
  x_beta <- 2 + x1 + x2 + x3
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

# 0.1 ou 0.3
alpha2_values <- c(0.1, 0.3)

sum_pik_2 <- function(alpha_1, alpha_2) {
  eta <- alpha_1 + 0.1*pop$x1 + 0.2*pop$x2 + alpha_2*pop$x3
  piA <- 1 / (1 + exp(-eta))
  sum(piA)
}

alpha1_values <- sapply(alpha2_values, function(alpha_2) {
  uniroot(function(alpha_1) sum_pik_2(alpha_1, alpha_2) - n_non_prob, 
          interval = c(-50, 50))$root
})

pop_long <- do.call(rbind, lapply(1:length(alpha1_values), function(i) {
  alpha1 <- alpha1_values[i]
  alpha2 <- alpha2_values[i]
  
  pop$Prob <- 1 / (1 + exp(-(alpha1 + 0.1*pop$x1 + 0.2*pop$x2 + alpha2*pop$x3)))
  pop$alpha2 <- factor(alpha2, levels = alpha2_values)  # Ajouter une colonne pour distinguer les valeurs
  
  return(pop)
}))
summary(pop_long %>% filter(alpha2==0.3))
library(gridExtra)
# Plage de données (limites de l'axe x)
x_limits <- c(0, 0.06)

# Calculer la largeur des bins en fonction de la plage des données
bin_width <- diff(x_limits) / 20  # Diviser l'intervalle par le nombre de bins souhaité

# Graphique pour alpha2 = 0.1
plot1 <- ggplot(pop_long[pop_long$alpha2 == 0.1, ], aes(x = Prob)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))),
                 binwidth = bin_width, fill = "#69b3a2", color = "white", alpha = 0.8) + 
  coord_cartesian(xlim = x_limits) + 
  labs(title = expression(alpha2 == 0.1), 
       x = "", 
       y = "Fréquence")+ 
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank())

# Graphique pour alpha2 = 0.3
plot2 <- ggplot(pop_long[pop_long$alpha2 == 0.3, ], aes(x = Prob)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))),
                 binwidth = bin_width, fill = "#69b3a2", color = "white", alpha = 0.8) + 
  coord_cartesian(xlim = x_limits) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(title = expression(alpha2 == 0.3), 
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


#########
library(tidyr)
library(tibble)
X <- pop %>% 
  select(-y) %>% 
  cor() %>% 
  as_tibble(rownames = "Variable") %>%  # Alternative plus propre
  pivot_longer(-Variable, names_to = "nom", values_to = "Cor") %>% 
  mutate(Cor = round(Cor,2),
         nom = factor(nom, levels = colnames(pop)),
         Variable = factor(Variable, levels = rev(colnames(pop))))

plot <- ggplot(X, aes(y=Variable)) + 
  geom_tile(aes(x = nom, fill = Cor), color = "white") + 
  geom_text(aes(x = nom, label = Cor)) +
  scale_fill_gradientn(colors = c("#E4003A","#286AC7"),
                       n.breaks = 5) +
  scale_x_discrete("", expand = c(0,0)) +
  scale_y_discrete("", expand = c(0,0)) +
  labs(fill = "Corrélation de Pearson") +
  theme(legend.position = "top", axis.text.x = element_text(size = 8))
plot
ggsave(paste0("png/correlogramme.png"), 
       plot = plot, width = 8, height = 8, dpi = 300, bg= "white")
