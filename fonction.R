# Fonction pour calculer la somme des πA pour un θ0 donné
sum_pik <- function(theta0) {
  eta <- theta0 + theta1*x1 + theta2*x2 + theta3*x3 + theta4*x4
  piA <- 1 / (1 + exp(-eta))
  sum(piA)
}


fonction_simulation <- function(n_prob, n_non_prob, theta1, theta2, theta3,
                                a, nb_GHR){
  
  
  # Pour l'échantillon non probabiliste
  # Résolution pour trouver θ0
  theta0 <- uniroot(function(theta) sum_pik(theta) - n_non_prob, interval = c(-50, 50))$root
  
  # Tirage de Poisson avec pik selon un modèle logistique
  pop$Prob <- 1 / (1 + exp(-(theta0 + theta1*x1 + theta2*x2 + theta3*x3 + theta4 * x4)))
  summary(pop$Prob)
  sum(pop$Prob)
  
  ech_non_prob <- sampling::UPpoisson(pop$Prob)
  ech_non_prob <- getdata(pop, ech_non_prob)  %>% 
    mutate(indic_participation = 1)
  summary(ech_non_prob$Prob)
  
  
  # Pour l'échantillon probabiliste
  
  #####
  # STSRS par rapport à x1
  Nh <- pop %>% group_by(x1) %>% summarise(Nh = n())
  Nh$alloc <- round(Nh$Nh*n_prob/N)
  sum(Nh$alloc)
  
  # Tirage de l'échantillon
  pop <- pop[order(pop$x1),]
  ech_prob <- sampling::strata(data = pop, stratanames = 'x1',
                               size = Nh$alloc, method = 'srswor')
  ech_prob <- getdata(pop, ech_prob) %>% 
    select(-Stratum) %>% 
    mutate(indic_participation = 0)
  summary(ech_prob$Prob)
  
  #####
  
  zi <- 0.45 + pop$x3 + 0.03*pop$y
  max(zi) / min(zi)
  
  pop$Prob <- n_prob*zi / sum(zi)
  summary(pop$Prob)
  
  ech_prob <- sampling::UPpoisson(pop$Prob)
  ech_prob <- getdata(pop, ech_prob)  %>% 
    mutate(indic_participation = 0)
  summary(ech_prob$Prob)
  
  ####
  
  data <- rbind(ech_prob, ech_non_prob)
  
  # model logistique
  modele_logistique_complet <- glm(indic_participation ~ x1 + x2 + x3, 
                                   data = data, 
                                   family = binomial)
  modele_logistique_incomplet <- glm(indic_participation ~ x1 + x2, 
                                     data = data, 
                                     family = binomial)
  
  data$prob_participation_complet <- predict(modele_logistique_complet, type = "response")
  data$prob_participation_incomplet <- predict(modele_logistique_incomplet, type = "response")
  
  ech_prob <- data %>% 
    filter(indic_participation == 0) %>% 
    mutate(rang_c = rank(prob_participation_complet),
           rang_inc = rank(prob_participation_incomplet))
  
  ech_non_prob <- data %>% 
    filter(indic_participation == 1) %>% 
    mutate(rang_c = rank(prob_participation_complet),
           rang_inc = rank(prob_participation_incomplet))
  
  # On suppose 3 GHR
  # Méthode de Frank : f(r_k) = log(1 + a*r_k / n_non_prob)
  log_a <- log(1 + a)
  
  ech_non_prob <- ech_non_prob %>% 
    mutate(frank_c = log(1 + a*rang_c/nrow(ech_non_prob)),
           frank_inc = log(1 + a*rang_inc/nrow(ech_non_prob))) %>% 
    mutate(
      GHR_c = case_when(
        frank_c <= log_a / nb_GHR ~ 1,
        frank_c <= 2 * log_a / nb_GHR ~ 2,
        frank_c <= 3 * log_a / nb_GHR ~ 3),
      GHR_inc = case_when(
        frank_inc <= log_a / nb_GHR ~ 1,
        frank_inc <= 2 * log_a / nb_GHR ~ 2,
        frank_inc <= 3 * log_a / nb_GHR ~ 3
      ))
  
  group_limits_c <- ech_non_prob %>%
    group_by(GHR_c) %>%
    summarise(
      min_prob = min(prob_participation_complet),
      max_prob = max(prob_participation_complet)
    ) %>%
    mutate(
      pt = lead(min_prob),
      midpoint_next = (max_prob + lead(min_prob)) / 2
    )
  
  group_limits_inc <- ech_non_prob %>%
    group_by(GHR_inc) %>%
    summarise(
      min_prob = min(prob_participation_incomplet),
      max_prob = max(prob_participation_incomplet)
    ) %>%
    mutate(
      pt = lead(min_prob),
      midpoint_next = (max_prob + lead(min_prob)) / 2
    )
  
  ech_prob <- ech_prob %>%
    mutate(
      GHR_c = case_when(
        prob_participation_complet < group_limits_c$midpoint_next[1] ~ 1,
        prob_participation_complet < group_limits_c$midpoint_next[2] ~ 2,
        prob_participation_complet >= group_limits_c$midpoint_next[2] ~ 3
      ),
      GHR_inc = case_when(
        prob_participation_incomplet < group_limits_inc$midpoint_next[1] ~ 1,
        prob_participation_incomplet < group_limits_inc$midpoint_next[2] ~ 2,
        prob_participation_incomplet >= group_limits_inc$midpoint_next[2] ~ 3
      )
    )
  
  Ng_c <- ech_prob %>% group_by(GHR_c) %>% summarise(Ng = sum(1/Prob))
  Ng_inc <- ech_prob %>% group_by(GHR_inc) %>% summarise(Ng = sum(1/Prob))
  
  ech_non_prob <- ech_non_prob %>% 
    mutate(
      poids_frank_c = case_when(
        GHR_c == 1 ~ Ng_c$Ng[1]/nrow(ech_non_prob[ech_non_prob$GHR_c==1, ]),
        GHR_c == 2 ~ Ng_c$Ng[2]/nrow(ech_non_prob[ech_non_prob$GHR_c==2, ]),
        GHR_c == 3 ~ Ng_c$Ng[3]/nrow(ech_non_prob[ech_non_prob$GHR_c==3, ]),
      ),
      poids_frank_inc = case_when(
        GHR_inc == 1 ~ Ng_inc$Ng[1]/nrow(ech_non_prob[ech_non_prob$GHR_inc==1, ]),
        GHR_inc == 2 ~ Ng_inc$Ng[2]/nrow(ech_non_prob[ech_non_prob$GHR_inc==2, ]),
        GHR_inc == 3 ~ Ng_inc$Ng[3]/nrow(ech_non_prob[ech_non_prob$GHR_inc==3, ]),
      )
    ) %>%
    mutate(produit_c = poids_frank_c * y,
           produit_inc = poids_frank_inc * y)
  
  # Calcul du total de y par la méthode de pondération inverse de Frank
  total_y_c <- round(sum(ech_non_prob$produit_c),1)
  total_y_inc <- round(sum(ech_non_prob$produit_inc),1)
  total_y_naif <- round(sum(ech_non_prob$y)*N/n_non_prob,1)
  
  biais_relatif_total_c <- round(100*(total_y_c - vrai_total)/vrai_total,3)
  biais_relatif_total_inc <- round(100*(total_y_inc - vrai_total)/vrai_total,3)
  biais_relatif_total_naif <- round(100*(total_y_naif - vrai_total)/vrai_total,3)
  
  moyenne_y_c <- round(sum(ech_non_prob$produit_c) / sum(ech_non_prob$poids_frank_c),3)
  moyenne_y_inc <- round(sum(ech_non_prob$produit_inc) / sum(ech_non_prob$poids_frank_inc),3)
  moyenne_y_naif <- mean(ech_non_prob$y)
  
  biais_relatif_moyenne_c <- round(100*(moyenne_y_c - vrai_moyenne)/vrai_moyenne,3)
  biais_relatif_moyenne_inc <- round(100*(moyenne_y_inc - vrai_moyenne)/vrai_moyenne,3)
  biais_relatif_moyenne_naif <- round(100*(moyenne_y_naif - vrai_moyenne)/vrai_moyenne,3)
  
  return(c(total_y_c , total_y_inc, total_y_naif, biais_relatif_total_c,
           biais_relatif_total_inc, biais_relatif_total_naif, 
           moyenne_y_c, moyenne_y_inc, moyenne_y_naif, biais_relatif_moyenne_c,
           biais_relatif_moyenne_inc, biais_relatif_moyenne_naif))
}