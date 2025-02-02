# Fonction pour calculer la somme des πA pour un θ0 donné
sum_pik <- function(theta0) {
  eta <- theta0 + theta1*x1 + theta2*x2 + theta3*x3
  piA <- 1 / (1 + exp(-eta))
  sum(piA)
}

tirage_non_proba <- function() {
  # Résolution pour trouver θ0
  theta0 <- uniroot(function(theta) sum_pik(theta) - n_non_prob, 
                    interval = c(-50, 50))$root
  
  # Tirage de Poisson avec pik selon un modèle logistique
  pop$Prob <- 1 / (1 + exp(-(theta0 + theta1*x1 + theta2*x2 + theta3*x3)))
  
  ech_non_prob <- sampling::UPpoisson(pop$Prob)
  ech_non_prob <- getdata(pop, ech_non_prob)  %>% 
    mutate(indic_participation = 1)
}

tirage_proba <- function(type) {
  if (type == 'Stratifié') {
    Nh <- pop %>% group_by(x1) %>% summarise(Nh = n())
    Nh$alloc <- round(Nh$Nh*n_prob/N)
    pop <- pop[order(pop$x1),]
    ech_prob <- sampling::strata(data = pop, stratanames = 'x1',
                                 size = Nh$alloc, method = 'srswor')
    ech_prob <- getdata(pop, ech_prob) %>% 
      select(-Stratum) %>% 
      mutate(indic_participation = 0)
    return(ech_prob)
  }
  if (type == 'Poisson') {
    zi <- 0.45 + pop$x3 + 0.03*pop$y
    max(zi) / min(zi)
    
    pop$Prob <- n_prob*zi / sum(zi)
    summary(pop$Prob)
    
    ech_prob <- sampling::UPpoisson(pop$Prob)
    ech_prob <- getdata(pop, ech_prob)  %>% 
      mutate(indic_participation = 0)
    return(ech_prob)
  }
}

fonction_simulation <- function(n_prob, n_non_prob, theta1, theta2, theta3,
                                a, nb_GHR){
  
  # Pour l'échantillon non probabiliste
  ech_non_prob <- tirage_non_proba()
  
  # Pour l'échantillon probabiliste
  # STSRS par rapport à x1
  type_tirage <- 'Stratifié'
  ech_prob <- tirage_proba(type = type_tirage)
  
  data <- rbind(ech_prob, ech_non_prob)
  
  modele_participation_complet <- glm(indic_participation ~ x1 + x2 + x3, 
                                      data = data, 
                                      family = binomial)
  
  modele_participation_incomplet <- glm(indic_participation ~ x1 + x2, 
                                        data = data, 
                                        family = binomial)
  
  data$prob_participation_complet <- predict(modele_participation_complet, type = "response")
  data$prob_participation_incomplet <- predict(modele_participation_incomplet, type = "response")
  
  ech_prob <- data %>% 
    filter(indic_participation == 0) %>% 
    mutate(rang_c = rank(prob_participation_complet),
           rang_inc = rank(prob_participation_incomplet))
  
  ech_non_prob <- data %>% 
    filter(indic_participation == 1) %>% 
    mutate(rang_c = rank(prob_participation_complet),
           rang_inc = rank(prob_participation_incomplet))
  
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
  
  # On suppose 10 GHR
  # Méthode de Frank : f(r_k) = log(1 + a*r_k / n_non_prob)
  log_a <- log(1 + a)
  
  ech_non_prob <- ech_non_prob %>% 
    mutate(
      frank_c = log(1 + a * rang_c / nrow(ech_non_prob)),
      frank_inc = log(1 + a * rang_inc / nrow(ech_non_prob))
    ) %>% 
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
        prob_participation_complet < group_limits_c$midpoint_next[3] ~ 3,
        prob_participation_complet < group_limits_c$midpoint_next[4] ~ 4,
        prob_participation_complet < group_limits_c$midpoint_next[5] ~ 5,
        prob_participation_complet < group_limits_c$midpoint_next[6] ~ 6,
        prob_participation_complet < group_limits_c$midpoint_next[7] ~ 7,
        prob_participation_complet < group_limits_c$midpoint_next[8] ~ 8,
        prob_participation_complet < group_limits_c$midpoint_next[9] ~ 9,
        prob_participation_complet >= group_limits_c$midpoint_next[9] ~ 10
      ),
      GHR_inc = case_when(
        prob_participation_incomplet < group_limits_inc$midpoint_next[1] ~ 1,
        prob_participation_incomplet < group_limits_inc$midpoint_next[2] ~ 2,
        prob_participation_incomplet < group_limits_inc$midpoint_next[3] ~ 3,
        prob_participation_incomplet < group_limits_inc$midpoint_next[4] ~ 4,
        prob_participation_incomplet < group_limits_inc$midpoint_next[5] ~ 5,
        prob_participation_incomplet < group_limits_inc$midpoint_next[6] ~ 6,
        prob_participation_incomplet < group_limits_inc$midpoint_next[7] ~ 7,
        prob_participation_incomplet < group_limits_inc$midpoint_next[8] ~ 8,
        prob_participation_incomplet < group_limits_inc$midpoint_next[9] ~ 9,
        prob_participation_incomplet >= group_limits_inc$midpoint_next[9] ~ 10
      ),
      produit = y / Prob
    )
  
  Ng_c <- ech_prob %>% group_by(GHR_c) %>% summarise(Ng = sum(1/Prob))
  Ng_inc <- ech_prob %>% group_by(GHR_inc) %>% summarise(Ng = sum(1/Prob))
  
  ech_non_prob <- ech_non_prob %>% 
    mutate(
      poids_frank_c = case_when(
        GHR_c == 1 ~ Ng_c$Ng[1] / nrow(ech_non_prob[ech_non_prob$GHR_c == 1, ]),
        GHR_c == 2 ~ Ng_c$Ng[2] / nrow(ech_non_prob[ech_non_prob$GHR_c == 2, ]),
        GHR_c == 3 ~ Ng_c$Ng[3] / nrow(ech_non_prob[ech_non_prob$GHR_c == 3, ]),
        GHR_c == 4 ~ Ng_c$Ng[4] / nrow(ech_non_prob[ech_non_prob$GHR_c == 4, ]),
        GHR_c == 5 ~ Ng_c$Ng[5] / nrow(ech_non_prob[ech_non_prob$GHR_c == 5, ]),
        GHR_c == 6 ~ Ng_c$Ng[6] / nrow(ech_non_prob[ech_non_prob$GHR_c == 6, ]),
        GHR_c == 7 ~ Ng_c$Ng[7] / nrow(ech_non_prob[ech_non_prob$GHR_c == 7, ]),
        GHR_c == 8 ~ Ng_c$Ng[8] / nrow(ech_non_prob[ech_non_prob$GHR_c == 8, ]),
        GHR_c == 9 ~ Ng_c$Ng[9] / nrow(ech_non_prob[ech_non_prob$GHR_c == 9, ]),
        GHR_c == 10 ~ Ng_c$Ng[10] / nrow(ech_non_prob[ech_non_prob$GHR_c == 10, ])
      ),
      
      poids_frank_inc = case_when(
        GHR_inc == 1 ~ Ng_inc$Ng[1] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 1, ]),
        GHR_inc == 2 ~ Ng_inc$Ng[2] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 2, ]),
        GHR_inc == 3 ~ Ng_inc$Ng[3] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 3, ]),
        GHR_inc == 4 ~ Ng_inc$Ng[4] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 4, ]),
        GHR_inc == 5 ~ Ng_inc$Ng[5] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 5, ]),
        GHR_inc == 6 ~ Ng_inc$Ng[6] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 6, ]),
        GHR_inc == 7 ~ Ng_inc$Ng[7] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 7, ]),
        GHR_inc == 8 ~ Ng_inc$Ng[8] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 8, ]),
        GHR_inc == 9 ~ Ng_inc$Ng[9] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 9, ]),
        GHR_inc == 10 ~ Ng_inc$Ng[10] / nrow(ech_non_prob[ech_non_prob$GHR_inc == 10, ])
      )
    ) %>% 
    mutate(
      produit_frank_c = poids_frank_c * y,
      produit_frank_inc = poids_frank_inc * y,
      produit_cart_c = y / propensity_c,
      produit_cart_inc = y / propensity_inc
    )
  
  tot_prob <- round(sum(ech_prob$produit),1)
  tot_naif <- round(sum(ech_non_prob$y)*N/n_non_prob,1)
  tot_frank_c <- round(sum(ech_non_prob$produit_frank_c),1)
  tot_frank_inc <- round(sum(ech_non_prob$produit_frank_inc),1)
  tot_cart_c <- round(sum(ech_non_prob$produit_cart_c),1)
  tot_cart_inc <- round(sum(ech_non_prob$produit_cart_inc),1)
  
  biais_r_tot_prob <- round(100*(tot_prob - vrai_tot)/vrai_tot,3)
  biais_r_tot_naif <- round(100*(tot_naif - vrai_tot)/vrai_tot,3)
  biais_r_tot_frank_c <- round(100*(tot_frank_c - vrai_tot)/vrai_tot,3)
  biais_r_tot_frank_inc <- round(100*(tot_frank_inc - vrai_tot)/vrai_tot,3)
  biais_r_tot_cart_c <- round(100*(tot_cart_c - vrai_tot)/vrai_tot,3)
  biais_r_tot_cart_inc <- round(100*(tot_cart_inc - vrai_tot)/vrai_tot,3)
  
  moy_prob <- round(sum(ech_prob$produit) / sum(1/ech_prob$Prob),3)
  moy_naif <- mean(ech_non_prob$y)
  moy_frank_c <- round(sum(ech_non_prob$produit_frank_c) / sum(ech_non_prob$poids_frank_c),3)
  moy_frank_inc <- round(sum(ech_non_prob$produit_frank_inc) / sum(ech_non_prob$poids_frank_inc),3)
  moy_cart_c <- round(sum(ech_non_prob$produit_cart_c) / sum(1/ech_non_prob$propensity_c),3)
  moy_cart_inc <- round(sum(ech_non_prob$produit_cart_inc) /  sum(1/ech_non_prob$propensity_inc),3)
  
  biais_r_moy_prob <- round(100*(moy_prob - vrai_moy)/vrai_moy,3)
  biais_r_moy_naif <- round(100*(moy_naif - vrai_moy)/vrai_moy,3)
  biais_r_moy_frank_c <- round(100*(moy_frank_c - vrai_moy)/vrai_moy,3)
  biais_r_moy_frank_inc <- round(100*(moy_frank_inc - vrai_moy)/vrai_moy,3)
  biais_r_moy_cart_c <- round(100*(moy_cart_c - vrai_moy)/vrai_moy,3)
  biais_r_moy_cart_inc <- round(100*(moy_cart_inc - vrai_moy)/vrai_moy,3)
  
  return(c(biais_r_tot_prob, biais_r_tot_naif,
           biais_r_tot_frank_c, biais_r_tot_frank_inc,
           biais_r_tot_cart_c, biais_r_tot_cart_inc, 
           biais_r_moy_prob, biais_r_moy_naif,
           biais_r_moy_frank_c, biais_r_moy_frank_inc,
           biais_r_moy_cart_c, biais_r_moy_cart_inc))
}