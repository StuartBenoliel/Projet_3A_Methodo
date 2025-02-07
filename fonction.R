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
  ech_non_prob <- getdata(pop, ech_non_prob) %>% 
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

U_complet <- function(alpha, data, tot_np){
  # Calcul de p_k(alpha) pour chaque k dans S_p
  p_k <- 1 / (1 + exp(-as.matrix(data[, c("x1", "x2", "x3")]) %*% alpha))
  
  # Calcul du second terme
  terme_2 <- colSums(data[, c("x1", "x2", "x3")] * 1/data$Prob * p_k)
  
  # Résidu de l'équation
  tot_np - terme_2
}

U_incomplet <- function(alpha, data, tot_np){
  # Calcul de p_k(alpha) pour chaque k dans S_p
  p_k <- 1 / (1 + exp(-as.matrix(data[, c("x1", "x2")]) %*% alpha))
  
  # Calcul du second terme
  terme_2 <- colSums(data[, c("x1", "x2")] * 1/data$Prob * p_k)
  
  # Résidu de l'équation
  tot_np[-3] - terme_2
}

fonction_simulation <- function(n_prob, n_non_prob, theta1, theta2, theta3,
                                a, nb_GHR){
  
  # Pour l'échantillon non probabiliste
  ech_non_prob <- tirage_non_proba()
  
  # Pour l'échantillon probabiliste
  # STSRS par rapport à x1
  type_tirage <- 'Stratifié'
  ech_prob <- tirage_proba(type = type_tirage)
  
  
  tot_np <- ech_non_prob %>% 
    summarise(tot_1 = sum(x1),
              tot_2 = sum(x2),
              tot_3 = sum(x3)) %>%
    as.numeric()
  
  alpha_complet <- newtonsys(U_complet, c(0, 0, 0), data = ech_prob, tot_np = tot_np)
  alpha_incomplet <- newtonsys(U_incomplet, c(0, 0), data = ech_prob, tot_np = tot_np)
  
  ech_prob <- ech_prob %>% 
    mutate(prob_participation_complet =
             1 / (1 + exp(-as.matrix(ech_prob[, c("x1", "x2", "x3")]) %*% 
                            alpha_complet$zero)),
           prob_participation_incomplet =
             1 / (1 + exp(-as.matrix(ech_prob[, c("x1", "x2")]) %*% 
                            alpha_incomplet$zero)),
           rang_c = rank(prob_participation_complet),
           rang_inc = rank(prob_participation_incomplet))
  
  ech_non_prob <- ech_non_prob %>% 
    mutate(prob_participation_complet =
             1 / (1 + exp(-as.matrix(ech_non_prob[, c("x1", "x2", "x3")]) %*% 
                            alpha_complet$zero)),
           prob_participation_incomplet =
             1 / (1 + exp(-as.matrix(ech_non_prob[, c("x1", "x2")]) %*% 
                            alpha_incomplet$zero)),
           rang_c = rank(prob_participation_complet),
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
  
  # On suppose 4 GHR
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
        prob_participation_complet >= group_limits_c$midpoint_next[3] ~ 4
      ),
      GHR_inc = case_when(
        prob_participation_incomplet < group_limits_inc$midpoint_next[1] ~ 1,
        prob_participation_incomplet < group_limits_inc$midpoint_next[2] ~ 2,
        prob_participation_incomplet < group_limits_inc$midpoint_next[3] ~ 3,
        prob_participation_incomplet >= group_limits_inc$midpoint_next[3] ~ 4
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
      produit_logit_c = y / prob_participation_complet,
      produit_logit_inc = y / prob_participation_incomplet,
      produit_frank_c = poids_frank_c * y,
      produit_frank_inc = poids_frank_inc * y,
      produit_cart_c = y / propensity_c,
      produit_cart_inc = y / propensity_inc
    )
  
  tot_prob <- sum(ech_prob$produit)
  tot_naif <- sum(ech_non_prob$y)*N/n_non_prob
  tot_logit_c <- sum(ech_non_prob$produit_logit_c)
  tot_logit_inc <- sum(ech_non_prob$produit_logit_inc)
  tot_frank_c <- sum(ech_non_prob$produit_frank_c)
  tot_frank_inc <- sum(ech_non_prob$produit_frank_inc)
  tot_cart_c <- sum(ech_non_prob$produit_cart_c)
  tot_cart_inc <- sum(ech_non_prob$produit_cart_inc)
  
  biais_r_tot_prob <- 100*(tot_prob - vrai_tot)/vrai_tot
  biais_r_tot_naif <- 100*(tot_naif - vrai_tot)/vrai_tot
  biais_r_tot_logit_c <- 100*(tot_logit_c - vrai_tot)/vrai_tot
  biais_r_tot_logit_inc <- 100*(tot_logit_inc - vrai_tot)/vrai_tot
  biais_r_tot_frank_c <- 100*(tot_frank_c - vrai_tot)/vrai_tot
  biais_r_tot_frank_inc <- 100*(tot_frank_inc - vrai_tot)/vrai_tot
  biais_r_tot_cart_c <- 100*(tot_cart_c - vrai_tot)/vrai_tot
  biais_r_tot_cart_inc <- 100*(tot_cart_inc - vrai_tot)/vrai_tot
  
  erreur_carre_tot_prob <- ((tot_prob - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_naif <- ((tot_naif - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_logit_c <- ((tot_logit_c - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_logit_inc <- ((tot_logit_inc - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_frank_c <- ((tot_frank_c - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_frank_inc <- ((tot_frank_inc - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_cart_c <- ((tot_cart_c - vrai_tot)/vrai_tot)^2
  erreur_carre_tot_cart_inc <- ((tot_cart_inc - vrai_tot)/vrai_tot)^2
  
  moy_prob <- sum(ech_prob$produit) / sum(1/ech_prob$Prob)
  moy_naif <- mean(ech_non_prob$y)
  moy_logit_c <- sum(ech_non_prob$produit_logit_c) / sum(ech_non_prob$prob_participation_complet)
  moy_logit_inc <- sum(ech_non_prob$produit_logit_inc) / sum(ech_non_prob$prob_participation_incomplet)
  moy_frank_c <- sum(ech_non_prob$produit_frank_c) / sum(ech_non_prob$poids_frank_c)
  moy_frank_inc <- sum(ech_non_prob$produit_frank_inc) / sum(ech_non_prob$poids_frank_inc)
  moy_cart_c <- sum(ech_non_prob$produit_cart_c) / sum(1/ech_non_prob$propensity_c)
  moy_cart_inc <- sum(ech_non_prob$produit_cart_inc) /  sum(1/ech_non_prob$propensity_inc)
  
  biais_r_moy_prob <- 100*(moy_prob - vrai_moy)/vrai_moy
  biais_r_moy_naif <- 100*(moy_naif - vrai_moy)/vrai_moy
  biais_r_moy_logit_c <- 100*(moy_logit_c - vrai_moy)/vrai_moy
  biais_r_moy_logit_inc <- 100*(moy_logit_inc - vrai_moy)/vrai_moy
  biais_r_moy_frank_c <- 100*(moy_frank_c - vrai_moy)/vrai_moy
  biais_r_moy_frank_inc <- 100*(moy_frank_inc - vrai_moy)/vrai_moy
  biais_r_moy_cart_c <- 100*(moy_cart_c - vrai_moy)/vrai_moy
  biais_r_moy_cart_inc <- 100*(moy_cart_inc - vrai_moy)/vrai_moy
  
  erreur_carre_moy_prob <- ((moy_prob - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_naif <- ((moy_naif - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_logit_c <- ((moy_logit_c - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_logit_inc <- ((moy_logit_inc - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_frank_c <- ((moy_frank_c - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_frank_inc <- ((moy_frank_inc - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_cart_c <- ((moy_cart_c - vrai_moy)/vrai_moy)^2
  erreur_carre_moy_cart_inc <- ((moy_cart_inc - vrai_moy)/vrai_moy)^2
  
  return(c(biais_r_tot_prob, biais_r_tot_naif,
           biais_r_tot_logit_c, biais_r_tot_logit_inc,
           biais_r_tot_frank_c, biais_r_tot_frank_inc,
           biais_r_tot_cart_c, biais_r_tot_cart_inc, 
           biais_r_moy_prob, biais_r_moy_naif,
           biais_r_moy_logit_c, biais_r_moy_logit_inc,
           biais_r_moy_frank_c, biais_r_moy_frank_inc,
           biais_r_moy_cart_c, biais_r_moy_cart_inc, 
           erreur_carre_tot_prob, erreur_carre_tot_naif,
           erreur_carre_tot_logit_c, erreur_carre_tot_logit_inc,
           erreur_carre_tot_frank_c, erreur_carre_tot_frank_inc,
           erreur_carre_tot_cart_c, erreur_carre_tot_cart_inc, 
           erreur_carre_moy_prob, erreur_carre_moy_naif,
           erreur_carre_moy_logit_c, erreur_carre_moy_logit_inc,
           erreur_carre_moy_frank_c, erreur_carre_moy_frank_inc,
           erreur_carre_moy_cart_c, erreur_carre_moy_cart_inc))
}