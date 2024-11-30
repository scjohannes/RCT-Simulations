## Johannes Schwenke - 2024-11-30

library(tidyverse)

nSims <- 2e4 # Number of simulations
intercept <- 3 # intercept of outcome (ie, no treatment, predictive variable not present)
beta_treat1 <- 0 # effect of treatment on outcome
beta_pred <- 3 # effect of predictive variable on outcome
lambda0 <- 20 #expected number of people in stratum without predictive variable
lambda1 <- 20 #expected number of people in stratum with predictive variable

nPatients_strat0 <- numeric()
nPatients_strat1 <- numeric()
nPatients <- numeric()
strat <- numeric()
mean_Y0 <- numeric(nSims)
mean_Y1 <- numeric(nSims)
obs_effect <- numeric(nSims)
p_value_adj <- numeric(nSims) 
p_value_unadj <- numeric(nSims)
beta_1_unadj <- numeric(nSims)
beta_1_adj <- numeric(nSims)

set.seed(1)


for (i in 1:nSims){
  
  nPatients_strat0 <- rpois(1, lambda0) #number of patients without predictive factor
  if (nPatients_strat0 %% 2 != 0){
    nPatients_strat0 <- nPatients_strat0 + 1
  }
  
  nPatients_strat1 <- rpois(1, lambda1) #number of patients with predictive factor
  if (nPatients_strat1 %% 2 != 0){
    nPatients_strat1 <- nPatients_strat1 + 1
  }
  
  nPatients <- sum(nPatients_strat0, nPatients_strat1)
  
  strat <- c(rep(0, nPatients_strat0), rep(1, nPatients_strat1)) #whole patient population
  
  #here we just permute the vector, because otherwise by using rep below, we'd create stratified randomization
  strat <- sample(strat, length(strat), replace = FALSE) 
  
  
  treatment <- rep(0:1, nPatients/2) #1:1 randomization
  
  
  outcome <- numeric(nPatients)
  outcome <- intercept + beta_treat1*treatment + beta_pred * strat 
  
  
  error <- rnorm(nPatients) #add error from normal dist with mean = 1 and SD = 1
  outcome <- outcome + error
  
  trialdat_non_strat <- tibble(outcome, treatment, strat)
  
  m1 <- lm(outcome ~ treatment + strat, data = trialdat_non_strat)
  m2 <- lm(outcome ~ treatment, data = trialdat_non_strat)
  
  mean_Y0[i] <- (trialdat_non_strat |>  group_by(treatment) |>  summarize(mean = mean(outcome)) |>  pull())[1]
  mean_Y1[i] <- (trialdat_non_strat |>  group_by(treatment) |>  summarize(mean = mean(outcome)) |>  pull())[2] 
  
  obs_effect[i] <- mean_Y1[i] - mean_Y0[i]
  
  
  p_value_adj[i] <- summary(m1)$coefficients[2,4]
  p_value_unadj[i] <- summary(m2)$coefficients[2,4]
  
  beta_1_adj[i] <- summary(m1)$coefficients[2,1]
  beta_1_unadj[i] <- summary(m2)$coefficients[2,1]
  
}

results_df <- tibble(p_value_adj, p_value_unadj, beta_1_adj, beta_1_unadj, obs_effect)

results_df |> 
  summarize(mean(p_value_adj <= 0.05), mean(p_value_unadj <= 0.05))

results_df |> 
  pivot_longer(cols = c(p_value_adj, p_value_unadj), names_to = "p_value_type", values_to = "p_value") |>
  ggplot(aes(x = p_value, group = p_value_type, color = p_value_type, fill = p_value_type)) +
    geom_histogram(alpha = 0.25, binwidth = 0.05, boundary = 0, position = "identity") +
    labs(x = "p-value", 
         y = "Frequency",
         title = paste0("Distribution of p-values under H0 (nSims = ", nSims, ")")) +
  theme_light()

results_df |> 
  pivot_longer(cols = c(beta_1_adj, beta_1_unadj), names_to = "beta_1_type", values_to = "beta_1_value") |>
  ggplot(aes(x = beta_1_value, group = beta_1_type, color = beta_1_type, fill = beta_1_type)) +
    geom_histogram(alpha = 0.25, position = "identity") +
    labs(x = "beta 1", 
         y = "Frequency",
         title = "Distribution of beta 1 estimates under H0") +
  theme_light()
