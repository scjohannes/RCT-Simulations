## Johannes Schwenke - 2024-04-12
#Simulation to show that stratification will induce correlation

library(tidyverse)


nSims <- 1e3 # Number of Simulation
intercept <- 3 # intercept of outcome (ie, no treatment, predicitive variable not present)
beta_treat1 <- 0 # effect of treatment on outcome
beta_strat <- 3 # effect of stratification factor (predictive variable) on outcome
lambda0 <- 100 #expected number of people in stratum without predictive variable
lambda1 <- 100 #expected number of people in stratum with predictive variable

# NOTE: stratification will only cause correlation if the strata are not always perfectly of equal size
# if always exactly 50% of the people have or haven't the predictive factor, then this will not result in a correlation
# That's why we vary the number of patients in each stratum for every iteration of the simulation by drawing from a poisson distribution

# Explanation in my own words:
# mean_Y1 and mean_Y0 are correlated, because the prognostic factor is always perfectly balanced if we stratify.
# If by chance overall slightly more people have the prognostic factor (which creates a higher Y), then both Y0 and Y1 will be higher, and vice versa
# In case of simple randomisation, when overall more people have the prognostic factors, some more of then end up in treatment == 1, some more in treatment == 0, 
# and overall this cancels out and we don't observe a correlation.

#NOTE: this is for a continuous outcome. I assume this would work similarly for a binary outcome, on a logodds scale

# Stratified randomization ------------------------------------------------

nPatients_strat0 <- numeric()
nPatients_strat1 <- numeric()
nPatients <- numeric()
strat <- numeric()
mean_Y0 <- numeric(nSims)
mean_Y1 <- numeric(nSims)
obs_effect <- numeric(nSims)
pvalue <- numeric(nSims)

set.seed(1)


for (i in 1:nSims){
  nPatients_strat0 <- rpois(1, lambda0) #number of patients stratum = 1
  if (nPatients_strat0 %% 2 != 0){
    nPatients_strat0 <- nPatients_strat0 + 1
  }
  
  nPatients_strat1 <- rpois(1, lambda1) #number of patients stratum = 0
  if (nPatients_strat1 %% 2 != 0){
    nPatients_strat1 <- nPatients_strat1 + 1
  }
  
  nPatients <- nPatients_strat0 + nPatients_strat1
  
  strat <- c(rep(0, nPatients_strat0), rep(1, nPatients_strat1)) #create vector of all patients who either have the strat factor or don't
  
  treatment <- rep(0:1, nPatients/2)
  # this creates a vector of "treatment allocations" which is actually just a sequence alternating between 0 and 1. (Block of 2)
  # because this is a simulation, and using rep we always have perfectly balanced allocation, no need to acutally randomize according to stratum
  # For non stratified randomization, we have to use rbinom or resample the strat vector with replacement. 
  
  
  outcome <- numeric(nPatients) # this creates an empty vector which we will use to assign the outcome for each patient
  outcome <- intercept + beta_treat1*treatment + beta_strat * strat  # this assigns the probability of death for patients receiving 'treatment 0' 
  
  error <- rnorm(nPatients) #create a vector of error terms 
  outcome <- outcome + error #add the error terms to the outcome
  
  trialdat_strat <- tibble(outcome, treatment)
  
  mean_Y0[i] <- (trialdat_strat %>% group_by(treatment) %>% summarize(mean = mean(outcome)) %>% pull())[1]
  mean_Y1[i] <- (trialdat_strat %>% group_by(treatment) %>% summarize(mean = mean(outcome)) %>% pull())[2]
  
  obs_effect[i] <- mean_Y1[i] - mean_Y0[i]
  
  t <- t.test((trialdat_strat %>% filter(treatment == 1))$outcome, (trialdat_strat %>% filter(treatment == 0))$outcome, var.equal = FALSE, alternative = "greater")
  pvalue[i] <- t$p.value
  
  
}

effect_size_strat <- obs_effect
hist(effect_size_strat)

var_strat_effect <- var(obs_effect) #variance in effect estimate distribution is smaller in stratified randomization
var_strat_effect

pvalues_strat <- pvalue
mean(pvalues_strat <= 0.05)
hist(pvalues_strat) #pvalues are no longer uniformly distributed --> model misspecified

plot_strat <- tibble(mean_Y0, mean_Y1) %>% 
  ggplot(aes(x = mean_Y1, y = mean_Y0)) +
  geom_point() +
  ggtitle(paste("Stratified Randomisation"))

plot_strat

rm(mean_Y0, mean_Y1, nPatients_strat0, nPatients_strat1, nPatients, strat, obs_effect, pvalue) #remove unused variables
  


# Non stratified randomization --------------------------------------------
nPatients_strat0 <- numeric()
nPatients_strat1 <- numeric()
nPatients <- numeric()
strat <- numeric()
mean_Y0 <- numeric(nSims)
mean_Y1 <- numeric(nSims)
obs_effect <- numeric(nSims)
pvalue <- numeric(nSims)

set.seed(1) # this sets the random seed for your results to be reproducible


for (i in 1:nSims){
  
  nPatients_strat0 <- rpois(1, lambda0) #number of patients stratum = 1
  if (nPatients_strat0 %% 2 != 0){
    nPatients_strat0 <- nPatients_strat0 + 1
  }
  
  nPatients_strat1 <- rpois(1, lambda1) #number of patients stratum = 0
  if (nPatients_strat1 %% 2 != 0){
    nPatients_strat1 <- nPatients_strat1 + 1
  }
  
  # NOTE: stratification will only cause correlation if the strata are not of equal size.
  # if exactly 50% of the people have or haven't the predictive factor, then this will not results in a correlation
  
  nPatients <- nPatients_strat0 + nPatients_strat1
  
  strat <- c(rep(0, nPatients_strat0), rep(1, nPatients_strat1))
  strat <- sample(strat, length(strat), replace = FALSE) #here we just permute the vector, because otherwise by using rep below, we'd create stratified randomization
  
  treatment <- rep(0:1, nPatients/2)
  
  
  outcome <- numeric(nPatients) # this creates an empty vector which we will use to assign the outcome for each patient
  outcome <- intercept + beta_treat1*treatment + beta_strat * strat  # this assigns the probability of death for patients receiving 'treatment 0' 
  
  
  error <- rnorm(nPatients) #create a vector of error terms 
  outcome <- outcome + error #add the error terms to the outcome
  
  trialdat_non_strat <- tibble(outcome, treatment)
  
  mean_Y0[i] <- (trialdat_non_strat %>% group_by(treatment) %>% summarize(mean = mean(outcome)) %>% pull())[1]
  mean_Y1[i] <- (trialdat_non_strat %>% group_by(treatment) %>% summarize(mean = mean(outcome)) %>% pull())[2] 
  
  obs_effect[i] <- mean_Y1[i] - mean_Y0[i]
  
  
  t <- t.test((trialdat_non_strat %>% filter(treatment == 1))$outcome, (trialdat_non_strat %>% filter(treatment == 0))$outcome, var.equal = FALSE, alternative = "greater")
  pvalue[i] <- t$p.value
  
}

effect_size_simple <- obs_effect
hist(effect_size_simple)

var_simple_effect <- var(obs_effect) #variance in effect estimate distribution is smaller in stratified randomization
var_simple_effect

pvalues_simple <- pvalue
mean(pvalues_simple <= 0.05) #expected number of type I errors, if we set effect size to 0
hist(pvalues_simple)  #p values are uniformely distributed --> model is not misspecified


plot_non_strat <- tibble(mean_Y0, mean_Y1) %>% 
  ggplot(aes(x = mean_Y1, y = mean_Y0)) +
  geom_point() +
  ggtitle(paste("Simple Randomization"))


plot_non_strat
