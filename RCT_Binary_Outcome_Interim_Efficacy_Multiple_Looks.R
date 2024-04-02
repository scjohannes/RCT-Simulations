# Trial Design Parameters - Part 1
# Here we will specify the basics: maximum total number of patients to enroll and event rate for each treatment arm
nPatients <- 1000 # here is where you specify the planned max number of patients you want included in each RCT 
death0 <- 0.4 # here is where you specify the event rate for patients receiving 'treatment 0' in these trials
death1 <- 0.4 # here is where you specify the event rate for patients receiving 'treatment 1' in these trials
# I have set this one up to test the power for a treatment that would reduce mortality from 40% in control group (0) to 30% in treatment group (1)
# If one wants to estimate the "type 1 error" under different interim approaches, simply make 'death0' and 'death1' the same (no treatment effect)

# Trial Design Parameters - Part 2
# Here we will define the interim analysis strategy and stopping rules
# For this trial we will include provisions for efficacy stopping only (no futility stopping boundaries)
# We will use the rpact package to compute the stopping/success thresholds at interim and final analysis 
# install.packages("rpact")
library(rpact)
library(tidyverse)

theme_set(theme_light())

nLooks<-5 # here is where you put the number of looks that will take place (INCLUDING the final analysis)
analyses_scheduled<-seq(from = 0.2, to = 1, by = 0.2) # here is where you list the information fraction (e.g. 50%, 75% and 100% information)
efficacy_thresholds<-numeric(nLooks)

design <- getDesignGroupSequential(sided=1, alpha=0.05, informationRates=analyses_scheduled, typeOfDesign = "P")

for(j in 1:nLooks){
efficacy_thresholds[j] = design$stageLevels[j]
}
analyses_nPatients <- analyses_scheduled*nPatients

analyses_scheduled
analyses_nPatients
efficacy_thresholds

# Simulation Parameters
nSims <- 1000
trialnum <- numeric(nSims)
or <- data.frame(matrix(ncol = nLooks, nrow = nSims))
lcl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
ucl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
pvalue <- data.frame(matrix(ncol = nLooks, nrow = nSims))
success <- data.frame(matrix(ncol = nLooks, nrow = nSims))

#provide column names
colnames(or) <- sprintf("or_%d", (1:nLooks))
colnames(lcl) <- sprintf("lcl_%d", (1:nLooks))
colnames(ucl) <- sprintf("ucl_%d", (1:nLooks))
colnames(pvalue) <- sprintf("pvalue_%d", (1:nLooks))
colnames(success) <- sprintf("success_%d", (1:nLooks))

overall_success <- numeric(nSims)

set.seed(2) # this sets the random seed for your results to be reproducible

for(i in 1:nSims) {
  trialnum[i] = i
  
  pid = seq(1, by = 1, len = nPatients)
  treatment = rep(0:1, nPatients / 2)
  deathprob <- numeric(nPatients)
  deathprob[treatment == 0] = death0
  deathprob[treatment == 1] = death1
  death = rbinom(nPatients, 1, deathprob)
  trialdata = data.frame(cbind(pid, treatment, death))
  
  for (j in 1:nLooks) {
    analysisdata <- subset(trialdata, pid <= analyses_nPatients[j])
    
    model <- glm(death ~ treatment,
              family = binomial(link = 'logit'),
              data = analysisdata)
    
    or[i, j] = round(exp(summary(model)$coefficients[2]), digits = 2)
    
    lcl[i, j] = round(exp(
      summary(model)$coefficients[2] - 1.96 * summary(model)$coefficients[4]
    ), digits = 2)
    
    ucl[i, j] = round(exp(
      summary(model)$coefficients[2] + 1.96 * summary(model)$coefficients[4]
    ), digits = 2)
    
    pvalue[i, j] = round(summary(model)$coefficients[8], digits = 4)
    
    success[i, j] = ifelse(or[i, j] < 1 &
                             pvalue[i, j] < efficacy_thresholds[j], 1, 0)
  }
  
  overall_success[i] <- 0
  
  for (j in 1:nLooks)
  {
    if (success[i, j] == 1)
    {
      overall_success[i] <- 1
    }
  }
  
}

simulation_results <- data.frame(cbind(trialnum,
or[1], lcl[1], ucl[1], pvalue[1], success[1],
or[2], lcl[2], ucl[2], pvalue[2], success[2],
or[3], lcl[3], ucl[3], pvalue[3], success[3],
or[4], lcl[4], ucl[4], pvalue[4], success[4],
or[5], lcl[5], ucl[5], pvalue[5], success[5],
overall_success))
head(simulation_results, n=10)

table(overall_success)
table(simulation_results$success_1, overall_success)
table(simulation_results$success_2, overall_success)
table(simulation_results$success_3, overall_success)
table(simulation_results$success_4, overall_success)
table(simulation_results$success_5, overall_success)


#visualize plots

#convert everything to long format and one df
or_long <- or %>% 
  mutate(trial = row_number()) %>%
  pivot_longer(
    cols = starts_with("or_"), 
    names_to = c("look"), 
    values_to = "OR"
  ) %>%
  mutate(look = as.numeric(gsub("or_", "", look)))

pvalue_long <- pvalue %>% 
  mutate(trial = row_number()) %>%
  pivot_longer(
    cols = starts_with("pvalue_"), 
    names_to = c("look"), 
    values_to = "pvalue"
  ) %>%
  mutate(look = as.numeric(gsub("pvalue_", "", look)))

success_long <- success %>% 
  mutate(trial = row_number()) %>%
  pivot_longer(
    cols = starts_with("success_"), 
    names_to = c("look"), 
    values_to = "success"
  ) %>%
  mutate(look = as.numeric(gsub("success_", "", look)))

overall_success_long <- tibble(overall_success = as.factor(rep(overall_success, each = nLooks)))

df_long <- right_join(or_long, pvalue_long)
df_long <- right_join(df_long, success_long)

df_long <- df_long %>% 
  bind_cols(overall_success_long) %>%
  mutate(nPat = case_when(
    look == 1 ~ 200, 
    look == 2 ~ 400,
    look == 3 ~ 600, 
    look == 4 ~ 800,
    look == 5 ~ 1000
  ))


#create dataframe which simulates trial stop after first sig results

df_stopped_interim <- df_long %>%
  group_by(trial) %>%
  mutate(first_success_index = min(which(success == 1))) %>%
  mutate(first_success_index = ifelse(first_success_index == Inf, 5, first_success_index))  %>%
  filter(row_number() <= first_success_index) %>%
  ungroup()

#find the ORs for each look, which will result in a significant p-value
#this will depend on the stopping rule
df_critical_or <- df_stopped_interim %>% 
  filter(overall_success == 1) %>% 
  group_by(look, first_success_index) %>% 
  summarize(OR = max(OR)) %>% 
  ungroup() %>% 
  filter(as.numeric(look) == first_success_index) %>%
  select(OR)

df_critical_or <- tibble(OR = df_critical_or$OR, nPat = analyses_nPatients)



#data viz  

df_stopped_interim %>%
  filter(overall_success == 0, trial < 100) %>%
  ggplot() +
  geom_point(aes(x = nPat, y = OR, group = trial), color = "black") +
  geom_line(aes(x = nPat, y = OR, group = trial), color = "black") +
  geom_point(data = df_stopped_interim %>% filter(overall_success == 1, trial < 1000), aes(x = nPat, y = OR, group = trial), color = "red") +
  geom_line(data = df_stopped_interim %>% filter(overall_success == 1, trial < 1000), aes(x = nPat, y = OR, group = trial), color = "red") +
  geom_line(data = df_critical_or, aes(x = nPat, y = OR), color = "blue", linewidth = 1, linetype = "dashed") +
  theme(legend.position = "bottom") +
  #ylim(0, 2) +
  ggtitle(paste("Design : ", design$typeOfDesign, ", event-rate = ", death0, "\nStopped at interim", sep = ""))


df_long %>%
  filter(overall_success == 0, trial < 100) %>%
  ggplot() +
  geom_point(aes(x = nPat, y = OR, group = trial), color = "black") +
  geom_line(aes(x = nPat, y = OR, group = trial), color = "black") +
  geom_point(data = df_long %>% filter(overall_success == 1, trial < 1000), aes(x = nPat, y = OR, group = trial), color = "red") +
  geom_line(data = df_long %>% filter(overall_success == 1, trial < 1000), aes(x = nPat, y = OR, group = trial), color = "red") +
  geom_line(data = df_critical_or, aes(x = nPat, y = OR), color = "blue", linewidth = 1, linetype = "dashed") +
  theme(legend.position = "bottom") +
  ggtitle(paste("Design : ", design$typeOfDesign, ", event-rate = ", death0, "\nNot stopped at interim", sep = ""))

#pvalue distribution

hist(c(pvalue$pvalue_1, pvalue$pvalue_2, pvalue$pvalue_3, pvalue$pvalue_4, pvalue$pvalue_5))
hist(pvalue$pvalue_5)
