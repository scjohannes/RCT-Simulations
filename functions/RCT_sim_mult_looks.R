#Johannes Schwenke - 2024-06-08

RCT_sim_mult_looks <- function(nPatients = 1000, nSims = 20, seed = NULL, death0 = 0.4, death1 = 0.3, nLooks = 5, alpha = 0.05, typeOfDesign = "OF"){
  library(rpact)
  library(tidyverse)
  library(glue)
  
  # Documentation
  #' RCT Simulation with Multiple Looks
  #' built on the rpact package
  #'
  #' This function performs a simulation of a randomized controlled trial (RCT) with multiple interim analyses.
  #'
  #' @param nPatients Integer. The planned maximum number of patients to be included in each trial. Default is 1000.
  #' @param nSims Integer. The number of simulations to run. Default is 20.
  #' @param seed Integer. A seed for the random number generator. Default is NULL.
  #' @param death0 Numeric. The event rate for patients receiving 'treatment 0'. Default is 0.4.
  #' @param death1 Numeric. The event rate for patients receiving 'treatment 1'. Default is 0.3.
  #' @param nLooks Integer. The number of interim looks (including final analysis). Default is 5.
  #' @param alpha Numeric. The significance level. Default is 0.05.
  #' @param typeOfDesign Character. The type of group sequential design. Can be "OF", "asOF", "P", or "asP". Default is "OF".
  #' @return A list containing three dataframes: df_long, df_stopped_interim, df_critical_z.
  #'
  #' @examples
  #' result <- RCT_sim_mult_looks(nPatients = 1000, nSims = 20, death0 = 0.4, death1 = 0.35, nLooks = 5, alpha = 0.05, typeOfDesign = "OF")
  
  trueOR <- (death1/(1-death1))/(death0/(1-death0))
  
  analyses_scheduled <- seq(from = 1/nLooks, to = 1, by = 1/nLooks) # here is where you list the information fraction (e.g. 50%, 75% and 100% information)
  efficacy_thresholds <- numeric(nLooks)
  
  design <- getDesignGroupSequential(sided=1, alpha=alpha, informationRates=analyses_scheduled, typeOfDesign = typeOfDesign)
  
  rule <- ifelse(design$typeOfDesign == "OF", "Standard O’Brien & Fleming boundaries", 
                 ifelse(design$typeOfDesign == "asOF", "O’Brien & Fleming type alpha-spending", 
                        ifelse(design$typeOfDesign == "P", "Standard Pocock and Haybittle & Peto boundaries", 
                               ifelse(design$typeOfDesign == "asP", "Pocock type alpha-spending", 
                                      ))))
  
  for(j in 1:nLooks){ #get the alpha for each look
    efficacy_thresholds[j] = design$stageLevels[j]
  }
  
  analyses_nPatients <- analyses_scheduled*nPatients
  
  analyses_scheduled
  analyses_nPatients
  efficacy_thresholds
  
  #NOTE: uncomment the line above, if you want to visualize type I error inflation if you don't use stopping rules
  #First number is the alpha used for all looks
  #efficacy_thresholds <- rep(0.05, 5)
  
  # Simulation Parameters
  trialnum <- numeric(nSims)
  or <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  lcl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  ucl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  pvalue <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  success <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  zvalue <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  
  #provide column names
  colnames(or) <- sprintf("or_%d", (1:nLooks))
  colnames(lcl) <- sprintf("lcl_%d", (1:nLooks))
  colnames(ucl) <- sprintf("ucl_%d", (1:nLooks))
  colnames(zvalue) <- sprintf("zvalue_%d", (1:nLooks))
  colnames(pvalue) <- sprintf("pvalue_%d", (1:nLooks))
  colnames(success) <- sprintf("success_%d", (1:nLooks))
  
  overall_success <- numeric(nSims)
  
  
  # Simulation --------------------------------------------------------------
  set.seed(seed) # this sets the random seed for your results to be reproducible
  
  for(i in 1:nSims) {
    trialnum[i] = i
    
    pid <- seq(1, by = 1, len = nPatients)
    treatment <- rep(0:1, nPatients / 2)
    deathprob <- numeric(nPatients)
    deathprob[treatment == 0] <- death0
    deathprob[treatment == 1] <- death1
    death = rbinom(nPatients, 1, deathprob)
    trialdata <- data.frame(cbind(pid, treatment, death))
    
    for (j in 1:nLooks) {
      analysisdata <- subset(trialdata, pid <= analyses_nPatients[j])
      
      model <- glm(death ~ treatment,
                   family = binomial(link = 'logit'),
                   data = analysisdata)
      
      or[i, j] <- round(exp(summary(model)$coefficients[2]), digits = 2)
      
      lcl[i, j] <- round(exp(
        summary(model)$coefficients[2] - 1.96 * summary(model)$coefficients[4]
      ), digits = 2)
      
      ucl[i, j] <- round(exp(
        summary(model)$coefficients[2] + 1.96 * summary(model)$coefficients[4]
      ), digits = 2)
      
      pvalue[i, j] <- round(summary(model)$coefficients[8], digits = 4)
      
      zvalue[i, j] <- round(summary(model)$coefficients[6], digits = 4)
      
      success[i, j] <- ifelse(or[i, j] < 1 &
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
  
  #Build dataframe (wide format)
  simulation_results <- tibble(trialnum)
  
  for(i in 1:nLooks){
    new_columns <- tibble(
      or[i],
      lcl[i],
      ucl[i],
      zvalue[i],
      pvalue[i],
      success[i]
    )
    
    simulation_results <- bind_cols(simulation_results, new_columns)
  }
  
  simulation_results <- bind_cols(simulation_results, tibble(overall_success))
  
  simulation_results <- simulation_results %>%
    mutate("or_{i}" := or[i],
           "lcl_{i}" := lcl[i],
           "ucl_{i}" := ucl[i],
           "zvalue_{i}" := zvalue[i],
           "pvalue_{i}" := pvalue[i],
           "sucess_{i}" := success[i]
    )
  
  # Data wrangling / conversion to long format ------------------------------
  # Custom function to convert to long format
  convert_to_long <- function(data, prefix) {
    data %>%
      mutate(trial = row_number()) %>%
      pivot_longer(
        cols = starts_with(prefix), 
        names_to = c("look"), 
        values_to = paste(prefix)
      ) %>%
      mutate(look = as.numeric(gsub(paste0(prefix, "_"), "", look)))
  }
  
  #apply function
  or_long <- convert_to_long(or, "or")
  zvalue_long <- convert_to_long(zvalue, "zvalue")
  pvalue_long <- convert_to_long(pvalue, "pvalue")
  success_long <- convert_to_long(success, "success")
  
  overall_success_long <- tibble(overall_success = as.factor(rep(overall_success, each = nLooks)))
  
  #create one df_long
  
  df_long <- suppressMessages(right_join(or_long, zvalue_long))
  df_long <- suppressMessages(right_join(df_long, pvalue_long))
  df_long <- suppressMessages(right_join(df_long, success_long))
  df_long <- df_long %>% bind_cols(overall_success_long)
  df_not_stopped_interim <- df_long %>% bind_cols(tibble(nPat = rep(NA, nLooks*nSims)))
  
  for(i in 1:nLooks){
    df_not_stopped_interim <- df_not_stopped_interim %>%
      mutate(
        nPat = case_when(
          look == i ~ analyses_nPatients[i], 
          look != i ~ nPat
        )
      )
  }
  
  suppressWarnings(
  df_not_stopped_interim <- df_not_stopped_interim %>%
    group_by(trial) %>%
    mutate(first_success_index = min(which(success == 1))) %>%
    mutate(first_success_index = ifelse(first_success_index == Inf, nLooks, first_success_index))  %>%
    ungroup() %>%
    arrange(first_success_index) %>%
    mutate(first_success_index = factor(first_success_index)) %>%
    arrange(trial, look)
  )
  
  suppressWarnings(
  #create dataframe which simulates trial stop if interim results sign
  df_stopped_interim <- df_not_stopped_interim %>%
    group_by(trial) %>%
    mutate(first_success_index = min(which(success == 1))) %>%
    mutate(first_success_index = ifelse(first_success_index == Inf, nLooks, first_success_index))  %>%
    ungroup() %>%
    filter(look <= first_success_index) %>%
    arrange(first_success_index) %>%
    mutate(first_success_index = factor(first_success_index)) %>%
    arrange(trial, look)
  )
  
  #create table with the critical z values
  df_critical_z <- tibble(zvalue = qnorm(efficacy_thresholds/2), nPat = analyses_nPatients)
  
  return(list(df_not_stopped_interim = df_not_stopped_interim, df_stopped_interim = df_stopped_interim, df_critical_z = df_critical_z, trueOR = trueOR, analyses_nPatients = analyses_nPatients, design = rule, death0 = death0, death1 = death1))
}