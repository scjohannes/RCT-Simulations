# Load custom functions ----
source("functions/plot_interim_or_FstScs.R")
source("functions/plot_interim_or_OvScs.R")
source("functions/plot_interim_z.R")
source("functions/RCT_sim_mult_looks.R")
source("functions/plot_interim_or_black.R")


# HA true ----
## Sim 1 ----
sim1 <- RCT_sim_mult_looks(seed = 2, death0 = 0.4, death1 = 0.3, nSims = 100, nPatients = 1e3, nLooks = 2, typeOfDesign = "P")

p1 <- plot_interim_or_FstScs(data = sim1$df_stopped_interim, true_effect = sim1$trueOR, x_breaks = sim1$analyses_nPatients, design = sim1$design, event_rate = sim1$death0, title_suffix = "Stopped at interim")

stopped_trials <- sim1$df_not_stopped_interim |> 
  filter(first_success_index == 1) |> 
  select(trial) |>
  unique() |> 
  pull()


## Sim 2 ----
sim2 <- RCT_sim_mult_looks(seed = 2, death0 = 0.4, death1 = 0.3, nSims = 100, nPatients = 1e3, nLooks = 10, typeOfDesign = "P")

One_trial_df <- sim2$df_not_stopped_interim |> filter(trial == 1)

p2.1 <- plot_interim_or_black(data = One_trial_df, true_effect = NULL, x_breaks = sim2$analyses_nPatients, event_rate = sim2$death0, title_suffix = ", not stopped at interim", ylim_range = c(0,2.5))
p2.2 <- p2.1 + geom_hline(yintercept = sim1$trueOR, color = "blue", linetype = "dashed", linewidth = 1)

p2.3 <- plot_interim_or_black(data = sim2$df_not_stopped_interim, true_effect = sim2$trueOR, x_breaks = sim2$analyses_nPatients, event_rate = sim2$death0, title_suffix = ", not stopped at interim", ylim_range = c(0,2.5))

p2.4 <- p2.3 + 
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", linewidth = 1)

sim2_stop <- sim2$df_not_stopped_interim |> filter(trial %in% stopped_trials) |> filter(look <= 5)
sim2_not_stop <- sim2$df_not_stopped_interim |> filter(!trial %in% stopped_trials)

p2.5 <- ggplot() +
  geom_point(data = sim2_not_stop, aes(x = nPat, y = or, group = trial), color = "black") +
  geom_line(data = sim2_not_stop, aes(x = nPat, y = or, group = trial), color = "black") +
  geom_point(data = sim2_stop, aes(x = nPat, y = or, group = trial), color = "red") +
  geom_line(data = sim2_stop, aes(x = nPat, y = or, group = trial), color = "red") +
  geom_hline(yintercept = sim1$trueOR, color = "blue", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = round(sim2$analyses_nPatients)) +
  ylim(0, 2.5) +
  theme(legend.position = "bottom") +
  ylab("Odds Ratio") +
  xlab("# Patients Analyzed") +
  ggtitle(paste("Event-rate = ", sim2$death0, ", stopped at interim"))


#H0 true ----

## Sim 3 ---- Null effect
sim3 <- RCT_sim_mult_looks(seed = 2, death0 = 0.4, death1 = 0.4, nSims = 100, nPatients = 1e3, nLooks = 2, typeOfDesign = "P")
sim4 <- RCT_sim_mult_looks(seed = 2, death0 = 0.4, death1 = 0.4, nSims = 100, nPatients = 1e3, nLooks = 10, typeOfDesign = "P")

stopped_trials_2 <- sim3$df_not_stopped_interim |> 
  filter(first_success_index == 1) |> 
  select(trial) |>
  unique() |> 
  pull()


p3.1 <- plot_interim_or_black(data = sim4$df_not_stopped_interim, true_effect = sim4$trueOR, x_breaks = sim4$analyses_nPatients, event_rate = sim4$death0, title_suffix = ", not stopped at interim", ylim_range = c(0,2.5))

p3.2 <- p3.1 + 
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", linewidth = 1)

sim4_stop <- sim4$df_not_stopped_interim |> filter(trial %in% stopped_trials_2) |> filter(look <= 5)
sim4_not_stop <- sim4$df_not_stopped_interim |> filter(!trial %in% stopped_trials_2)

p3.3 <- ggplot() +
  geom_point(data = sim4_not_stop, aes(x = nPat, y = or, group = trial), color = "black") +
  geom_line(data = sim4_not_stop, aes(x = nPat, y = or, group = trial), color = "black") +
  geom_point(data = sim4_stop, aes(x = nPat, y = or, group = trial), color = "red") +
  geom_line(data = sim4_stop, aes(x = nPat, y = or, group = trial), color = "red") +
  geom_hline(yintercept = sim3$trueOR, color = "blue", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = round(sim4$analyses_nPatients)) +
  ylim(0, 2.5) +
  theme(legend.position = "bottom") +
  ylab("Odds Ratio") +
  xlab("# Patients Analyzed") +
  ggtitle(paste("Event-rate = ", sim3$death0, ", stopped at interim"))

#sim3$df_not_stopped_interim |> filter(overall_success == 1) |> select(trial) |> unique() |> count()
sim3_Error_1_look <- sim3$df_not_stopped_interim |> filter(look == 2) |> select(pvalue) |> pull()
sim3_Error_1_look <- mean(sim3_Error_1_look < 0.05)*100
