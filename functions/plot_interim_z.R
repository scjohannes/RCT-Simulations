library(tidyverse)

plot_interim_z <- function(data_trials, data_crit_z, filter_trial = NULL, x_breaks = NULL, design = NULL, event_rate = NULL, title_suffix = "", ylim_range = NULL) {
  p <- data_trials %>%
    {if (!is.null(filter_trial)) filter(., trial == filter_trial) else .} %>%
    ggplot() +
    geom_point(aes(x = nPat, y = zvalue, group = trial, color = overall_success)) +
    geom_line(aes(x = nPat, y = zvalue, group = trial, color = overall_success)) +
    geom_line(data = data_crit_z, aes(x = nPat, y = zvalue), color = "blue", linewidth = 1, linetype = "dashed") +
    {if(!is.null(x_breaks)) scale_x_continuous(breaks = round(x_breaks)) else NULL} +
    theme(legend.position = "bottom") +
    {if(!is.null(ylim_range)) ylim(ylim_range) else NULL} +
    ylab("Z-value") +
    xlab("Number of Patients") +
    ggtitle(paste("Design : ", design, ", event-rate = ", event_rate, "\n", title_suffix, sep = ""))
  
  return(p)
}