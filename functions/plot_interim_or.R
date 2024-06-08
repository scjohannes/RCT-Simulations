library(tidyverse)

plot_interim_or <- function(data, filter_trial = NULL, true_effect = NULL, x_breaks = NULL, design = NULL, event_rate = NULL, title_suffix = "", ylim_range = NULL) {
  p <- data %>%
    {if(!is.null(filter_trial)) filter(., trial == filter_trial) else .} %>% #. at very end returns data as is
    ggplot(aes(x = nPat, y = or, group = trial, color = overall_success)) +
    geom_point() +
    geom_line() +
    {if(!is.null(true_effect)) geom_hline(yintercept = true_effect, color = "blue", linewidth = 1) else NULL} +
    {if(!is.null(x_breaks)) scale_x_continuous(breaks = round(x_breaks)) else NULL} +
    theme(legend.position = "bottom") +
    {if(!is.null(ylim_range)) ylim(ylim_range) else NULL} +
    ylab("Odds ratio") +
    xlab("Number of Patients") +
    ggtitle(paste("Design : ", design, ", event-rate = ", event_rate, "\n", title_suffix, sep = ""))
  
  return(p)
}