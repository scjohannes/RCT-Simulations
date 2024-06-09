plot_interim_or_black <- function(data, 
                                   filter_trial = NULL, 
                                   true_effect = NULL,
                                   event_rate = NULL, 
                                   title_suffix = "", 
                                   ylim_range = NULL,
                                   x_breaks = NULL) {
  library(tidyverse)
  library(ggplot2)
  
  p <- data %>%
    {if(!is.null(filter_trial)) filter(., trial == filter_trial) else .} %>% #. at very end returns data as is
    ggplot(aes(x = nPat, y = or, group = trial)) +
    geom_point() +
    geom_line() +
    {if(!is.null(true_effect)) geom_hline(yintercept = true_effect, color = "blue", linewidth = 1, linetype = "dashed") else NULL} +
    {if(!is.null(x_breaks)) scale_x_continuous(breaks = round(x_breaks)) else NULL} +
    theme(legend.position = "bottom") +
    {if(!is.null(ylim_range)) ylim(ylim_range) else NULL} +
    ylab("Odds Ratio") +
    xlab("# Patients Analyzed") +
    ggtitle(paste("Event-rate = ", event_rate, title_suffix, sep = ""))
  
  return(p)
}