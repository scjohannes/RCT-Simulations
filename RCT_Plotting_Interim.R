# Plotting ------------------------------------------------------------
## Set up palette ----
library(ggokabeito)

options(
  # set default colors in ggplot2 to colorblind-friendly
  # Okabe-Ito and Viridis palettes
  ggplot2.discrete.colour = ggokabeito::palette_okabe_ito(),
  ggplot2.discrete.fill = ggokabeito::palette_okabe_ito(),
  ggplot2.continuous.colour = "viridis",
  ggplot2.continuous.fill = "viridis",
  # set theme font and size
  book.base_family = "sans",
  book.base_size = 14
)

# set default theme
theme_set(
  theme_minimal(
    base_size = getOption("book.base_size"),
    base_family = getOption("book.base_family")
  ) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
)

## Load simulation ----
source("RCT_Binary_Outcome_Interim_Efficacy_Multiple_Looks.r")

## Load custom functions ----
source("functions/plot_interim_or.R")
source("functions/plot_interim_z.R")

## Plot ORs ----
#stopped at interim
plot1 <- plot_interim_or(data = df_stopped_interim, true_effect = trueOR, x_breaks = analyses_nPatients, design = design$typeOfDesign, event_rate = death0, title_suffix = "Stopped at interim")
#not stopped at interim
plot2 <- plot_interim_or(data = df_long, true_effect = trueOR, x_breaks = analyses_nPatients, design = design$typeOfDesign, event_rate = death0, title_suffix = "Not stopped at interim")

plot1
plot2

## Plot z values -------------------------------------------------------

plot3 <- plot_interim_z(data = df_stopped_interim, data_crit_z = df_critical_z, x_breaks = analyses_nPatients, design = design$typeOfDesign, event_rate = death0, title_suffix = "Stopped at interim")
plot4 <- plot_interim_z(data = df_long, data_crit_z = df_critical_z, x_breaks = analyses_nPatients, design = design$typeOfDesign, event_rate = death0, title_suffix = "Not stopped at interim")

plot3
plot4

# Histograms for various statistics ---------------------------------------

#Histogram of ORs according to look (and therefore increase in sample size per look)
#or_long %>% ggplot(aes(x = or)) + geom_histogram(binwidth = 0.1) + facet_wrap(~look, ncol = 2)

#Histogram of zvalues according to look (and therefore increase in sample size per look)
#zvalue_long %>% ggplot(aes(x = zvalue)) + geom_histogram(binwidth = 0.5) + facet_wrap(~look, ncol = 2)
