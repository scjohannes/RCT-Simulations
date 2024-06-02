library(ggplot2)
library(gifski)
library(tidyverse)

# set default theme
theme_set(theme_light())

set.seed(929)

n <- 200 # total number of datapoints (per condition) after initial 10
d <- 0.0 # effect size d

p <- numeric(n) # store p-values
x <- numeric(n) # store x-values
y <- numeric(n) # store y-values

n <- n + 10 # add 10 to number of datapoints

for (i in 11:n) { # for each simulated participants after the first 10
  x[i] <- rnorm(n = 1, mean = 0, sd = 1)
  y[i] <- rnorm(n = 1, mean = d, sd = 1)
  p[i] <- t.test(x[1:i], y[1:i], var.equal = TRUE)$p.value
}

p <- tibble(p = p[11:n], n = as.numeric(1:(n-10))) # Remove first 10 empty p-values

filenames <- character(n-10)

# Create the plot
for(i in 1:(n-10)){
  filenames[i] <- tempfile(fileext = ".png")
  png(filenames[i], width = 800, height = 500, units = "px")
  plot <- p %>%
    filter(n <= i) %>%
    ggplot() +
    geom_line(mapping = aes(x = n, y = p)) +
    geom_hline(yintercept = 0.05, linewidth = 1, linetype = "dashed", color = "darkgrey") +
    scale_x_continuous(limits = c(0, n-10)) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = "Sample Size", y = "P-value")
  print(plot)
  dev.off()
}

# save last image for PDF
file.copy(filenames[i], "images/animatep.png", overwrite = TRUE)

# create gif
gifski::gifski(png_files = filenames,  
               gif_file = "images/animatep.gif", 
               delay = 0.01)

min(p$p) # Return lowest p-value from all looks
cat("The lowest p-value was observed at sample size", which.min(p$p) + 10) 
cat("The p-value dropped below 0.05 for the first time at sample size:", 
    ifelse(is.na(which(p < 0.05)[1] + 10), "NEVER", which(p < 0.05)[1] + 10)) 
