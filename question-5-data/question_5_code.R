# Install and load packages
install.packages(c("janitor", "ggplot2", "dplyr", "magrittr"))
library(janitor)
library(ggplot2)
library(dplyr)
library(magrittr)

# Load in the data
dsvirus <- read.csv("Cui_etal2014.csv")

# View the data
names(dsvirus)
head(dsvirus)

# Clean the data
clean_dsvirus <- dsvirus %>%
  clean_names() %>%
  select(c("genome_length_kb", "virion_volume_nm_nm_nm"))

# Log transform the data
clean_dsvirus$log_genome_length <- log(clean_dsvirus$genome_length_kb)
clean_dsvirus$log_virion_volume <- log(clean_dsvirus$virion_volume_nm_nm_nm)

# Plot the data
plot_dsvirus <- ggplot(clean_dsvirus, aes(x = log_genome_length, y = log_virion_volume)) +
    geom_point(color = "black",
               size = 2,
               alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(x = "log[Genome length(kb)])",
         y = "log[Virion volume (nm³)]",) +
    theme_bw()
print(plot_dsvirus)

# Make a linear model
dsvirus_lm <- lm(log_virion_volume ~ log_genome_length, data = clean_dsvirus)

# Print the results
summary(dsvirus_lm)
coef(dsvirus_lm)

exp(7.074800)
1181.807 * 300^(1.515228)
