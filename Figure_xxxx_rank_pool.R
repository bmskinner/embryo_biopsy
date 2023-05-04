# Create the rank pool figure

# Aim is to shows that biopsies are useful for ranking embryos from a pool,
# not just two embryos at a time.

# Data generated from rankOrder.R

source("parameters.R")
source("functions.R")

out.file <- "data/Rank_pool_results.csv"
pool.results <- read.csv(out.file, header = T)

# Summarise the raw data
summ <- pool.results %>%
  dplyr::group_by(disp) %>%
  dplyr::summarise(
    MeanPct = mean(pct.correct),
    MedianPct = median(pct.correct),
    SdPct = sd(pct.correct)
  )

# Make the plot
ggplot(summ, aes(x = disp, y = MeanPct, group = disp)) +
  geom_hline(yintercept = 50) +
  geom_ribbon(aes(ymin = MeanPct - SdPct, ymax = MeanPct + SdPct)) +
  geom_errorbar(aes(ymin = MeanPct - SdPct, ymax = MeanPct + SdPct)) +
  geom_point() +
  labs(x = "Dispersal", y = "Mean percent correctly selected embryos") +
  coord_cartesian(ylim = c(0, 110)) +
  scale_y_continuous(breaks = seq(0, 110, 10)) +
  theme_bw()
