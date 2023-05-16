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
  dplyr::group_by(disp, pool.size, best.size) %>%
  dplyr::summarise(
    MeanPct = mean(pct.correct),
    MedianPct = median(pct.correct),
    SdPct = sd(pct.correct)
  )


full.plot.data <- summ %>% dplyr::filter(best.size %in% c(2, 4, 6, 8))

# Make the full plot
p1 <- ggplot(full.plot.data, aes(x = disp, y = MeanPct)) +
  # geom_hline(yintercept = 50) +
  geom_ribbon(aes(x = disp, ymin = MeanPct - SdPct, ymax = MeanPct + SdPct), fill = "darkgrey", alpha = 0.5) +
  geom_point() +
  labs(x = "Dispersal", y = "Mean percent correctly selected embryos") +
  coord_cartesian(ylim = c(0, 110)) +
  scale_y_continuous(
    breaks = seq(0, 150, 20),
    sec.axis = sec_axis(~., name = "Numner of selected embryos", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo pool size", breaks = NULL, labels = NULL)) +
  facet_grid(best.size ~ pool.size) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

save.double.width(p1, height = 170, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_xxxx_rank_pool"))


# Make a plot with just a single pool size
filt <- summ %>% dplyr::filter(pool.size == 6, best.size == 3)
p2 <- ggplot(filt, aes(x = disp, y = MeanPct)) +
  # geom_hline(yintercept = 50) +
  geom_ribbon(aes(x = disp, ymin = MeanPct - SdPct, ymax = MeanPct + SdPct), fill = "darkgrey", alpha = 0.5) +
  geom_point() +
  labs(x = "Dispersal", y = "Mean percent correctly selected embryos") +
  coord_cartesian(ylim = c(0, 110)) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  theme_bw()
save.single.width(p2, height = 85, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_xxxx_rank_pool_single"))
