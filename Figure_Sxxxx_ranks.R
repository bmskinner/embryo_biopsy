# Create supplementary rank plots
source("parameters.R")
source("functions.R")

# # Read the saved raw values
output <- read.csv("Rank_results.csv", header = T)

# Aggregate results by difference in aneuploidy and calculate means, SDs
filt <- output %>%
  # dplyr::filter(Low.disp == 0 & High.disp==0) %>%
  dplyr::group_by(Aneu.diff, Low.disp, High.disp) %>%
  dplyr::summarise(
    Mean.correct = mean(Correct.rank),
    SD.correct = sd(Correct.rank),
    Mean.incorrect = mean(Incorrect.rank),
    SD.incorrect = sd(Incorrect.rank),
    Mean.no.rank = mean(No.rank),
    SD.no.rank = sd(No.rank),
    Mean.adj.correct = mean(Adj.correct.rank),
    SD.adj.correct = sd(Adj.correct.rank),
    Mean.adj.incorrect = mean(Adj.incorrect.rank),
    SD.adj.incorrect = sd(Adj.incorrect.rank)
  )


# Show correct, incorrect and no ranks on single plot
out.combined.plot <- ggplot(filt, aes(x = Aneu.diff * 100)) +
  annotate("rect", xmin = 0, xmax = 20, ymin = 0, ymax = Inf, fill = "lightgray") +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct, col = "Correct", shape = "Correct")) +
  geom_errorbar(aes(
    ymin = Mean.correct - SD.correct,
    ymax = Mean.correct + SD.correct, col = "Correct"
  ), size = 0.5) +
  geom_point(aes(y = Mean.incorrect, col = "Incorrect", shape = "Incorrect")) +
  geom_errorbar(aes(
    ymin = Mean.incorrect - SD.incorrect,
    ymax = Mean.incorrect + SD.incorrect, col = "Incorrect"
  ), size = 0.5) +
  geom_point(aes(y = Mean.no.rank, col = "Tied", shape = "Tied")) +
  geom_errorbar(aes(
    ymin = Mean.no.rank - SD.no.rank,
    ymax = Mean.no.rank + SD.no.rank, col = "Tied"
  ), size = 0.5) +
  scale_colour_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red")
  ) +
  scale_shape_manual(
    name = element_blank(),
    values = c(Correct = 16, Tied = 15, Incorrect = 17)
  ) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Aneuploidy difference (%)", y = "Percent embryos with rank order...") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = c(0.9, 0.83),
    legend.background = element_blank()
  )
save.double.width(out.combined.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure xxxx - Ranks_all"), 170)

# Split the values for no.rank between the correct and incorrect
# Can't get error bars direct from this - recalculate from original values
out.split.plot <- ggplot(filt, aes(x = Aneu.diff * 100)) +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.adj.incorrect, col = "Incorrect", shape = "Incorrect")) +
  geom_errorbar(aes(
    ymin = Mean.adj.incorrect - SD.adj.incorrect,
    ymax = Mean.adj.incorrect + SD.adj.incorrect, col = "Incorrect"
  ), size = 0.5) +
  geom_point(aes(y = Mean.adj.correct, col = "Correct", shape = "Correct")) +
  geom_errorbar(aes(
    ymin = Mean.adj.correct - SD.adj.correct,
    ymax = Mean.adj.correct + SD.adj.correct, col = "Correct"
  ), size = 0.5) +
  scale_colour_manual(
    name = element_blank(),
    values = c(Correct = "blue", Incorrect = "red")
  ) +
  scale_shape_manual(
    name = element_blank(),
    values = c(Correct = 16, Incorrect = 17)
  ) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Aneuploidy difference (%)", y = "Effective percent embryos with rank order...") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = c(0.9, 0.85),
    legend.background = element_blank()
  )
save.double.width(out.split.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure xxxx - Ranks_split"), 170)
