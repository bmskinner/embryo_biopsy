# Figure 7: effect of aneuploidy differences on rank ordering

source("parameters.R")
source("functions.R")

# # Read the saved raw values
output <- read.csv("data/Rank_results.csv", header = T)

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

################################################################################


# Show correct, incorrect and no ranks for clustered embryos
fig7A.plot <- ggplot(filt[filt$Low.disp == 0 & filt$High.disp == 0, ], aes(x = Aneu.diff * 100)) +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct, col = "Correct", shape = "Correct"), size = 0.5) +
  geom_line(aes(y = Mean.correct, col = "Correct")) +
  geom_ribbon(aes(
    ymin = Mean.correct - SD.correct,
    ymax = Mean.correct + SD.correct,
    fill = "Correct"
  ), alpha = 0.5, col = NA) +
  geom_point(aes(y = Mean.incorrect, col = "Incorrect", shape = "Incorrect"), size = 0.5) +
  geom_line(aes(y = Mean.incorrect, col = "Incorrect")) +
  geom_ribbon(aes(
    ymin = Mean.incorrect - SD.incorrect,
    ymax = Mean.incorrect + SD.incorrect,
    fill = "Incorrect"
  ), alpha = 0.5, col = NA) +
  geom_point(aes(y = Mean.no.rank, col = "Tied", shape = "Tied"), size = 0.5) +
  geom_line(aes(y = Mean.no.rank, col = "Tied")) +
  geom_ribbon(aes(
    ymin = Mean.no.rank - SD.no.rank,
    ymax = Mean.no.rank + SD.no.rank,
    fill = "Tied"
  ), alpha = 0.5, col = NA) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
  ) +
  scale_colour_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red")
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red"),
    guide = "none"
  ) +
  scale_shape_manual(
    name = element_blank(),
    values = c(Correct = 16, Tied = 15, Incorrect = 17)
  ) +
  labs(
    x = "Aneuploidy difference (%)", y = "Percent embryos with rank order...",
    color = "Z", shape = "Z"
  ) +
  theme_bw() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = c(0.8, 0.47),
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_blank()
  )

# Split the values for no.rank between the correct and incorrect
# Can't get error bars direct from this - recalculate from original values
fig7B.plot <- ggplot(filt[filt$Low.disp == 0 & filt$High.disp == 0, ], aes(x = Aneu.diff * 100)) +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.adj.incorrect, col = "Incorrect", shape = "Incorrect"), size = 0.5) +
  geom_line(aes(y = Mean.adj.incorrect, col = "Incorrect")) +
  geom_ribbon(aes(
    ymin = Mean.adj.incorrect - SD.adj.incorrect,
    ymax = Mean.adj.incorrect + SD.adj.incorrect,
    fill = "Incorrect"
  ), alpha = 0.5, col = NA) +
  geom_point(aes(y = Mean.adj.correct, col = "Correct", shape = "Correct"), size = 0.5) +
  geom_line(aes(, y = Mean.adj.correct, col = "Correct")) +
  geom_ribbon(aes(
    ymin = Mean.adj.correct - SD.adj.correct,
    ymax = Mean.adj.correct + SD.adj.correct,
    fill = "Correct"
  ), alpha = 0.5, col = NA) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
  ) +
  scale_colour_manual(
    name = element_blank(),
    values = c(Correct = "blue", Incorrect = "red")
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = c(Correct = "blue", Incorrect = "red"),
    guide = "none"
  ) +
  scale_shape_manual(
    name = element_blank(),
    values = c(Correct = 16, Incorrect = 17)
  ) +
  labs(
    x = "Aneuploidy difference (%)", y = "Effective percent embryos\nwith rank order...",
    color = "Z", shape = "Z"
  ) +
  theme_bw() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = c(0.8, 0.52),
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_blank()
  )


fig7 <- fig7A.plot + fig7B.plot + patchwork::plot_annotation(tag_levels = c("A"))
save.double.width(fig7, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7_Ranks_zero_dispersal"), 85)


# Make plots for talks - build each plot up in turn

fig7A.part1.plot <- ggplot(filt[filt$Low.disp == 0 & filt$High.disp == 0, ], aes(x = Aneu.diff * 100)) +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct, col = "Correct", shape = "Correct"), size = 0.5) +
  geom_line(aes(y = Mean.correct, col = "Correct")) +
  geom_ribbon(aes(
    ymin = Mean.correct - SD.correct,
    ymax = Mean.correct + SD.correct,
    fill = "Correct"
  ), alpha = 0.5, col = NA) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
  ) +
  scale_colour_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red")
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red"),
    guide = "none"
  ) +
  scale_shape_manual(
    name = element_blank(),
    values = c(Correct = 16, Tied = 15, Incorrect = 17)
  ) +
  labs(
    x = "Aneuploidy difference (%)", y = "Percent embryos with rank order...",
    color = "Z", shape = "Z"
  ) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = c(0.8, 0.47),
    legend.background = element_blank()
  )

fig7A.part2.plot <- ggplot(filt[filt$Low.disp == 0 & filt$High.disp == 0, ], aes(x = Aneu.diff * 100)) +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct, col = "Correct", shape = "Correct"), size = 0.5) +
  geom_line(aes(y = Mean.correct, col = "Correct")) +
  geom_ribbon(aes(
    ymin = Mean.correct - SD.correct,
    ymax = Mean.correct + SD.correct,
    fill = "Correct"
  ), alpha = 0.5, col = NA) +
  geom_point(aes(y = Mean.incorrect, col = "Incorrect", shape = "Incorrect"), size = 0.5) +
  geom_line(aes(y = Mean.incorrect, col = "Incorrect")) +
  geom_ribbon(aes(
    ymin = Mean.incorrect - SD.incorrect,
    ymax = Mean.incorrect + SD.incorrect,
    fill = "Incorrect"
  ), alpha = 0.5, col = NA) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
  ) +
  scale_colour_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red")
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = c(Correct = "blue", Tied = "black", Incorrect = "red"),
    guide = "none"
  ) +
  scale_shape_manual(
    name = element_blank(),
    values = c(Correct = 16, Tied = 15, Incorrect = 17)
  ) +
  labs(
    x = "Aneuploidy difference (%)", y = "Percent embryos with rank order...",
    color = "Z", shape = "Z"
  ) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = c(0.8, 0.52),
    legend.background = element_blank()
  )

save.single.width(fig7A.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7_part3"), 85)
save.single.width(fig7A.part2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7_part2"), 85)
save.single.width(fig7A.part1.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7_part1"), 85)
save.single.width(fig7B.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7_part4"), 85)
