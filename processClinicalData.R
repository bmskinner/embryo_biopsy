# Read Manuel's clinical data

library(readxl)
library(FSA)
library(ggpubr)

source("parameters.R")
source("functions.R")

in.file <- "data/Mosaic Embryo Transfers Sent.xlsx"

# Read the spreadsheet and clean the data
# Assign the embryos to appropriate categories
data <- readxl::read_xlsx(in.file,
  sheet = "Sheet2", range = "A2:N1734",
  col_names = c(
    "embryo", "mosaic_result", "level_aneuploidy", "type", "reported_outcome",
    "implantation_sac", "ongoing_birth", "unknown", "n_embryos_in_transfer",
    "notes_on_transfer", "n_total_blastocysts_in_cycle", "n_euploid_embryos_in_cycle",
    "n_mosaic_embryos_in_cycle", "n_aneuploid_embryos_in_cycle"
  ),
  col_types = c(
    "numeric", "text", "numeric", "text", "text", "numeric", "numeric", "numeric",
    "numeric", "text", "numeric", "numeric", "numeric", "numeric"
  )
) %>%
  dplyr::mutate(
    implantation_sac = factor(implantation_sac, levels = c(0, 1)),
    ongoing_birth = factor(ongoing_birth, levels = c(0, 1)),
    bin = factor(cut(level_aneuploidy, breaks = seq(0, 100, 10), include.lowest = T, right = F)),
    pgdis_class = factor(to.pgdis.class(level_aneuploidy / 100), levels = c("Euploid", "Low level", "High level", "Aneuploid")),
    split_class = factor(case_when(
      level_aneuploidy > 80 ~ "Aneuploid",
      level_aneuploidy >= 50 ~ "High level",
      level_aneuploidy >= 20 ~ "Low level",
      T ~ "Euploid"
    ), levels = c("Euploid", "Low level", "High level", "Aneuploid")),
    f_euploid_in_cycle = n_euploid_embryos_in_cycle / n_total_blastocysts_in_cycle,
    f_mosaic_in_cycle = n_mosaic_embryos_in_cycle / n_total_blastocysts_in_cycle,
    f_aneuploid_in_cycle = n_aneuploid_embryos_in_cycle / n_total_blastocysts_in_cycle,
    generic_outcome = case_when(
      reported_outcome == "Birth" ~ "Birth",
      reported_outcome == "Terminated" ~ "Terminated",
      str_detect(reported_outcome, "Abortion") ~ "Abortion",
      T ~ reported_outcome
    ),
    has.segmental = str_detect(mosaic_result, "[S|s]egmental") | str_detect(type, "[S|s]egmental") | str_detect(mosaic_result, "partial"),
    has.whole = str_detect(mosaic_result, "([W|w]ho[el|le])+") | str_detect(type, "(Monosomy|Trisomy)"),
    has.pq = str_detect(mosaic_result, "[\\d+|X|Y][p|q]+"),
    has.gain.loss = str_detect(mosaic_result, "(gain|loss)+"),
    has.mono.tri = str_detect(mosaic_result, "([M|m]ono(somy)?|[T|t]ri(somy)?)"),
    just.complex = str_detect(mosaic_result, "mos\\([C|c]omplex\\)"),
    affects.sex.chr = str_detect(mosaic_result, "X|Y"),
    not.X = str_detect(mosaic_result, "^(?:(?!X).)*$"),
    not.Y = str_detect(mosaic_result, "^(?:(?!Y).)*$"),
    is.autosomal.only = !affects.sex.chr,
    has.autosome = str_detect(mosaic_result, "[^pq]\\d+"),
    is.auto.and.sex = affects.sex.chr & has.autosome,
    seg_type = case_when(
      (has.segmental | has.pq) & has.whole ~ "Whole chromosome", # match definition in Viotti 2021 -
      has.segmental | has.pq ~ "Segmental only",
      has.whole | has.mono.tri ~ "Whole chromosome",
      just.complex ~ "Whole chromosome", # match definition in Viotti 2021 -When mosaicism
      # was present in more than two chromosomes, the mosaicism
      # type was considered <U+2018><U+2018>complex,<U+2019><U+2019> including combinations of
      # whole chromosomes and segmental regions.
      T ~ "Whole chromosome"
    ), # remainder that cannot be assigned segmental
    chr_type = case_when(
      is.autosomal.only ~ "Autosomal only",
      is.auto.and.sex ~ "Auto and sex chrs",
      T ~ "Sex chrs only"
    ),
    isOBP = generic_outcome == "Birth" | generic_outcome == "Ongoing Pregnancy" # ongoing pregnacy or birth
  )

################################################################################

# Implantation results grouped by segmental/whole chromosome
imp.data <- data %>%
  dplyr::group_by(bin, seg_type) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, implantation_sac, seg_type) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos,
    p_implanted_embryos = f_implanted_embryos * 100
  ) %>%
  dplyr::select(bin, seg_type, total_embryos, implanated_embryos, p_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1) %>%
  dplyr::mutate(bin_data = substr(bin, start = 2, stop = 6)) %>%
  tidyr::separate(bin_data, c("bin_min", "bin_max"), sep = ",", convert = T) %>%
  dplyr::mutate(label = paste0(
    bin_min, "-", bin_max - 1,
    "% n=",
    total_embryos
  ))

# Birth outcomes, as for the implantation data
birth.data <- data %>%
  dplyr::group_by(bin, seg_type) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, isOBP, seg_type) %>%
  dplyr::mutate(
    n_outcomes = n(),
    f_outcomes = n_outcomes / total_embryos,
    p_outcomes = f_outcomes * 100
  ) %>%
  dplyr::select(bin, seg_type, total_embryos, n_outcomes, p_outcomes) %>%
  dplyr::distinct() %>%
  dplyr::filter(isOBP == T) %>%
  rbind(., list(
    "isOBP" = T, "bin" = factor("[70,80)"), "seg_type" = "Whole chromosome",
    "total_embryos" = 23, "n_outcomes" = 0, "p_outcomes" = 0.000
  )) %>%
  dplyr::mutate(bin_data = substr(bin, start = 2, stop = 6)) %>%
  tidyr::separate(bin_data, c("bin_min", "bin_max"), sep = ",", convert = T) %>%
  dplyr::mutate(label = paste0(
    bin_min, "-", bin_max - 1,
    "% n=",
    total_embryos
  ))


# Implantation results with whole and segmental combined
imp.combined.data <- data %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos,
    p_implanted_embryos = f_implanted_embryos * 100
  ) %>%
  dplyr::select(bin, total_embryos, implanated_embryos, p_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1) %>%
  dplyr::mutate(bin_data = substr(bin, start = 2, stop = 6)) %>%
  tidyr::separate(bin_data, c("bin_min", "bin_max"), sep = ",", convert = T) %>%
  dplyr::mutate(label = paste0(
    bin_min, "-", bin_max - 1,
    "% n=",
    total_embryos
  ))

birth.combined.data <- data %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, isOBP) %>%
  dplyr::mutate(
    n_outcomes = n(),
    f_outcomes = n_outcomes / total_embryos,
    p_outcomes = f_outcomes * 100
  ) %>%
  dplyr::select(bin, total_embryos, n_outcomes, p_outcomes) %>%
  dplyr::distinct() %>%
  dplyr::filter(isOBP == T) %>%
  dplyr::mutate(bin_data = substr(bin, start = 2, stop = 6)) %>% # Add labels for x-axis
  tidyr::separate(bin_data, c("bin_min", "bin_max"), sep = ",", convert = T) %>%
  dplyr::mutate(label = paste0(
    bin_min, "-", bin_max - 1,
    "% n=",
    total_embryos
  ))

################################################################################

# make p-vals for linear regressions to display on each chart
imp.whole.glm <- glm(p_implanted_embryos ~ as.integer(bin), data = imp.data %>% dplyr::filter(seg_type == "Whole chromosome"))
imp.whole.pval <- round(coef(summary(imp.whole.glm))[2, 4], digits = 4)

imp.seg.glm <- glm(p_implanted_embryos ~ as.integer(bin), data = imp.data %>% dplyr::filter(seg_type == "Segmental only"))
imp.seg.pval <- round(coef(summary(imp.seg.glm))[2, 4], digits = 4)

imp.all.glm <- glm(p_implanted_embryos ~ as.integer(bin), data = imp.data)
imp.all.pval <- round(coef(summary(imp.all.glm))[2, 4], digits = 4)

birth.whole.glm <- glm(p_outcomes ~ as.integer(bin), data = birth.data %>% dplyr::filter(seg_type == "Whole chromosome"))
birth.whole.pval <- round(coef(summary(birth.whole.glm))[2, 4], digits = 4)

birth.seg.glm <- glm(p_outcomes ~ as.integer(bin), data = birth.data %>% dplyr::filter(seg_type == "Segmental only"))
birth.seg.pval <- round(coef(summary(birth.seg.glm))[2, 4], digits = 4)

birth.all.glm <- glm(p_outcomes ~ as.integer(bin), data = birth.data)
birth.all.pval <- round(coef(summary(birth.all.glm))[2, 4], digits = 4)


################################################################################

# Create Viotti et al 2021 Fig2A row 1 equivalent

YMAX.PCT <- 70

whole.xlabels <- imp.data %>%
  dplyr::filter(seg_type == "Whole chromosome") %>%
  dplyr::select(bin, label) %>%
  dplyr::arrange(bin)

# Whole chromosome (and whole with segmental)
whole.plot <- ggplot(
  imp.data %>% dplyr::filter(seg_type == "Whole chromosome"),
  aes(x = as.integer(bin), y = p_implanted_embryos)
) +
  geom_smooth(method = "glm", col = "black") +
  geom_smooth(
    data = birth.data %>% dplyr::filter(seg_type == "Whole chromosome"),
    aes(x = as.integer(bin), y = p_outcomes),
    method = "glm", col = "red", fill = "red"
  ) +
  geom_point(shape = 17) +
  geom_point(
    data = birth.data %>% dplyr::filter(seg_type == "Whole chromosome"),
    aes(x = as.integer(bin), y = p_outcomes), col = "red"
  ) +
  annotate("text",
    label = paste0("p=", imp.whole.pval, " (ns)"),
    x = 4.4, y = 50, size = 2, colour = "black"
  ) +
  annotate("text",
    label = paste0("p=", birth.whole.pval, " (ns)"),
    x = 3.4, y = 10, size = 2, colour = "red"
  ) +
  # geom_text(aes(label = total_embryos, y = 0.05), size = 2) +
  coord_cartesian(ylim = c(0, YMAX.PCT)) +
  scale_x_continuous(
    labels = function(x) str_wrap(whole.xlabels$label[x], width = 7),
    breaks = seq(0, 7, 1)
  ) +
  labs(x = "Mosaic level", y = "Positive outcome (%)", title = "Whole chromosome mosaics") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.7, size = 5),
    plot.title = element_text(hjust = 0.5, size = 8)
  )

seg.xlabels <- imp.data %>%
  dplyr::filter(seg_type == "Segmental only") %>%
  dplyr::select(bin, label) %>%
  dplyr::arrange(bin)

# Segmental only plot (no whole chromosome aneuploidies)
seg.plot <- ggplot(
  imp.data %>% dplyr::filter(seg_type == "Segmental only"),
  aes(x = as.integer(bin), y = p_implanted_embryos)
) +
  geom_smooth(method = "glm", col = "black") +
  geom_smooth(
    data = birth.data %>% dplyr::filter(seg_type == "Segmental only"),
    aes(x = as.integer(bin), y = p_outcomes),
    method = "glm", col = "red", fill = "red"
  ) +
  geom_point(shape = 17) +
  geom_point(data = birth.data %>% dplyr::filter(seg_type == "Segmental only"), aes(x = as.integer(bin), y = p_outcomes), col = "red") +
  annotate("text",
    label = paste0("p=", imp.seg.pval, " (ns)"),
    x = 4.4, y = 60, size = 2, colour = "black"
  ) +
  annotate("text",
    label = paste0("p=", birth.seg.pval, " (ns)"),
    x = 4.4, y = 20, size = 2, colour = "red"
  ) +
  # geom_text(aes(label = total_embryos, y = 0.05), size = 2) +
  coord_cartesian(ylim = c(0, YMAX.PCT)) +
  scale_x_continuous(
    labels = function(x) str_wrap(seg.xlabels$label[x], width = 7),
    breaks = seq(0, 7, 1)
  ) +
  labs(x = "Mosaic level", y = "Positive outcome (%)", title = "Segmental mosaics") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.7, size = 5),
    plot.title = element_text(hjust = 0.5, size = 8)
  )


combined.xlabels <- imp.combined.data %>%
  dplyr::select(bin, label) %>%
  dplyr::arrange(bin)

# Both whole and segmental data combined
combined.plot <- ggplot(imp.combined.data, aes(x = as.integer(bin), y = p_implanted_embryos)) +
  geom_smooth(method = "glm", col = "black") +
  geom_smooth(
    data = birth.combined.data,
    aes(x = as.integer(bin), y = p_outcomes),
    method = "glm", col = "red", fill = "red"
  ) +
  geom_point(aes(col = "Implantation", shape = "Implantation")) +
  geom_point(data = birth.combined.data, aes(
    x = as.integer(bin), y = p_outcomes,
    col = "Ongoing Pregnancy / Birth",
    shape = "Ongoing Pregnancy / Birth"
  )) +
  annotate("text",
    label = paste0("p=", imp.all.pval, " (ns)"),
    x = 4.4, y = 50, size = 2, colour = "black"
  ) +
  annotate("text",
    label = paste0("p=", birth.all.pval, " (ns)"),
    x = 4.4, y = 15, size = 2, colour = "red"
  ) +
  # geom_text(aes(label = total_embryos, y = 0.05), size = 2) +
  coord_cartesian(ylim = c(0, YMAX.PCT)) +
  scale_x_continuous(
    labels = function(x) str_wrap(combined.xlabels$label[x], width = 7),
    breaks = seq(0, 7, 1)
  ) +
  scale_shape_manual(
    name = "State",
    labels = c("Implantation", "Ongoing Pregnancy / Birth"),
    values = c("Ongoing Pregnancy / Birth" = 16, "Implantation" = 17)
  ) +
  scale_color_manual(
    name = "State",
    labels = c("Implantation", "Ongoing Pregnancy / Birth"),
    values = c("Ongoing Pregnancy / Birth" = "red", "Implantation" = "black")
  ) +
  labs(x = "Mosaic level", y = "Positive outcome (%)", title = "All mosaics") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.7, size = 5),
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.title = element_blank(),
    legend.position = c(0.5, 0.9),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "mm")
  )



################################################################################

# Recreate Viotti et al Fig 2A row 2: split into two groups, aneuploidy <50 or >=50%
# Compare rates of implantation or OP/B in each group


imp.halves.data <- data %>%
  dplyr::mutate(group = case_when(
    level_aneuploidy < 50 ~ "<50%",
    T ~ "\u226550%"
  )) %>%
  dplyr::group_by(group, seg_type) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(group, implantation_sac, seg_type) %>%
  dplyr::mutate(
    n_embryos = n(),
    f_implanted_embryos = n_embryos / total_embryos,
    p_outcome = f_implanted_embryos * 100
  ) %>%
  dplyr::select(group, seg_type, total_embryos, n_embryos, p_outcome) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1) %>%
  dplyr::mutate(Type = "Implantation") %>%
  dplyr::ungroup() %>%
  dplyr::select(Type, group, seg_type, total_embryos, n_embryos, p_outcome) %>%
  dplyr::mutate(x_embryo = total_embryos - n_embryos)

birth.halves.data <- data %>%
  dplyr::mutate(group = case_when(
    level_aneuploidy < 50 ~ "<50%",
    T ~ "\u226550%"
  )) %>%
  dplyr::group_by(group, seg_type) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(group, isOBP, seg_type) %>%
  dplyr::mutate(
    n_embryos = n(),
    f_implanted_embryos = n_embryos / total_embryos,
    p_outcome = f_implanted_embryos * 100
  ) %>%
  dplyr::select(group, seg_type, total_embryos, n_embryos, p_outcome) %>%
  dplyr::distinct() %>%
  dplyr::filter(isOBP == T) %>%
  dplyr::mutate(Type = "OP/B") %>%
  dplyr::ungroup() %>%
  dplyr::select(Type, group, seg_type, total_embryos, n_embryos, p_outcome) %>%
  dplyr::mutate(x_embryo = total_embryos - n_embryos)

combined.halves.data <- rbind(imp.halves.data, birth.halves.data) %>%
  dplyr::mutate(x_embryo = total_embryos - n_embryos) %>%
  dplyr::mutate(
    seg_type = factor(seg_type, levels = c("Whole chromosome", "Segmental only")),
    cgroup = paste0(Type, group)
  )


# combine the whole and segmental results together and bind
all.halves.data <- combined.halves.data %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group, Type) %>%
  dplyr::summarise(
    seg_type = "All",
    total_embryos = sum(total_embryos),
    n_embryos = sum(n_embryos),
    p_outcome = n_embryos / total_embryos * 100,
    x_embryo = total_embryos - n_embryos,
    cgroup = paste0(Type, group)
  ) %>%
  dplyr::distinct() %>%
  rbind(., combined.halves.data)


# Test associations
tests <- expand.grid(
  seg_type = c("Whole chromosome", "Segmental only", "All"),
  Type = c("Implantation", "OP/B")
)

run.test <- function(s, t) {
  filt <- all.halves.data %>%
    dplyr::ungroup() %>%
    dplyr::filter(seg_type == s & Type == t) %>%
    dplyr::select(n_embryos, x_embryo)
  ch <- chisq.test(filt)
  ch$p.value
}

tests$test.results <- mapply(run.test, tests$seg_type, tests$Type)
tests$p.adj <- p.adjust(tests$test.results, method = "bonferroni")
tests$x <- ifelse(tests$Type == "Implantation", 1.5, 3.5)


# Make the plots
# Plot for whole chromosomes only
whole.halves.data <- combined.halves.data %>%
  dplyr::filter(seg_type == "Whole chromosome")

whole.lines.data <- data.frame(
  "seg_type" = c(
    "Whole chromosome", "Whole chromosome", "Whole chromosome",
    "Whole chromosome", "Whole chromosome", "Whole chromosome"
  ),
  "Type" = c(
    "Implantation", "Implantation", "Implantation",
    "OP/B", "OP/B", "OP/B"
  ),
  "xstart" = c(
    1, 1, 2,
    3, 3, 4
  ),
  "xend" = c(
    1, 2, 2,
    3, 4, 4
  ),
  "ystart" = c(
    44.25, 46.25, 46.25,
    36, 38, 38
  ),
  "yend" = c(
    46.25, 46.25, 34,
    38, 38, 22
  )
)

opb.colour <- "#ff6d6dff" # eyedropper from transparent region of upper panels, but with more saturation because there is full red in the upper panels

whole.halves.fig <- ggplot(whole.halves.data, aes(x = cgroup, y = p_outcome)) +
  geom_col(aes(fill = Type)) +
  geom_segment(data = whole.lines.data, aes(x = xstart, y = ystart, xend = xend, yend = yend)) +
  geom_text(aes(y = p_outcome - 5, label = round(p_outcome, digits = 2)), size = 2) +
  geom_text(aes(y = 5, label = paste0("n=", round(total_embryos, digits = 2))), size = 2) +
  geom_text(
    data = tests[tests$seg_type == "Whole chromosome", ], aes(y = 53, x = x, label = paste0("p=", round(p.adj, digits = 4))),
    size = 2
  ) +
  labs(y = "Positive outcome (%)", title = "Whole chromosome mosaics") +
  scale_x_discrete(labels = c("<50%", "\u226550%", "<50%", "\u226550%")) +
  coord_cartesian(ylim = c(0, YMAX.PCT)) +
  scale_fill_manual(values = c("grey", opb.colour)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "mm")
  )


seg.halves.data <- combined.halves.data %>%
  dplyr::filter(seg_type == "Segmental only")

seg.lines.data <- data.frame(
  "seg_type" = c(
    "Segmental only", "Segmental only", "Segmental only",
    "Segmental only", "Segmental only", "Segmental only"
  ),
  "Type" = c(
    "Implantation", "Implantation", "Implantation",
    "OP/B", "OP/B", "OP/B"
  ),
  "xstart" = c(
    1, 1, 2,
    3, 3, 4
  ),
  "xend" = c(
    1, 2, 2,
    3, 4, 4
  ),
  "ystart" = c(
    54, 55, 55,
    45, 47, 47
  ),
  "yend" = c(
    55, 55, 49,
    47, 47, 39
  )
)

# Segmental only, 50/50 split figure
seg.halves.fig <- ggplot(seg.halves.data, aes(x = cgroup, y = p_outcome)) +
  geom_col(aes(fill = Type)) +
  geom_segment(data = seg.lines.data, aes(x = xstart, y = ystart, xend = xend, yend = yend)) +
  geom_text(aes(y = p_outcome - 5, label = round(p_outcome, digits = 2)), size = 2) +
  geom_text(aes(y = 5, label = paste0("n=", round(total_embryos, digits = 2))), size = 2) +
  geom_text(
    data = tests[tests$seg_type == "Segmental only", ], aes(y = 58, x = x, label = paste0("p=", round(p.adj, digits = 4))),
    size = 2
  ) +
  labs(y = "Positive outcome (%)", title = "Segmental mosaics") +
  scale_x_discrete(labels = c("<50%", "\u226550%", "<50%", "\u226550%")) +
  coord_cartesian(ylim = c(0, YMAX.PCT)) +
  scale_fill_manual(values = c("grey", opb.colour)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "mm")
  )

both.halves.data <- combined.halves.data %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Type, group) %>%
  summarise(
    total_embryos = sum(total_embryos),
    n_embryos = sum(n_embryos),
    p_outcome = n_embryos / total_embryos * 100,
    cgroup = paste0(Type, group)
  ) %>%
  dplyr::distinct()


both.lines.data <- data.frame(
  "Type" = c(
    "Implantation", "Implantation", "Implantation",
    "OP/B", "OP/B", "OP/B"
  ),
  "xstart" = c(
    1, 1, 2,
    3, 3, 4
  ),
  "xend" = c(
    1, 2, 2,
    3, 4, 4
  ),
  "ystart" = c(
    50, 52, 52,
    40, 42, 42
  ),
  "yend" = c(
    52, 52, 42,
    42, 42, 30
  )
)

both.halves.fig <- ggplot(both.halves.data, aes(x = cgroup, y = p_outcome)) +
  geom_col(aes(fill = Type)) +
  geom_segment(data = both.lines.data, aes(x = xstart, y = ystart, xend = xend, yend = yend)) +
  geom_text(aes(y = p_outcome - 5, label = round(p_outcome, digits = 2)), size = 2) +
  geom_text(aes(y = 5, label = paste0("n=", round(total_embryos, digits = 2))), size = 2) +
  geom_text(
    data = tests[tests$seg_type == "All", ], aes(y = 58, x = x, label = paste0("p=", round(p.adj, digits = 4))),
    size = 2
  ) +
  labs(y = "Positive outcome (%)", title = "All mosaics") +
  scale_x_discrete(labels = c("<50%", "\u226550%", "<50%", "\u226550%")) +
  coord_cartesian(ylim = c(0, YMAX.PCT)) +
  scale_fill_manual(values = c("grey", opb.colour)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "mm")
  )


fig9.full <- (combined.plot + whole.plot + seg.plot) / (both.halves.fig + whole.halves.fig + seg.halves.fig) +
  plot_annotation(tag_levels = list(c("A", "", "", "B", "", "")))
save.double.width(fig9.full, filename = "figure/Figure_9_full", height = 140)




################################################################################

# Create glms to test associations between aneuploidy, segmentation type and outcomes

# Check aneuploidy and segmental type in glm
imp.glm <- glm(implantation_sac ~ level_aneuploidy + seg_type, family = binomial(link = "logit"), data = data)
imp.sum <- broom::tidy(imp.glm) %>%
  dplyr::mutate(
    odds.ratio = exp(1)^estimate, # e^<U+03B2>
    odds.ci.min = exp(1)^(estimate - 1.96 * std.error),
    odds.ci.max = exp(1)^(estimate + 1.96 * std.error),
    odds.dec = (odds.ratio - 1) * 100,
    odds.ci.min.pct = (odds.ci.min - 1) * 100,
    odds.ci.max.pct = (odds.ci.max - 1) * 100,
  ) # e(<U+03B2> +/- 1.96*std error).

# Is the model useful? Strong lack of support for the null hypothesis
pchisq(imp.glm$null.deviance - imp.glm$deviance, imp.glm$df.null - imp.glm$df.residual, lower.tail = F)

# What does the model explain?
imp.deviance.diff <- imp.glm$null.deviance - imp.glm$deviance
imp.deviance.diff / imp.glm$null.deviance * 100
# Very little of the variation (~1%) is explained by this model


opb.glm <- glm(ongoing_birth ~ level_aneuploidy + seg_type, family = binomial(link = "logit"), data = data)
summary(opb.glm)

opb.sum <- broom::tidy(opb.glm) %>%
  dplyr::mutate(
    odds.ratio = exp(1)^estimate, # e^<U+03B2>
    odds.ci.min = exp(1)^(estimate - 1.96 * std.error),
    odds.ci.max = exp(1)^(estimate + 1.96 * std.error),
    odds.dec = (odds.ratio - 1) * 100,
    odds.ci.min.pct = (odds.ci.min - 1) * 100,
    odds.ci.max.pct = (odds.ci.max - 1) * 100,
  ) # e(<U+03B2> +/- 1.96*std error).

# Is the model useful? Strong lack of support for the null hypothesis
pchisq(opb.glm$null.deviance - opb.glm$deviance, opb.glm$df.null - opb.glm$df.residual, lower.tail = F)


# What does the model explain?
opb.deviance.diff <- opb.glm$null.deviance - opb.glm$deviance
opb.deviance.diff / opb.glm$null.deviance * 100
# Very little of the variation (~1%) is explained by this model



################################################################################


# Make logistic plots to show the glm data

# Predict values using the model
new.data <- expand.grid(
  level_aneuploidy = seq(20, 80, 5),
  seg_type = c("Whole chromosome", "Segmental only")
)

imp.preds <- predict(imp.glm, newdata = new.data, type = "response", se = T)
new.data$imp.pred.full <- imp.preds$fit
new.data$imp.ymin <- new.data$imp.pred.full - 2 * imp.preds$se.fit
new.data$imp.ymax <- new.data$imp.pred.full + 2 * imp.preds$se.fit


opb.preds <- predict(opb.glm, newdata = new.data, type = "response", se = T)
new.data$opb.pred.full <- opb.preds$fit
new.data$opb.ymin <- new.data$opb.pred.full - 2 * opb.preds$se.fit
new.data$opb.ymax <- new.data$opb.pred.full + 2 * opb.preds$se.fit


# Plot the real data and add the predicted values & error

glm.fig <- ggplot() +
  geom_jitter(data = data, aes(x = level_aneuploidy, y = as.integer(implantation_sac) - 1), height = 0.1) +
  geom_ribbon(data = new.data, aes(x = level_aneuploidy, y = imp.pred.full, ymin = imp.ymin, ymax = imp.ymax), alpha = 0.25) +
  geom_line(data = new.data, aes(x = level_aneuploidy, y = imp.pred.full, colour = "Implantation"), linewidth = 2) +
  geom_ribbon(data = new.data, aes(x = level_aneuploidy, y = opb.pred.full, ymin = opb.ymin, ymax = opb.ymax), fill = "red", alpha = 0.25) +
  geom_line(data = new.data, aes(x = level_aneuploidy, y = opb.pred.full, colour = "Ongoing pregnancy / birth"), linewidth = 2) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~seg_type) +
  labs(x = "Mosaic level (%)", y = "Probability of successful outcome") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.15, 0.665),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(1.5, "mm")
  )

save.double.width(glm.fig, filename = "figure/Figure_9_glm", height = 85)
