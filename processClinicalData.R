# Read Manuel's clinical data

library(tessera)
library(tidyverse)
library(parallel)
library(patchwork)
library(svglite)
library(readxl)
library(FSA)

source("parameters.R")
source("functions.R")

in.file <- "data/Mosaic Embryo Transfers Sent.xlsx"

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
    )
  )


################################################################################

# Implantation results

imp.data <- data %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(bin, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1)


imp.plot <- ggplot(imp.data, aes(x = bin, y = f_implanted_embryos)) +
  geom_col() +
  labs(x = "Aneuploidy (%)", y = "Fraction of embryos implanted") +
  theme_classic() +
  theme(legend.position = "none")

# Remove segmental results

implantation.no.segmentatals <- data %>%
  dplyr::filter(!(str_detect(type, "Segmental"))) %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(bin, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1)

imp.no.seg.plot <- ggplot(implantation.no.segmentatals, aes(x = bin, y = f_implanted_embryos)) +
  geom_col() +
  labs(x = "Aneuploidy (%)", y = "Fraction of embryos implanted") +
  theme_classic() +
  theme(legend.position = "none")

# Implantation by PGDIS classes

imp.pgdis.data <- data %>%
  dplyr::group_by(pgdis_class) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(pgdis_class, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(pgdis_class, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct()

imp.pgdis.plot <- ggplot(imp.pgdis.data, aes(x = pgdis_class, y = f_implanted_embryos)) +
  geom_col() +
  labs(x = "PGDIS class", y = "Fraction of embryos implanted") +
  theme_classic() +
  theme(legend.position = "none")

imp.pgdis.no.seg.data <- data %>%
  dplyr::group_by(pgdis_class) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(pgdis_class, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(pgdis_class, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1)

imp.pgdis.no.seg.plot <- ggplot(imp.pgdis.no.seg.data, aes(x = pgdis_class, y = f_implanted_embryos)) +
  geom_col() +
  labs(x = "PGDIS class", y = "Fraction of embryos implanted") +
  theme_classic() +
  theme(legend.position = "none")



# Implantation by split classes
imp.split.data <- data %>%
  dplyr::group_by(split_class) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(split_class, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(split_class, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1)

imp.split.plot <- ggplot(imp.split.data, aes(x = split_class, y = f_implanted_embryos)) +
  geom_col() +
  labs(x = "50% split class", y = "Fraction of embryos implanted") +
  theme_classic() +
  theme(legend.position = "none")

imp.split.no.seg.data <- data %>%
  dplyr::group_by(split_class) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(split_class, implantation_sac) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(split_class, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1)

imp.split.no.seg.plot <- ggplot(imp.split.no.seg.data, aes(x = split_class, y = f_implanted_embryos)) +
  geom_col() +
  labs(x = "50% split class", y = "Fraction of embryos implanted") +
  theme_classic() +
  theme(legend.position = "none")

(imp.plot + imp.no.seg.plot) / (imp.pgdis.plot + imp.pgdis.no.seg.plot) / (imp.split.plot + imp.split.no.seg.plot) + plot_annotation(tag_levels = c("A"))




################################################################################

# Outcome results
outcome.data <- data %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, reported_outcome) %>%
  dplyr::mutate(
    f_outcome_embryos = n(),
    f_outcome_embryos = f_outcome_embryos / total_embryos
  ) %>%
  dplyr::select(bin, reported_outcome, f_outcome_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(reported_outcome != "n/a")

ggplot(outcome.data, aes(x = bin, y = f_outcome_embryos)) +
  geom_col() +
  labs(x = "Aneuploidy (%)", y = "Fraction of embryos implanted") +
  facet_wrap(~reported_outcome) +
  theme_classic()

# Outcome generic results
outcome.data <- data %>%
  dplyr::group_by(bin) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, generic_outcome) %>%
  dplyr::mutate(
    f_outcome_embryos = n(),
    f_outcome_embryos = f_outcome_embryos / total_embryos
  ) %>%
  dplyr::select(bin, generic_outcome, f_outcome_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(generic_outcome != "n/a")

outcome.plot <- ggplot(outcome.data, aes(x = bin, y = f_outcome_embryos)) +
  geom_col() +
  labs(x = "Aneuploidy (%)", y = "Fraction of embryos with outcome") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~generic_outcome) +
  theme_classic()

# Outcome by PGDIS class
# Outcome generic results
pgdis.outcome.data <- data %>%
  dplyr::group_by(pgdis_class) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(pgdis_class, generic_outcome) %>%
  dplyr::mutate(
    f_outcome_embryos = n(),
    f_outcome_embryos = f_outcome_embryos / total_embryos
  ) %>%
  dplyr::select(pgdis_class, generic_outcome, f_outcome_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(generic_outcome != "n/a")


pgdis.outcome.plot <- ggplot(pgdis.outcome.data, aes(x = pgdis_class, y = f_outcome_embryos)) +
  geom_col() +
  labs(x = "PGDIS class", y = "Fraction of embryos with outcome") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~generic_outcome) +
  theme_classic()

# Outcome by split class
# Outcome generic results
split.outcome.data <- data %>%
  dplyr::group_by(split_class) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(split_class, generic_outcome) %>%
  dplyr::mutate(
    f_outcome_embryos = n(),
    f_outcome_embryos = f_outcome_embryos / total_embryos
  ) %>%
  dplyr::select(split_class, generic_outcome, f_outcome_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(generic_outcome != "n/a")


split.outcome.plot <- ggplot(split.outcome.data, aes(x = split_class, y = f_outcome_embryos)) +
  geom_col() +
  labs(x = "50% split class", y = "Fraction of embryos with outcome") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~generic_outcome) +
  theme_classic()


outcome.plot + pgdis.outcome.plot + split.outcome.plot + plot_annotation(tag_levels = c("A"))


################################################################################


# Measure the contribution of class/bin to outcomes

kruskal.test(generic_outcome ~ split_class, data = data)
kruskal.test(reported_outcome ~ split_class, data = data)
kruskal.test(generic_outcome ~ pgdis_class, data = data)
kruskal.test(generic_outcome ~ bin, data = data)


kruskal.test(implantation_sac ~ bin, data = data)
kruskal.test(implantation_sac ~ pgdis_class, data = data)
kruskal.test(implantation_sac ~ split_class, data = data)


m <- glm(implantation_sac ~ level_aneuploidy, family = binomial(link = "logit"), data = data)
summary(m)
# Is the model useful to explain variation in measurements? Calc chi sq and convert to pvalue
pchisq(m$null.deviance - m$deviance, 1, lower.tail = T)

# Test using a multinomial model for reported outcome data
glm.fit <- multinom(reported_outcome ~ level_aneuploidy, data = data)
summary(glm.fit)

# Separate test and training data
alpha <- 0.7
d <- sort(sample(nrow(data), nrow(data) * alpha))
train <- data[d, ]
test <- data[-d, ]
glm.fit <- multinom(reported_outcome ~ level_aneuploidy, data = train)
test$predict <- predict(glm.fit, test)

test$correct <- test$predict == test$reported_outcome
sum(test$correct) / length(test$correct)


################################################################################


# Look at fraction of euploid vs mosaic embryos vs outcomes


ggplot(data, aes(x = reported_outcome, y = f_mosaic_in_cycle)) +
  geom_violin() +
  geom_jitter() +
  facet_wrap(~type)

ggplot(data, aes(x = reported_outcome, y = f_euploid_in_cycle)) +
  geom_violin() +
  geom_jitter() +
  facet_wrap(~type)
ggplot(data, aes(x = reported_outcome, y = f_aneuploid_in_cycle)) +
  geom_violin() +
  geom_jitter() +
  facet_wrap(~type)


# Look at how total faction of mosaic embryos vary with aneuploidy levels

# You don't have more mosaic embryos when aneuploidy is higher
ggplot(data, aes(x = level_aneuploidy, y = f_mosaic_in_cycle)) +
  geom_point() +
  geom_smooth()


ggplot(data, aes(x = level_aneuploidy, y = f_euploid_in_cycle)) +
  geom_point() +
  geom_smooth()
ggplot(data, aes(x = level_aneuploidy, y = f_aneuploid_in_cycle)) +
  geom_point() +
  geom_smooth()

################################################################################

# How many double transfers?


filt <- data %>%
  filter(n_embryos_in_transfer > 1) %>%
  dplyr::mutate(identity_deduced = str_detect(notes_on_transfer, "identity deduced")) %>%
  filter(identity_deduced)

ggplot(data, aes(x = reported_outcome, y = level_aneuploidy)) +
  geom_violin() +
  geom_jitter()

ggplot(data, aes(x = type, y = level_aneuploidy)) +
  geom_violin()

ggplot(data, aes(x = type, y = level_aneuploidy)) +
  geom_violin()
