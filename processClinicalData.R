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
    ),
    has.segmental = str_detect(mosaic_result, "[S|s]egmental") | str_detect(type, "[S|s]egmental")| str_detect(mosaic_result, "partial"),
    has.whole     = str_detect(mosaic_result, "([W|w]ho[el|le])+") | str_detect(type, "(Monosomy|Trisomy)"),
    has.pq        = str_detect(mosaic_result, "[\\d+|X|Y][p|q]+"),
    has.gain.loss = str_detect(mosaic_result, "(gain|loss)+"),
    has.mono.tri  = str_detect(mosaic_result, "([M|m]ono(somy)?|[T|t]ri(somy)?)"),
    just.complex  = str_detect(mosaic_result, "mos\\([C|c]omplex\\)"),
    affects.sex.chr = str_detect(mosaic_result, "X|Y"),
    not.X         = str_detect(mosaic_result, "^(?:(?!X).)*$"),
    not.Y         = str_detect(mosaic_result, "^(?:(?!Y).)*$"),
    is.autosomal.only  = !affects.sex.chr,
    has.autosome   =  str_detect(mosaic_result, "[^pq]\\d+"),
    is.auto.and.sex    = affects.sex.chr & has.autosome,
    seg_type = case_when(
      (has.segmental| has.pq) & has.whole ~ "Segmental and whole",
      has.segmental | has.pq ~ "Segmental",
      has.whole | has.mono.tri ~ "Whole chromosome",
      just.complex ~ "Unclear",
      T ~ "Whole chromosome"), # remainder that cannot be assigned segmental
    chr_type = case_when(
      is.autosomal.only ~ "Autosomal only",
      is.auto.and.sex ~ "Auto and sex chrs",
      T~ "Sex chrs only"
    )
    )
  


################################################################################

# Implantation results

imp.data <- data %>%
  dplyr::group_by(bin, seg_type, chr_type) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, implantation_sac, seg_type, chr_type) %>%
  dplyr::mutate(
    implanated_embryos = n(),
    f_implanted_embryos = implanated_embryos / total_embryos
  ) %>%
  dplyr::select(bin,seg_type, total_embryos, implanated_embryos, f_implanted_embryos) %>%
  dplyr::distinct() %>%
  dplyr::filter(implantation_sac == 1)

# Check aneuploidy and segmntal type in glm
imp.glm <- glm(implantation_sac ~ level_aneuploidy+seg_type+chr_type, family = binomial(link = "logit"), data = data)
summary(imp.glm)
# Is the model useful? Strong lack of support for the null hypothesis
pchisq(imp.glm$null.deviance - imp.glm$deviance, imp.glm$df.null-imp.glm$df.residual, lower.tail = F)

# What does the model explain?
deviance.diff = imp.glm$null.deviance-imp.glm$deviance
deviance.diff/imp.glm$null.deviance*100
# Very little of the variation (~1%) is explained by this model

imp.plot <- ggplot(imp.data, aes(x = as.integer(bin), y = f_implanted_embryos)) +
  geom_hline(yintercept = 0.5, size=1)+
  geom_point() +
  # geom_line()+
  geom_smooth(method = "lm")+
  geom_text(aes(label=total_embryos, y =0.05))+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(labels = function(x) levels(data$bin)[x], 
                     breaks = seq(0, 7, 1))+
  labs(x = "Aneuploidy bin", y = "Fraction of embryos implanted") +
  theme_bw()+
  facet_grid(chr_type~seg_type)

save.double.width(imp.plot, filename = "figure/Figure_8_implantation", height = 170)





# Birth outcomes, as for the implantation data
birth.data <- data %>%
  dplyr::group_by(bin, seg_type, chr_type) %>%
  dplyr::mutate(total_embryos = n()) %>%
  dplyr::group_by(bin, ongoing_birth, seg_type, chr_type) %>%
  dplyr::mutate(
    n_outcomes = n(),
    f_outcomes = n_outcomes / total_embryos
  ) %>%
  dplyr::select(bin,seg_type, total_embryos, n_outcomes, f_outcomes) %>%
  dplyr::distinct() %>%
  dplyr::filter(ongoing_birth == 1)

# Check aneuploidy and segmntal type in glm
birth.glm <- glm(ongoing_birth ~ level_aneuploidy+seg_type+chr_type, family = binomial(link = "logit"), data = data)
summary(birth.glm)
# Is the model useful? Strong lack of support for the null hypothesis
pchisq(birth.glm$null.deviance - birth.glm$deviance,  birth.glm$df.null-birth.glm$df.residual, lower.tail = F)

# What does the model explain?
deviance.diff = birth.glm$null.deviance-birth.glm$deviance
deviance.diff/birth.glm$null.deviance*100
# Very little of the variation (~1%) is explained by this model




birth.plot <- ggplot(birth.data, aes(x = as.integer(bin), y = f_outcomes)) +
  geom_hline(yintercept = 0.5, size=1)+
  geom_point() +
  # geom_line()+
  geom_smooth(method = "lm")+
  geom_text(aes(label=total_embryos, y =0.05))+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(labels = function(x) levels(data$bin)[x], 
                     breaks = seq(0, 7, 1))+
  labs(x = "Aneuploidy bin", y = "Fraction of births") +
  theme_bw()+
  facet_grid(chr_type~seg_type)

save.double.width(birth.plot, filename = "figure/Figure_9_births", height = 170)



















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


imp.plot + imp.no.seg.plot

# Split to all, whole chr and segmental mosaics

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
