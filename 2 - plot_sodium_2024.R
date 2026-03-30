# TITLE: Sodium plots - Sodium sampling 2024
# AUTHOR: Ashley Darst
# DATA INPUT: 2024 Nectar, leaf, and soil sodium quantification
# DATA OUTPUT: Analyses and main figure 1A-B
# DATE: 2024-09-25
# DESCRIPTION: Clean nectar, leaf, and soil sampling data, conduct analysis, and create figure panels

# ***MUST RUN THIS SCRIPT BEFORE plot_sodium_2025.R***

# Load libraries
library(tidyverse)
library(janitor)
library(ggforce)
library(sjPlot)
library(lme4)
library(lmerTest)
library(bbmle)
library(emmeans)
library(car)
library(RColorBrewer)

# Data cleaning ----
# Load data
soil <- read.csv("data/soil_icp.csv", na.strings=c("","NA"))
leaf <- read.csv("data/leaf_icp.csv")
nectar <- read.csv("data/nectar_icp.csv")


# Clean column names
soil <- clean_names(soil)
leaf <- clean_names(leaf)
nectar <- clean_names(nectar)

# Replace c with control and na with sodium
leaf <- leaf %>% 
  mutate(trt = ifelse(trt == "c", "control", "sodium"))
nectar <- nectar %>% 
  mutate(trt = ifelse(trt == "c", "control", "sodium"))
 
# Date to date format
soil$date <- as.Date(soil$date)
leaf$date <- as.Date(leaf$date)
nectar$date <- as.Date(nectar$date)

# trt and species to factor 
soil$trt <- as.factor(soil$trt)
leaf$trt <- as.factor(leaf$trt)
leaf$species <- as.factor(leaf$species)
nectar$trt <- as.factor(nectar$trt)

# Do average for repeats
(110.62 + 124.29) / 2
soil$ug_na_g_soil[which(soil$soil_id == 70)] <- 117.455
(99.44 + 113.09) / 2
soil$ug_na_g_soil[which(soil$soil_id == 71)] <- 106.265

# Remove unwanted rows
soil <- soil[-c(41:42), ]
leaf <- leaf[-c(20:23), ]

# Add pair_id to soil
soil$pair_id <- sub("_[^_]+$", "", soil$plot_id)
leaf$pair_id <- sub("_[^_]+$", "", leaf$plot_id)

# Change date to factor
soil$date_factor <- as.factor(soil$date)

# Plots ----
## Soil plot ----
soil.plot <- soil %>%
  ggplot(aes(x = trt, y = ug_na_g_soil)) +
  geom_sina(alpha = 0.2, jitter_y = F) +
  stat_summary(fun.data = "mean_cl_boot") +
  labs(y = "Soil sodium (ppm)", x = "Treatment") +
  scale_x_discrete(labels=c("control", "sodium")) +
  theme_bw() +
  theme(text = element_text(size = 15))

## Leaf plot ----
leaf.plot <- leaf %>%
  ggplot(aes(x = trt, y = ug_na_g_leaf_powder, col = species)) +
  geom_sina(alpha = 0.3) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge2(0.7)) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  labs(y = "Leaf sodium (ppm)", x = "Treatment", col = NULL) +
  scale_x_discrete(labels=c("control", "sodium")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.text = element_text(face = "italic"))

## Figure S2 ----
# Including values below LOD
# Have not standardized by gram
nectar %>%
  ggplot(aes(x = trt, y = total_ng_na_in_icp_sample)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=0,ymax=100),
            fill="salmon") +
  geom_sina(alpha = 0.2) +
  geom_hline(yintercept=100, linetype="dashed", color = "red") +
  labs(y = "total sample Na (ng) for nectar", x = "treatment") +
  scale_x_discrete(labels=c("control", "sodium")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  stat_summary(fun.data = "mean_cl_boot")

# Analyses ----
## Soil ----
# Simple soil model
soil.lm.rm <- lm(ug_na_g_soil ~ trt, data = soil)
summary(soil.lm.rm)
anova(soil.lm.rm)
plot(soil.lm.rm)

# Ln-transformed to help normality
soil.lm.rm.ln <- lm(log(ug_na_g_soil) ~ trt, data = soil)
summary(soil.lm.rm.ln)
anova(soil.lm.rm.ln)
plot(soil.lm.rm.ln)

# Soil model with plot nested in pair as random effect (for repeated sampling) and time point
soil.lmer <- lmer(log(ug_na_g_soil) ~ trt + (1|pair_id) + (1|plot_id:pair_id), data = soil)
summary(soil.lmer) # random effects explain zero variance

### Final model stats ----
# Add date to soil model
soil.lmer.date <- lmer(log(ug_na_g_soil) ~ date_factor + trt + (1|pair_id) + (1|plot_id:pair_id), data = soil)
summary(soil.lmer.date)
Anova(soil.lmer.date, type = "II")
plot(soil.lmer.date)
qqnorm(resid(soil.lmer.date))
qqline(resid(soil.lmer.date))
(exp(2.40158) - 1) * 100 # about 1004% more Na in sodium plots

soil.lm.date <- lm(log(ug_na_g_soil) ~ date_factor + trt, data = soil)
summary(soil.lm.date)
Anova(soil.lm.date, type = "II")
plot(soil.lm.date)

soil %>%
  group_by(trt) %>%
  summarize(mean = mean(ug_na_g_soil))

## Leaf ----
# Simple leaf model
leaf.lm <- lm(ug_na_g_leaf_powder ~ trt, data = leaf)
summary(leaf.lm)
plot(leaf.lm)

# Ln-transformed to help normality and residuals
leaf.lm.ln <- lm(log(ug_na_g_leaf_powder) ~ trt, data = leaf)
summary(leaf.lm.ln)
plot(leaf.lm.ln) # still not perfect

# Leaf model with species and ln-transformed
leaf.lm.ln2 <- lm(log(ug_na_g_leaf_powder) ~ species + trt, data = leaf)
summary(leaf.lm.ln2)
anova(leaf.lm.ln2)
plot(leaf.lm.ln2) # bad

# Leaf model with species and trt interaction
leaf.lm.ln.int <- lm(log(ug_na_g_leaf_powder) ~ species*trt, data = leaf)
leaf.aov.ln.int.aov <- aov(log(ug_na_g_leaf_powder) ~ species*trt, data = leaf)
summary(leaf.lm.ln.int)
Anova(leaf.lm.ln.int, type = "II")
emmeans(leaf.lm.ln.int, revpairwise ~ trt|species, adjust = "tukey")

### Final model stats ----
# With plot pair random effect
leaf.lmer.ln.int <- lmer(log(ug_na_g_leaf_powder) ~ species*trt + (1|pair_id), data = leaf)
summary(leaf.lmer.ln.int)
Anova(leaf.lmer.ln.int, type = "III")
emmeans(leaf.lmer.ln.int, revpairwise ~ trt|species, adjust = "tukey")
(exp(3.092) - 1) * 100 # about 2102% more Na in C. stoebe in sodium plots
(exp(0.872) - 1) * 100 # about 139% more Na in M. fistulosa in sodium plots
(exp(3.123) - 1) * 100 # about 2171% more Na in P. digitalis in sodium plots

## Nectar ----
# Can't do stats for 2024 nectar :(

