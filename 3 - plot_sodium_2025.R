# TITLE: Sodium plots - Sodium sampling 2025
# AUTHOR: Ashley Darst
# DATA INPUT: 2025 Nectar and leaf sodium quantification
# DATA OUTPUT: Analyses and main figure 1
# DATE: 2025-12-03
# DESCRIPTION: Clean nectar and leaf sampling data, conduct analysis, and create figure

# ***MUST RUN plot_sodium_2024.R FIRST***

# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(car)
library(sjPlot)
library(bbmle)
library(ggforce)
library(patchwork)
library(cowplot)
library(emmeans)

# Data cleaning ----
# Load data
leaf_meta <- read.csv("data/leaf_2025.csv")
nectar_meta <- read.csv("data/nectar_2025.csv")
icp <- read.csv("data/2025_sodium.csv")

# Clean datasets
leaf_meta$sample <- paste(leaf_meta$plot_id, leaf_meta$leaf_id,sep="_")
leaf_meta$type <- "leaf"

sodium <- icp %>%
  select(c(sample, SampleWeight, Na23))

nectar_meta$nectar_id <- as.character(nectar_meta$nectar_id)
nectar_meta$type <- "nectar"
nectar_meta <- rename(nectar_meta, sample = nectar_id)
nectar_meta <- nectar_meta %>%
  select(-c(empty_weight_mg, full_weight_mg, time))

# Merge datasets
meta <- bind_rows(nectar_meta, leaf_meta)
data <- full_join(sodium, meta)

# Add pair id
data$pair_id <- sub("_[^_]+$", "", data$plot_id)

# Merge a different way for correlations
data2 <- data %>%
  select(c(Na23, date, plot_id, pair_id, treatment, species, type))

data2 <- data2 %>%
  pivot_wider(names_from = type, values_from = Na23)

# Average the leaf and nectar that were sampled twice
data_av <- data2 %>%
  group_by(plot_id, pair_id, treatment, species) %>%
  summarize(nectar = mean(nectar),
            leaf = mean(leaf))

# Add log to leaf column
data_av$leaf_ln <- log(data_av$leaf)

# Add log to nectar column
data_av$nectar_ln <- log(data_av$nectar)
  
# Plots ----
## Nectar ----
nectar.plot <- data_av %>%
  ggplot(aes(x = treatment, y = nectar, col = species)) +
  geom_sina(alpha = 0.5) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge2(0.7)) +
  labs(y = "Nectar sodium (ppm)", x = "Treatment") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Monarda fistulosa" = "#d95f02", "Penstemon digitalis" = "#7570b3"))

# Summary ----
data2 %>%
  group_by(species) %>%
  summarize(N_L_ratio = mean(nectar/leaf)) # Ellen Welti paper to compare Penstemon

# Analysis ----

## Nectar ----
# Simple model
nectar.lmer <- lmer(nectar ~ treatment + (1|pair_id), data = data_av)
summary(nectar.lmer)
check_model(nectar.lmer) # not great

### Final model stats ----
nectar_sp.lmer <- lmer(nectar ~ treatment + species + (1|pair_id), data = data_av)
summary(nectar_sp.lmer) # pair_id explains no variance
Anova(nectar_sp.lmer, type = "II")
emmeans(nectar_sp.lmer, revpairwise ~ treatment|species, adjust = "tukey", infer = T)
check_model(nectar_sp.lmer) # better
D <- cooks.distance(nectar_sp.lmer)
which(D > 0.5) # None are highly influential

# Try interaction
nectar_int.lmer <- lmer(nectar ~ treatment*species + (1|pair_id), data = data_av)
summary(nectar_int.lmer) # pair_id explains little variance
Anova(nectar_int.lmer, type = "III")
check_model(nectar_sp.lmer) # same

# Compare models # No evidence for interaction
compare_performance(nectar.lmer, nectar_sp.lmer, nectar_int.lmer) # using AICc, simple model
anova(nectar.lmer, nectar_sp.lmer) # nectar_sp.lmer marginally better
anova(nectar_sp.lmer, nectar_int.lmer) # int not better

## Leaf ----
# Simple leaf model
leaf.lmer <- lmer(leaf ~ treatment + (1|pair_id), data = data_av)
summary(leaf.lmer) # pair explains no variance
check_model(leaf.lmer) # bad
data_av$leaf[which(cooks.distance(leaf.lmer)>0.5)] # 39.54425 influential
data_av$leaf[which(cooks.distance(leaf.lmer)>1)] # none highly influential

# Add species
leaf_sp.lmer <- lmer(leaf ~ treatment + species + (1|pair_id), data = data_av)
summary(leaf_sp.lmer) # pair explains no variance
check_model(leaf_sp.lmer) #bad

# Interaction
leaf_int.lmer <- lmer(leaf ~ treatment*species + (1|pair_id), data = data_av)
summary(leaf_int.lmer) # pair explains no variance
check_model(leaf_int.lmer) # bad

# Compare models # No evidence for interaction
compare_performance(leaf.lmer, leaf_sp.lmer, leaf_int.lmer) # using AICc, simple model
anova(leaf.lmer, leaf_sp.lmer) # species not better
anova(leaf.lmer, leaf_int.lmer) # int not better

# Log scale
# Try to improve diagnostics plot
leaf.lmer.ln <- lmer(leaf_ln ~ treatment + (1|pair_id), data = data_av)
summary(leaf.lmer.ln)
check_model(leaf.lmer.ln) # better

### Final model stats ----
leaf_sp.lmer.ln <- lmer(leaf_ln ~ treatment + species + (1|pair_id), data = data_av)
summary(leaf_sp.lmer.ln) # pair explains no variance
Anova(leaf_sp.lmer.ln, type = "II")
check_model(leaf_sp.lmer.ln) # ok
(exp(0.8754) - 1) * 100 # about 139%

leaf_int.lmer.ln <- lmer(log(leaf) ~ treatment*species + (1|pair_id), data = data_av)
summary(leaf_int.lmer.ln) # pair explains no variance
check_model(leaf_int.lmer.ln) # okay

## Correlation ----
cor.lmer <- lmer(nectar ~ leaf + (1|pair_id), data = data_av)
cor.lmer.reml <- lmer(nectar ~ leaf + (1|pair_id), data = data_av, REML = F)
summary(cor.lmer) # pair explains no variance
check_model(cor.lmer) # not great

cor_sp.lmer <- lmer(nectar ~ leaf + species + (1|pair_id), data = data_av)
cor_sp.lmer.reml <- lmer(nectar ~ leaf + species + (1|pair_id), data = data_av, REML = F)
summary(cor_sp.lmer) # pair explains no variance
check_model(cor_sp.lmer) # better

cor_int.lmer <- lmer(nectar ~ leaf*species + (1|pair_id), data = data_av)
cor_int.lmer.reml <- lmer(nectar ~ leaf*species + (1|pair_id), data = data_av, REML = F)
summary(cor_int.lmer) # pair explains no variance
check_model(cor_int.lmer)

# Try log scale to improve model diagnostics
cor.lmer.ln <- lmer(nectar_ln ~ leaf_ln + (1|pair_id), data = data_av)
cor.lmer.ln.reml <- lmer(nectar_ln ~ leaf_ln + (1|pair_id), data = data_av, REML = F)
summary(cor.lmer.ln) # pair explains no variance
check_model(cor.lmer.ln) # ok

### Final model stats ----
cor_sp.lmer.ln <- lmer(nectar_ln ~ leaf_ln + species + (1|pair_id), data = data_av)
cor_sp.lmer.ln.reml <- lmer(nectar_ln ~ leaf_ln + species + (1|pair_id), data = data_av, REML = F)
summary(cor_sp.lmer.ln) # pair explains no variance
Anova(cor_sp.lmer.ln, type = "II")
emmeans(cor_sp.lmer.ln, ~ leaf_ln|species, adjust = "tukey", infer = T)
check_model(cor_sp.lmer.ln) # better

cor_int.lmer.ln <- lmer(nectar_ln ~ leaf_ln*species + (1|pair_id), data = data_av)
cor_int.lmer.ln.reml <- lmer(nectar_ln ~ leaf_ln*species + (1|pair_id), data = data_av, REML = F)
summary(cor_int.lmer.ln) # pair explains no variance
Anova(cor_int.lmer.ln, type = "III")
check_model(cor_int.lmer.ln)

AICctab(cor.lmer.ln.reml, cor_sp.lmer.ln.reml, cor_int.lmer.ln.reml) # simple model by AICc
anova(cor.lmer.ln.reml, cor_sp.lmer.ln.reml) # species better
anova(cor_sp.lmer.ln.reml, cor_int.lmer.ln.reml) # int marginally better

# Plot model ----
## Correlation ----
cor.plot <- plot_model(cor_sp.lmer.ln,
           type = "pred",
           terms = c("leaf_ln", "species"),
           show.data = TRUE,
           title = "",
           jitter = c(0.1, 0),
           colors = c("#d95f02", "#7570b3")) +
  labs(x = "ln(Leaf sodium (ppm))", y = "ln(Nectar sodium (ppm))") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(legend.position = "none")

# Figure 1 ----
legend <- get_legend(leaf.plot)
leaf.plot <- leaf.plot + theme(legend.position = "none")

plot.nl <- soil.plot + leaf.plot + nectar.plot + cor.plot & plot_annotation(tag_levels = 'A')

plot.wl <- cowplot::plot_grid(
  legend,
  plot.nl,
  ncol = 1,
  rel_heights = c(0.05, 1))

pdf("plot_fig_1.pdf", width = 7, height = 7)
plot.wl
dev.off()


