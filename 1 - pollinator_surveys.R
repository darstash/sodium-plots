# TITLE: Sodium plots - Pollinator Surveys
# AUTHOR: Ashley Darst
# DATA INPUT: Pollinator and floral survey data
# DATA OUTPUT: Analyses and main figures 2-4
# DATE: 2024-06-26
# DESCRIPTION: Clean pollinator and floral survey data, conduct main analyses, and create figures

# Load libraries
library(tidyverse)
library(vegan)
library(lme4)
library(glmmTMB)
library(sjPlot)
library(lmerTest)
library(car)
library(emmeans)
library(bbmle)
library(bipartite)
library(DHARMa)
library(MuMIn)
library(lubridate)
library(ggpubr)
library(ggforce)
library(patchwork)
library(ggeffects)

# Data cleaning ----
# Load data
pollinator <- read.csv("data/pollinator_surveys.csv")
floral <- read.csv("data/floral_surveys.csv")

# Look at datasets
glimpse(pollinator)
glimpse(floral)

# Fix date
pollinator$date <- as.Date(pollinator$date)
floral$date <- as.Date(floral$date)

# Change trt to factor
pollinator$trt <- as.factor(pollinator$trt)
floral$trt <- as.factor(floral$trt)

# Add underscore to species
floral$plant_sp <- sub(" ", "_", floral$plant_sp)
pollinator$pollinator_sp <- gsub(" ", "_", pollinator$pollinator_sp)

# Add NA to time columns
pollinator$start_time[pollinator$start_time==""] <- NA
pollinator$end_time[pollinator$end_time==""] <- NA

# Check time range of surveys
unique(pollinator$start_time, na.rm=TRUE)
unique(pollinator$end_time, na.rm=TRUE)

# Check unique floral species
unique(floral$plant_sp)

# Floral wide
floral_wide <- floral %>%
  select(-c(photo_id, notes)) %>%
  group_by(plot_id, pair_id, trt, date, week, surveyor) %>%
  pivot_wider(names_from = plant_sp, values_from = flower_num)

# Replace NAs with zeros
floral_wide <- floral_wide %>%
  replace(is.na(.), 0)

# drop NA column
floral_wide <- floral_wide %>%
  select(-"NA")

# Sum rows for floral abundance
floral_wide <- floral_wide %>%
  ungroup() %>%
  rowwise() %>%
  mutate(floral_abundance = sum(c_across(Penstemon_digitalis:Verbena_urticifolia), na.rm = TRUE))

# Floral species richness
floral_wide <- floral_wide %>%
  rowwise() %>%
  mutate(floral_richness = specnumber(c_across(Penstemon_digitalis:Verbena_urticifolia)))

# Floral diversity using shannons diversity index
floral_wide <- floral_wide %>%
  rowwise() %>%
  mutate(floral_shannon = diversity(c_across(Penstemon_digitalis:Verbena_urticifolia)))

# Make subset of floral surveys by plot and week
floral_summary <- floral_wide %>%
  select(-c(date, surveyor:Verbena_urticifolia))

# Find unique pollinator groups
unique(pollinator$pollinator_sp)

# Make a column in pollinator with broader classes
# fly, wasp, green_sweat_bee, tiny_dark_bee, medium_dark_bee, striped_sweat_bee, hairy_belly_bee, cuckoo bee, carpenter_bee, honey_bee, metallic_hairy_belly_bee, and chap_leg_bee the same
# Lepidoptera = Polites_peckius, Phyciodes_selenis_tharos, Thymelicus_lineola, Pieris_rapae, skipper_sp., Epargyreus_clarus, Vanessa_virginiensis, Vanessa_atalanta, moth, Anatrytone_logan, Erynnis_baptisiae, and Vernia_verna
# Bombus = Bombus_impatiens, Bombus_sp., Bombus_bimaculatus, Bombus_perplexus, Bombus_vagans, Bombus_griseocollis, and Bombus_auricomus
# hairy belly bee : hairy belly bee and metallic hairy belly bee

pol_groups <- pollinator %>%
  mutate(pollinator_groups = ifelse((str_detect(pollinator_sp, "Bombus")), "Bombus", pollinator_sp)) %>%
  mutate(pollinator_groups = ifelse(pollinator_sp == "Polites_peckius" | pollinator_sp == "Phyciodes_selenis_tharos" | pollinator_sp == "Anatrytone_logan" | pollinator_sp == "Thymelicus_lineola" | pollinator_sp == "Pieris_rapae" | pollinator_sp == "skipper_sp." | pollinator_sp == "Epargyreus_clarus" | pollinator_sp == "Vanessa_virginiensis" | pollinator_sp == "Vanessa_atalanta" | pollinator_sp == "moth" | pollinator_sp == "Erynnis_baptisiae" | pollinator_sp == "Vernia_verna","Lepidoptera", pollinator_groups)) %>%
  mutate(pollinator_groups = ifelse(pollinator_sp == "hairy_belly_bee" | pollinator_sp == "metallic_hairy_belly_bee", "hairy_belly_bee", pollinator_groups))

# Check unique pollinator groups in pol_groups
unique(pol_groups$pollinator_groups)

# Summarize pollinator abundance for each plot/date by group 
pol_groups  <- pol_groups %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, pollinator_groups) %>%
  summarize(pol_abundance = n()) %>%
  mutate(pol_abundance = ifelse(is.na(pollinator_groups), 0, pol_abundance))

# Add zero when a group wasn't seen at a plot-date combo
# First add an id column for plot-date
pol_groups$unique_id <- paste(pol_groups$plot_id, pol_groups$date)

# Add rows for all the plot-date combos with zeros for each pollinator group
pol_groups_zero <- pol_groups %>%
   ungroup() %>%
   complete(nesting(unique_id, plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud), pollinator_groups, fill = list(pol_abundance = 0))

# Remove NA from being a pollinator group
pol_groups_zero <- pol_groups_zero %>%
  drop_na(pollinator_groups)

# Add floral metrics to pol_groups and pol_groups_zero
pol_groups_floral <- left_join(pol_groups, floral_summary)
pol_groups_zero_floral <- left_join(pol_groups_zero, floral_summary)

# Look at number of pollinators in each group
pol_groups_zero_floral %>%
  group_by(pollinator_groups) %>%
  summarize(sum = sum(pol_abundance))

# Remove groups with less than 20 observations
# bee, chap_leg_bee, cuckoo_bee
pol_groups_filter <- pol_groups_zero_floral %>%
  filter(pollinator_groups != "bee" & pollinator_groups != "chap_leg_bee" & pollinator_groups != "cuckoo_bee")

# Check filter
pol_groups_filter %>%
  group_by(pollinator_groups) %>%
  summarize(sum = sum(pol_abundance))

# Pollinator wide
pol_wide <- pollinator %>%
  select(-c(flower, flower2, photo_id, notes)) %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, pollinator_sp) %>%
  summarize(sum = n()) %>%
  pivot_wider(names_from = pollinator_sp, values_from = sum)

# Duplicated date-plot combinations # None
sum(duplicated(pol_wide[,c('plot_id','date')]))

# Replace NAs with zeros
pol_wide <- pol_wide %>%
  mutate_at(vars("green_sweat_bee":"Bombus_vagans"), ~replace_na(., 0))

# drop NA column
pol_wide <- pol_wide %>%
  select(-"NA")

# Sum rows for pollinator abundance
pol_wide <- pol_wide %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pol_abundance = sum(c_across(green_sweat_bee:Bombus_vagans)))

# Pollinator species richness # Not species richness, group richness
pol_wide <- pol_wide %>%
  rowwise() %>%
  mutate(pol_richness = specnumber(c_across(green_sweat_bee:Bombus_vagans)))

# Pollinator diversity using shannons diversity index
pol_wide <- pol_wide %>%
  rowwise() %>%
  mutate(pol_shannon = diversity(c_across(green_sweat_bee:Bombus_vagans)))

# Join floral metrics with pollinator dataset
pol_floral <- left_join(pol_wide, floral_summary)

# Get total surveying time per plot
pol_wide %>%
  group_by(plot_id) %>%
  summarize(total = sum(active_time)/60)

# Get mean surveying time per plot
pol_wide %>%
  group_by(plot_id) %>%
  summarize(total = sum(active_time)/60) %>%
  ungroup() %>%
  summarize(mean = mean(total))

# Convert date to format for ou temporal autocorrelatiom
pol_floral <- (pol_floral
        %>% arrange(date)
        %>% mutate(times = lubridate::decimal_date(date) %% 1)
        %>% ungroup())

# Convert date to format for ou temporal autocorrelatiom
pol_groups_filter <- (pol_groups_filter
               %>% arrange(date)
               %>% mutate(times = lubridate::decimal_date(date) %% 1)
               %>% ungroup())

# Make pollinator groups a factor
pol_groups_filter$pollinator_groups <- as.factor(pol_groups_filter$pollinator_groups)

# Calculate polliantor group richness
# Pivot pol_groups wider (NOT filtered < 20 individuals)
pol_groups_wide <- pol_groups_zero_floral %>%
  group_by(plot_id, pair_id, trt, date, week, surveyor) %>%
  pivot_wider(names_from = pollinator_groups, values_from = pol_abundance)

# Pollinator species richness # Not species richness, group richness
pol_groups_wide <- pol_groups_wide %>%
  rowwise() %>%
  mutate(pol_group_richness = specnumber(c_across(bee:wasp)))


# Exploratory plots ---- 
# Pollinator abundance vs floral abundance
pol_floral %>%
  ggplot(aes(x = floral_abundance, y = pol_abundance, col = trt)) +
  geom_point() +
  geom_smooth(method = lm)

# Pollinator diversity vs floral diversity
pol_floral %>%
  ggplot(aes(x = floral_shannon, y = pol_shannon, col = trt)) +
  geom_point() +
  geom_smooth(method = lm)

# Pollinator abundance salty vs control
pol_wide %>%
  ggplot(aes(x = week, y = pol_abundance, col = trt)) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2))

# Pollinator diversity salty vs control
pol_wide %>%
  ggplot(aes(x = week, y = pol_shannon, col = trt)) +
  geom_jitter() +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2))

# Pollinator abundance in salty vs control by group
pol_groups %>%
  ggplot(aes(x = pollinator_groups, y = pol_abundance, col = trt)) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(0.2)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2))

pol_groups_floral %>%
  ggplot(aes(x = floral_abundance, y = pol_abundance, col = pollinator_groups)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~trt)

# Pollinator abundance in salty vs control by group with zeros added
pol_groups_filter %>%
  ggplot(aes(x = pollinator_groups, y = pol_abundance, col = trt)) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(0.2)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2))

pol_groups_filter %>%
  ggplot(aes(x = pollinator_groups, y = pol_abundance, col = trt)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Models ----
## Abundance model ----
# Variance >> mean, so we should use negative binomial
pol_floral %>%
  group_by(trt) %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

# Not zero inflated, probably
barplot(table(cut(pol_floral$pol_abundance, c(0, seq(1, 70, 1)), right = FALSE)), space = 0)

pol_abun.pois <- glmmTMB(pol_abundance ~ trt*week + scale(floral_abundance) + (1|pair_id) + (1|plot_id:pair_id), data = pol_floral, family = "poisson")
pol_abun.nb <- glmmTMB(pol_abundance ~ trt*week + scale(floral_abundance) + (1|pair_id) + (1|plot_id:pair_id), data = pol_floral, family = "nbinom2")
pol_abun.nb.zi3 <- glmmTMB(pol_abundance ~ trt*week + scale(floral_abundance) + (1|pair_id) + (1|plot_id:pair_id), data = pol_floral, family = "nbinom2", ziformula = ~ floral_abundance + week)
pol_abun.nb.add.zi3 <- glmmTMB(pol_abundance ~ trt + week + scale(floral_abundance) + (1|pair_id) + (1|plot_id:pair_id), data = pol_floral, family = "nbinom2", ziformula = ~ floral_abundance + week)

# Add floral richness
pol_abun.div <- glmmTMB(pol_abundance ~ trt*week + scale(floral_abundance) + floral_richness + (1|pair_id) + (1|plot_id:pair_id), data = pol_floral, family = "nbinom2", ziformula = ~ floral_abundance + week)
pol_abun.div.add <- glmmTMB(pol_abundance ~ trt + week + scale(floral_abundance) + floral_richness + (1|pair_id) + (1|plot_id:pair_id), data = pol_floral, family = "nbinom2", ziformula = ~ floral_abundance + week)

anova(pol_abun.pois, pol_abun.nb) # nb better
anova(pol_abun.nb.add.zi3, pol_abun.nb.zi3) # additive better
anova(pol_abun.div, pol_abun.nb.zi3) # diversity model better!
anova(pol_abun.div.add, pol_abun.div) # additive diversity model better

AICctab(pol_abun.pois, pol_abun.nb, pol_abun.nb.zi3, pol_abun.nb.add.zi3, pol_abun.div.add, pol_abun.div)

# Interaction model with diversity
simres <- simulateResiduals(pol_abun.div)
plot(simres) # Residual vs predicted not great

# Model statistics
summary(pol_abun.div) # Interaction is not significant, drop for better additive model
Anova(pol_abun.div, type = "III")

### Final model stats ----
summary(pol_abun.div.add)
Anova(pol_abun.div.add, type = "II")
confint(pol_abun.div.add)
simres <- simulateResiduals(pol_abun.div.add)
plot(simres) # not perfect

## Richness model ----
# Variance ~= mean, so we should use poisson
pol_groups_wide %>%
  group_by(trt) %>%
  summarize(mean.count = mean(pol_group_richness), 
            var.count = var(pol_group_richness))

# Not zero inflated, probably
barplot(table(cut(pol_groups_wide$pol_group_richness, c(0, seq(1, 8, 1)), right = FALSE)), space = 0)

# Add floral richness
pol_group_rich.add <- glmmTMB(pol_group_richness ~ trt + week + scale(floral_abundance) + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_wide, family = "poisson")
pol_group_rich.add.div <- glmmTMB(pol_group_richness ~ trt + week + scale(floral_abundance) + floral_richness + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_wide, family = "poisson")
pol_group_rich.add.div.zi <- glmmTMB(pol_group_richness ~ trt + week + scale(floral_abundance) + floral_richness + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_wide, family = "poisson", ziformula = ~ floral_abundance + week)
pol_group_rich.div.zi <- glmmTMB(pol_group_richness ~ trt*week + scale(floral_abundance) + floral_richness + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_wide, family = "poisson", ziformula = ~ floral_abundance + week)

anova(pol_group_rich.add, pol_group_rich.add.div) # richness better
anova(pol_group_rich.add.div, pol_group_rich.add.div.zi) # zi better

AICctab(pol_group_rich.add, pol_group_rich.add.div, pol_group_rich.add.div.zi, pol_group_rich.div.zi)

# Look at interaction model
summary(pol_group_rich.div.zi) # drop interaction for better fitting model

# Model check
simres4 <- simulateResiduals(pol_group_rich.add.div.zi)
plot(simres4) # Bad
testOutliers(simres4, type = "bootstrap") # looks fine

### Final model stats ----
summary(pol_group_rich.add.div.zi)
Anova(pol_group_rich.add.div.zi, type = "II")


## Pollinator group model ----
# Using filtered dataset removing groups where n < 20
# Variance > mean generally
pol_groups_filter %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

# Variance vs mean depends on the group...
pol_groups_filter %>%
  group_by(pollinator_groups) %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

# Zero inflated when including all groups
barplot(table(cut(pol_groups_filter$pol_abundance, c(0, seq(1, 50, 4)), right = FALSE)), space = 0)

pol_filter.pois <- glmmTMB(pol_abundance ~ trt + pollinator_groups + scale(floral_abundance) + week + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_filter, family = "poisson")
pol_filter.nb <- glmmTMB(pol_abundance ~ trt + pollinator_groups + scale(floral_abundance) + week + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_filter, family = "nbinom2")
pol_filter.nb.int.zi3 <- glmmTMB(pol_abundance ~ trt*pollinator_groups + scale(floral_abundance) + week + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_filter, family = "nbinom2", ziformula =  ~ floral_abundance + week)
pol_filter.nb.int2.zi3 <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(floral_abundance) + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_filter, family = "nbinom2", ziformula =  ~ floral_abundance + week)

# Try model with floral diversity
pol_filter.div <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(floral_abundance) + floral_richness + week + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_filter, family = "nbinom2", ziformula =  ~ floral_abundance + week)
pol_filter.div2 <- glmmTMB(pol_abundance ~ trt*pollinator_groups + pollinator_groups*week + scale(floral_abundance) + floral_richness + week + (1|pair_id) + (1|plot_id:pair_id), data = pol_groups_filter, family = "nbinom2", ziformula =  ~ floral_abundance + week)

anova(pol_filter.pois, pol_filter.nb) # nb better
anova(pol_filter.nb.int.zi3, pol_filter.nb.int2.zi3) # three way interaction better, makes biological sense
anova(pol_filter.div, pol_filter.nb.int2.zi3) # adding richness is way better
anova(pol_filter.div2, pol_filter.div) # three way diversity model better

AICctab(pol_filter.pois, pol_filter.nb, pol_filter.nb.int.zi3, pol_filter.nb.int2.zi3, pol_filter.div, pol_filter.div2)
# Anova and AICc conflict

# Model checks
simres2 <- simulateResiduals(pol_filter.div)
plot(simres2) # Looks ok, some outliers in residual vs predicted
testOutliers(simres2, type = "bootstrap") # not significant here

### Final model stats ----
# We care about the three way interaction and have biological reason to include it, so we will include it since on border
summary(pol_filter.div)
Anova(pol_filter.div, type = "III") # Three way interaction significant
results <- Anova(pol_filter.div, type = "III")
emmip(pol_filter.div, trt ~ week | pollinator_groups, cov.reduce = range)

joint_tests(pol_filter.div, by = "pollinator_groups")
emmeans(pol_filter.div, revpairwise ~ trt|pollinator_groups, adjust = "tukey")
(exp(0.6800) - 1) * 100 # about 97% more honey bees in Na plots
(exp(-0.5812) - 1) * 100 # about 44% less wasps in control plots
emtrends(pol_filter.div, revpairwise ~ trt|pollinator_groups, adjust = "tukey", var = "week")

# Check model without trt*week interaction
# Model checks
simres2 <- simulateResiduals(pol_filter.div2)
plot(simres2) # Looks ok

# Model statistics
summary(pol_filter.div2)
Anova(pol_filter.div2, type = "III") # Two way interactions significant
emmip(pol_filter.div2, trt ~ week | pollinator_groups, cov.reduce = range)
joint_tests(pol_filter.div2, by = "pollinator_groups") # similar to three way interaction model

## Specific floral species ----
# Pivot pollinator dataset longer to make each flower visit a row
pol_network2 <- pollinator %>% 
  pivot_longer(cols = starts_with("flower"),
               names_to = "flower_num",
               names_prefix = "flower",
               values_to = "flower_sp",
               values_drop_na = F)
pol_network2$flower_sp[pol_network2$flower_sp == ""] <- NA
pol_network2 <- pol_network2 %>%
  filter(!(flower_num  %in% c("2", "3") & is.na(flower_sp)))

# Convert to pol groups
pol_network_group <- pol_network2 %>%
  mutate(pollinator_groups = ifelse((str_detect(pollinator_sp, "Bombus")), "Bombus", pollinator_sp)) %>%
  mutate(pollinator_groups = ifelse(pollinator_sp == "Polites_peckius" | pollinator_sp == "Phyciodes_selenis_tharos" | pollinator_sp == "Anatrytone_logan" | pollinator_sp == "Thymelicus_lineola" | pollinator_sp == "Pieris_rapae" | pollinator_sp == "skipper_sp." | pollinator_sp == "Epargyreus_clarus" | pollinator_sp == "Vanessa_virginiensis" | pollinator_sp == "Vanessa_atalanta" | pollinator_sp == "moth" | pollinator_sp == "Erynnis_baptisiae" | pollinator_sp == "Vernia_verna","Lepidoptera", pollinator_groups)) %>%
  mutate(pollinator_groups = ifelse(pollinator_sp == "hairy_belly_bee" | pollinator_sp == "metallic_hairy_belly_bee", "hairy_belly_bee", pollinator_groups))

# Filter out species with > 20 
pol_network_group <- pol_network_group %>%
  filter((pollinator_groups != "bee" & pollinator_groups != "chap_leg_bee" & pollinator_groups != "cuckoo_bee") %>% replace_na(TRUE))

pol_network_group <- pol_network_group %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, pollinator_groups, flower_sp) %>%
  summarize(pol_abundance = n()) %>%
  mutate(pol_abundance = ifelse(is.na(pollinator_groups), 0, pol_abundance))

# Calculate pollinator group richness
pol_network_group <- pol_network_group %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp) %>%
  mutate(pol_group_richness = n_distinct(pollinator_groups))
pol_network_group$pol_group_richness[pol_network_group$pol_abundance == 0] <- 0

# Format floral info
floral_summary2 <- floral_wide %>%
  dplyr::select(-c(date, surveyor))

# Add floral information # Dropping times we never encountered a pollinator # Add zeros for every pol/flower combo?
pol_network_floral <- left_join(pol_network_group, floral_summary2)

# What are the most abundant flowers?
colSums(floral_summary2[,c(5:32)])
# Securigera_varia: 8370
# Erigeron_annuus: 8012
# Monarda_fistulosa: 6374
# Penstemon_digitalis: 2710
# Hypericum_perforatum: 1438
# Hesperis_matronalis: 680

# What are most interacted with flowers?
pol_network_floral %>%
  group_by(flower_sp) %>%
  summarize(n = n()) %>%
  arrange(-n)
# Monarda fistulosa         656
# Securigera varia          227
# Hypericum perforatum      226
# Penstemon digitalis       222
# Erigeron annuus           208
# Achillea millefolium      189

# Subset by top five flowers
pol_network_floral_filter <- pol_network_floral %>%
  filter(flower_sp == "Securigera varia" | flower_sp == "Erigeron annuus" | flower_sp == "Monarda fistulosa" | flower_sp == "Penstemon_digitalis" | flower_sp == "Hypericum perforatum")

# Model top five flowers # Rank deficient
# pol_abun_sp_top <- glmmTMB(pol_abundance ~ trt*pollinator_groups*flower_sp + scale(floral_abundance) + floral_richness + (1 + week|plot_id/pair_id), data = pol_network_floral_filter, family = "nbinom2")

# Need to make a model for each flower species

### Securigera varia ----
pol_network_floral_sec <- pol_network_floral %>%
  filter(flower_sp == "Securigera varia" | (Securigera_varia > 0 & is.na(flower_sp)))

# Look at number of pollinators in each group
pol_network_floral_sec %>%
  group_by(pollinator_groups) %>%
  summarize(sum = sum(pol_abundance))

# Remove pol groups with less than 20 interactions
pol_network_floral_sec_filter <- pol_network_floral_sec %>%
  filter((pollinator_groups == "Bombus" | pollinator_groups == "fly" | pollinator_groups == "honey_bee" | pollinator_groups == "tiny_dark_bee") %>% replace_na(TRUE))

# Fill with zeros
pol_network_floral_sec_filter$unique_id <- paste(pol_network_floral_sec_filter$plot_id, pol_network_floral_sec_filter$date)

pol_network_floral_sec_filter_zero <- pol_network_floral_sec_filter %>%
  ungroup() %>%
  complete(nesting(unique_id, plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, Securigera_varia, floral_abundance), pollinator_groups, fill = list(pol_abundance = 0)) %>%
  drop_na(pollinator_groups)
  
# Variance >> mean for Bombus and honey bee, NB?
pol_network_floral_sec_filter_zero %>%
  group_by(pollinator_groups) %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

barplot(table(cut(pol_network_floral_sec_filter_zero$pol_abundance, c(0, seq(1, 50, 4)), right = FALSE)), space = 0)

# Final model
pol_abun_sec_var_z.int <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(Securigera_varia) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_sec_filter_zero, family = "nbinom2")

# Model checks
simres2 <- simulateResiduals(pol_abun_sec_var_z.int)
plot(simres2) # Looks good!

summary(pol_abun_sec_var_z.int)
Anova(pol_abun_sec_var_z.int, type = "III")

# Richness
# Drop NA values to get richness 
pol_network_floral_sec_rich <- pol_network_floral_sec_filter %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, pol_group_richness, Securigera_varia, floral_abundance) %>%
  summarize(pol_group_richness = first(pol_group_richness))

# Mean > variance
pol_network_floral_sec_rich %>%
  group_by(trt) %>%
  summarize(mean.count = mean(pol_group_richness), 
            var.count = var(pol_group_richness))

barplot(table(cut(pol_network_floral_sec_rich$pol_group_richness, c(0, seq(1, 7, 1)), right = FALSE)), space = 0)

pol_abun_sec_var.rich <- glmmTMB(pol_group_richness ~ trt*week + scale(Securigera_varia) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_sec_rich, family = "poisson")
pol_abun_sec_var.rich.add <- glmmTMB(pol_group_richness ~ trt + week + scale(Securigera_varia) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_sec_rich, family = "poisson")

AICctab(pol_abun_sec_var.rich, pol_abun_sec_var.rich.add)
anova(pol_abun_sec_var.rich.add, pol_abun_sec_var.rich)

summary(pol_abun_sec_var.rich.add)
Anova(pol_abun_sec_var.rich.add, type = "II")

### Monarda fistulosa ----
pol_network_floral_mon <- pol_network_floral %>%
  filter(flower_sp == "Monarda fistulosa" | (Monarda_fistulosa > 0 & is.na(flower_sp)))

# Look at number of pollinators in each group
pol_network_floral_mon %>%
  group_by(pollinator_groups) %>%
  summarize(sum = sum(pol_abundance))

# Remove pol groups with less than 20 interactions
pol_network_floral_mon_filter <- pol_network_floral_mon %>%
  filter((pollinator_groups == "Bombus" | pollinator_groups == "fly" | pollinator_groups == "green_sweat_bee" | pollinator_groups == "honey_bee" | pollinator_groups == "striped_sweat_bee" | pollinator_groups == "tiny_dark_bee" | pollinator_groups == "wasp") %>% replace_na(TRUE))

# Fill with zeros
pol_network_floral_mon_filter$unique_id <- paste(pol_network_floral_mon_filter$plot_id, pol_network_floral_mon_filter$date)

pol_network_floral_mon_filter_zero <- pol_network_floral_mon_filter %>%
  ungroup() %>%
  complete(nesting(unique_id, plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, Monarda_fistulosa, floral_abundance), pollinator_groups, fill = list(pol_abundance = 0)) %>%
  drop_na(pollinator_groups)

# Variance >> mean for most groups
pol_network_floral_mon_filter_zero %>%
  group_by(pollinator_groups) %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

barplot(table(cut(pol_network_floral_mon_filter_zero$pol_abundance, c(0, seq(1, 50, 4)), right = FALSE)), space = 0)

pol_abun_mon <- glmmTMB(pol_abundance ~ trt*pollinator_groups + scale(Monarda_fistulosa) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_mon_filter_zero, family = "nbinom2")
pol_abun_mon.pois <- glmmTMB(pol_abundance ~ trt*pollinator_groups + scale(Monarda_fistulosa) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_mon_filter_zero, family = "poisson")
pol_abun_mon.int <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(Monarda_fistulosa) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_mon_filter_zero, family = "nbinom2")

AICctab(pol_abun_mon, pol_abun_mon.pois, pol_abun_mon.int)

# Model checks
simres2 <- simulateResiduals(pol_abun_mon.int)
plot(simres2) # Ok
testOutliers(simres2, type = "bootstrap")

summary(pol_abun_mon.int)
Anova(pol_abun_mon.int, type = "III")

# Richness
# Drop NA values to get richness 
pol_network_floral_mon_rich <- pol_network_floral_mon_filter %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, pol_group_richness, Monarda_fistulosa, floral_abundance) %>%
  summarize(pol_group_richness = first(pol_group_richness))

# Mean > variance
pol_network_floral_mon_rich %>%
  group_by(trt) %>%
  summarize(mean.count = mean(pol_group_richness), 
            var.count = var(pol_group_richness))

barplot(table(cut(pol_network_floral_mon_rich$pol_group_richness, c(0, seq(1, 7, 1)), right = FALSE)), space = 0)

pol_abun_mon.rich <- glmmTMB(pol_group_richness ~ trt*week + scale(Monarda_fistulosa) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_mon_rich, family = "poisson")
pol_abun_mon.rich.add <- glmmTMB(pol_group_richness ~ trt + week + scale(Monarda_fistulosa) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_mon_rich, family = "poisson")

AICctab(pol_abun_mon.rich, pol_abun_mon.rich.add)
anova(pol_abun_mon.rich.add, pol_abun_mon.rich)

summary(pol_abun_mon.rich.add)
Anova(pol_abun_mon.rich.add, type = "II")

### Penstemon digitalis ----
pol_network_floral_pen <- pol_network_floral %>%
  filter(flower_sp == "Penstemon digitalis" | (Penstemon_digitalis > 0 & is.na(flower_sp)))

# Look at number of pollinators in each group
pol_network_floral_pen %>%
  group_by(pollinator_groups) %>%
  summarize(sum = sum(pol_abundance))

# Remove pol groups with less than 20 interactions
pol_network_floral_pen_filter <- pol_network_floral_pen %>%
  filter((pollinator_groups == "Bombus" | pollinator_groups == "fly" | pollinator_groups == "green_sweat_bee" | pollinator_groups == "tiny_dark_bee") %>% replace_na(TRUE))

# Fill with zeros
pol_network_floral_pen_filter$unique_id <- paste(pol_network_floral_pen_filter$plot_id, pol_network_floral_pen_filter$date)

pol_network_floral_pen_filter_zero <- pol_network_floral_pen_filter %>%
  ungroup() %>%
  complete(nesting(unique_id, plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, Penstemon_digitalis, floral_abundance), pollinator_groups, fill = list(pol_abundance = 0)) %>%
  drop_na(pollinator_groups)

# Variance >> mean for some groups
pol_network_floral_pen_filter_zero %>%
  group_by(pollinator_groups) %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

barplot(table(cut(pol_network_floral_pen_filter_zero$pol_abundance, c(0, seq(1, 50, 4)), right = FALSE)), space = 0)

pol_abun_pen.int <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(Penstemon_digitalis) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_pen_filter_zero, family = "nbinom2")
pol_abun_pen.int.pois <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(Penstemon_digitalis) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_pen_filter_zero, family = "poisson")
pol_abun_pen.add <- glmmTMB(pol_abundance ~ trt + pollinator_groups + week + scale(Penstemon_digitalis) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_pen_filter_zero, family = "nbinom2")

AICctab(pol_abun_pen.int, pol_abun_pen.int.pois, pol_abun_pen.add)
anova(pol_abun_pen.add, pol_abun_pen.int) # int better

# Model checks
simres2 <- simulateResiduals(pol_abun_pen.int)
plot(simres2) # Good!

simres2 <- simulateResiduals(pol_abun_pen.add)
plot(simres2) # Bad, go with interactive model (that's our hypothesis anyways)

summary(pol_abun_pen.int)
Anova(pol_abun_pen.int, type = "III")
emmeans(pol_abun_pen.int, revpairwise ~ trt)
(exp(0.609) - 1) * 100 # about 84% more in Na plots

# Richness
# Drop NA values to get richness 
pol_network_floral_pen_rich <- pol_network_floral_pen_filter %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, pol_group_richness, Penstemon_digitalis, floral_abundance) %>%
  summarize(pol_group_richness = first(pol_group_richness))

# Mean > variance
pol_network_floral_pen_rich %>%
  group_by(trt) %>%
  summarize(mean.count = mean(pol_group_richness), 
            var.count = var(pol_group_richness))

barplot(table(cut(pol_network_floral_pen_rich$pol_group_richness, c(0, seq(1, 7, 1)), right = FALSE)), space = 0)

pol_abun_pen.rich <- glmmTMB(pol_group_richness ~ trt*week + scale(Penstemon_digitalis) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_pen_rich, family = "poisson")
pol_abun_pen.rich.add <- glmmTMB(pol_group_richness ~ trt + week + scale(Penstemon_digitalis) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_pen_rich, family = "poisson")

AICctab(pol_abun_pen.rich, pol_abun_pen.rich.add)
anova(pol_abun_pen.rich.add, pol_abun_pen.rich)

summary(pol_abun_pen.rich.add)
Anova(pol_abun_pen.rich.add, type = "II")
(exp(0.33043) - 1) * 100 # about 39% more in Na plots

### Hypericum perforatum ----
pol_network_floral_hyp <- pol_network_floral %>%
  filter(flower_sp == "Hypericum perforatum" | (Hypericum_perforatum > 0 & is.na(flower_sp)))

# Look at number of pollinators in each group
pol_network_floral_hyp %>%
  group_by(pollinator_groups) %>%
  summarize(sum = sum(pol_abundance))

# Remove pol groups with less than 20 interactions
pol_network_floral_hyp_filter <- pol_network_floral_hyp %>%
  filter((pollinator_groups == "Bombus" | pollinator_groups == "fly" | pollinator_groups == "tiny_dark_bee") %>% replace_na(TRUE))

# Fill with zeros
pol_network_floral_hyp_filter$unique_id <- paste(pol_network_floral_hyp_filter$plot_id, pol_network_floral_hyp_filter$date)

pol_network_floral_hyp_filter_zero <- pol_network_floral_hyp_filter %>%
  ungroup() %>%
  complete(nesting(unique_id, plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, Hypericum_perforatum, floral_abundance), pollinator_groups, fill = list(pol_abundance = 0)) %>%
  drop_na(pollinator_groups)

# Variance >> mean for some groups
pol_network_floral_hyp_filter_zero %>%
  group_by(pollinator_groups) %>%
  summarize(mean.count = mean(pol_abundance), 
            var.count = var(pol_abundance))

barplot(table(cut(pol_network_floral_hyp_filter_zero$pol_abundance, c(0, seq(1, 50, 4)), right = FALSE)), space = 0)

pol_abun_hyp.int <- glmmTMB(pol_abundance ~ trt*pollinator_groups*week + scale(Hypericum_perforatum) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_hyp_filter_zero, family = "nbinom2")
pol_abun_hyp <- glmmTMB(pol_abundance ~ trt*pollinator_groups + week + scale(Hypericum_perforatum) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_hyp_filter_zero, family = "nbinom2")

AICctab(pol_abun_hyp.int, pol_abun_hyp)

# Model checks
simres2 <- simulateResiduals(pol_abun_hyp.int)
plot(simres2) # Good!

summary(pol_abun_hyp.int)
Anova(pol_abun_hyp.int, type = "III")

# Richness
# Drop NA values to get richness 
pol_network_floral_hyp_rich <- pol_network_floral_hyp_filter %>%
  group_by(plot_id, pair_id, trt, date, week, start_time, end_time, active_time, surveyor, temp, wind, cloud, flower_sp, pol_group_richness, Hypericum_perforatum, floral_abundance) %>%
  summarize(pol_group_richness = first(pol_group_richness))

# Mean > variance
pol_network_floral_hyp_rich %>%
  group_by(trt) %>%
  summarize(mean.count = mean(pol_group_richness), 
            var.count = var(pol_group_richness))

barplot(table(cut(pol_network_floral_hyp_rich$pol_group_richness, c(0, seq(1, 7, 1)), right = FALSE)), space = 0)

pol_abun_hyp.rich <- glmmTMB(pol_group_richness ~ trt*week + scale(Hypericum_perforatum) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_hyp_rich, family = "poisson")
pol_abun_hyp.rich.add <- glmmTMB(pol_group_richness ~ trt + week + scale(Hypericum_perforatum) + scale(floral_abundance) + (1|plot_id/pair_id), data = pol_network_floral_hyp_rich, family = "poisson")

AICctab(pol_abun_hyp.rich, pol_abun_hyp.rich.add)
anova(pol_abun_hyp.rich.add, pol_abun_hyp.rich)

# Model checks
simres2 <- simulateResiduals(pol_abun_hyp.rich.add)
plot(simres2) # ok

summary(pol_abun_hyp.rich.add)
Anova(pol_abun_hyp.rich.add, type = "II")


# Plot models ----
# Abundance model
pred.abun <- ggpredict(pol_abun.div.add, terms = "trt")
abun.plot <- ggplot() +
  geom_sina(data = model.frame(pol_abun.div.add), aes(x = trt, y = pol_abundance, color = trt, shape = trt), alpha = 0.1, orientation = "x") +
  geom_point(data = pred.abun, aes(x = x, y = predicted, color = x, shape = x), size = 3) +
  geom_errorbar(data = pred.abun, aes(x = x, ymin = conf.low, ymax = conf.high, color = x), width = 0) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Treatment", y = "Pollinator abundance", color = "Treatment") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15))

# Richness model
pred.rich <- ggpredict(pol_group_rich.add.div.zi, terms = "trt")
rich.plot <- ggplot() +
  geom_sina(data = model.frame(pol_group_rich.add.div.zi), aes(x = trt, y = pol_group_richness, color = trt, shape = trt), alpha = 0.1, orientation = "x") +
  geom_point(data = pred.rich, aes(x = x, y = predicted, color = x, shape = x), size = 3) +
  geom_errorbar(data = pred.rich, aes(x = x, ymin = conf.low, ymax = conf.high, color = x), width = 0) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Treatment", y = "Pollinator group richness", color = "Treatment") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15))

pred.groups <- ggpredict(pol_filter.div, terms = c("pollinator_groups", "trt"))
group.plot <- ggplot() +
  geom_point(data = pred.groups, aes(x = x, y = predicted, color = group, shape = group), size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = pred.groups, aes(x = x, ymin = conf.low, ymax = conf.high, color = group), width = 0, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Pollinator groups", y = "Pollinator abundance", color = "Treatment", shape = "Treatment") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("bumble bee", "carpenter bee", "fly", "green sweat bee", "hairy belly bee", "honey bee", "lepidopteran", "medium dark bee", "striped sweat bee", "tiny dark bee", "wasp"))

## Figure 2 ----
pdf("plot_fig_2_report.pdf", width = 7, height = 7)
(abun.plot + rich.plot) / group.plot + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A') 
dev.off()

## Figure S3 ----
emmip.plot <- emmip(pol_filter.div, trt ~ week | pollinator_groups,
            cov.reduce = range, plotit = FALSE)
ggplot(emmip.plot, aes(x = week, y = yvar, color = trt, linetype = trt)) +
  geom_line() +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  scale_linetype_manual(values = c("control" = 1,"sodium" = 2)) +
  theme_bw() +
  labs(color = "Treatment", linetype = "Treatment", y = "Linear prediction") +
  facet_wrap(~ pollinator_groups, labeller = as_labeller(
    c(Bombus = "bumble bee", carpenter_bee = "carpenter bee", fly = "fly", green_sweat_bee = "green sweat bee", hairy_belly_bee = "hairy belly bee", honey_bee = "honey bee", Lepidoptera = "lepidopteran", medium_dark_bee = "medium dark bee", striped_sweat_bee = "striped sweat bee", tiny_dark_bee = "tiny dark bee", wasp = "wasp"))) +
  theme(text = element_text(size = 15))

## Specific flower models ----
## Securigera varia # Nectar on outside of calyx, pollen
pred.sec.abun <- ggpredict(pol_abun_sec_var_z.int, terms = c("pollinator_groups", "trt"))
sec_var.plot <- ggplot() +
  geom_point(data = pred.sec.abun, aes(x = x, y = predicted, color = group, shape = group), size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = pred.sec.abun, aes(x = x, ymin = conf.low, ymax = conf.high, color = group), width = 0, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  theme_bw() +
  labs(y = "Pollinator abundance", color = "Treatment", x = "", shape = "Treatment", title = "Securigera varia") +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 12, face = "italic")) +
  scale_x_discrete(labels = c("bumble bee", "fly", "honey bee", "tiny dark bee"))

pred.sec <- ggpredict(pol_abun_sec_var.rich.add, terms = "trt")
sec_var_rich.plot <- ggplot() +
  geom_sina(data = model.frame(pol_abun_hyp.rich.add), aes(x = trt, y = pol_group_richness, color = trt, shape = trt), alpha = 0.2, orientation = "x") +
  geom_point(data = pred.sec, aes(x = x, y = predicted, color = x, shape = x), size = 3) +
  geom_errorbar(data = pred.sec, aes(x = x, ymin = conf.low, ymax = conf.high, color = x), width = 0) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Treatment", y = "Pollinator richness", color = "Treatment", shape = "Treatment") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15))

# Monarda fistulosa # Nectar and pollen
mon.abun <- ggpredict(pol_abun_mon.int, terms = c("pollinator_groups", "trt"))
mon.plot <- ggplot() +
  geom_point(data = mon.abun, aes(x = x, y = predicted, color = group, shape = group), size = 3, position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mon.abun, aes(x = x, ymin = conf.low, ymax = conf.high, color = group), width = 0, position = position_dodge(width = 0.8)) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  theme_bw() +
  labs(y = "Pollinator abundance", color = "Treatment", x = "", shape = "Treatment", title = "Monarda fistulosa") +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 12, face = "italic")) +
  scale_x_discrete(labels = c("bumble bee", "fly", "green sweat bee", "honey bee", "striped sweat bee", "tiny dark bee", "wasp"))

pred.mon <- ggpredict(pol_abun_mon.rich.add, terms = "trt")
mon_rich.plot <- ggplot() +
  geom_sina(data = model.frame(pol_abun_hyp.rich.add), aes(x = trt, y = pol_group_richness, color = trt, shape = trt), alpha = 0.2, orientation = "x") +
  geom_point(data = pred.mon, aes(x = x, y = predicted, color = x, shape = x), size = 3) +
  geom_errorbar(data = pred.mon, aes(x = x, ymin = conf.low, ymax = conf.high, color = x), width = 0) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Treatment", y = "Pollinator richness", color = "Treatment", shape = "Treatment") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15))

# Penstemon digitalis # Nectar and pollen
pen.abun <- ggpredict(pol_abun_pen.int, terms = c("pollinator_groups", "trt"))
pen.plot <- ggplot() +
  geom_point(data = pen.abun, aes(x = x, y = predicted, color = group, shape = group), size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(data = pen.abun, aes(x = x, ymin = conf.low, ymax = conf.high, color = group), width = 0, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  theme_bw() +
  labs(y = "Pollinator abundance", color = "Treatment", x = "", shape = "Treatment", title = "Penstemon digitalis") +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 12, face = "italic")) +
  scale_x_discrete(labels = c("bumble bee", "fly", "green sweat bee", "tiny dark bee"))

pred.pen <- ggpredict(pol_abun_pen.rich.add, terms = "trt")
pen_rich.plot <- ggplot() +
  geom_sina(data = model.frame(pol_abun_hyp.rich.add), aes(x = trt, y = pol_group_richness, color = trt, shape = trt), alpha = 0.2, orientation = "x") +
  geom_point(data = pred.pen, aes(x = x, y = predicted, color = x, shape = x), size = 3) +
  geom_errorbar(data = pred.pen, aes(x = x, ymin = conf.low, ymax = conf.high, color = x), width = 0) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Treatment", y = "Pollinator richness", color = "Treatment", shape = "Treatment") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15))

# Hypericum perforatum # Only pollen
hyp.abun <- ggpredict(pol_abun_hyp.int, terms = c("pollinator_groups", "trt"))
hyp.plot <- ggplot() +
  geom_point(data = hyp.abun, aes(x = x, y = predicted, color = group, shape = group), size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = hyp.abun, aes(x = x, ymin = conf.low, ymax = conf.high, color = group), width = 0, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  theme_bw() +
  labs(y = "Pollinator abundance", color = "Treatment", x = "", shape = "Treatment", title = "Hypericum perforatum") +
  theme(text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 12, face = "italic")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("bumble bee", "fly", "tiny dark bee"))

pred.hyp <- ggpredict(pol_abun_hyp.rich.add, terms = "trt")
hyp_rich.plot <- ggplot() +
  geom_sina(data = model.frame(pol_abun_hyp.rich.add), aes(x = trt, y = pol_group_richness, color = trt, shape = trt), alpha = 0.2, orientation = "x") +
  geom_point(data = pred.hyp, aes(x = x, y = predicted, color = x, shape = x), size = 3) +
  geom_errorbar(data = pred.hyp, aes(x = x, ymin = conf.low, ymax = conf.high, color = x), width = 0) +
  scale_color_manual(values = c("control" = "#e41a1c","sodium" = "#377eb8")) +
  labs(x = "Treatment", y = "Pollinator richness", color = "Treatment", shape = "Treatment") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15))

## Figure 3 ----
pdf("plot_fig_4.1.pdf", width = 11, height = 6)
mon.plot + sec_var.plot + hyp.plot + pen.plot + mon_rich.plot + sec_var_rich.plot + hyp_rich.plot + pen_rich.plot + plot_layout(ncol = 4) + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A') 
dev.off()


# NMDS analysis ----

set.seed(100)

# Subset community data from pol_floral
comm.data <- pol_floral[rowSums(pol_floral[13:44])>0,]

# Run NMDS
pol.mds <- metaMDS(comm.data[13:44])
pol.mds  # Stress = 0.1438, ok!
plot(pol.mds, type="t")

# Convert to data frame (only extracting site NMDSs)
mds.scores <- as.data.frame(scores(pol.mds, display="site"))  

# Add plot info back in
mds.scores$plot_id <- comm.data$plot_id
mds.scores$pair_id <- comm.data$pair_id
mds.scores$trt <- comm.data$trt
mds.scores$date <- comm.data$date
mds.scores$floral_abundance <- comm.data$floral_abundance
mds.scores$floral_richness <- comm.data$floral_richness

# Plot NMDS
# Weird outlier...
ggplot(data=mds.scores) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=trt),level = 0.50) +
  geom_point(aes(x=NMDS1,y=NMDS2,shape=trt,colour=trt),size=4, alpha = 0.5)

# Stress plot # Bad
stressplot(pol.mds, lwd = 0.5)

# Remove Cabbage White outlier
comm.data2 <- comm.data %>%
  filter(!(plot_id == "bird_3_c" & date == "2024-06-09"))

# Run NMDS
pol.mds2 <- metaMDS(comm.data2[13:44], k = 3) # need to up dimensions due to stress
pol.mds2  # Stress = 0.175, ok!
plot(pol.mds2, type="t")

# Convert to data frame (only extracting site NMDSs)
mds.scores2 <- as.data.frame(scores(pol.mds2, display="site"))  

# Add plot info back in
mds.scores2$plot_id <- comm.data2$plot_id
mds.scores2$pair_id <- comm.data2$pair_id
mds.scores2$trt <- comm.data2$trt
mds.scores2$date <- comm.data2$date
mds.scores2$floral_abundance <- comm.data2$floral_abundance
mds.scores2$floral_richness <- comm.data2$floral_richness

# Plot NMDS
ggplot(data=mds.scores2) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=trt)) +
  geom_point(aes(x=NMDS1,y=NMDS2,color=trt), alpha = 0.5) +
  theme_bw() +
  labs(color = "treatment") +
  scale_color_manual(values = c("#e41a1c", "#377eb8"))

# Plot mean values by plot id
mds.scores.mean2 <- mds.scores2 %>%
  group_by(plot_id, trt) %>%
  summarize(mean.NMDS1 = mean(NMDS1),
            mean.NMDS2 = mean(NMDS2))

pol.nmds.mean <- ggplot(data=mds.scores.mean2) + 
  stat_ellipse(aes(x=mean.NMDS1,y=mean.NMDS2,color=trt, linetype = trt)) +
  geom_point(aes(x=mean.NMDS1,y=mean.NMDS2,color=trt, shape = trt), size = 3, alpha = 0.5) +
  theme_bw() +
  labs(x = "NMDS1", y = "NMDS2", color = "Treatment", shape = "Treatment", linetype = "Treatment") +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  scale_linetype_manual(values = c(1, 2)) +
  coord_equal() +
  scale_y_continuous(breaks=c(-0.4, 0, 0.4)) +
  scale_x_continuous(breaks=c(-0.4, 0, 0.4))

# Stress plot # Better
stressplot(pol.mds2, lwd = 0.5)

# Model
# Create permutation structure to account for random effects
h.out <- with(comm.data2,
           how(plots = Plots(strata = plot_id),
               blocks = pair_id,
               nperm = 4999)
)

pol.adon.rmout <- adonis2(comm.data2[13:44] ~ floral_abundance + trt, method = "bray", perm = h.out, by = "margin", data = comm.data2)
pol.adon.rmout # trt and floral abundance significant

# Try adding floral richness
pol.adon.rmout2 <- adonis2(comm.data2[13:44] ~ comm.data2$floral_abundance + comm.data2$trt + comm.data2$floral_richness, method = "bray", perm = h.out, by = "margin")
pol.adon.rmout2 # trt and floral abundance still significant, floral richness not, drop?

# Try different permutation
pol.adon.rmout3 <- adonis2(comm.data2[13:44] ~ floral_abundance + trt, method = "bray", strata = comm.data2$pair_id, by = "margin", data = comm.data2, permutations = 4999)
pol.adon.rmout3 # trt and floral abundance significant

## Floral communities ----
# Subset community data from pol_floral
floral.comm.data <- floral_wide[rowSums(floral_wide[7:34])>0,]

# Run NMDS
floral.mds <- metaMDS(floral.comm.data[7:34]) # warning that stress is near zero
floral.mds  # Stress ~ 0
plot(floral.mds, type="t") # weird outlier says B. vulgaris

# Convert to data frame (only extracting site NMDSs)
floral.mds.scores <- as.data.frame(scores(floral.mds, display="site"))  

# Add plot info back in
floral.mds.scores$plot_id <- floral.comm.data$plot_id
floral.mds.scores$pair_id <- floral.comm.data$pair_id
floral.mds.scores$trt <- floral.comm.data$trt
floral.mds.scores$date <- floral.comm.data$date

# Plot NMDS
# Weird outlier...
ggplot(data=floral.mds.scores) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=trt),level = 0.50) +
  geom_point(aes(x=NMDS1,y=NMDS2,shape=trt,colour=trt),size=4, alpha = 0.5)

# Stress plot # Bad
stressplot(floral.mds, lwd = 0.5)

# Try to remove B. vulgaris outlier
floral.comm.data2 <- floral.comm.data[!floral.comm.data$Barbarea_vulgaris > 0,]
floral.comm.data2 <- floral.comm.data2 %>%
  select(-Barbarea_vulgaris)

# Run NMDS
floral.mds2 <- metaMDS(floral.comm.data2[7:33], k = 3, trymax = 1000) # no warning
floral.mds2  # Stress ~ 0.1
plot(floral.mds2, type="t") # No major outliers!

# Convert to data frame (only extracting site NMDSs)
floral.mds.scores2 <- as.data.frame(scores(floral.mds2, display="site"))  

# Add plot info back in
floral.mds.scores2$plot_id <- floral.comm.data2$plot_id
floral.mds.scores2$pair_id <- floral.comm.data2$pair_id
floral.mds.scores2$trt <- floral.comm.data2$trt
floral.mds.scores2$date <- floral.comm.data2$date

# Plot NMDS
ggplot(data=floral.mds.scores2) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=trt),level = 0.50) +
  geom_point(aes(x=NMDS1,y=NMDS2,shape=trt,colour=trt),size=4, alpha = 0.5)

ggplot(data=floral.mds.scores2) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=trt)) +
  geom_point(aes(x=NMDS1,y=NMDS2,color=trt), alpha = 0.5) +
  theme_bw() +
  labs(color = "treatment") +
  scale_color_manual(values = c("#e41a1c", "#377eb8"))

# Plot mean values by plot id
floral.mds.scores.mean2 <- floral.mds.scores2 %>%
  group_by(plot_id, trt) %>%
  summarize(mean.NMDS1 = mean(NMDS1),
            mean.NMDS2 = mean(NMDS2))

floral.nmds.mean <- ggplot(data=floral.mds.scores.mean2) + 
  stat_ellipse(aes(x=mean.NMDS1,y=mean.NMDS2,color=trt, linetype = trt)) +
  geom_point(aes(x=mean.NMDS1,y=mean.NMDS2,color=trt, shape = trt), size = 3, alpha = 0.5) +
  theme_bw() +
  labs(x = "NMDS1", y = "NMDS2", color = "Treatment", shape = "Treatment", linetype = "Treatment") +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  scale_linetype_manual(values = c(1, 2)) +
  coord_equal() +
  scale_y_continuous(breaks=c(-1, 0, 1)) +
  scale_x_continuous(breaks=c(-1, 0, 1))

# Stress plot # Still has weird uptick at the end but better
stressplot(floral.mds2, lwd = 0.5)

# Try to add random effects by using permutation structure
h.plant <- with(floral.comm.data2,
  how(plots = Plots(strata = plot_id),
      blocks = pair_id,
      nperm = 4999)
)

floral.adon6 <- adonis2(floral.comm.data2[7:33] ~ floral.comm.data2$trt, method = "bray", perm = h.plant)
floral.adon6

# Try different permutation structure
floral.adon7 <- adonis2(floral.comm.data2[7:33] ~ trt, method = "bray", data = floral.comm.data2, strata = floral.comm.data2$pair_id, permutations = 4999)
floral.adon7 # marginally significant, but now not accounting for repeated sampling

## Figure S4 ----
# Figure 5
jpeg("plot_fig_s4.jpeg", width = 7, height = 4, units = "in", res = 300)
pol.nmds.mean + floral.nmds.mean + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A')
dev.off()
