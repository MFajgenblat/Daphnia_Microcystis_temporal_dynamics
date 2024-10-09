################################################################################
# DAPHNIA - MICROCYSTIS TEMPORAL DYNAMICS ANALYSIS
################################################################################

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(patchwork)
library(ggtext)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(lme4)

#-------------------------------------------------------------------------------
# Reading and preprocessing data
#-------------------------------------------------------------------------------

summary_data <- read.csv("summary_data.csv", sep=";") %>%
  drop_na() %>%
  rowid_to_column(var = "jar_ID") %>%
  mutate(Daph_season = factor(Daph_season, levels = c("early", "late")),
         MC_season = factor(MC_season, levels = c("early", "late")),
         Synchronicity = factor(c("Allotemporal", "Contemporal")[(Daph_season == MC_season) + 1], levels = c("Allotemporal", "Contemporal")),
         Daph_year = factor(Daph_year, levels = c(2018, 2019)),
         Daph_clone = factor(Daph_clone),
         MC_strain = factor(MC_strain))

mean(summary_data$D12_alive/summary_data$N_individuals)
sd(summary_data$D12_alive/summary_data$N_individuals)/sqrt(sum(summary_data$N_individuals))*100

#-------------------------------------------------------------------------------
# Exploratory visualization
#-------------------------------------------------------------------------------

summary_data %>%
  select(jar_ID, MC_season, Daph_season, Daph_year, N_individuals, D0_alive, D5_alive, D9_alive, D12_alive) %>%
  mutate(Daph_year = factor(Daph_year),
         Daph_season = case_when(Daph_season == "early" ~ "Early <i>Daphnia</i>",
                                 Daph_season == "late" ~ "Late <i>Daphnia</i>"),
         MC_season = case_when(MC_season == "early" ~ "Early <i>Microcystis</i>",
                               MC_season == "late" ~ "Late <i>Microcystis</i>")) %>%
  melt(id.vars = c("jar_ID", "MC_season", "Daph_season", "Daph_year", "N_individuals")) %>%
  mutate(day = as.numeric(gsub("_alive", "", gsub("D", "", variable)))) %>%
  ggplot(aes(x = day, y = value/N_individuals, group = jar_ID, color = Daph_year)) +
  geom_line(linewidth = 0.2, position = position_dodge(width = 0.4)) +
  geom_point(size = 0.05, position = position_dodge(width = 0.4)) +
  scale_x_continuous("Day", breaks = seq(0, 12, by = 2)) +
  scale_y_continuous("Fraction of individuals surviving", breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  scale_fill_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  coord_cartesian(clip = "on", xlim = c(0,12), ylim = c(0,1)) +
  facet_grid(MC_season ~ Daph_season) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_markdown(face = "bold", size = 8),
        strip.text.y = element_markdown(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave("Exploratory_visualization.png", width = 16, height = 10, units = "cm", dpi = 600)

summary(summary_data$D12_alive / summary_data$N_individuals)

#-------------------------------------------------------------------------------
# Binomial GLMM - Data preparation
#-------------------------------------------------------------------------------

datalist_GpG <- list(N = nrow(summary_data),
                     N_clones = length(levels(summary_data$Daph_clone)),
                     N_strains = length(levels(summary_data$MC_strain)),
                     N_clonestraincombinations = length(unique(paste0(summary_data$Daph_clone, "_", summary_data$MC_strain))),
                     Survived_individuals = summary_data$D12_alive,
                     Total_individuals = summary_data$N_individuals,
                     Daph_season = as.numeric(summary_data$Daph_season) - 1,
                     MC_season = as.numeric(summary_data$MC_season) - 1,
                     Synchronicity = as.numeric(summary_data$Daph_season == summary_data$MC_season),
                     Daph_year = as.numeric(summary_data$Daph_year),
                     Daph_clone = as.numeric(summary_data$Daph_clone),
                     MC_strain = as.numeric(summary_data$MC_strain),
                     Clonestraincombination = as.numeric(factor(paste0(summary_data$Daph_clone, "_", summary_data$MC_strain))),
                     GxG = 0,
                     beta_prior = c(0,3),
                     re_sd_prior = c(0,3))

#-------------------------------------------------------------------------------
# Binomial GLMM - Model fitting
#-------------------------------------------------------------------------------

# Reading and compiling Stan model
model <- cmdstan_model("Binomial_GLMM.stan", cpp_options = list(stan_threads = F))

# Original model
stanfit <- model$sample(data = datalist_GpG, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_12_GpG <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_12_GpG, "fit_12_GpG.rds")

# Model with Genotype x Genotype interaction
datalist_12_GxG <- datalist_GpG
datalist_12_GxG$GxG <- 1
stanfit <- model$sample(data = datalist_12_GxG, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_12_GxG <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_12_GxG, "fit_12_GxG.rds")

# Model at day 9
datalist_9_GpG <- datalist_GpG
datalist_9_GpG$Survived_individuals <- summary_data$D9_alive
stanfit <- model$sample(data = datalist_9_GpG, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_9_GpG <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_9_GpG, "fit_9_GpG.rds")

# Model at day 5
datalist_5_GpG <- datalist_GpG
datalist_5_GpG$Survived_individuals <- summary_data$D5_alive
stanfit <- model$sample(data = datalist_5_GpG, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_5_GpG <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_5_GpG, "fit_5_GpG.rds")

# Original model with tighter priors
datalist_GpG_tighterpriors <- datalist_GpG
datalist_GpG_tighterpriors$beta_prior <- c(0, 1)
datalist_GpG_tighterpriors$re_sd_prior <- c(0, 1)
stanfit <- model$sample(data = datalist_GpG_tighterpriors, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_12_GpG_tighterpriors <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_12_GpG_tighterpriors, "fit_12_GpG_tighterpriors.rds")

# Original model with wider priors
datalist_GpG_widerpriors <- datalist_GpG
datalist_GpG_widerpriors$beta_prior <- c(0, 9)
datalist_GpG_widerpriors$re_sd_prior <- c(0, 9)
stanfit <- model$sample(data = datalist_GpG_widerpriors, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_12_GpG_widerpriors <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_12_GpG_widerpriors, "fit_12_GpG_widerpriors.rds")

#-------------------------------------------------------------------------------
# Binomial GLMM - Reading the saved posterior files
#-------------------------------------------------------------------------------

fit_12_GpG <- readRDS("fit_12_GpG.rds")
fit_12_GxG <- readRDS("fit_12_GxG.rds")
fit_9_GpG <- readRDS("fit_9_GpG.rds")
fit_5_GpG <- readRDS("fit_5_GpG.rds")
fit_12_GpG_tighterpriors <- readRDS("fit_12_GpG_tighterpriors.rds")
fit_12_GpG_widepriors <- readRDS("fit_12_GpG_widerpriors.rds")

#-------------------------------------------------------------------------------
# Binomial GLMM - Metadata and helper functions
#-------------------------------------------------------------------------------

covariate_metadata <- data.frame(name = c("beta_0", "beta_1", "beta_2", "beta_3"),
                                 Covariate = factor(c("Intercept", "<i>Daphnia</i> timing", "<i>Microcystis</i> timing", "Contemporality"),
                                                    levels = c("Intercept", "<i>Daphnia</i> timing", "<i>Microcystis</i> timing", "Contemporality")))

expit <- function(x) {exp(x)/(1+exp(x))}

#-------------------------------------------------------------------------------
# Binomial GLMM - Convergence checks
#-------------------------------------------------------------------------------

summary(apply(fit_12_GpG, 2, rhat))
summary(apply(fit_12_GpG, 2, ess_bulk))
summary(apply(fit_12_GpG, 2, ess_tail))

fit_12_GpG[,colnames(fit_12_GpG)[!(substr(colnames(fit_12_GpG), 1, 3) %in% c("mu[", "Sur", "ave") | substr(colnames(fit_12_GpG), 1, 8) %in% c("clone_ef", "strain_e"))]] %>%
  melt() %>%
  left_join(data.frame(draw = 1:8000, chain = rep(1:8, each = 1000), iteration = rep(1:1000, 8))) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = value, color = factor(chain)), linewidth = 0.05, alpha = 0.5) +
  scale_color_brewer("Chain", palette = "Set1", direction = -1,
                     guide = guide_legend(nrow = 1, override.aes = list(linewidth = 1, alpha = 1))) +
  scale_x_continuous("Iteration") +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ variable, scales = "free_y", strip.position = "left", ncol = 4) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside",
        axis.title.x = element_text(face = "bold", size = 8),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 8),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave("BinomialGLMM_Traceplots.png", width = 16, height = 20, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Posterior predictive checks
#-------------------------------------------------------------------------------

PPC1 <- fit_12_GpG %>%
  spread_draws(Survived_individuals_rep[jar_ID]) %>%
  left_join(summary_data) %>%
  group_by(Daph_season, .draw) %>%
  summarize(Mean = mean(Survived_individuals_rep),
            Median = median(Survived_individuals_rep),
            SD = sd(Survived_individuals_rep)) %>%
  pivot_longer(c(Mean, Median, SD)) %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = Daph_season), alpha = 0.5, position = "identity") +
  scale_color_manual("Daphnia time", values = c("#8cad11","#f55600"), labels = c("Early <i>Daphnia</i>", "Late <i>Daphnia</i>")) +
  scale_fill_manual("Daphnia time", values = c("#8cad11","#f55600"), labels = c("Early <i>Daphnia</i>", "Late <i>Daphnia</i>")) +
  geom_vline(data = pivot_longer(summarize(group_by(summary_data, Daph_season),
                                           Mean = mean(D12_alive),
                                           Median = median(D12_alive),
                                           SD = sd(D12_alive)),
                                 c(Mean, Median, SD)),
             aes(xintercept = value, color = Daph_season)) +
  facet_wrap(~ name, scales = "free", switch = "x", ncol = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_markdown(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
PPC2 <- fit_12_GpG %>%
  spread_draws(Survived_individuals_rep[jar_ID]) %>%
  left_join(summary_data) %>%
  group_by(MC_season, .draw) %>%
  summarize(Mean = mean(Survived_individuals_rep),
            Median = median(Survived_individuals_rep),
            SD = sd(Survived_individuals_rep)) %>%
  pivot_longer(c(Mean, Median, SD)) %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = MC_season), alpha = 0.5, position = "identity") +
  geom_vline(data = pivot_longer(summarize(group_by(summary_data, MC_season),
                                           Mean = mean(D12_alive),
                                           Median = median(D12_alive),
                                           SD = sd(D12_alive)),
                                 c(Mean, Median, SD)),
             aes(xintercept = value, color = MC_season)) +
  scale_color_manual("Microcystis time", values = c("#8cad11","#f55600"), labels = c("Early <i>Microcystis</i>", "Late <i>Microcystis</i>")) +
  scale_fill_manual("Microcystis time", values = c("#8cad11","#f55600"), labels = c("Early <i>Microcystis</i>", "Late <i>Microcystis</i>")) +
  facet_wrap(~ name, scales = "free", switch = "x", ncol = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_markdown(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
PPC3 <- fit_12_GpG %>%
  spread_draws(Survived_individuals_rep[jar_ID]) %>%
  left_join(summary_data) %>%
  group_by(Synchronicity, .draw) %>%
  summarize(Mean = mean(Survived_individuals_rep),
            Median = median(Survived_individuals_rep),
            SD = sd(Survived_individuals_rep)) %>%
  pivot_longer(c(Mean, Median, SD)) %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = Synchronicity), alpha = 0.5, position = "identity") +
  geom_vline(data = pivot_longer(summarize(group_by(summary_data, Synchronicity),
                                           Mean = mean(D12_alive),
                                           Median = median(D12_alive),
                                           SD = sd(D12_alive)),
                                 c(Mean, Median, SD)),
             aes(xintercept = value, color = Synchronicity)) +
  scale_color_manual("Temporality", values = c("#8cad11","#f55600"), labels = c("Allotemporal", "Contemporal")) +
  scale_fill_manual("Temporality", values = c("#8cad11","#f55600"), labels = c("Allotemporal", "Contemporal")) +
  facet_wrap(~ name, scales = "free", switch = "x", ncol = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_markdown(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
PPC4 <- fit_12_GpG %>%
  spread_draws(Survived_individuals_rep[jar_ID]) %>%
  left_join(summary_data) %>%
  group_by(Daph_year, .draw) %>%
  summarize(Mean = mean(Survived_individuals_rep),
            Median = median(Survived_individuals_rep),
            SD = sd(Survived_individuals_rep)) %>%
  pivot_longer(c(Mean, Median, SD)) %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = Daph_year), alpha = 0.5, position = "identity") +
  geom_vline(data = pivot_longer(summarize(group_by(summary_data, Daph_year),
                                           Mean = mean(D12_alive),
                                           Median = median(D12_alive),
                                           SD = sd(D12_alive)),
                                 c(Mean, Median, SD)),
             aes(xintercept = value, color = Daph_year)) +
  scale_color_manual("Year", values = c("#8cad11","#f55600")) +
  scale_fill_manual("Year", values = c("#8cad11","#f55600")) +
  facet_wrap(~ name, scales = "free", switch = "x", ncol = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_markdown(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
PPC5 <- fit_12_GpG %>%
  spread_draws(Survived_individuals_rep[jar_ID]) %>%
  left_join(summary_data) %>%
  group_by(Daph_clone, .draw) %>%
  summarize(Mean = mean(Survived_individuals_rep),
            Median = median(Survived_individuals_rep),
            SD = sd(Survived_individuals_rep)) %>%
  pivot_longer(c(Mean, Median, SD)) %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = Daph_clone), alpha = 0.75, position = "identity") +
  geom_vline(data = pivot_longer(summarize(group_by(summary_data, Daph_clone),
                                           Mean = mean(D12_alive),
                                           Median = median(D12_alive),
                                           SD = sd(D12_alive)),
                                 c(Mean, Median, SD)),
             aes(xintercept = value, color = Daph_clone)) +
  scale_color_brewer("<i>Daphnia</i> clone", palette = "Set3", guide = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 12)) +
  scale_fill_brewer("<i>Daphnia</i> clone", palette = "Set3") +
  facet_wrap(~ name, scales = "free", switch = "x", ncol = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_markdown(face = "bold", size = 8),
        legend.text = element_markdown(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
PPC6 <- fit_12_GpG %>%
  spread_draws(Survived_individuals_rep[jar_ID]) %>%
  left_join(summary_data) %>%
  group_by(MC_strain, .draw) %>%
  summarize(Mean = mean(Survived_individuals_rep),
            Median = median(Survived_individuals_rep),
            SD = sd(Survived_individuals_rep)) %>%
  pivot_longer(c(Mean, Median, SD)) %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = MC_strain), alpha = 0.75, position = "identity") +
  geom_vline(data = pivot_longer(summarize(group_by(summary_data, MC_strain),
                                           Mean = mean(D12_alive),
                                           Median = median(D12_alive),
                                           SD = sd(D12_alive)),
                                 c(Mean, Median, SD)),
             aes(xintercept = value, color = MC_strain)) +
  scale_color_brewer("<i>Microcystis</i> strain", palette = "Set3", guide = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 12)) +
  scale_fill_brewer("<i>Microcystis</i> strain", palette = "Set3") +
  facet_wrap(~ name, scales = "free", switch = "x", ncol = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_markdown(face = "bold", size = 8),
        legend.text = element_markdown(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
PPC1 / PPC2 / PPC3 / PPC4 / PPC5 / PPC6
ggsave("BinomialGLMM_PPC.png", width = 16, height = 20, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Prior sensitivity analysis
#-------------------------------------------------------------------------------

rbind(data.frame(Prior = "Tighter (N(0,1))", spread_draws(fit_12_GpG_tighterpriors, beta_0[year], beta_1[year], beta_2[year], beta_3[year])),
      data.frame(Prior = "Regular (N(0,3))", spread_draws(fit_12_GpG, beta_0[year], beta_1[year], beta_2[year], beta_3[year])),
      data.frame(Prior = "Wider (N(0,9))", spread_draws(fit_12_GpG_widerpriors, beta_0[year], beta_1[year], beta_2[year], beta_3[year]))) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year, Prior)) %>%
  left_join(covariate_metadata, by = "name") %>%
  mutate(name = factor(paste0(Covariate, " - ", c(2018,2019)[year]),
                       levels = rev(c("Intercept - 2018", "Intercept - 2019",
                                      "<i>Daphnia</i> timing - 2018", "<i>Daphnia</i> timing - 2019",
                                      "<i>Microcystis</i> timing - 2018", "<i>Microcystis</i> timing - 2019",
                                      "Contemporality - 2018", "Contemporality - 2019"))),
         Prior = factor(Prior, levels = c("Tighter (N(0,1))", "Regular (N(0,3))", "Wider (N(0,9))"))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "grey50", size = 0.4) +
  stat_interval(aes(y = name, x = value, color = factor(Prior)), size = 2, .width = c(0.5, 0.8, 0.95, 0.99), alpha = 1/4, position = position_dodge(width = 0.5)) +
  scale_y_discrete("Variable") +
  scale_x_continuous("Effect on survival (logit scale)", expand = c(0,0)) +
  scale_color_brewer("Prior\nspecification", palette = "Dark2", guide = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_rect(fill = NA, color = NA),
        panel.grid.minor.x = element_line(color = "grey95"),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.4),
        axis.text.y = element_markdown(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank())
ggsave("BinomialGLMM_Prior_sensitivity_analysis.png", width = 16, height = 10, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Numerical results
#-------------------------------------------------------------------------------

# Inference on parameters average across both years
fit_12_GpG %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  mutate(year = c(2018,2019)[year]) %>%
  left_join(covariate_metadata, by = "name") %>%
  filter(Covariate != "Intercept") %>%
  group_by(Covariate) %>%
  summarize(Posterior_mean = mean(value),
            Posterior_median = median(value),
            Posterior_sd = sd(value),
            Posterior_CrI_0025 = quantile(value, 0.025),
            Posterior_CrI_0975 = quantile(value, 0.975),
            Posterior_probability_positive = mean(value > 0)*100) %>%
  arrange(Covariate)

# Reduction in survival percentage points average across both years on the natural scale
pp_reduction <- filter(spread_draws(fit_12_GpG, average_survival[year,treatment]), treatment %in% c(1,4))$average_survival - filter(spread_draws(fit_12_GpG, average_survival[year,treatment]), treatment %in% c(2,3))$average_survival
mean(pp_reduction)
quantile(pp_reduction, 0.025)
quantile(pp_reduction, 0.975)

# Inference on parameters per year
fit_12_GpG %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  mutate(year = c(2018,2019)[year]) %>%
  left_join(covariate_metadata, by = "name") %>%
  filter(Covariate != "Intercept") %>%
  group_by(year, Covariate) %>%
  summarize(Posterior_mean = mean(value),
            Posterior_median = median(value),
            Posterior_sd = sd(value),
            Posterior_CrI_0025 = quantile(value, 0.025),
            Posterior_CrI_0975 = quantile(value, 0.975),
            Posterior_probability_positive = mean(value > 0)*100) %>%
  arrange(Covariate)

# Inference on parameters being higher in 2019 compared to 2018
fit_12_GpG %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  mutate(year = c(2018,2019)[year]) %>%
  left_join(covariate_metadata, by = "name") %>%
  pivot_wider(names_from = year, values_from = value) %>%
  mutate(Difference = `2019` - `2018`) %>%
  group_by(Covariate) %>%
  summarize(Posterior_mean = mean(Difference),
            Posterior_median = median(Difference),
            Posterior_sd = sd(Difference),
            Posterior_CrI_0025 = quantile(Difference, 0.025),
            Posterior_CrI_0975 = quantile(Difference, 0.975),
            Posterior_probability_positive = mean(Difference > 0)*100) %>%
  arrange(Covariate)

# Comparison of worst and best performing Daphnia clone
mean((filter(spread_draws(fit_12_GpG, clone_effect[clone]), clone == 5)$clone_effect - filter(spread_draws(fit_12_GpG, clone_effect[clone]), clone == 4)$clone_effect) > 0)

# Comparison of most and least toxic Microcystis strain
mean((filter(spread_draws(fit_12_GpG, strain_effect[strain]), strain == 12)$strain_effect - filter(spread_draws(fit_12_GpG, strain_effect[strain]), strain == 10)$strain_effect) > 0)

#-------------------------------------------------------------------------------
# Binomial GLMM - Visualizing the fixed effects on the logit scale
#-------------------------------------------------------------------------------

fit_12_GpG %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  left_join(covariate_metadata, by = "name") %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50", size = 0.4) +
  stat_eye(aes(y = value, x = Covariate, color = factor(year), fill = factor(year)), position = position_dodge(width = 0.8), alpha = 0.4, point_alpha = 1, interval_alpha = 1, interval_size = 0.05, adjust = 1.5) +
  scale_x_discrete("Variable") +
  scale_y_continuous("Effect on survival (logit scale)", expand = c(0,0)) +
  scale_color_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  scale_fill_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  coord_cartesian(ylim = c(-3.9,2.6)) +
  theme(panel.background = element_rect(fill = NA, color = NA),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.major.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.4),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_markdown(size = 7),
        axis.title = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_blank())
ggsave("BinomialGLMM_Fixed_effects.png", width = 16, height = 7, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Visualizing the estimated survival per category
#-------------------------------------------------------------------------------

fit_12_GpG %>%
  spread_draws(average_survival[year,treatment]) %>%
  mutate(year = factor(2017 + year),
         combination = c("Early <i>Daphnia</i><br>Early <i>Microcystis</i>",
                         "Early <i>Daphnia</i><br>Late <i>Microcystis</i>",
                         "Late <i>Daphnia</i><br>Early <i>Microcystis</i>",
                         "Late <i>Daphnia</i><br>Late <i>Microcystis</i>")[treatment]) %>%
  group_by(year, combination) %>%
  summarise(Mean = mean(average_survival),
            LI = quantile(average_survival, 0.25),
            UI = quantile(average_survival, 0.75)) %>%
  ggplot(aes(x = combination, y = Mean, fill = year, color = year)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.65), width = 0.5) +
  geom_errorbar(aes(ymin = LI, ymax = UI), color = "black", width = 0, position = position_dodge(width = 0.65), size = 0.4) +
  scale_y_continuous("Average survival probability", expand = c(0,0), limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  scale_fill_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "grey93"),
        axis.text.x = element_markdown(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 8),
        axis.text.y = element_text(size = 7),
        axis.line.x = element_line(color = "black", size = 0.4),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave("BinomialGLMM_Estimated_survival.png", width = 16, height = 8, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Estimating the clone and strain random effects
#-------------------------------------------------------------------------------

clone_effects <- fit_12_GpG %>%
  spread_draws(clone_effect[clone]) %>%
  left_join(data.frame(clone = 1:11,
                       Clone = levels(summary_data$Daph_clone))) %>%
  ggplot() +
  stat_eye(aes(y = Clone, x = clone_effect), alpha = 0.4, point_alpha = 1, interval_alpha = 1, interval_size = 0.05, color = "#2b6a70", fill = "#2b6a70") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_color_brewer("Credible\ninterval") +
  scale_x_continuous("Effect on <i>Daphnia</i> survival (logit scale)", expand = c(0,0)) +
  scale_y_discrete("<i>Daphnia</i> clone", limits = rev(levels(summary_data$Daph_clone))) +
  coord_cartesian(xlim = c(-4,4)) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey95"),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(size = 8, face = "bold"),
        axis.title.x = element_markdown(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(colour = "black", size = 0.4),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.justification = "top",
        legend.key.width = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(size = 9, face = "bold", angle = 0),
        strip.placement = "outside")
strain_effects <- fit_12_GpG %>%
  spread_draws(strain_effect[strain]) %>%
  left_join(data.frame(strain = 1:12,
                       Strain = levels(summary_data$MC_strain))) %>%
  ggplot() +
  stat_eye(aes(y = Strain, x = strain_effect), alpha = 0.4, point_alpha = 1, interval_alpha = 1, interval_size = 0.05, color = "#99435c", fill = "#99435c") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_color_brewer("Credible\ninterval") +
  scale_x_continuous("Effect on <i>Daphnia</i> survival (logit scale)", expand = c(0,0)) +
  scale_y_discrete("<i>Microcystis</i> strain", limits = rev(levels(summary_data$MC_strain))) +
  coord_cartesian(xlim = c(-4,4)) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey95"),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(size = 8, face = "bold"),
        axis.title.x = element_markdown(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(colour = "black", size = 0.4),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.justification = "top",
        legend.key.width = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(size = 9, face = "bold", angle = 0),
        strip.placement = "outside")
clone_effects + strain_effects
ggsave("BinomialGLMM_Random_effects.png", width = 16, height = 8, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Visualizing the Genotype plus Genotype effects using a heatmap
#-------------------------------------------------------------------------------

fit_12_GpG %>%
  spread_draws(mu[jar_ID]) %>%
  group_by(jar_ID) %>%
  summarise(mu = mean(expit(mu))) %>%
  right_join(summary_data) %>%
  ggplot() +
  geom_tile(aes(x = Daph_clone, y = MC_strain, fill = mu)) +
  geom_vline(xintercept = 3.5, color = "white", size = 1) +
  geom_hline(yintercept = 3.5, color = "white", size = 1) +
  geom_text(aes(x = Daph_clone, y = MC_strain, label = substr(Synchronicity, 1, 1)), color = "white", size = 2) +
  scale_x_discrete("<i>Daphnia</i> clone ID", expand = c(0,0)) +
  scale_y_discrete("<i>Microcystis</i> strain ID", expand = c(0,0)) +
  scale_fill_viridis_c("Survival\nprobability\n", option = "B", limits = c(0,1), breaks = seq(0, 1, by = 0.2),
                       guide = guide_colorbar(barheight = unit(4.5, "cm"))) +
  facet_wrap(~ Daph_year, scales = "free") +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 8),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        axis.text = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7))
ggsave("BinomialGLMM_GpG_heatmap.png", width = 16, height = 8, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Comparison against saturated model
#-------------------------------------------------------------------------------

rbind(data.frame(Type = "Additive", spread_draws(fit_12_GpG, mu[jar_ID])),
      data.frame(Type = "Interaction", spread_draws(fit_12_GxG, mu[jar_ID]))) %>%
  pivot_wider(names_from = Type, values_from = mu) %>%
  group_by(jar_ID) %>%
  summarise(Additive_mean = mean(Additive),
            Additive_L = quantile(Additive, 0.025),
            Additive_U = quantile(Additive, 0.975),
            Interaction_mean = mean(Interaction),
            Interaction_L = quantile(Interaction, 0.025),
            Interaction_U = quantile(Interaction, 0.975)) %>%
  ggplot() +
  geom_errorbar(aes(x = Additive_mean, ymin = Interaction_L, ymax = Interaction_U), linewidth = 0.2) +
  geom_errorbarh(aes(xmin = Additive_L, xmax = Additive_U, y = Interaction_mean), linewidth = 0.2) +
  geom_point(aes(x = Additive_mean, y = Interaction_mean), size = 1) +
  scale_x_continuous("Linear predictor in the main, additive model") +
  scale_y_continuous("Linear predictor in the saturated model") +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7))
ggsave("BinomialGLMM_GxG.png", width = 12, height = 16, units = "cm", dpi = 600)

posterior_correlation <- sapply(seq(1, 8000, by = 10), function(i) cor(filter(spread_draws(fit_12_GpG, mu[jar_ID]), .draw == i)$mu, filter(spread_draws(fit_12_GxG, mu[jar_ID]), .draw == i)$mu))
hist(posterior_correlation)
1-mean(posterior_correlation)
1-median(posterior_correlation)
1-quantile(posterior_correlation, 0.025)
1-quantile(posterior_correlation, 0.975)

#-------------------------------------------------------------------------------
# Binomial GLMM - Frequentist analysis
#-------------------------------------------------------------------------------

summary_data$D12_dead <- summary_data$N_individuals - summary_data$D12_alive
summary_data$MC_strain <- factor(summary_data$MC_strain)

fit <- glmer(cbind(D12_alive, D12_dead) ~ Daph_season + MC_season + Synchronicity + (1 | Daph_clone) + (1 | MC_strain), family = binomial("logit"), data = subset(summary_data, Daph_year == 2018))
summary(fit)
fit <- glmer(cbind(D12_alive, D12_dead) ~ Daph_season + MC_season + Synchronicity + (1 | Daph_clone) + (1 | MC_strain), family = binomial("logit"), data = subset(summary_data, Daph_year == 2019))
summary(fit)

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Data preparation
#-------------------------------------------------------------------------------

individual_data <- summary_data[unlist(sapply(1:nrow(summary_data), function(i) rep(i, each = summary_data$N_individuals[i]))), c("jar_ID", "Daph_season", "MC_season", "Daph_year", "Daph_clone", "MC_strain")]
individual_data$Interval <- "[0,5)"
for (i in summary_data$jar_ID) {
  if (summary_data$D5_alive[summary_data$jar_ID == i] > 0) {
    individual_data$Interval[individual_data$jar_ID == i][1:summary_data$D5_alive[summary_data$jar_ID == i]] <- "[5,9)"
  }
  if (summary_data$D9_alive[summary_data$jar_ID == i] > 0) {
    individual_data$Interval[individual_data$jar_ID == i][1:summary_data$D9_alive[summary_data$jar_ID == i]] <- "[9,12)"
  }
  if (summary_data$D12_alive[summary_data$jar_ID == i] > 0) {
    individual_data$Interval[individual_data$jar_ID == i][1:summary_data$D12_alive[summary_data$jar_ID == i]] <- "[12,Inf)"
  }
}
individual_data$Interval <- factor(individual_data$Interval, levels = c("[0,5)", "[5,9)", "[9,12)", "[12,Inf)"))

datalist_survival <- list(N = nrow(individual_data),
                          N_clones = length(levels(individual_data$Daph_clone)),
                          N_strains = length(levels(individual_data$MC_strain)),
                          N_clonestraincombinations = length(unique(paste0(individual_data$Daph_clone, "_", individual_data$MC_strain))),
                          Interval = as.numeric(individual_data$Interval),
                          Daph_season = as.numeric(individual_data$Daph_season) - 1.5,
                          MC_season = as.numeric(individual_data$MC_season) - 1.5,
                          Synchronicity = as.numeric(individual_data$Daph_season == individual_data$MC_season),
                          Daph_year = as.numeric(individual_data$Daph_year),
                          Daph_clone = as.numeric(individual_data$Daph_clone),
                          MC_strain = as.numeric(individual_data$MC_strain),
                          Clonestraincombination = as.numeric(factor(paste0(individual_data$Daph_clone, "_", individual_data$MC_strain))),
                          GxG = 0,
                          N_preds = 200,
                          pred_grid = seq(0, 7.2, length.out = 200),
                          beta_prior = c(0,3),
                          re_sd_prior = c(0,3))

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Model fitting
#-------------------------------------------------------------------------------

model_survival <- cmdstan_model("Interval_censored_survival.stan", cpp_options = list(stan_threads = F))

# Main model
stanfit <- model_survival$sample(datalist_survival, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_survival <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")
fit_survival <- fit_survival$post_warmup_draws
saveRDS(fit_survival, "fit_survival.rds")

# Original model with tighter priors
datalist_survival_tighterpriors <- datalist_survival
datalist_survival_tighterpriors$beta_prior <- c(0, 1)
datalist_survival_tighterpriors$re_sd_prior <- c(0, 1)
stanfit <- model_survival$sample(data = datalist_survival_tighterpriors, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_survival_tighterpriors <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_survival_tighterpriors, "fit_survival_tighterpriors.rds")

# Original model with wider priors
datalist_survival_widerpriors <- datalist_survival
datalist_survival_widerpriors$beta_prior <- c(0, 9)
datalist_survival_widerpriors$re_sd_prior <- c(0, 9)
stanfit <- model_survival$sample(data = datalist_survival_widerpriors, chains = 8, iter_warmup = 1000, iter_sampling = 1000, parallel_chains = 8)
fit_survival_widerpriors <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit_survival_widerpriors, "fit_survival_widerpriors.rds")

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Reading the saved posterior files
#-------------------------------------------------------------------------------

fit_survival <- readRDS("fit_survival.rds")
fit_survival_tighterpriors <- readRDS("fit_survival_tighterpriors.rds")
fit_survival_widerpriors <- readRDS("fit_survival_widerpriors.rds")

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Convergence checks
#-------------------------------------------------------------------------------

summary(apply(fit_survival, 2, rhat))
summary(apply(fit_survival, 2, ess_bulk))
summary(apply(fit_survival, 2, ess_tail))


fit_survival[,colnames(fit_survival)[!(substr(colnames(fit_survival), 1, 3) %in% c("mu[", "Sur", "ave") | substr(colnames(fit_survival), 1, 8) %in% c("clone_ef", "strain_e", "S_averag"))]] %>%
  melt() %>%
  left_join(data.frame(draw = 1:8000, chain = rep(1:8, each = 1000), iteration = rep(1:1000, 8))) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = value, color = factor(chain)), linewidth = 0.05, alpha = 0.5) +
  scale_color_brewer("Chain", palette = "Set1", direction = -1,
                     guide = guide_legend(nrow = 1, override.aes = list(linewidth = 1, alpha = 1))) +
  scale_x_continuous("Iteration") +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ variable, scales = "free_y", strip.position = "left", ncol = 4) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside",
        axis.title.x = element_text(face = "bold", size = 8),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 8),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave("Survival_Traceplots.png", width = 16, height = 20, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Model criticism
#-------------------------------------------------------------------------------

fit_survival %>%
  spread_draws(Survival_times[rowid]) %>%
  left_join(rowid_to_column(individual_data)) %>%
  group_by(jar_ID, .draw) %>%
  summarize(Survival_times = mean(Survival_times)) %>%
  ggplot() +
  stat_interval(aes(y = Survival_times*10, x = NA), .width = c(0.5,0.8,0.95,0.99), size = 4) +
  geom_linerange(data = mutate(group_by(individual_data, jar_ID), indid = row_number()), aes(x = NA, ymin = c(0,5,9,12)[as.numeric(Interval)], ymax = c(5,9,12,Inf)[as.numeric(Interval)], group = indid),
                 position = position_dodge(width = 0.5), color = "#ba164a", alpha = 1, linewidth = 0.3) +
  scale_x_discrete("Experimental unit") +
  scale_y_continuous("Survival time (days)", expand = c(0,0)) +
  scale_color_brewer("Posterior\npredictive\ncredible\ninterval") +
  coord_cartesian(ylim = c(0,24)) +
  facet_wrap(~ jar_ID, strip.position = "bottom", ncol = 22, scales = "free_x") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "grey93"),
        panel.spacing.x = unit(0, "cm"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(face = "bold", size = 8),
        axis.line.x = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave("Survival_PPC.png", width = 16, height = 10, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Binomial GLMM - Prior sensitivity analysis
#-------------------------------------------------------------------------------

rbind(data.frame(Prior = "Tighter (N(0,1))", spread_draws(fit_survival_tighterpriors, beta_0[year], beta_1[year], beta_2[year], beta_3[year])),
      data.frame(Prior = "Regular (N(0,3))", spread_draws(fit_survival, beta_0[year], beta_1[year], beta_2[year], beta_3[year])),
      data.frame(Prior = "Wider (N(0,9))", spread_draws(fit_survival_widerpriors, beta_0[year], beta_1[year], beta_2[year], beta_3[year]))) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year, Prior)) %>%
  left_join(covariate_metadata, by = "name") %>%
  mutate(name = factor(paste0(Covariate, " - ", c(2018,2019)[year]),
                       levels = rev(c("Intercept - 2018", "Intercept - 2019",
                                      "<i>Daphnia</i> timing - 2018", "<i>Daphnia</i> timing - 2019",
                                      "<i>Microcystis</i> timing - 2018", "<i>Microcystis</i> timing - 2019",
                                      "Contemporality - 2018", "Contemporality - 2019"))),
         Prior = factor(Prior, levels = c("Tighter (N(0,1))", "Regular (N(0,3))", "Wider (N(0,9))"))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "grey50", size = 0.4) +
  stat_interval(aes(y = name, x = value, color = factor(Prior)), size = 2, .width = c(0.5, 0.8, 0.95, 0.99), alpha = 1/4, position = position_dodge(width = 0.5)) +
  scale_y_discrete("Variable") +
  scale_x_continuous("Effect on survival", expand = c(0,0)) +
  scale_color_brewer("Prior\nspecification", palette = "Dark2", guide = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_rect(fill = NA, color = NA),
        panel.grid.minor.x = element_line(color = "grey95"),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.4),
        axis.text.y = element_markdown(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank())
ggsave("Survival_Prior_sensitivity_analysis.png", width = 16, height = 10, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Numerical results
#-------------------------------------------------------------------------------

# Inference on parameters per year
fit_survival %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  mutate(year = c(2018,2019)[year]) %>%
  left_join(covariate_metadata, by = "name") %>%
  filter(Covariate != "Intercept") %>%
  group_by(year, Covariate) %>%
  summarize(Posterior_mean = mean(value),
            Posterior_median = median(value),
            Posterior_sd = sd(value),
            Posterior_CrI_0025 = quantile(value, 0.025),
            Posterior_CrI_0975 = quantile(value, 0.975),
            Posterior_probability_positive = mean(value > 0)*100) %>%
  arrange(Covariate)

# Inference on parameters being higher in 2019 compared to 2018
fit_survival %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  mutate(year = c(2018,2019)[year]) %>%
  left_join(covariate_metadata, by = "name") %>%
  pivot_wider(names_from = year, values_from = value) %>%
  mutate(Difference = `2019` - `2018`) %>%
  group_by(Covariate) %>%
  summarize(Posterior_mean = mean(Difference),
            Posterior_median = median(Difference),
            Posterior_sd = sd(Difference),
            Posterior_CrI_0025 = quantile(Difference, 0.025),
            Posterior_CrI_0975 = quantile(Difference, 0.975),
            Posterior_probability_positive = mean(Difference > 0)*100) %>%
  arrange(Covariate)

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Visualizing the fixed effects
#-------------------------------------------------------------------------------

fit_survival %>%
  spread_draws(beta_0[year], beta_1[year], beta_2[year], beta_3[year]) %>%
  pivot_longer(!c(.chain, .iteration, .draw, year)) %>%
  left_join(covariate_metadata, by = "name") %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50", size = 0.4) +
  stat_eye(aes(y = value, x = Covariate, color = factor(year), fill = factor(year)), position = position_dodge(width = 0.8), alpha = 0.4, point_alpha = 1, interval_alpha = 1, interval_size = 0.05, adjust = 1.5) +
  scale_x_discrete("Variable") +
  scale_y_continuous("Effect on survival", expand = c(0,0)) +
  scale_color_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  scale_fill_manual("Year", values = c("#8cad11","#f55600"), labels = c(2018, 2019)) +
  coord_cartesian(ylim = c(-1.8,1.3)) +
  theme(panel.background = element_rect(fill = NA, color = NA),
        panel.grid.minor.y = element_line(color = "grey95"),
        panel.grid.major.y = element_line(color = "grey95"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.4),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_markdown(size = 7),
        axis.title = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_blank())
ggsave("Survival_Fixed_effects.png", width = 16, height = 7, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Interval-censored survival analysis - Visualizing the estimated survival curves
#-------------------------------------------------------------------------------

a <- fit_survival %>%
  spread_draws(S_averaged[treatment, level, time]) %>%
  filter(treatment == 1) %>%
  mutate(days = 72*(time-1)/199) %>%
  ggplot() +
  stat_lineribbon(aes(x = days, y = S_averaged, color = factor(level), fill = factor(level)), alpha = 1/5, .width = c(0.5, 0.8, 0.95), size = 0.5) +
  stat_lineribbon(aes(x = days, y = S_averaged, color = factor(level), fill = factor(level)), .width = NA, size = 0.5) +
  scale_color_manual(values = c("#0090C8", "#CE178C"), labels = c("Late <i>Daphnia</i>", "Early <i>Daphnia</i>"),
                     guide = guide_legend(override.aes = list(linetype = "solid", alpha = 1, fill = NA))) +
  scale_fill_manual(values = c("#0090C8", "#CE178C"), labels = c("Late <i>Daphnia</i>", "Early <i>Daphnia</i>")) +
  scale_x_continuous("", breaks = seq(0, 72, by = 4), expand = c(0,0)) +
  scale_y_continuous("Survival probability", breaks = seq(0, 1, by = 0.2), expand = c(0,0), limits = c(0,1)) +
  coord_cartesian(xlim = c(0, 24)) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95", size = 0.3),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_markdown(size = 7),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_line(size = 0.4),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 7))
b <- fit_survival %>%
  spread_draws(S_averaged[treatment, level, time]) %>%
  filter(treatment == 2) %>%
  mutate(days = 72*(time-1)/199) %>%
  ggplot() +
  stat_lineribbon(aes(x = days, y = S_averaged, color = factor(level), fill = factor(level)), alpha = 1/5, .width = c(0.5, 0.8, 0.95), size = 0.5) +
  stat_lineribbon(aes(x = days, y = S_averaged, color = factor(level), fill = factor(level)), .width = NA, size = 0.5) +
  scale_color_manual(values = c("#0090C8", "#CE178C"), labels = c("Late <i>Microcystis</i>", "Early <i>Microcystis</i>"),
                     guide = guide_legend(override.aes = list(linetype = "solid", alpha = 1, fill = NA))) +
  scale_fill_manual(values = c("#0090C8", "#CE178C"), labels = c("Late <i>Microcystis</i>", "Early <i>Microcystis</i>")) +
  scale_x_continuous("Days", breaks = seq(0, 72, by = 4), expand = c(0,0)) +
  scale_y_continuous("", breaks = seq(0, 1, by = 0.2), expand = c(0,0), limits = c(0,1)) +
  coord_cartesian(xlim = c(0, 24)) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95", size = 0.3),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_markdown(size = 7),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_line(size = 0.4),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 7))
c <- fit_survival %>%
  spread_draws(S_averaged[treatment, level, time]) %>%
  filter(treatment == 3) %>%
  mutate(days = 72*(time-1)/199) %>%
  ggplot() +
  stat_lineribbon(aes(x = days, y = S_averaged, color = factor(level), fill = factor(level)), alpha = 1/5, .width = c(0.5, 0.8, 0.95), size = 0.5) +
  stat_lineribbon(aes(x = days, y = S_averaged, color = factor(level), fill = factor(level)), .width = NA, size = 0.5) +
  scale_color_manual(values = c("#0090C8", "#CE178C"), labels = c("Contemporal", "Allotemporal"),
                     guide = guide_legend(override.aes = list(linetype = "solid", alpha = 1, fill = NA))) +
  scale_fill_manual(values = c("#0090C8", "#CE178C"), labels = c("Contemporal", "Allotemporal")) +
  scale_x_continuous("", breaks = seq(0, 72, by = 4), expand = c(0,0)) +
  scale_y_continuous("", breaks = seq(0, 1, by = 0.2), expand = c(0,0), limits = c(0,1)) +
  coord_cartesian(xlim = c(0, 24)) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95", size = 0.3),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_markdown(size = 7),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_line(size = 0.4),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key = element_blank(),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 7))
a+b+c
ggsave("Survival_Curves.png", width = 16, height = 6, units = "cm", dpi = 600)
