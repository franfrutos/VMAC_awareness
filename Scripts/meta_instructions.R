# Author: Francisco Garre-Frutos

# Meta-analysis for studies with and without instructions

source("scripts/functions.R") # This will load functions and data into the work space

# Loading packages----

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(
  dplyr,# Data wrangling
  tidyr,# Data wrangling
  ggplot2,# Plotting
  stringr,
  metafor,
  effectsize,
  pwr,
  brms
)

# Formate data ----
d_meta <- read_excel(path="Input/meta_analysis/VMAC_meta_search.xlsx", sheet = 2) %>%
  filter(included == 1)

tmp <- read_excel("Input/data/Raw_data_Exp1_Laurent_.xlsx", sheet = 2)
d_meta[which(d_meta$Study == "Gregoire et al. (2022) - Experiment 1"),"statistic"] <- as.character(t.test(tmp$Reward, tmp$Neutral, paired = T)$statistic)

tmp <- read_excel("Input/data/Raw_data_Exp2_Laurent_.xlsx", sheet = 2)
d_meta[which(d_meta$Study == "Gregoire et al. (2022) - Experiment 2"),"statistic"] <- as.character(t.test(tmp$Reward, tmp$Neutral, paired = T)$statistic)

nrow(d_meta) # Number of valid contrasts

# Compute effect sizes:
d_meta$statistic <- as.numeric(d_meta$statistic)
d_meta$t_statistic <- ifelse(d_meta$test == "F", sqrt(d_meta$statistic), d_meta$statistic)
d_meta$dz <- d_meta$t_statistic/sqrt(d_meta$N)

# Hedges correction
J <- 1 - (3 / (4*(d_meta$N - 1) - 1))
d_meta$g <- d_meta$g_yi * J
d_meta$vg <- d_meta$g_vi * J^2

# Format moderators:
d_meta$Measure <- factor(d_meta$Measure)
contrasts(d_meta$Measure) <- contr.sum(2)/2
d_meta$Instructions <- factor(d_meta$Instructions)
contrasts(d_meta$Instructions) <- contr.sum(2)/2
d_meta$Training_length_c <- log(d_meta$Training_length) - log(288)


# Meta-analysis ----
# Overall VMAC effect:
rma(g_yi, g_vi, data = d_meta, method = "REML")

# Moderator analysis:
(fit_meta <- rma(g_yi, g_vi,  mods = ~ Measure + Instructions + Training_length_c,
    data = d_meta, method = "REML"))

# Only for studies with no instructions: 
(rma(g_yi, g_vi, mods = ~ Training_length_c + Measure,
    data = d_meta %>% filter(Instructions == 0) %>%
      mutate(Measure = ifelse(Measure == 0, .5, -.5)), method = "REML"))

# Power to find this effect? 
pwr.t.test(d = 0.1740, power = .80, type = "paired")

# Power to detect the same effect on instructed studies?
(rma(g_yi, g_vi, mods = ~ Training_length_c + Measure,
     data = d_meta %>% filter(Instructions == 1) %>%
       mutate(Measure = ifelse(Measure == 0, .5, -.5)), method = "REML"))

# Power to find this effect? 
pwr.t.test(d = 0.5729, power = .80, type = "paired")

# Bayesian analysis:
fit_bf <- 
  brm(data = d_meta %>% filter(Instructions == 0) %>%
        mutate(Measure = ifelse(Measure == 0, .5, -.5)), family = gaussian,
      g_yi | se(g_vi) ~ Measure + Training_length_c + (1 | Study),
      iter = 20000, warmup = 2000, cores = 4, chains = 4,
      save_pars = save_pars(all = T),
      seed = 14)

fit_bf2 <- 
  brm(data = d_meta %>% filter(Instructions == 0) %>%
        mutate(Measure = ifelse(Measure == 0, .5, -.5)), family = gaussian,
      g_yi | se(g_vi) ~ 0 + Measure + Training_length_c + (1 | Study),
      iter = 20000, warmup = 2000, cores = 4, chains = 4,
      save_pars = save_pars(all = T),
      seed = 14)

# Bayes Factor for the null
1/bayes_factor(fit_bf, fit_bf2)$bf

# Moderator plot ----
Training_length_c <- seq(log(30) - log(288), max(d_meta$Training_length_c), by = .005)

Instructions <- c(0.5, -.5)
Measure <- c(0.5, -.5)

# Create newdata
newdata <- expand.grid(
  Instructions = Instructions,
  Measure = Measure,
  Training_length_c = Training_length_c
)

pred <- predict(fit_meta, newmods = model.matrix(~ Measure + Instructions + Training_length_c, newdata)[, -1])

# Back-transform predictions
newdata$pred_z <- pred$pred
newdata$se_z <- pred$se
newdata$ci.lb_z <- pred$ci.lb
newdata$ci.ub_z <- pred$ci.ub
newdata$pi.lb_z <- pred$pi.lb
newdata$pi.ub_z <- pred$pi.ub

newdata <- newdata %>%
  mutate(Measure = ifelse(Measure == .5, "RTs", "Eye-tracking"),
         Instructions = ifelse(Instructions == .5, "No instructions", "Instructions"),
         Training_length = exp(Training_length_c + log(288)))

d_meta2 <- d_meta %>%
  mutate(Measure = ifelse(Measure == 0, "RTs", "Eye-tracking"),
         Instructions = ifelse(Instructions == 0, "No instructions", "Instructions"))

# Summary taking into account moderators
avg_data <- cbind(Rmisc::summarySE(data=d_meta2,
                 measurevar = "g_yi",
                 groupvars = c("Measure","Instructions")) %>%
  dplyr::rename("Measure1" = "g_yi",
                "se1" = "se") %>%
  select(Measure, Instructions, Measure1, se1),
  Rmisc::summarySE(data=d_meta2,
                 measurevar = "Training_length",
                 groupvars = c("Measure", "Instructions"))%>%
    dplyr::rename("Measure2" = "Training_length",
                  "se2" = "se") %>%
    select(Measure2, se2))

# Summary across all studies
avg_data2 <- cbind(Rmisc::summarySE(data=d_meta2,
                                   measurevar = "g_yi") %>%
                    dplyr::rename("Measure1" = "g_yi",
                                  "se1" = "se") %>%
                    select(,Measure1, se1),
                  Rmisc::summarySE(data=d_meta2,
                                   measurevar = "Training_length")%>%
                    dplyr::rename("Measure2" = "Training_length",
                                  "se2" = "se") %>%
                    select(Measure2, se2))

# Differences in training length?
t.test(Training_length ~ Instructions, data = d_meta2, var.equal = F)

# Plot moderators
ggplot(newdata, aes(Training_length, pred_z, color = Instructions)) +
  geom_point(data=d_meta2, aes(Training_length, g_yi, shape = Measure, size = N, fill = Instructions), alpha = .2) +
  geom_point(data=avg_data, aes(Measure2, Measure1, shape = Measure, fill = Instructions), size = 3) +
  geom_point(data=avg_data2, aes(Measure2, Measure1), size = 3, color = "black") + 
  geom_line(aes(linetype = Measure), size = 1) +
  geom_errorbar(data=avg_data, aes(ymin = Measure1 - se1, ymax = Measure1 + se1, y = Measure1, x = Measure2),
                width = 0) +
  geom_errorbar(data=avg_data2, aes(ymin = Measure1 - se1, ymax = Measure1 + se1, y = Measure1, x = Measure2),
                width = 0, color = "black") +
  geom_errorbarh(data=avg_data, aes(xmin = Measure2 - se2, xmax = Measure2 + se2, y = Measure1, x = Measure2),
                 height = 0) +
  geom_errorbarh(data=avg_data2, aes(xmin = Measure2 - se2, xmax = Measure2 + se2, y = Measure1, x = Measure2),
                 height = 0, color = "black") +
  coord_cartesian(ylim = c(-.25, 1.5), xlim = c(50, 2000)) +
  scale_y_continuous(breaks = seq(-.5, 1.500, .25)) +
  scale_x_continuous(breaks = seq(0, 2000, 250)) +
  scale_size(range = c(2,15), limits = c(min(d_meta$N), max(d_meta$N)),
             guide = "none") +
  scale_fill_manual(values=c("#252373", "#d16606")) +
  scale_color_manual(values=c("#252373", "#d16606")) +
  geom_ribbon(aes(ymin = ci.lb_z, ymax = ci.ub_z, fill = Instructions,
                  group = interaction(Measure, Instructions)),
              alpha = .15, color = NA) +
  scale_shape_manual(values = c(21, 24, 25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Effect size (Hedge's g)",
       x = "Training length") +
  theme_Publication(base_size = 12, text_size = 12) +
  theme(legend.key.size = unit(.5, 'cm'))

ggsave(filename = "Output/plots/figure5.png",
       height = 17,
       width = 20,
       units = "cm",
       dpi = 900)
