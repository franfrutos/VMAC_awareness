# Author: Francisco Garre-Frutos
# Date: 20/06/2024

# Analysis for all the analyses related to experiment 1

source("scripts/functions.R") # This will load functions and data into the workspace

# Loading packages----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(
  dplyr,# Data wrangling
  tidyr,# Data wrangling
  lme4,# LMMs
  lmerTest,# P-values in summary
  marginaleffects,# Model predictions and conditional effects
  hypr,# Set the contrast matrix
  sjPlot,# Tables
  ggplot2,# Plotting
  ggExtra,# fancier plotting
  ggpubr # Combined plots
)


# RT analysis ----
d_RT <-
  filter_data(
    raw_e1,
    f.absent = F,
    sd_filter = 2,
    experiment = "e1"
  )
d_RT$log_RT <- log(d_RT$rt)
d_RT$Singleton <-
  factor(d_RT$Singleton, levels = c("High", "Low", "Absent"))

HcRep <- hypr(
  VMAC = High ~ Low,
  AC = Low ~ Absent,
  levels = c("High", "Low", "Absent")
) # Contrast for VMAC and AC effects

contrasts(d_RT$Singleton) <- contr.hypothesis(HcRep)

fit_pow <-
  lmer(
    log_RT ~ Singleton * scale(log(Block)) + (Singleton * scale(log(Block)) |
                                                ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT
  )

summary(fit_pow) # Singular fit

fit.2_pow <-
  lmer(
    log_RT ~ Singleton * scale(log(Block)) + (Singleton + scale(log(Block)) |
                                                ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT
  ) # Maximal feasible model

summary(fit.2_pow) # model summary


# Selected model predictions for each type of singleton
predictions(
  fit.2_pow,
  newdata = datagrid(
    Singleton = unique,
    Block = seq(1, 12, 1),
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(fit.2_pow) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects comparisons:
avg_comparisons(
  fit.2_pow,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block = seq(1, 12, 1), ID = NA),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(fit.2_pow) ^ 2) /
                               2) - exp(lo + (sigma(fit.2_pow) ^ 2) / 2)
)

# ACC analysis ----

d_acc <-
  filter_data(raw_e1,
              f.absent = F,
              acc = F,
              experiment = "e1")
d_acc$Singleton <-
  factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))

contrasts(d_acc$Singleton) <- contr.hypothesis(HcRep)

fit.acc <-
  glmer(
    acc ~ Singleton * Block + (Singleton * scale(Block) | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    family = binomial()
  )

summary(fit.acc) # Singular fit

fit.acc2 <-
  glmer(
    acc ~ Singleton * Block + (Singleton + scale(Block) | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    family = binomial()
  )

summary(fit.acc2) # Singular fit, low SingletonAC random slope variance


fit.acc3 <-
  glmer(
    acc ~ Singleton * scale(Block) + (scale(Block) | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    family = binomial()
  )

summary(fit.acc3) # Maximal model

# Predictions for each Singleton:
predictions(
  fit.acc3,
  newdata = datagrid(
    Singleton = unique,
    Block = seq(1, 12, 1),
    ID = NA
  ),
  re.form = NA,
  by = "Singleton"
) # Predicted accuracy by Singleton averaged over block

# Contrasts for high-low and low-absent:
avg_comparisons(
  fit.acc3,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block = seq(1, 12, 1), ID = NA),
  re.form = NA
)


# Tables from both models ----
sjPlot::tab_model(fit.2_pow, digits = 3, digits.re = 3) # table 1 (left)
sjPlot::tab_model(fit.acc3, digits = 3, digits.re = 3) # table 1 (right)

# Plots ----
# Model predictions:
# Returns a list with raw averaged data and model preds
preds <- get_predictions(d_RT, fit.2_pow)

# Plot predictions
ppreds <- ggplot(data = preds[["mod"]],
                 aes(
                   y = estimate,
                   x = Block,
                   color = Singleton,
                   fill = Singleton
                 )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  geom_point(data = preds[["raw"]],
             aes(x = as.numeric(Block)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = preds[["raw"]],
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(600, 900, 50)) +
  scale_color_brewer(palette = "Set1", name = "Singleton type") +
  scale_fill_brewer(palette = "Set1", name = "Singleton type") +
  labs(y = "Response time (ms)") +
  guides(colour = guide_legend(position = "inside"),
         fill = guide_legend(position = "inside")) +
  theme_Publication(text_size = 12) +
  theme(
    legend.position = "top",
    legend.spacing.y = unit(.2, 'cm'),
    plot.margin = unit(c(2, 2, 2, 2), "mm")
  )

# Analysis of contingency rating:

raw_diff <- raw_e1[!raw_e1$ID %in% acc_ex$e1 &
                     raw_e1$awareness_estimate_High != "None", ] %>%
  mutate(
    points_h = as.numeric(awareness_estimate_High),
    points_l = as.numeric(awareness_estimate_Low)
  ) %>%
  dplyr::summarise(points = max(points_h) - max(points_l), .by = "ID")
effs <- d_RT %>%
  dplyr::summarise(RT = mean(rt), .by = c("ID", "Singleton")) %>%
  tidyr::spread(Singleton, RT) %>%
  mutate(VMAC = High - Low) %>%
  left_join(., raw_diff, by = "ID")


car::qqPlot(effs$points) # Non-normally distributed
car::qqPlot(effs$VMAC) # Almost normally distributed

wilcox.test(effs$points)# Significant differences

effs$sd_diff <- scale(effs$points)[, 1]
effs$sd_vmac <- scale(effs$VMAC)[, 1]

# Descriptive statistics
mean(effs$points) # 3,478.9 points
median(effs$points) # 0
sd(effs$points) # 11,861.12

# 48.78% of the participants estimated the contingency in the correct direction (i.e. High > Low)
with(effs, (length(points[points > 0]) /
              length(points)) * 100)
length(effs$points[effs$points > 0]) # 40 participants

# Correlation between VMAC and points difference
with(effs, cor.test(sd_diff, sd_vmac, method = "spearman"))

plot_points <- ggplot(effs, aes(sd_diff, sd_vmac)) +
  geom_point(alpha = .5,
             color = "darkorange",
             fill = "darkorange") +
  geom_smooth(method = "lm",
              color = "darkorange",
              fill = "darkorange") +
  labs(x = "Standarized contingency rating", y = "Standarized VMAC effect") +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  theme_Publication()

plot_points <- ggMarginal(
  plot_points,
  fill = "darkorange",
  color = "darkorange",
  alpha = .3,
  type = "densigram",
  size = 8
)

ggarrange(ppreds, plot_points, labels = c("a", "b"))
ggsave(
  "Output/plots/figure2.png",
  height = 15,
  width = 30,
  unit = "cm",
  dpi = 900
)

# Complementary analysis: Is contingency awareness associated with the temporal dynamics of VMAC? ----
# Merge data frames:
d_RTc <- d_RT %>% left_join(., raw_diff, by = "ID") %>% dplyr::rename("Contingency" = "points.y")

fit_c <-
  lmer(
    log_RT ~ Singleton * scale(log(Block)) * scale(Contingency) + (Singleton + scale(log(Block)) |
                                                                     ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RTc
  ) # Same model structure than the previous model

summary(fit_c) # Significant VMAC x Block x Contingency interaction

# Creating VMAC and AC conditional effects as a function of contingency
comps_ni <- comparisons(
  fit_c,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Block = seq(1, 12, .01),
    Contingency = function(x)
      quantile(x, probs = c(.30, .80)),
    # Selecting predictions for different quantiles
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(fit_c) ^ 2) / 2) - exp(lo +
                                                                  (sigma(fit_c) ^ 2) / 2)
) %>%
  mutate(Contingency = case_when(
    Contingency < quantile(raw_diff$points, probs = .6)[1] ~ "Low Contingency",
    T ~ "High Contingency",
  ))

comps_ni$Contingency <- factor(comps_ni$Contingency,
                               levels = c("Low Contingency", "High Contingency"))

comps_ni$Effect <- ifelse(comps_ni$contrast == "High, Low", "VMAC", "AC")

comps_ni$Effect <- factor(comps_ni$Effect, levels = c("VMAC", "AC"))

# Get raw data summary statistics splitting the sample
rawCont <- get_effectCont(d_RTc, probs = .6)

# Plot both model conditional effects and raw data together
ggplot(data = comps_ni, aes(
  y = estimate,
  x = Block,
  color = Effect,
  fill = Effect,
)) +
  geom_line() +
  facet_wrap( ~ Contingency) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  #facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(data = rawCont, aes(x = as.numeric(Block)), position = position_dodge(.5)) +
  geom_errorbar(
    data = rawCont,
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.6),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_Publication(text_size = 12) + # text size is adjusted for DPI
  theme(
    legend.spacing.x = unit(.2, 'cm'),
    plot.margin = margin(t = ".2", l = ".5", unit = "cm"),
    legend.position = "bottom"
  )



# Complementary analysis: Same as before but on accuracy
# Format data
d_accC <- d_acc %>% left_join(., raw_diff, by = "ID") %>% dplyr::rename("Contingency" = "points.y")

fit.acc_cont <-
  glmer(
    acc ~ Singleton * Block * scale(Contingency) + (Block | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_accC,
    family = binomial()
  )

summary(fit.acc_cont)

sjPlot::tab_model(fit_c, digits = 3, digits.re = 3) # table S2 (left)
sjPlot::tab_model(fit.acc_cont, digits = 3, digits.re = 3) # table S2 (right)
