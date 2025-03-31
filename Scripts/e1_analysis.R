# Author: Francisco Garre-Frutos

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
  ggpubr, # Combined plots
  afex,
  BayesFactor, # BF
  effectsize, # cohen's d
  Rmisc, # Within-subject summary statistics 
  car # Anova tables for (G)LMMs
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

Anova(fit.2_pow) # ANOVA-like table of the selected model

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

# Bayes factor for the VMAC effect
dBF <- d_RT %>%
  filter(Singleton != "Absent") %>%
  dplyr::summarise(rt = mean(rt),
                   .by = c(ID, Singleton)) %>%
  spread(Singleton, rt)

t.test(dBF$High, dBF$Low, paired = T)
t_to_d(0.39793, 81, paired = T)

1/BayesFactor::ttestBF(dBF$High, dBF$Low, paired = T)

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
Anova(fit.acc3) # Anova-like table

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
(ppreds <- ggplot(data = preds[["mod"]],
                 aes(
                   y = estimate,
                   x = Block,
                   color = Singleton,
                   fill = Singleton
                 )) +
  geom_line(size = 1, aes(linetype = Singleton)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .15,
              color = NA) +
  geom_point(data = preds[["raw"]],
             aes(x = as.numeric(Block), shape = Singleton),
             position = position_dodge(.5),
             size = 2.5) +
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
  scale_y_continuous(breaks = seq(600, 900, 25)) +
  scale_color_manual(values = list_colors$singleton, name = "Singleton type") +
  scale_fill_manual(values = list_colors$singleton, name = "Singleton type") +
  scale_shape(name = "Singleton type") +
  scale_linetype(name = "Singleton type") +
  labs(y = "Response time (ms)") +
  theme_Publication(text_size = 14) +
  theme(
    legend.position = "top",
    legend.spacing.y = unit(.2, 'cm'),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    legend.key.size = unit(1, 'cm')
  ))

# Analysis of contingency rating:
raw_diff <- raw_e1[!raw_e1$ID %in% acc_ex$e1 &
                     raw_e1$awareness_estimate_High != "None", ] %>%
  mutate(
    points_h = as.numeric(awareness_estimate_High),
    points_l = as.numeric(awareness_estimate_Low)
  ) %>%
  dplyr::summarise(points = max(points_h) - max(points_l), .by = "ID") %>%
  mutate(points_norm = points/max(points))
effs <- d_RT %>%
  dplyr::summarise(RT = mean(rt), .by = c("ID", "Singleton")) %>%
  tidyr::spread(Singleton, RT) %>%
  mutate(VMAC = High - Low) %>%
  left_join(., raw_diff, by = "ID")


car::qqPlot(effs$points) # Non-normally distributed
car::qqPlot(effs$VMAC) # Almost normally distributed

wilcox.test(effs$points) # Significant differences

# Descriptive statistics
mean(effs$points) # 3,478.9 points
median(effs$points) # 0
sd(effs$points) # 11,861.12

# 48.78% of the participants estimated the contingency in the correct direction (i.e. High > Low)
with(effs, (length(points[points > 0]) /
              length(points)) * 100)
length(effs$points[effs$points > 0]) # 40 participants

# Correlation between VMAC and points difference
with(effs, cor.test(points, VMAC, method = "spearman"))
sp_correction(.12, .75, .7)

avg_data <- cbind(Rmisc::summarySE(data=effs,
                                    measurevar = "VMAC") %>%
                     dplyr::rename("Measure1" = "VMAC",
                                   "se1" = "se") %>%
                     select(Measure1, se1),
                   Rmisc::summarySE(data=effs,
                                    measurevar = "points_norm")%>%
                     dplyr::rename("Measure2" = "points_norm",
                                   "se2" = "se") %>%
                     select(Measure2, se2))

(plot_points <- ggplot(effs, aes(points_norm, VMAC)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(alpha = .4,
             color = "#d16606",
             fill = "#d16606") +
  geom_smooth(method = "lm",
              color = "#d16606",
              fill = "#d16606",
              size = .9, alpha = .15,
              linetype = "dashed") +
  labs(x = "Normalized contingency rating", y = "VMAC scores (ms)") +
  geom_point(data=avg_data, aes(Measure2, Measure1), size = 3, color = "#D2691E") +  
  geom_errorbar(data=avg_data, aes(ymin = Measure1 - se1, ymax = Measure1 + se1, y = Measure1, x = Measure2),
                width = 0, color = "#D2691E", size = 1) +
  geom_errorbarh(data=avg_data, aes(xmin = Measure2 - se2, xmax = Measure2 + se2, y = Measure1, x = Measure2),
                 height = 0, color = "#D2691E", size = 1) +
  scale_y_continuous(breaks = seq(-40, 80, 40)) +
  theme_Publication(text_size = 14))

plot_points <- ggMarginal(
  plot_points,
  fill = "#d16606",
  color = "#d16606",
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
(pcontingency1 <- ggplot(data = comps_ni,
                  aes(
                    y = estimate,
                    x = Block,
                    color = Effect,
                    fill = Effect,
                    shape = Effect
                  )) +
    geom_line(size = 1, aes(linetype = Effect)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = .15,
                color = NA) +
    geom_point(data = rawCont,
               aes(x = as.numeric(Block), shape = Effect),
               position = position_dodge(.5),
               size = 2.5) +
    geom_errorbar(
      data = rawCont,
      aes(
        x = as.numeric(Block),
        ymin = estimate - se,
        ymax = estimate + se
      ),
      position = position_dodge(.5),
      width = .01
    ) +
    facet_wrap( ~ Contingency) +
    scale_x_continuous(breaks = seq(1, 24, 1)) +
    scale_y_continuous(breaks = seq(-20, 60, 20)) +
    scale_color_manual(values = list_colors$Effects, name = "Contrast") +
    scale_fill_manual(values = list_colors$Effects, name = "Contrast") +
    scale_shape(name = "Contrast") +
    scale_linetype(name = "Contrast") +
    labs(y = "Contrast (ms)") +
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme_Publication(text_size = 14) + # text size is adjusted for DPI
    theme(legend.position = "none",
          plot.margin = unit(c(8, 2, .5, 2), "mm"),
          strip.text.x = element_text(size = 14)
    ))


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
