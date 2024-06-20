# Author: Francisco Garre-Frutos
# Date: 20/06/2024

# Analysis for all the analyses related to experiment 2

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
# Rewarded phase:
d_RT <-
  filter_data(raw_e2, f.absent = F,
              sd_filter = 2,
              experiment = "e2")
d_RT$log_RT <- log(d_RT$rt)
d_RT$Singleton <-
  factor(d_RT$Singleton, levels = c("High", "Low", "Absent"))
d_RT$Group <- ifelse(d_RT$group =="A", .5, -.5) # Deviation coding fo Group predictor

HcRep <- hypr(
  VMAC = High ~ Low,
  AC = Low ~ Absent,
  levels = c("High", "Low", "Absent")
) # Contrast for VMAC and AC effects

contrasts(d_RT$Singleton) <- contr.hypothesis(HcRep)

fit_pow <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * Group + (Singleton * scale(log(Block_num)) |
                                                       ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT
  )

summary(fit_pow) # Singular fit

fit.2_pow <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * Group + (Singleton + scale(log(Block_num)) |
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
    Block_num = seq(1, 12, 1),
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(fit.2_pow) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects comparisons:
comparisons(
  fit.2_pow,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block_num = seq(1, 12, 1),
                     ID = NA,
                     Group = c(-.5, .5)),
  re.form = NA,
  by = c("Singleton", "Group"),
  comparison = \(hi, lo) exp(hi + (sigma(fit.2_pow) ^ 2) /
                               2) - exp(lo + (sigma(fit.2_pow) ^ 2) / 2)
)

# ACC analysis ----
d_acc <-
  filter_data(raw_e2, f.absent = F, acc = F, sd_filter = 2, experiment = "e2")
d_acc$Singleton <-
  factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))
d_acc$Group <- ifelse(d_acc$group =="A", .5, -.5) # Deviation coding for Group predictor

contrasts(d_acc$Singleton) <- contr.hypothesis(HcRep)

fit.acc <-
  glmer(
    acc ~ Singleton * scale(Block_num) * Group + (Singleton * scale(Block_num)|ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    family = binomial()
  ) 

summary(fit.acc)

fit.acc2 <-
  glmer(
    acc ~ Singleton * scale(Block_num) * Group + (Singleton + scale(Block_num)|ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    family = binomial()
  )

summary(fit.acc2) 

fit.acc3 <-
  glmer(
    acc ~ Singleton * scale(Block_num) * Group + (scale(Block_num)|ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    family = binomial()
  ) # Maximal model

summary(fit.acc3) # Summary 

# Predictions for each Singleton:
predictions(
  fit.acc3,
  newdata = datagrid(
    Singleton = unique,
    Block_num = seq(1, 12, 1),
    ID = NA,
    Group = c(.5, -.5)
  ),
  re.form = NA,
  by = c("Group", "Singleton")
) # Predicted accuracy by Group averaged over block

# Contrasts for high-low and low-absent:

avg_comparisons(
  fit.acc3,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block_num = seq(1, 12, 1),
                     ID = NA,
                     Group = c(.5, -.5)),
  re.form = NA,
  by = c("Group")
)  

# Tables for RT and ACC analysis ----
sjPlot::tab_model(fit.2_pow, digits = 3, digits.re = 3) # table 2 (left)
sjPlot::tab_model(fit.acc3, digits = 3, digits.re = 3) # table 2 (right)


# Plots RT----
# Model predictions:
# Returns a list with raw averaged data and model preds
preds <- get_predictions(d_RT, fit.2_pow, epoch = T, exp="e2")

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
  facet_wrap(.~Group) +
  theme(legend.position = "top",
        legend.spacing.y = unit(.2, 'cm'),
        plot.margin = unit(c(2, 2, 2, 2), "mm")
  )

ppreds

# Conditional effect of Singleton:
# Returns a list with raw averaged data and model conditional effects
comps <- get_comparisons(d_RT, fit.2_pow, exp = "e2", epoch = T)

comps[["mod"]] <- comps[["mod"]] %>% mutate(Block = Block_num, 
                                            Group = ifelse(Group == .5, "Instructions", "No instructions"))

comps[["raw"]] <- comps[["raw"]] %>% mutate(Group = ifelse(Group == .5, "Instructions", "No instructions"))

pcomps <- ggplot(data = comps[["mod"]],
                 aes(
                   y = estimate,
                   x = Block,
                   color = Effect,
                   fill = Effect
                 )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ Group, scales = "free_x") +
  geom_point(data = comps[["raw"]],
             aes(x = as.numeric(Block)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]],
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 12) + # text size is adjusted for DPI
  theme(legend.position = "bottom",
        legend.spacing.y = unit(.2, 'cm'),
        plot.margin = unit(c(2, 2, 2, 2), "mm"),
        strip.text = element_blank()
  )

ggarrange(ppreds, pcomps, labels = c("a", "b"), nrow = 2)

ggsave(
  "Output/plots/figure3.png",
  height = 20,
  width = 22,
  unit = "cm",
  dpi = 1200
)

# Analysis of contingency awareness ----
df_cont <- d_RT %>%
  dplyr::summarise(Contingency = max(Contingency),
                   Confidence = max(Confidence),
                   .by = c("ID", "group")) %>%
  drop_na() %>%
  left_join(., get_effects(d_RT)[,c("ID", "VMAC")], by = "ID")

car::qqPlot(df_cont$Contingency[df_cont$group == "A"]) # Normally distributed
car::qqPlot(df_cont$Contingency[df_cont$group == "B"]) # Normally distributed
car::qqPlot(df_cont$Confidence[df_cont$group == "A"]) # Normally distributed
car::qqPlot(df_cont$Confidence[df_cont$group == "B"]) # Normally distributed
car::qqPlot(df_cont$VMAC[df_cont$group == "A"]) # Normally distributed
car::qqPlot(df_cont$VMAC[df_cont$group == "B"]) # Normally distributed

# Correlation between contingency and Confidence
cor.test(~Contingency + Confidence, data = df_cont[df_cont$group == "A",])
cor.test(~Contingency + Confidence, data = df_cont[df_cont$group == "B",])


# There is an effect of information on Contingency?
t.test(Contingency~group, data = df_cont)


# There is an effect of information on Confidence?
t.test(Confidence~group, data = df_cont)

# Correlation between contingency and VMAC
cor.test(~Contingency + VMAC, data = df_cont[df_cont$group == "A",])
cor.test(~Contingency + VMAC, data = df_cont[df_cont$group == "B",])

# Plots contingency awareness
df_cont$sVMAC <- c(scale(df_cont$VMAC)[,1])
df_cont$sCont <- c(scale(df_cont$Contingency)[,1])
df_cont$sConf <- c(scale(df_cont$Confidence)[,1])

df_cont$Group <- ifelse(df_cont$group == "A", "Instructions", "No instructions")

pp <- ggplot(df_cont, aes(sConf, sCont, color = Group, fill = Group)) +
  geom_point(alpha = .3) +
  geom_smooth(method = "lm") +
  scale_fill_manual(values=c("darkblue", "darkorange")) +
  scale_color_manual(values=c("darkblue", "darkorange")) +
  labs(x = "Scaled Confidence rating", y = "Scaled Contingency rating") +
  theme_Publication() +
  theme(legend.key.size = unit(.7, "cm"),
        legend.position = "top")
pp2 <- ggplot(df_cont, aes(sCont, sVMAC, color = Group, fill = Group)) +
  geom_point(alpha = .3) +
  geom_smooth(method = "lm") +
  scale_fill_manual(values=c("darkblue", "darkorange")) +
  scale_color_manual(values=c("darkblue", "darkorange")) +
  labs(x = "Scaled Contingency rating", y = "Scaled VMAC") +
  theme_Publication() +
  theme(legend.key.size = unit(.7, "cm"))

leg <- get_legend(pp)

plot_points <- ggMarginal(rm_legend(pp), alpha = .1,
                          type = "densigram", groupColour = T, groupFill = T)

plot_points2 <- ggMarginal(rm_legend(pp2), alpha = .1,
                           type = "densigram", groupColour = T, groupFill = T)

legend_p <- ggarrange(leg)
top <- ggarrange(plot_points, plot_points2,
                 labels = c("a", "b"), ncol = 2, common.legend = F)

# Complementary analysis: Is contingency awareness associated with the temporal dynamics of VMAC? ----
d <- d_RT[d_RT$group == "B",] # No instructions group

fit_ni <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * scale(Contingency) + (Singleton + scale(log(Block_num)) |
                                                                         ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d
  ) # Previous maximal model

summary(fit_ni) # Significant VMAC x Block x Contingency interaction

# Creating VMAC and AC conditional effects as a function of contingency
comps_ni <- comparisons(
  fit_ni,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Block_num = seq(1, 12, .01),
    Contingency = function(x) quantile(x, probs = c(.30, .80)), # Selecting predictions for different quantiles
    ID = NA),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(fit_ni) ^ 2) / 2) - exp(lo +
                                                                   (sigma(fit_ni) ^ 2) / 2)
) %>%
  mutate(Contingency = case_when(
    Contingency < quantile(df_cont$Contingency, probs = .6)[1]~"Low Contingency",
    T~"High Contingency",
  ))

comps_ni$Contingency <- factor(comps_ni$Contingency, levels = c("Low Contingency", "High Contingency"))

comps_ni$Effect <- ifelse(comps_ni$contrast == "High, Low", "VMAC", "AC")

comps_ni$Effect <- factor(comps_ni$Effect, levels = c("VMAC", "AC"))

# Get raw data summary statistics splitting the sample
rawCont <- get_effectCont(d, probs = .6)

# Plot both model conditional effects and raw data together
pcomp_ni <- ggplot(data = comps_ni,
                   aes(
                     y = estimate,
                     x = Block_num,
                     color = Effect,
                     fill = Effect,
                   )) +
  geom_line() +
  facet_wrap(~Contingency)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA)+
  #facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(data = rawCont,
             aes(x = as.numeric(Block)),
             position = position_dodge(.5)) +
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
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x="Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 12) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'),
        plot.margin = margin(t=".2", l = ".5", unit = "cm"),
        legend.position = "bottom")


# Complementary analysis: Same as before but on accuracy
# Format data

dacc <- d_acc[d_acc$group == "B",]

fitacc_ni <-
  glmer(
    acc ~ Singleton * scale(Block_num) * scale(Contingency) + (scale(log(Block_num)) |ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = dacc, family = binomial()
  ) # Same random effects structure as previous maximal model

summary(fitacc_ni) # No significant Contingency predictor nor interaction

sjPlot::tab_model(fit_ni, digits = 3, digits.re = 3) # table S2 (left)
sjPlot::tab_model(fitacc_ni, digits = 3, digits.re = 3) # table S2 (right)

