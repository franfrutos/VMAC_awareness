# Author: Francisco Garre-Frutos

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
  ggpubr, # Combined plots
  Rmisc,
  effectsize
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
Anova(fit.2_pow) # Anova-like output

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

# BF for instructed group
dBF <- d_RT %>%
  filter(Singleton != "Absent") %>%
  dplyr::summarise(rt = mean(rt),
                   .by = c(group, ID, Singleton)) %>%
  filter(group == "A") %>%
  spread(Singleton, rt)
  
t.test(dBF$High, dBF$Low, paired = T)
t_to_d(4.3052, 80, paired = T)
BayesFactor::ttestBF(dBF$High, dBF$Low, paired = T)

# BF for uninstructed group
dBF <- d_RT %>%
  filter(Singleton != "Absent") %>%
  dplyr::summarise(rt = mean(rt),
                   .by = c(group, ID, Singleton)) %>%
  filter(group == "B") %>%
  spread(Singleton, rt)

t.test(dBF$High, dBF$Low, paired = T)
t_to_d(1.4278, 80, paired = T)
1/BayesFactor::ttestBF(dBF$High, dBF$Low, paired = T)

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
car::Anova(fit.acc3) # ANOVA-like table

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
  theme_Publication(text_size = 12) + 
  facet_wrap(.~Group) +
  theme(legend.position = "top",
        legend.key.size = unit(1, "cm"),
        legend.spacing.y = unit(.2, 'cm'),
        plot.margin = unit(c(2, 2, 2, 2), "mm")
  ))


# Conditional effect of Singleton:
# Returns a list with raw averaged data and model conditional effects
comps <- get_comparisons(d_RT, fit.2_pow, exp = "e2", epoch = T)

comps[["mod"]] <- comps[["mod"]] %>% mutate(Block = Block_num, 
                                            Group = ifelse(Group == .5, "Instructions", "No instructions"))

comps[["raw"]] <- comps[["raw"]] %>% mutate(Group = ifelse(Group == .5, "Instructions", "No instructions"))

(pcomps <- ggplot(data = comps[["mod"]],
                 aes(
                   y = estimate,
                   x = Block,
                   color = Effect,
                   fill = Effect
                 )) +
  geom_line(size = 1, aes(linetype = Effect)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .15,
              color = NA) +
  geom_point(data = comps[["raw"]],
             aes(x = as.numeric(Block), shape = Effect),
             position = position_dodge(.5),
             size = 2.5) +
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
  facet_wrap(.~Group) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_manual(values = list_colors$Effects, name = "Contrast") +
  scale_fill_manual(values = list_colors$Effects, name = "Contrast") +
  scale_shape(name = "Contrast") +
  scale_linetype(name = "Contrast") +
  labs(y = "Contrast (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 12) + # text size is adjusted for DPI
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.spacing.y = unit(.2, 'cm'),
        plot.margin = unit(c(2, 2, 2, 2), "mm"),
        strip.text = element_blank()
  ))

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


# Correlation between contingency and VMAC
cor.test(~Contingency + VMAC, data = df_cont[df_cont$group == "A",])
sp_correction(.016, .58, .7)

cor.test(~Contingency + VMAC, data = df_cont[df_cont$group == "B",])
sp_correction(.13, .42, .7)


# Plots contingency awareness
df_cont$Group <- ifelse(df_cont$group == "A", "Instructions", "No instructions")

avg_data_CC <- cbind(Rmisc::summarySE(data=df_cont,
                                   measurevar = "Confidence",
                                   groupvars = c("Group")) %>%
                    dplyr::rename("Measure1" = "Confidence",
                                  "se1" = "se") %>%
                    select(Group, Measure1, se1),
                  Rmisc::summarySE(data=df_cont,
                                   measurevar = "Contingency",
                                   groupvars = c("Group"))%>%
                    dplyr::rename("Measure2" = "Contingency",
                                  "se2" = "se") %>%
                    select(Measure2, se2))

avg_data_CVMAC <- cbind(Rmisc::summarySE(data=df_cont,
                                      measurevar = "VMAC",
                                      groupvars = c("Group")) %>%
                       dplyr::rename("Measure1" = "VMAC",
                                     "se1" = "se") %>%
                       select(Group, Measure1, se1),
                     Rmisc::summarySE(data=df_cont,
                                      measurevar = "Contingency",
                                      groupvars = c("Group"))%>%
                       dplyr::rename("Measure2" = "Contingency",
                                     "se2" = "se") %>%
                       select(Measure2, se2))


(pp <- ggplot(df_cont, aes(Confidence/100, Contingency/100, color = Group, fill = Group, shape = Group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_point(alpha = .3) +
  geom_point(data=avg_data_CC, aes(Measure1/100, Measure2/100), size = 3) +  
  geom_errorbar(data=avg_data_CC, aes(ymin = Measure2/100 - se2/100, ymax = Measure2/100 + se2/100, y = Measure2/100, x = Measure1/100),
                width = 0, size = 1) +
  geom_errorbarh(data=avg_data_CC, aes(xmin = Measure1/100 - se1/100, xmax = Measure1/100 + se1/100, y = Measure2/100, x = Measure1/100),
                 height = 0, size = 1) +
  geom_smooth(method = "lm", alpha = .15, aes(linetype = Group)) +
  scale_fill_manual(values=c("#252373", "#d16606")) +
  scale_color_manual(values=c("#252373", "#d16606")) +
  labs(x = "Confidence rating", y = "Contingency rating") +
  theme_Publication(text_size = 14) +
  theme(legend.key.size = unit(1, "cm"),
        legend.position = "top"))
(pp2 <- ggplot(df_cont, aes(Contingency/100, VMAC, color = Group, fill = Group, shape = Group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_point(alpha = .3) +
  geom_point(data=avg_data_CVMAC, aes(y=Measure1, x=Measure2/100), size = 3) + 
  geom_errorbarh(data=avg_data_CVMAC, aes(xmin = Measure2/100 - se2/100, xmax = Measure2/100 + se2/100, x = Measure2/100, y = Measure1),
                height = 0, size = 1) +
  geom_errorbar(data=avg_data_CVMAC, aes(ymin = Measure1 - se1, ymax = Measure1 + se1, x = Measure2/100, y = Measure1),
                 width = 0, size = 1) +
  geom_smooth(method = "lm", alpha = .15, aes(linetype = Group)) +
  scale_fill_manual(values=c("#252373", "#d16606")) +
  scale_color_manual(values=c("#252373", "#d16606")) +
  labs(x = "Contingency rating", y = "VMAC scores (ms)") +
  theme_Publication(text_size = 14) +
  theme(legend.key.size = unit(1, "cm")))

leg <- get_legend(pp)

plot_points <- ggMarginal(rm_legend(pp), alpha = .1,
                          type = "densigram", groupColour = T, groupFill = T)

plot_points2 <- ggMarginal(rm_legend(pp2), alpha = .1,
                           type = "densigram", groupColour = T, groupFill = T)

legend_p <- ggarrange(leg)
top <- ggarrange(plot_points, plot_points2,
                 labels = c("a", "b"), ncol = 2, common.legend = F)
top
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
(pcontingency2 <- ggplot(data = comps_ni,
                         aes(
                           y = estimate,
                           x = Block_num,
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
    labs(y = "Contrast (ms)", x = "Block") +
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme_Publication(text_size = 14) + # text size is adjusted for DPI
    theme(legend.position = "bottom",
          legend.key.size = unit(1, "cm"),
          legend.spacing.y = unit(.2, 'cm'),
          plot.margin = unit(c(2, 2, 2, 2), "mm"),
          strip.text = element_blank()
    ))


ggarrange(leg, top, pcontingency1, pcontingency2, nrow= 4, heights= c(.07, rep((1-.07)/3, times = 3)),
          labels = c("", "", "c", "d"))

ggsave(
  "Output/plots/figure4.png",
  height = 32,
  width = 22,
  dpi = 900,
  units = "cm"
) 

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

