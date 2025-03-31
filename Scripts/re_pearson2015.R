# Author: Francisco Garre-Frutos

# re-analysis of Pearson et al. (2015)

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
  stringr
)

contingencyP <- read_xlsx("Input/data/contingencyPearson.xlsx") %>%
  dplyr::rename("ID" = "subNum") %>%
  mutate(Correct = ifelse(Awareness > 0, 1, 0))

Rmisc::summarySE(
  data = contingencyP,
  measurevar = "Awareness"
)

ggplot(contingencyP, aes(Awareness)) +
  geom_histogram(aes(y = after_stat(..density..)), alpha = .2, color = "black", bins = 11) +
  geom_density(alpha = .2, fill = "gray") +
  geom_jitter(aes(y = -.02), height = .01, shape = 1, width = 0) +
  theme_Publication() +
  scale_x_continuous(breaks = seq(-5,5,1)) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(y = "Denisity", x = "Contingency Awareness") +
  theme(legend.position = "none")

ggsave(filename = "Output/plots/figureS2_2.png",
       height = 12,
       width = 15,
       units = "cm",
       dpi = 900)

# High overall contingency
t.test(contingencyP$Awareness)

# High probability of reporting the correct contingency
binom.test(sum(contingencyP$Correct), nrow(contingencyP), p = .5)

# Relationship between contingency test and condifence level
glm_fit <- glm(Correct ~ Awareness, data = contingencyP %>% mutate(Awareness  = abs(Awareness)),
               family = binomial)

summary(glm_fit) # Significant effect of awareness

plot_predictions(glm_fit, condition = "Awareness") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_jitter(data=contingencyP %>% mutate(Awareness  = abs(Awareness)),
             aes(y = ifelse(Correct == 1, 1, .2), x = Awareness), height = .01, width = 0,
             shape = 1) +
  labs(y = "Pr(Contingency Awareness = 1)", x = "Confidence level") +
  coord_cartesian(ylim = c(0.2, 1)) +
  theme_Publication()

ggsave(filename = "Output/plots/figureS3_2.png",
       height = 12,
       width = 15,
       units = "cm",
       dpi = 900)


# Re-analysis
d <- read.csv(file="Input/data/pearson2015.csv")[,-1]

# Filter data
d <- d[which((d$timeout != 1 & d$trialPropGoodSamples > .25) &
               !d$trial %in% c(seq(1, 480, 96),seq(2, 480, 96))),]

d_agg <- d %>%
  dplyr::summarise(Prop = mean(omissionTrial),
                   Awareness = max(Awareness),
                   N = n(), .by = c(ID, block, Singleton)) %>%
  mutate(Correct = ifelse(Awareness > 0, 1, 0)) %>%
  left_join(., contingencyP[, c("ID", "Awareness", "Correct")], by = "ID") %>%
  select(ID, block, Singleton, Prop, N, Awareness.x, Correct.x) %>%
  dplyr::rename("Awareness" = "Awareness.x",
         "Correct" = "Correct.x")
HcRep <- hypr(
  VMAC = High ~ Low,
  AC = Low ~ Absent,
  levels = c("High", "Low", "Absent")
) 

d_agg$Singleton <- factor(d_agg$S, levels = c("High", "Low", "Absent"))
contrasts(d_agg$Singleton) <- contr.hypothesis(HcRep)

# Abalysis with the same awareness measure
fit_pearson <- glmer(Prop ~ scale(block)*Singleton*scale(Awareness) + (scale(block)*Singleton | ID),
      data = d_agg, family = binomial, weights = N,
      control = glmerControl(optimizer = 'bobyqa'),
) # Singlular fit

fit_pearson <- glmer(Prop ~ scale(block)*Singleton*scale(Awareness) + (scale(block)+Singleton | ID),
                     data = d_agg, family = binomial, weights = N,
                     control = glmerControl(optimizer = 'bobyqa'),
) 

summary(fit_pearson)

sjPlot::tab_model(fit_pearson)

plot_predictions(fit_pearson,
                 condition = c("block", "Singleton", "Awareness"),
                 re.form = NA)

# Analysis as a function of whether participants selected the correct contingency or not. 
fit_pearson2 <- glmer(Prop ~ scale(block)*Singleton*Correct + (scale(block)+Singleton | ID),
                     data = d_agg, family = binomial, weights = N,
                     control = glmerControl(optimizer = 'bobyqa'),
) 

summary(fit_pearson2)
sjPlot::tab_model(fit_pearson2)


d_agg$Block <- create_epochs(d_agg$block)
agg_pearson <- summarySEwithin2(
  d_agg, measurevar = "Prop",
  betweenvars = "Correct",
  withinvars = c("Singleton", "Block")
)
agg_pearson$Block <- as.numeric(agg_pearson$Block) * 2 - .5

labels <- list(
  "0"="Unaware",
  "1"="Aware"
)
phase_labeller <- function(variable,value){
  return(labels[value])
}

predictions(
  fit_pearson2,
  newdata = datagrid(
    ID = unique(d_agg$ID),
    Correct = c(0, 1),
    Singleton = c("High", "Low", "Absent"),
    block = 1:20
  ),
  by = c("Correct", "Singleton", "block")
) %>%
  ggplot(aes(x = block, y = estimate)) +
  geom_line(aes(color = Singleton), size = 1) +
  geom_point(data = agg_pearson,
             aes(x = as.numeric(Block), color = Singleton, y = Prop),
             position = position_dodge(.5),
             size = 2.5) +
  geom_ribbon(
    aes(
      x = as.numeric(block),
      ymin = conf.low,
      ymax = conf.high,
      y = estimate,
      fill = Singleton
    ),
    alpha = .15,
    position = position_dodge(.5),
  )  +
  geom_errorbar(
    data = agg_pearson,
    aes(
      x = as.numeric(Block),
      ymin = Prop - se,
      ymax = Prop + se,
      y = Prop,
      color = Singleton
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  facet_wrap(.~Correct, labeller = phase_labeller) +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  coord_cartesian(ylim = c(0, .25), xlim = c(1, 20)) +
  scale_color_manual(values = list_colors$singleton, name = "Singleton") +
  scale_fill_manual(values = list_colors$singleton, name = "Singleton") +
  labs(x = "Blocks", y = "Propotion of omissions") +
  theme_Publication(text_size = 14)
  
ggsave(
  "Output/plots/figureS11.png",
  height = 15,
  width = 25,
  unit = "cm",
  dpi = 900
)


