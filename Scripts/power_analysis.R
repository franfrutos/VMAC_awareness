# Author: Francisco Garre-Frutos

# Power analysis from Garre-Frutos et al. (2024)

# Load packages ----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr, lme4, lmerTest, hypr, simr, pbapply, future, parallel, ggplot2, pwr)

source("scripts/functions.R") # Load data and all relevant functions

# If the file with the power analysis result is detected, do not run power analysis
if (!file.exists("Output/power.rds")) {
  
  # Setting model contrasts with hypr package
  HcRep <- hypr(
    VMAC = High ~ Low,
    AC = Low ~ Absent,
    levels = c("High","Low", "Absent")
  )
  
  d <- filter_data(raw_fm, f.absent = F, experiment = "original") %>% filter(Phase == "Rewarded")
  d$rt <- log(d$RT)
  d$Singleton <- factor(d$Singleton, levels= c("High","Low", "Absent"))
  contrasts(d$Singleton) <- HcRep
  d$c_block <- scale(log(d$Block), scale = F)[,1]
  
  d_test <- d %>%
    filter(Singleton != "Absent") %>%
    dplyr::summarise(rt = mean(rt),
                     .by = c(ID, Singleton)) %>%
    spread(Singleton, rt)
  
  t.test(d_test$High, d_test$Low, paired = T)

  fit.2_power <- lmer(rt~Singleton*c_block+(Singleton+c_block|ID),control=lmerControl(optimizer='bobyqa'),
                      data = d)
  
  
  summary(fit.2_power)
  
  ncls <- 20
  steps <- seq(20, 120, 20)
  nsim <- 1e3
  seed <- 123
  li <- list()
  for (i in 1:length(steps)) {
    print(paste("N:", steps[i], "Participants"))
    set.seed(seed)
    model <- extend(fit.2_power, along = "ID", n = steps[i]) #extend the model to include n participants
    dat <- getData(model)
    cl <- makeCluster(ncls)
    clusterExport(cl, c("sim_lme4", "doSim", "lmer", "dat", "model", "fit.2_power"))
    power <- pbreplicate(n = nsim, expr = sim_lme4(dat, model), simplify = T, cl = cl)
    stopCluster(cl)
    li[[paste(i)]]$VMAC <- binom::binom.confint(sum(unlist(power[1,])), nsim, method = "exact")
    li[[paste(i)]]$Int <- binom::binom.confint(sum(unlist(power[2,])), nsim, method = "exact")
    print(li[[i]])
  }
  
  
  power <- do.call(rbind,do.call(rbind, li))
  
  power$N <- rep(steps, length.out = nrow(power_1b))
  
  power$Effect <- rep(c("VMAC effect", "VMAC x Block interaction"), each = 6)
  
  saveRDS(power, file = "Output/power.rds")
  
} else {
  power <- readRDS("Output/power.rds")
}

power$Effect <- ifelse(power$Effect == "VMAC effect", "VMAC", "VMAC x Block interaction")

plot_power <- 
  ggplot(power, aes(x = N, y = mean, color = Effect)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = unique(power$N), labels = unique(power$N)) +
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(labels = paste0(seq(0, 1, .1)*100, "%"), breaks = seq(0, 1, .1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Effect), alpha = 1/3, color = NA)+
  geom_hline(yintercept = .80, linetype = "dashed") +
  labs(y = "Power", x = "Number of participants")+
  scale_fill_manual(values=c("darkblue", "darkorange")) +
  scale_color_manual(values=c("darkblue", "darkorange")) +
  theme_Publication() +
  theme(text = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0))

ggsave(plot = plot_power,
       filename = "Output/plots/figureS1.png",
       height = 8,
       width = 10,
       units = "cm",
       dpi = 900)

# Sensitivity analysis:
# Experiment 1
pwr.r.test(n = 82, power = .80)
pwr.r.test(n = 81, power = .80)


