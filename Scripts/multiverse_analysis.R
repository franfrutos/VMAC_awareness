# Author: Francisco Garre-Frutos
# Date: 20/06/2024

# Multiverse analysis reported as supplementary material

source("scripts/functions.R") # This will load functions and data into the workspace

# Load and install packages ----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

# Install a custom function apadted from Parsons (2021) splith package
# In the main text is cited as Parsons (2021).
if (!require(multi.s)) {
  if (!require(devtools))
    install.packages("devtools")
  devtools::install_github("franfrutos/multi.s") 
  library(multi.s)
}

p_load(here, # Function to easily control paths
       dplyr,# Data wrangling
       tidyr,# Data wrangling
       rstatix, # Convert as factor
       doParallel,# Parallel computing
       snow, # If you want to run the multiverse using parallel computing
       ggplot2, # Plotting
       patchwork, # Functions to create multiplots
       ggrain # Rainclould plots
)

source("scripts/functions.R") # Load functions and raw data

# Multiverse analysis ----

# creating speciffcations:
specifications <- list(
  RT_sd_cutoff      = c(0, 2, 2.5, 3),
  RT_fix_cutof = c(T, F),
  averaging_method  = c("mean", "median"),
  logtransform     = c("yes", "no"),
  NBlocks          = c(6, 12),
  Two_Trials       = c(T, F)
)

# Constructing a data frame with all possible combinations
spec <- expand.grid(specifications) %>%
  mutate(Universe = 1:nrow(.)) %>%
  as_tibble() %>%
  group_by(Universe) %>%
  nest()

# Multiverse e1 ----

if (!file.exists("Output/data_for_multi1.rds")) {
  # remove the file in Output/ to generate datasets again
  dlist = list()
  iter = 0
  
  v <-
    c(seq(1, 577, 24), seq(2, 577, 24)) # vector to filter the 2 first trials of each block
  
  #constructing data frames of all universes:
  for (U in spec$data) {
    acc_ex_m <- raw_e1 %>%
      group_by(ID) %>%
      dplyr::summarise(mean_ACC = mean(acc)) %>% filter(mean_ACC < .7) %>% dplyr::pull(ID)
    
    
    if (U$Two_Trials == TRUE) {
      rawFPar <- raw_e1 %>%
        filter(!(Trials %in% v))
    }
    
    rawfiltered <- rawFPar %>%
      filter(!ID %in% acc_ex_m)
    
    
    rawRT <- rawfiltered %>%
      filter(acc == 1, !Singleton %in% c("Absent", "None","Singleton"))
    
    
    if (U$RT_fix_cutof == T)
      rawRT <- rawRT[which(rawRT$rt > 150 & rawRT$rt < 1800), ]
    
    if (U$RT_sd_cutoff != 0) {
      rawRT <- rawRT %>%
        group_by(ID) %>%
        dplyr::mutate(
          high  = mean(rt) + (U$RT_sd_cutoff * sd(rt)),
          low   = mean(rt) - (U$RT_sd_cutoff * sd(rt))
        ) %>%
        dplyr::filter(rt >= low & rt <= high) %>%
        ungroup()
    }
    
    if (U$NBlocks == 6) {
      rawRT <- rawRT %>%
        filter(Block > 6)
    }
    
    if (U$logtransform == "yes") {
      rawRT$rt <- log(rawRT$rt)
    }
    iter = iter + 1
    Sys.sleep(.01)
    cat("\r", iter, "of", nrow(spec), "Specifications")
    dlist[[iter]] <- rawRT
  }
  
  
  # Save data for later use:
  saveRDS(dlist, "output/data_for_multi1.rds")
}

if (!exists("dlist"))
  dlist <- readRDS(here::here("output/data_for_multi1.rds"))
if (!file.exists("output/multiverse_estimates_rel1.rds")) {
  # If you want to run the multiverse, delete multiverse_estimates.rds from Output/
  permutations <- 5e3
  iter <- nrow(spec)
  
  #cl <- makeCluster(2, type = "SOCK") # Uncomment this for parallel computing
  #registerDoParallel(cl) # Uncomment this for parallel computing
  
  # This may take some time...
  multiverse_analisysR <-
    foreach(
      i = icount(iter),
      .combine = f(iter),
      .packages = c("multi.s", "dplyr"),
      .errorhandling = 'stop'
    ) %do% {
      # change %do% for %dopar% for parallel computing
      multi.s::splith(
        outcome = "rt",
        data = dlist[[i]],
        permutations = permutations,
        variable = "Singleton",
        subject = "ID",
        include_block = T,
        block = "Block",
        average = spec$data[[i]][[3]]
      )$final_estimates
    }
  #stopCluster(cl) # Uncomment this for parallel computing
  
  # Prepare multiverse data to be saved
  multiverse_analisysR <- multiverse_analisysR %>%
    mutate(Universe = 1:iter) %>%
    group_by(Universe) %>%
    nest()
  
  # Data with specification information and reliability estimates (only final estimates)
  multi_data <- tibble(spec,
                       Estimates = multiverse_analisysR$data)
  
  # Saving data
  saveRDS(multi_data, "output/multiverse_estimates_rel1.rds")
} else {
  estimates <- readRDS(here::here("output/multiverse_estimates_rel1.rds"))
}

if (!exists("multi_data")) {
  multi_data <- readRDS("output/multiverse_estimates_rel1.rds")
}


estimates <- multi_data[, 3] %>% unnest(cols = c("Estimates"))

# Descriptives:
estimates %>% arrange(spearmanbrown) %>% .[nrow(.)/2,]
range(estimates$spearmanbrown)
mean(estimates$spearmanbrown > .7) *100

# Create plot for reliability multiverse:

# Create plot for reliability multiverse:
n_spec <- cbind(estimates,multi_data[,2] %>% unnest()) %>%
  arrange(spearmanbrown) %>%
  mutate(Universe = 1:128) %>%
  filter(RT_sd_cutoff == 2, RT_fix_cutof == T,
         averaging_method == "mean", logtransform == "no",
         NBlocks == 12, Two_Trials == F) %>% pull(Universe)


mr_plot_e1 <- multiverse_plot(estimates, multi_data, spec_seg = n_spec)

ggsave(
  plot = mr_plot_e1,
  height = 20,
  width = 17,
  dpi = 900,
  filename = "Output/plots/figureS2.png",
  units = "cm"
)

# Multiverse e2 ----

if (!file.exists("Output/data_for_multi2.rds")) {
  # remove the file in Output/ to generate datasets again
  dlist = list()
  iter = 0
  
  v <-
    c(seq(1, 577, 24), seq(2, 577, 24)) # vector to filter the 2 first trials of each block
  
  #constructing data frames of all universes:
  for (U in spec$data) {
    acc_ex_m <- raw_e2 %>%
      group_by(ID) %>%
      dplyr::summarise(mean_ACC = mean(acc)) %>% filter(mean_ACC < .7) %>% dplyr::pull(ID)
    
    
    if (U$Two_Trials == TRUE) {
      rawFPar <- raw_e2 %>%
        filter(!(Trials %in% v))
    }
    
    rawfiltered <- rawFPar %>%
      filter(!ID %in% acc_ex_m)
    
    
    rawRT <- rawfiltered %>%
      filter(acc == 1, Singleton != "Absent")
    
    if (U$RT_fix_cutof == T)
      rawRT <- rawRT[which(rawRT$rt > 150 & rawRT$rt < 1800), ]
    
    if (U$RT_sd_cutoff != 0) {
      rawRT <- rawRT %>%
        group_by(ID) %>%
        dplyr::mutate(
          high  = mean(rt) + (U$RT_sd_cutoff * sd(rt)),
          low   = mean(rt) - (U$RT_sd_cutoff * sd(rt))
        ) %>%
        dplyr::filter(rt >= low & rt <= high) %>%
        ungroup()
    }
    
    
    if (U$NBlocks == 6) {
      rawRT <- rawRT %>%
        filter(Block_num > 6)
    }
    
    if (U$logtransform == "yes") {
      rawRT$rt <- log(rawRT$rt)
    }
    iter = iter + 1
    Sys.sleep(.01)
    cat("\r", iter, "of", nrow(spec), "Specifications")
    dlist[[iter]] <- rawRT
  }
  
  # Save data for later use:
  saveRDS(dlist, "output/data_for_multi2.rds")
}

if (!exists("dlist"))
  dlist <- readRDS(here::here("output/data_for_multi2.rds"))

if (!file.exists("output/multiverse_estimates_rel2.rds")) {
  # If you want to run the multiverse, delete multiverse_estimates.rds from Output/
  permutations <- 5e3
  iter <- nrow(spec)
  
  #cl <- makeCluster(2, type = "SOCK") # Uncomment this for parallel computing
  #registerDoParallel(cl) # Uncomment this for parallel computing
  
  # This may take some time...
  multiverse_analisysR <- list()
  
  for (group_i in c("A", "B")) {
    multiverse_analisysR[[group_i]] <-
      foreach(
        i = icount(iter),
        .packages = c("multi.s", "dplyr"),
        .errorhandling = 'stop'
      ) %do% {
        # change %do% for %dopar% for parallel computing
        multi.s::splith(
          outcome = "rt",
          data = dlist[[i]] %>% filter(group == group_i),
          permutations = permutations,
          variable = "Singleton",
          #condition = "group",
          subject = "ID",
          include_block = T,
          block = "Block_num",
          average = spec$data[[i]][[3]]
        )$final_estimates
        #stopCluster(cl) # Uncomment this for parallel computing
      }
  }
  
  # Prepare multiverse data to be saved
  multiverse_analisys <- rbind(do.call(rbind, multiverse_analisysR[["A"]]),
                               do.call(rbind, multiverse_analisysR[["B"]])) %>%
    mutate(Universe = rep(1:iter, length.out = 128*2),
           Group = rep(c("A", "B"), each = 128)) %>%
    group_by(Universe) %>%
    nest()
  
  # Data with specification information and reliability estimates (only final estimates)
  multi_data <- tibble(spec,
                       Estimates = multiverse_analisys$data)
  
  # Saving data
  saveRDS(multi_data, "output/multiverse_estimates_rel2.rds")
} else {
  multi_data <- readRDS(here::here("output/multiverse_estimates_rel2.rds"))
}

if (!exists("multi_data")) {
  multi_data <- readRDS("output/multiverse_estimates_rel2.rds")
}

estimates <- multi_data[, 3] %>% unnest(cols = c("Estimates"))

# Descriptives:
# Instructions group:
estimates[estimates$Group == "A",] %>% arrange(spearmanbrown) %>% .[nrow(.)/2,]
range(estimates[estimates$Group == "A",]$spearmanbrown)
mean(estimates[estimates$Group == "A",]$spearmanbrown > .7) *100

# Create plot for reliability multiverse:
n_specA <- cbind(estimates[estimates$Group == "A",],spec[,2] %>% unnest()) %>%
  arrange(spearmanbrown) %>%
  mutate(Universe = 1:128) %>%
  filter(RT_sd_cutoff == 2, RT_fix_cutof == T,
         averaging_method == "mean", logtransform == "no",
         NBlocks == 12, Two_Trials == F) %>%pull(Universe)

mr_plot_e2_i <- multiverse_plot(estimates[estimates$Group == "A",], multi_data,
                                spec_seg = n_specA)

ggsave(
  plot = mr_plot_e2_i,
  height = 20,
  width = 17,
  dpi = 900,
  filename = "Output/plots/figureS3.png",
  units = "cm"
)

# No instructions group:
estimates[estimates$Group == "B",] %>% arrange(spearmanbrown) %>% .[nrow(.)/2,]
range(estimates[estimates$Group == "B",]$spearmanbrown)
mean(estimates[estimates$Group == "B",]$spearmanbrown > .7) *100

# Create plot for reliability multiverse:
n_specB <- cbind(estimates[estimates$Group == "B",],spec[,2] %>% unnest()) %>%
  arrange(spearmanbrown) %>%
  mutate(Universe = 1:128) %>%
  filter(RT_sd_cutoff == 2, RT_fix_cutof == T,
         averaging_method == "mean", logtransform == "no",
         NBlocks == 12, Two_Trials == F) %>% pull(Universe)

mr_plot_e2_in <- multiverse_plot(estimates[estimates$Group == "B",], multi_data,
                                 spec_seg = n_specB)
ggsave(
  plot = mr_plot_e2_in,
  height = 20,
  width = 17,
  dpi = 900,
  filename = "Output/plots/figureS4.png",
  units = "cm"
)

# Multiverse on the correlation between VMAC an the contigency rating for experiment 2 ----
# Compute correlation for each specification
dlist <- readRDS("output/data_for_multi2.rds")

cors <- lapply(dlist, \(x) {
  df <- x %>%
    dplyr::summarise(Contingency = max(Contingency),
                     .by = c("group", "ID")) %>%
    drop_na() %>%
    left_join(., get_effects(x)[,c("group", "ID", "VMAC")], by = c("group","ID"))
  
  A_cor <- with(df[df$group == "A",],
                cor.test(VMAC, Contingency))
  B_cor <- with(df[df$group == "B",],
                cor.test(VMAC, Contingency))
  
  return(
    data.frame(
      cor = c(A_cor$estimate, B_cor$estimate),
      lwr.ci = c(A_cor$conf.int[1],B_cor$conf.int[1]),
      ipr.ci = c(A_cor$conf.int[2],B_cor$conf.int[2]),
      Group = c("Instructions", "No instructions")
    )
  )
})

cors <- do.call(rbind, cors) 

multi_data <- readRDS("output/multiverse_estimates_rel2.rds")

estimates <- multi_data[, 3] %>% unnest(cols = c("Estimates"))


cor_multi <- rbind(
  cbind(cors[cors$Group == "Instructions",], 
        estimates[estimates$Group == "A", -ncol(estimates)]),
  cbind(cors[cors$Group == "No instructions",],
        estimates[estimates$Group == "B", -ncol(estimates)])
) %>% ungroup() 


# Plot multiverse for the instructions group:
corM_e2i <- multiverse_plot(cor_multi[cor_multi$Group == "Instructions",], multi_data, cor = T,
                            minimun_t = F, spec_seg = n_specA, y_name = "Correlation coeficient")

ggsave(
  plot = corM_e2i,
  height = 20,
  width = 17,
  dpi = 900,
  filename = "Output/plots/figureS5.png",
  units = "cm"
)

# Plot multiverse for the no instructions group:
corM_e2ni <- multiverse_plot(cor_multi[cor_multi$Group == "No instructions",], multi_data, cor = T,
                             minimun_t = F, spec_seg = n_specB, y_name = "Correlation coeficient")

ggsave(
  plot = corM_e2ni,
  height = 20,
  width = 17,
  dpi = 900,
  filename = "Output/plots/figureS6.png",
  units = "cm"
)

# Check the change in correlation as a function of the number of Blocks employed:
cor_multi <- cbind(cor_multi, spec %>% unnest()) %>% group_by(NBlocks)%>% mutate(ID = 1:128) %>% ungroup()

cor_multi$NBlocks <- factor(cor_multi$NBlocks, levels = c("12", "6"))

# Raincloud Plot:
rain <- ggplot(cor_multi, aes(as.factor(NBlocks), cor, color = Group, fill = Group), alpha = .3)+
  geom_rain(rain.side = 'f2x2', id.long.var = "ID", alpha = .3,
            violin.args = list(adjust = 2, trim = F, alpha = .3))+ 
  theme_Publication(text_size = 12) +
  scale_y_continuous(limits = c(-.3, .5), breaks = round(seq(-.3, .5, .1), 2)) +
  scale_fill_manual(values=c("darkblue", "darkorange")) +
  scale_color_manual(values=c("darkblue", "darkorange")) +
  labs(x = "Number of Blocks", y = "Correlation coefficient") +
  scale_x_discrete(labels = c("All blocks", "Last six blocks")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(color = "none", fill = "none") +
  theme(plot.margin = margin(r = "5", l = "5", b = "1", t = "1", unit = "cm"))

# To create this plot you need to run e2_analysis scipt to generate the top and bottom plots
ggarrange(leg, top, rain, pcomp_ni, nrow= 4, heights= c(.07, rep((1-.07)/3, times = 3)),
          labels = c("", "", "c", "d"))

ggsave(
  "Output/plots/figure4.png",
  height = 30,
  width = 22,
  dpi = 900,
  units = "cm"
) 

# Multiverse on the correlation between VMAC an the contigency rating for Experiment 1 ----
dlist <- readRDS("output/data_for_multi1.rds")
corse1 <- lapply(dlist, \(x) {
  df <- get_effects(x, exp = "e1")[,c("ID", "VMAC")] %>%
    left_join(raw_diff, ., by = c("ID"))
  
  cor <- with(df,
              cor.test(VMAC, points, method = "spearman"))
  return(
    tibble(
      cor = cor$estimate
    )
  )
})

corse1 <- do.call(rbind, corse1)

core1_multi <- cbind(corse1, spec %>% unnest()) %>% group_by(NBlocks) %>% mutate(ID = 1:64) %>% ungroup()

core1_multi$NBlocks <- factor(core1_multi$NBlocks, levels = c("12", "6"))

ggplot(core1_multi, aes(as.factor(NBlocks), cor))+
  geom_rain(rain.side = 'f1x1', id.long.var = "ID",
            violin.args = list(adjust = 2, trim = F, color = "darkorange",
                               fill = "darkorange", alpha = .3),
            boxplot.args = list(color = "darkorange",
                                fill = "darkorange", alpha = .3),
            point.args = list(color = "darkorange",
                              fill = "darkorange", alpha = .3),
            line.args = list(color = "darkorange", alpha = .3))+ 
  theme_Publication(text_size = 12) +
  scale_y_continuous(limits = c(-.1, .2), breaks = round(seq(-.1, .2, .1), 2)) +
  labs(x = "Number of Blocks", y = "Correlation coefficient") +
  scale_x_discrete(labels = c("All blocks", "Last six blocks")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(color = "none", fill = "none")

ggsave(
  "Output/plots/figureS7.png",
  height = 12,
  width = 15,
  unit = "cm",
  dpi = 1200
)

