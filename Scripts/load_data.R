# Author: Francisco Garre-Frutos
# Date: 20/06/2024

# Script to load data;

# Installing (if needed) and loading the packages ----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(here, dplyr, tidyr, readxl)

setwd(here::here())

acc_ex <- list() # List to exclude participants based on accuracy

# Loading data from Garre-Frutos et al. (2024): ----
# Loading raw data of the previous experiment
raw_m <- read.csv("Input/raw_data_multi.csv")

# There are participants with less observations than the max number of observations?
check <- raw_m %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Rewarded", "Unrewarded")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 23 participants with missing data

# There are some participants with missing data. We will exclude those participants

raw_fm <- raw_m %>% filter(!ID %in% check$ID)

length(unique(raw_fm$ID)) # N = 193 participants

# Participants with less than .7 of accuracy?

acc_ex[["original"]] <- raw_fm %>%
  filter(Phase %in% c("Rewarded", "Unrewarded")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(Accuracy)) %>% filter(mean_ACC < .7) %>% pull(ID)

length(unique(acc_ex[["original"]])) # 11 participants with less than .7 of accuracy that will be excluded from the analysis

# How many participants?
raw_fm %>%
  filter(!ID %in% acc_ex[["original"]], Phase %in% c("Rewarded", "Unrewarded")) %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  group_by(Phase) %>%
  dplyr::count() # 182 participants with complete observations and ACC > .7

# Loading data from Experiment 1: ----

raw_e1 <- read.csv("Input/experiment_1.csv")

check <- raw_e1 %>% filter(Phase == "Rewarded") %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID) %>%
  filter(n < 288 |
           is.na(n)) # All participants had complete observations

check

# Participants with less than .7 of accuracy?

acc_ex[["e1"]] <- raw_e1 %>%
  filter(Phase %in% c("Rewarded")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(acc)) %>% filter(mean_ACC < .7) %>% pull(ID)


length(unique(acc_ex[["e1"]])) # 1 participants with less than .7 of accuracy that will be excluded from the analysis

# Number of Participants per group
length(unique(raw_e1$ID))

# Descriptions of data:
# % of excluded trials
(1 - (nrow(
  filter_data(
    raw_e1,
    f.absent = F,
    acc = T,
    experiment = "e1",
    fixed = F,
    sd_filter = NULL
  )
) /
  nrow(
    filter_data(
      raw_e1,
      f.absent = F,
      acc = F,
      experiment = "e1",
      fixed = F,
      sd_filter = NULL
    )
  )))*100
# 5.08% incorrect responses

(1 - (nrow(
  filter_data(
    raw_e1,
    f.absent = F,
    acc = T,
    experiment = "e1",
    fixed = T,
    sd_filter = NULL
  )
) /
  nrow(
    filter_data(
      raw_e1,
      f.absent = F,
      acc = T,
      experiment = "e1",
      fixed = F,
      sd_filter = NULL
    )
  )))*100

# 0.08 % outliers (RT < 1800 & > 150 ms)


(1 - (nrow(
  filter_data(
    raw_e1,
    f.absent = F,
    acc = T,
    experiment = "e1",
    fixed = T,
    sd_filter = 2
  )
) /
  nrow(
    filter_data(
      raw_e1,
      f.absent = F,
      acc = T,
      experiment = "e1",
      fixed = T,
      sd_filter = NULL
    )
  )))*100

# 4.25 % outliers (RT <> 2 SDs)

# Loading data from Experiment 2: ----

raw_e2 <- read.csv("Input/experiment_2.csv")

raw_e2 <- raw_e2[raw_e2$Phase == "Rewarded", ] %>%
  dplyr::rename("Trials" = "trial_num", "acc" = "correct")

check <- raw_e2 %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  filter(n < 288 |
           is.na(n))# All participants had complete observations

check

# Participants with less than .7 of accuracy?

acc_ex[["e2"]] <- raw_e2 %>%
  filter(Phase %in% c("Rewarded")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(acc)) %>% filter(mean_ACC < .7) %>% pull(ID)

length(unique(acc_ex[["e2"]])) # 3 participants with less than .7 of accuracy that will be excluded from the analysis

# Number of Participants per group
raw_e2 %>%
  group_by(ID, Phase, group) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID) %>%
  filter(!ID %in% acc_ex[["e2"]], !is.na(ID)) %>%
  dplyr::summarise(n = n(), .by = "group")

# Descriptions of data:
# % of excluded trials
(1 - (nrow(
  filter_data(
    raw_e2,
    f.absent = F,
    acc = T,
    experiment = "e2",
    fixed = F,
    sd_filter = NULL
  )
) /
  nrow(
    filter_data(
      raw_e2,
      f.absent = F,
      acc = F,
      experiment = "e2",
      fixed = F,
      sd_filter = NULL
    )
  )))*100

# 6.54% incorrect responses

(1 - (nrow(
  filter_data(
    raw_e2,
    f.absent = F,
    acc = T,
    experiment = "e2",
    fixed = T,
    sd_filter = NULL
  )
) /
  nrow(
    filter_data(
      raw_e2,
      f.absent = F,
      acc = T,
      experiment = "e2",
      fixed = F,
      sd_filter = NULL
    )
  )))*100

# 0.45 % outliers (RT < 1800 & > 150 ms)


(1 - (nrow(
  filter_data(
    raw_e2,
    f.absent = F,
    acc = T,
    experiment = "e2",
    fixed = T,
    sd_filter = 2
  )
) /
  nrow(
    filter_data(
      raw_e2,
      f.absent = F,
      acc = T,
      experiment = "e2",
      fixed = T,
      sd_filter = NULL
    )
  )))*100

# 4.41 % outliers (RT <> 2 SDs)
