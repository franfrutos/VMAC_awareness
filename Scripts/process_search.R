# Author: Francisco Garre-Frutos

# Process search for meta-analysis

# Load packages
if(!require(pacman)){
  install.packages("pacman")
  library(pacman)
}

p_load(readxl, dplyr, here)

# WOS search ----

cols_f <- c("Publication Type",
            "Authors",
            "Publication Year",
            "Article Title",
            "Abstract",
            "DOI")

wos_l <- list()
files <- list.files(path = "Input/search", pattern = "wos")

for (file in files) {
  wos_l[[file]] <- read_excel(here(paste0("Input/search/", file)))[, cols_f] %>% filter(`Publication Type` == "J") %>% select(-`Publication Type`) %>%
    rename("Title" = "Article Title",
           "Year" = "Publication Year")
}

wos <- do.call(rbind, wos_l)
wos <- wos[which(!duplicated(wos$DOI) & !is.na(wos$DOI)),]

# WOS research papers
nrow(wos) # N = 301

# Scopus search ----
cols_f <- c("Document.Type",
            "Authors",
            "Title",
            "Year",
            "Abstract",
            "DOI")


scopus_l <- list()
files <- list.files(path = "Input/search", pattern = "scopus")
for (file in files) {
  scopus_l[[file]] <- read.csv(here(paste0("Input/search/", file)))[, cols_f] %>%
    filter(`Document.Type` == "Article") %>%
    select(-Document.Type)
}

scopus <- do.call(rbind, scopus_l)
scopus <- scopus[which(!duplicated(scopus$DOI) & !is.na(scopus$DOI)),]

# WOS research papers
nrow(scopus) # N = 164

# Joint WOS and SCOPUS ----

all <- bind_rows(wos, scopus)

nrow(all) # N = 593

# Removing duplicates
all_distinct <- all[which(!duplicated(all$DOI) & !is.na(all$DOI)),]
nrow(all_distinct) # N = 231

# Save processed search
if (exists(here("Output/processed_Search.csv"))) {
  all_current <- read.csv(here("Output/processed_Search.csv"))
} else {
  write.csv(all_distinct, here("Output/processed_search.csv"))
}



