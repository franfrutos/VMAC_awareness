# Author: Francisco Garre-Frutos
# Date: 20/06/2024

# Script to create the project structure. Only run at the beginning of the project. 

# Load relevant packages ----
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(
  here # Package to manage paths relative to the working directory
  )

# Creating project structure ----
# Folders for data and analysis scripts of the main analyses:
dir.create("Input")
dir.create("Output")
dir.create(here("Output/plots"))
dir.create("Scripts")

# Aditional folders for materials and the re-analysis of Rusz et al. (2020)
dir.create("Materials")
dir.create("Re_rusz_2020")

# Creating script files to perform analyses
file.create(here("Scripts/functions.R"))
file.create(here("Scripts/load_data.R"))
file.create(here("Scripts/power_analysis.R"))
file.create(here("Scripts/e1_analysis.R"))
file.create(here("Scripts/e2_analysis.R"))
file.create(here("Scripts/multiverse_analysis.R"))