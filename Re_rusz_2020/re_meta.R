#------------------------------------------------------------#
#         META-ANALYSIS ON REWARD-DRIVEN DISTRACTION 
#------------------------------------------------------------#

# This code is a re-analysis of Rusz et al. (2020). 
# Here I check for the effect of explicit instructions about stimulus-reward contingency
# This script is modified version of the one provided by Rusz et al. (2020): https://osf.io/rgeb6/

# load packages
library(metafor) #for meta-analysis
library(plyr) # for frequencies
library(dplyr) #for data manipulation
library(readxl) #for reading excel files

# read main meta-analysis dataset 
MA_final <- read_excel("MA_RDD_final_data.xlsx")
MA_final


#check variables
sapply(MA_final, class)
# it seems that many are characters so they have to be transformed, but only after creating a subset of the data


# in this part of the script, I only include studies 
# that measured high distractor reward against low or neutral distractor rewards

# I create subset for analysis on HIGH VS LOW 
# only those studies (rows) should be in the new file that contain values for high and low)
MA_high_low <- subset(MA_final, high_low == "y", select = -c(RT_M_no, RT_SD_no, ACC_M_high, ACC_M_low, ACC_SD_high, ACC_SD_low, ACC_M_no, ACC_SD_no))

# check if data file contains NRs or NAs, if not (=0) then go on
sum(MA_high_low$RT_M_high == 'NR', na.rm = T) 
sum(is.na(MA_high_low$RT_M_high))

#change Ms and SDs to numeric
MA_high_low$RT_M_high <- as.numeric(MA_high_low$RT_M_high)
MA_high_low$RT_M_low <- as.numeric(MA_high_low$RT_M_low)
MA_high_low$RT_SD_high <- as.numeric(MA_high_low$RT_SD_high)
MA_high_low$RT_SD_low <- as.numeric(MA_high_low$RT_SD_low)

# Compute effect sizes for high vs low RT #
# note, standardized mean change is the best measure, as recommended in R documentation
# https://cran.r-project.org/web/packages/metafor/metafor.pdf < page 100-101
SMD_high_low <- escalc(measure = "SMCC", # choose effect size measure: standardized mean change (change score standardization)
                       m1 = RT_M_high,            # high reward Mean
                       m2 = RT_M_low,             # low reward Mean
                       sd1 = RT_SD_high,          # high reward Standard Deviation
                       sd2 = RT_SD_low,           # low reward Standrad Devation
                       ni = nr_pp,                # sample size(N)
                       ri = ri,                   # correlation coefficient
                       data = MA_high_low)        # data set to reference

#check the new datasat > should contain 2 new columns, yi = effect size, vi = variance

# some of the effect sizes should be flipped:
# when RT -> high condition is higher (people take longer)
# when ACC -> high condition is lower (meaning people make more mistakes)
# so now, we flip the sign in front of the effect sizes that were calculated from accuracy measures
# I marked these studies in the data file with RT_Rev == reverse

# so, first we create a new variable
SMD_high_low$yi2 <- NA

# then, we take yi (effect size, but flip some of them)
SMD_high_low$yi2 <- ifelse(SMD_high_low$RT_Rev == "reverse", -SMD_high_low$yi, SMD_high_low$yi)


#-----------------------------------------------------------------------
# Do explicit information moderate the reward driven distraction effect?
#-----------------------------------------------------------------------

# Here I perform the moderator analysis taking into account information regarding stimulus-reward contingnecies:

count(SMD_high_low, feature_reward_inf) # 73 studies with no instructions, 15 studies with instructions and 3 studies where is not clear in the method section

# no = intercept; coeficients: no vs not stated and no vs yes
moderator_INF <- rma(yi2, vi, mods = ~ factor(feature_reward_inf), data = SMD_high_low[SMD_high_low$feature_reward_inf != "not stated",]) 

moderator_INF # significant increase in the effect size for studies that report using instructions vs no reported

# Predictions for instructed and no instructed learning:
instructions_predict <- rma(yi2, vi,
                            data = SMD_high_low[SMD_high_low$feature_reward_inf == "yes",])

no_instructions_predict <- rma(yi2, vi,
                            data = SMD_high_low[SMD_high_low$feature_reward_inf == "no",])


predict(instructions_predict) # Predictions for instructions
predict(no_instructions_predict) # Predictions for no instructions

# forest plot for instructions only studies:
# d <- SMD_high_low[SMD_high_low$feature_reward_inf == "yes",]
# 
# ES_INF <- rma(yi2, vi, data = d)
# 
# tiff("forest_high_low.tiff", width = 5500, height = 8000, units = "px", res = 600) # in this line of code, i tell R to start saving the file and give the limits
# par(mar=c(4,4,1,3))
# forest(ES_INF, 
#        slab = paste(d$Study), 
#        xlim = c(-5,5),
#        ylim=c(-1, nrow(d)+3),
#        addfit = T,
#        alim=c(-3,3),
#        cex.lab = .9,
#        digits=c(2,2), cex=.70, xlab = "Standardized Mean Change", order = "obs")
# text(-5.0, 93, "Author(s), Year, Study", cex=.9, pos=4, font = 2)
# text(5, 93, "SMCC [95% CI]", cex=.9, pos=2, font = 2)
# text(0, 93.5, "High vs Low", cex = 1.2, font = 2)
# dev.off()


# Check: Assuming that not stated are indeed no information provided
SMD_high_low$feature_reward_inf_2 <- ifelse(SMD_high_low$feature_reward_inf == "not stated" | SMD_high_low$feature_reward_inf == "no", "no", "yes")

moderator_INF_check <- rma(yi2, vi, mods = ~ factor(feature_reward_inf_2), data = SMD_high_low) 

# no = intercept; 1 coeficient yes vs no
moderator_INF_check # significant increase in the effect size for studies that report using instructions vs no instructions


# Check 2: Assuming that not stated are indeed information provided
SMD_high_low$feature_reward_inf_3 <- ifelse(SMD_high_low$feature_reward_inf == "not stated" | SMD_high_low$feature_reward_inf == "yes", "yes", "no")

moderator_INF_check2 <- rma(yi2, vi, mods = ~ factor(feature_reward_inf_3), data = SMD_high_low) 

# no = intercept; 1 coeficient yes vs no
moderator_INF_check2 # significant increase in the effect size for studies that report using instructions vs no instructions

