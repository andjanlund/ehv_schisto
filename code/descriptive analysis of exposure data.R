###################################################
#
#           descriptive analysis of 
#       water contact and exposure indices
#               Code by: Andrea Lund
#         Last updated: April 23, 2020

# packages
library(dplyr)

# set working directory 
setwd("C:/Users/andre/Box Sync/Papers/Ch3Prawn-CEA")

# source script that calculates exposure indices
source("./Code/constructing exposure indices.R")

###################################################
#
#   visualize distributions of exposure indices
#
#

# weekly exposure
ggplot(data = kidsHH) + 
  geom_histogram(aes(F_WkAvg), col = "black", fill = 1, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-60
  geom_histogram(aes(FD_WkAvg), col = "black", fill = 2, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-940
  geom_histogram(aes(FDB_WkAvg), col = "black", fill = 3, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-300
  xlab("Exposure") + 
  ylab("Count") +
  #scale_fill_manual(values = c(1:3), labels = c("F", "FD", "FDB")) +
  theme_bw()

# monthly exposure
ggplot(data = kidsHH) + 
  geom_histogram(aes(F_MoAvg), col = "black", fill = 1, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-400
  geom_histogram(aes(FD_MoAvg), col = "black", fill = 2, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-6000
  geom_histogram(aes(FDB_MoAvg), col = "black", fill = 3, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-2000
  xlab("Exposure") + 
  ylab("Count") +
  #scale_fill_manual(values = c(1:3), labels = c("F", "FD", "FDB")) +
  theme_bw()

###################################################
#
#     calculate values for bsaXkids descriptive table

# summarize the frequency of each activity
activityFreq_all <- kidsHH %>% 
  dplyr::summarise(laundryFreqWk_avg = mean(F_laundryWk, na.rm = T), # average frequency
            dishesFreqWk_avg = mean(F_dishesWk, na.rm = T),
            collectFreqWk_avg = mean(F_collectWk, na.rm = T), 
            irrigFreqWk_avg = mean(F_irrigWk, na.rm = T),
            animalsFreqWk_avg = mean(F_livestockWk, na.rm = T),
            fishingFreqWk_avg = mean(F_fishingWk, na.rm = T),
            bathingFreqWk_avg = mean(F_bathingWk, na.rm = T),
            # minimum weekly frequency
            laundryFreqWk_min = mean(F_laundryWk_min, na.rm = T),
            dishesFreqWk_min = mean(F_dishesWk_min, na.rm = T),
            collectFreqWk_min = mean(F_collectWk_min, na.rm = T), 
            irrigFreqWk_min = mean(F_irrigWk_min, na.rm = T),
            animalsFreqWk_min = mean(F_livestockWk_min, na.rm = T),
            fishingFreqWk_min = mean(F_fishingWk_min, na.rm = T),
            bathingFreqWk_min = mean(F_fishingWk_min, na.rm = T),
            # maximum weekly frequency
            laundryFreqWk_max = mean(F_laundryWk_max, na.rm = T),
            dishesFreqWk_max = mean(F_dishesWk_max, na.rm = T),
            collectFreqWk_max = mean(F_collectWk_max, na.rm = T), 
            irrigFreqWk_max = mean(F_irrigWk_max, na.rm = T),
            animalsFreqWk_max = mean(F_livestockWk_max, na.rm = T),
            fishingFreqWk_max = mean(F_fishingWk_max, na.rm = T),
            bathingFreqWk_max = mean(F_bathingWk_max, na.rm = T)) %>% 
  gather() %>% 
  separate(key, into = c("activity", "measure"), sep = "_") %>% 
  pivot_wider(names_from = "measure", values_from = "value") %>% 
  dplyr::mutate(DEM01.x = NA)

activityFreq_sex <- kidsHH %>% 
  dplyr::group_by(DEM01.x) %>% 
  dplyr::summarise(laundryFreqWk_avg = mean(F_laundryWk, na.rm = T), # average frequency
            dishesFreqWk_avg = mean(F_dishesWk, na.rm = T),
            collectFreqWk_avg = mean(F_collectWk, na.rm = T), 
            irrigFreqWk_avg = mean(F_irrigWk, na.rm = T),
            animalsFreqWk_avg = mean(F_livestockWk, na.rm = T),
            fishingFreqWk_avg = mean(F_fishingWk, na.rm = T),
            bathingFreqWk_avg = mean(F_bathingWk, na.rm = T),
            # minimum weekly frequency
            laundryFreqWk_min = mean(F_laundryWk_min, na.rm = T),
            dishesFreqWk_min = mean(F_dishesWk_min, na.rm = T),
            collectFreqWk_min = mean(F_collectWk_min, na.rm = T), 
            irrigFreqWk_min = mean(F_irrigWk_min, na.rm = T),
            animalsFreqWk_min = mean(F_livestockWk_min, na.rm = T),
            fishingFreqWk_min = mean(F_fishingWk_min, na.rm = T),
            bathingFreqWk_min = mean(F_fishingWk_min, na.rm = T),
            # maximum weekly frequency
            laundryFreqWk_max = mean(F_laundryWk_max, na.rm = T),
            dishesFreqWk_max = mean(F_dishesWk_max, na.rm = T),
            collectFreqWk_max = mean(F_collectWk_max, na.rm = T), 
            irrigFreqWk_max = mean(F_irrigWk_max, na.rm = T),
            animalsFreqWk_max = mean(F_livestockWk_max, na.rm = T),
            fishingFreqWk_max = mean(F_fishingWk_max, na.rm = T),
            bathingFreqWk_max = mean(F_bathingWk_max, na.rm = T)) %>% 
  #filter(!is.na(DEM01.x)) %>% 
  pivot_longer(cols = -DEM01.x, names_to = c("activity", "measure"), names_sep = "_") %>% 
  pivot_wider(names_from = "measure", values_from = "value") 

# export descriptive frequency tables, overall and by sex
activityFreq_combo <- rbind(activityFreq_sex, activityFreq_all)
activityFreq_combo$DEM01.x <- ifelse(is.na(activityFreq_combo$DEM01.x), "all", activityFreq_combo$DEM01.x)
write.csv(activityFreq_combo, file = "./Figures & Tables/activityFrequencies_wBathing.csv")

activityFreqPlot <- kidsHH %>%
  dplyr::select(village, DEM01.x, OD10p, concession, F_laundryWk:F_fishingWk) %>% 
  pivot_longer(cols = starts_with("F_"),
               names_to = "activity", 
               values_to = "frequency")

ggplot(data = activityFreqPlot) + 
  geom_histogram(aes(x = frequency, fill = activity),
                 position = "dodge", bins = 15, na.rm = T) + 
  xlab("Weekly Frequency") + 
  ylab("Count") +
  facet_grid(DEM01.x ~ . ) +
  theme_bw()

# how many children are closest to a water point where snail sampling wasn't performed?
table(kidsHH$nearestWatPt)
# Bountou Yaye Tabara - 28
# Canal Mayoro - 29
# Gae 1 - 3
# Gae 2 - 3
# Mboltogne - 1
# Medina 1 - 7
# Medina 2 - 3
# Medina 3 - 19
# Pont Ndiawdoune - 2
# Pont Ndiol - 16
# Tamakh - 3
# --------------------
# Total - 118 (9.7%)

table(kidsHH$nearestSnail)

# do frequencies of activities differ by village?
ggplot(data = kidsHH) + 
  geom_histogram(aes(F_WkAvg), col = "black", fill = 1, bins = 100,
                 position = "dodge", alpha = 0.2, na.rm = T) + # range 0-60
  facet_grid(LakeYN ~ .) +
  xlab("Exposure") + 
  ylab("Count") +
  #scale_fill_manual(values = c(1:3), labels = c("F", "FD", "FDB")) +
  theme_bw()

# how are environmental variables distributed? 
hist(kidsHH$pctHabitat_peak, breaks = 40)
hist(kidsHH$habitatArea_peak, breaks = 40)
hist(kidsHH$habitatAreaPeak_d, breaks = 40)
hist(kidsHH$habitatArea_yr, breaks = 40)
hist(kidsHH$habitatAreaYr_d, breaks = 40)

# summarize numeric data for tables 1 and 2
demosNum <- kidsHH %>% 
  summarise(age_mean = mean(DEM04), age_sd = sd(DEM04),
            meanDistance = mean(distanceSnail), 
            sdDistance = sd(distanceSnail),
            meanHabitat = mean(habitatArea_peak),
            sdHabitat = sd(habitatArea_peak),
            visitsMean = mean(ENF02),
            visitsSD = sd(ENF02),
            FMean = mean(F_WkAvg),
            FSD = sd(F_WkAvg),
            FDMean = mean(FD_WkAvg),
            FDSD = sd(FD_WkAvg),
            FDBMean = mean(FDB_WkAvg),
            FDBSD = sd(FDB_WkAvg))

# range for peak habitat
min(kidsHH$habitatArea_peak)
max(kidsHH$habitatArea_peak)

freqMeans <- kidsHH %>% 
  summarise(laundryFreqWk_avg = mean(F_laundryWk, na.rm = T), # average frequency
            laundry_SD = sd(F_laundryWk),
            dishesFreqWk_avg = mean(F_dishesWk, na.rm = T),
            dishes_SD = sd(F_dishesWk),
            collectFreqWk_avg = mean(F_collectWk, na.rm = T), 
            collect_SD = sd(F_collectWk),
            irrigFreqWk_avg = mean(F_irrigWk, na.rm = T),
            irrig_SD = sd(F_irrigWk, na.rm = T),
            animalsFreqWk_avg = mean(F_livestockWk, na.rm = T),
            animals_SD = sd(F_livestockWk, na.rm = T),
            fishingFreqWk_avg = mean(F_fishingWk, na.rm = T),
            fishing_SD = sd(F_fishingWk, na.rm = T),
            bathingFreqWk_avg = mean(F_bathingWk, na.rm = T),
            bathing_SD = sd(F_bathingWk))

demosNumSex <- kidsHH %>% 
  group_by(DEM01.x) %>% 
  summarise(age_mean = mean(DEM04), age_sd = sd(DEM04),
            meanDistance = mean(distanceSnail), 
            sdDistance = sd(distanceSnail),
            meanHabitat = mean(habitatArea_peak),
            sdHabitat = sd(habitatArea_peak),
            visitsMean = mean(ENF02),
            visitsSD = sd(ENF02),
            FMean = mean(F_WkAvg),
            FSD = sd(F_WkAvg),
            FDMean = mean(FD_WkAvg),
            FDSD = sd(FD_WkAvg),
            FDBMean = mean(FDB_WkAvg),
            FDBSD = sd(FDB_WkAvg))

habitatRangeSex <- kidsHH %>% 
  group_by(DEM01.x) %>% 
  summarise(minHabitat = min(habitatArea_peak),
            maxHabitat = max(habitatArea_peak))

freqMeansSex <- kidsHH %>% 
  group_by(DEM01.x) %>% 
  summarise(laundryFreqWk_avg = mean(F_laundryWk, na.rm = T), # average frequency
            laundry_SD = sd(F_laundryWk),
            dishesFreqWk_avg = mean(F_dishesWk, na.rm = T),
            dishes_SD = sd(F_dishesWk),
            collectFreqWk_avg = mean(F_collectWk, na.rm = T), 
            collect_SD = sd(F_collectWk),
            irrigFreqWk_avg = mean(F_irrigWk, na.rm = T),
            irrig_SD = sd(F_irrigWk, na.rm = T),
            animalsFreqWk_avg = mean(F_livestockWk, na.rm = T),
            animals_SD = sd(F_livestockWk, na.rm = T),
            fishingFreqWk_avg = mean(F_fishingWk, na.rm = T),
            fishing_SD = sd(F_fishingWk, na.rm = T),
            bathingFreqWk_avg = mean(F_bathingWk, na.rm = T),
            bathing_SD = sd(F_bathingWk))

# summarize categorical data for tables 1 and 2
table(kidsHH$ENF02)
table(kidsHH$surface) 
table(kidsHH$surface, kidsHH$DEM01.x) 

table(kidsHH$sanitation)
table(kidsHH$sanitation, kidsHH$DEM01.x)

table(kidsHH$wealthQuintile1)

# summarize infection data
table(kidsHH$Sh_presence_18)
table(kidsHH$Sh_presence_18, kidsHH$DEM01.x)

gMean_all <- kidsHH %>% 
  mutate(Sh_trans = log(Sh_mean_18 + 1)) %>% 
  summarise(Sh_gMean = (exp(mean(Sh_trans, na.rm = T))-1),
            Sh_sd  = sd(Sh_mean_18, na.rm = T)) 

gMean_sex <- kidsHH %>% 
  group_by(DEM01.x) %>% 
  mutate(Sh_trans = log(Sh_mean_18 + 1)) %>% 
  summarise(Sh_gMean = (exp(mean(Sh_trans, na.rm = T))-1),
            Sh_sd  = sd(Sh_mean_18, na.rm = T)) 

# how many children respond for themselves with raw frequency variable?
table(kidsHH$ENF00, kidsHH$ENF02)

# do survey data support that households use surface water regardless of SES?
table(kidsHH$surface, kidsHH$wealthQuintile1)
kruskal.test(surface ~ wealthQuintile1, data = kidsHH)
# H0: distribution is the same in each group; Chisq = 19.117, df = 4, p < 0.001

summary(kidsHH$habitatArea_peak)

# correlation for all variables considered in models
library(corrplot)
corVars <- select(kidsHH, Sh_presence_18, Sh_median_18,
                          ENF02, F_WkAvgC, FD_WkAvgC, FDB_WkAvgC,
                          habitatArea_peakC, habitatArea_yrC,
                          habitatAreaPeak_dC, habitatAreaYr_dC,
                          habitatAreaV_peakC, habitatAreaV_yrC,
                          surface, surfaceYN, sanitation, 
                          privateSan, wealthQuintile1)
colnames(corVars) <- c("Sh_presence", "Sh_intensity",
                       "F(raw)", "F(sum)", "FD", "FDB",
                       "areaPeak", "areaYear", 
                       "areaPeakD", "areaYearD",
                       "areaPeakV", "areaYearV",
                       "surface", "surfaceYN",
                       "sanitation", "privateSan", "SES")
corMatrix <- cor(corVars)
corrplot(corMatrix, method = "color", 
         type = "lower",
         tl.col = "black", tl.srt = 45)

corrplot.mixed(corMatrix)
