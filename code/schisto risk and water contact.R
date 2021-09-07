#############################################################
#
#                 Creation of exposure indices 
#       Schistosomiasis risk by water contact activity
#                   Dissertation chapter 3
#                   Code by: Andrea Lund
#                   Created: 14 June 2019
#                Last updated: 6 April 2020
#
#############################################################

library(lme4)
library(glmmTMB)
library(GGally)  

# set working directory 
setwd("C:/Users/andre/Box Sync/Papers/Ch3Prawn-CEA")

# import data
kidsHH <- read.csv("./kidsHH.csv")

###################################################
#
#         create new frequency variables
#

#http://www.cookbook-r.com/Manipulating_data/Recoding_data/
oldvalues <- c(0:5)
wkavgvals <- c(0, 0, 1, 4, 7, 14)
wkminvals <- c(0, 0, 1, 3, 7, 10)
wkmaxvals <- c(0, 0, 1, 5, 7, 21)
moavgvals <- c(0, 2, 4, 14, 30, 90)
mominvals <- c(0, 1, 4, 8, 30, 60)
momaxvals <- c(1,3, 4, 21, 90, 120)

### weekly frequencies by activity
# average
kidsHH$WC08wavg <- wkavgvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09wavg <- wkavgvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10wavg <- wkavgvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11wavg <- wkavgvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12wavg <- wkavgvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13wavg <- wkavgvals[match(kidsHH$WC13, oldvalues)]

# minimum
kidsHH$WC08wmin <- wkminvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09wmin <- wkminvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10wmin <- wkminvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11wmin <- wkminvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12wmin <- wkminvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13wmin <- wkminvals[match(kidsHH$WC13, oldvalues)]

# maximum
kidsHH$WC08wmax <- wkmaxvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09wmax <- wkmaxvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10wmax <- wkmaxvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11wmax <- wkmaxvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12wmax <- wkmaxvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13wmax <- wkmaxvals[match(kidsHH$WC13, oldvalues)]

### monthly frequencies by activity
# average
kidsHH$WC08mavg <- moavgvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09mavg <- moavgvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10mavg <- moavgvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11mavg <- moavgvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12mavg <- moavgvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13mavg <- moavgvals[match(kidsHH$WC13, oldvalues)]

# minimum
kidsHH$WC08mmin <- mominvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09mmin <- mominvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10mmin <- mominvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11mmin <- mominvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12mmin <- mominvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13mmin <- mominvals[match(kidsHH$WC13, oldvalues)]

# maximum
kidsHH$WC08mmax <- momaxvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09mmax <- momaxvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10mmax <- momaxvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11mmax <- momaxvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12mmax <- momaxvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13mmax <- momaxvals[match(kidsHH$WC13, oldvalues)]

###################################################
#
#         create duration and bsa variables
#

kidsBSA = 1.130 # derived from Mosteller 1987 and cited in Sudat et al 2010
duration = c(13.3, 13.3, 2.0, 7.9, 9.3, 23.8) # derived from Sow et al 2011
pctBSA = c(0.358, 0.318, 0.453, 0.225, 0.415, 0.476) # derived from BSA interviews
BSA <- kidsBSA*pctBSA # calculate estimated m2 exposed per activity for children

###################################################
#
#         calculate exposure indices
#

kidsHH <- kidsHH %>% 
  mutate(F_MoAvg = WC01*WC08mavg +  # monthly average frequency
                   WC02*WC09mavg +
                   WC03*WC10mavg + 
                   WC04*WC11mavg + 
                   WC05*WC12mavg +
                   WC06*WC13mavg,
         F_MoMin = WC01*WC08mmin + # monthly minimum frequency
                   WC02*WC09mmin +
                   WC03*WC10mmin + 
                   WC04*WC11mmin + 
                   WC05*WC12mmin +
                   WC06*WC13mmin,
         F_MoMax = WC01*WC08mmax + # monthly maximum frequency
                   WC02*WC09mmax +
                   WC03*WC10mmax + 
                   WC04*WC11mmax + 
                   WC05*WC12mmax +
                   WC06*WC13mmax,
         FD_MoAvg = WC01*WC08mavg*duration[1] +  # monthly average time exposed
                    WC02*WC09mavg*duration[2] +
                    WC03*WC10mavg*duration[3] + 
                    WC04*WC11mavg*duration[4] + 
                    WC05*WC12mavg*duration[5] +
                    WC06*WC13mavg*duration[6],
         FD_MoMin = WC01*WC08mmin*duration[1] + # monthly minimum time exposed
                    WC02*WC09mmin*duration[2] +
                    WC03*WC10mmin*duration[3] + 
                    WC04*WC11mmin*duration[4] + 
                    WC05*WC12mmin*duration[5]+
                    WC06*WC13mmin*duration[6],
         FD_MoMax = WC01*WC08mmax*duration[1] + # monthly maximum time exposed
                    WC02*WC09mmax*duration[2] +
                    WC03*WC10mmax*duration[3] + 
                    WC04*WC11mmax*duration[4] + 
                    WC05*WC12mmax*duration[5] +
                    WC06*WC13mmax*duration[6],
         FDB_MoAvg = WC01*WC08mavg*duration[1]*BSA[1] + # monthly average bsa-adjusted time exposed
                       WC02*WC09mavg*duration[2]*BSA[2] +
                       WC03*WC10mavg*duration[3]*BSA[3] + 
                       WC04*WC11mavg*duration[4]*BSA[4] + 
                       WC05*WC12mavg*duration[5]*BSA[5] +
                       WC06*WC13mavg*duration[6]*BSA[6], 
         FDB_MoMin = WC01*WC08mmin*duration[1]*BSA[1] + # monthly average bsa-adjusted time exposed
                       WC02*WC09mmin*duration[2]*BSA[2] +
                       WC03*WC10mmin*duration[3]*BSA[3] + 
                       WC04*WC11mmin*duration[4]*BSA[4] + 
                       WC05*WC12mmin*duration[5]*BSA[5] +
                       WC06*WC13mmin*duration[6]*BSA[6],
         FDB_MoMax = WC01*WC08mmax*duration[1]*BSA[1] + # monthly average bsa-adjusted time exposed
                       WC02*WC09mmax*duration[2]*BSA[2] +
                       WC03*WC10mmax*duration[3]*BSA[3] + 
                       WC04*WC11mmax*duration[4]*BSA[4] + 
                       WC05*WC12mmax*duration[5]*BSA[5] +
                       WC06*WC13mmax*duration[6]*BSA[6],
         F_WkAvg = WC01*WC08wavg +  # weekly average frequency
                   WC02*WC09wavg +
                   WC03*WC10wavg + 
                   WC04*WC11wavg + 
                   WC05*WC12wavg +
                   WC06*WC13wavg,
         F_WkMin = WC01*WC08wmin + # weekly minimum frequency
                   WC02*WC09wmin +
                   WC03*WC10wmin + 
                   WC04*WC11wmin + 
                   WC05*WC12wmin +
                   WC06*WC13wmin,
         F_WkMax = WC01*WC08wmax + # weekly maximum frequency
                   WC02*WC09wmax +
                   WC03*WC10wmax + 
                   WC04*WC11wmax + 
                   WC05*WC12wmax +
                   WC06*WC13wmax,
         FD_WkAvg = WC01*WC08wavg*duration[1] +  # weekly average time exposed
                    WC02*WC09wavg*duration[2] +
                    WC03*WC10wavg*duration[3] + 
                    WC04*WC11wavg*duration[4] + 
                    WC05*WC12wavg*duration[5] +
                    WC06*WC13wavg*duration[6],
         FD_WkMin = WC01*WC08wmin*duration[1] + # weekly minimum time exposed
                    WC02*WC09wmin*duration[2] +
                    WC03*WC10wmin*duration[3] + 
                    WC04*WC11wmin*duration[4] + 
                    WC05*WC12wmin*duration[5] +
                    WC06*WC13wmin*duration[6],
         FD_WkMax = WC01*WC08wmax*duration[1] + # weekly maximum time exposed
                     WC02*WC09wmax*duration[2] +
                     WC03*WC10wmax*duration[3] + 
                     WC04*WC11wmax*duration[4] + 
                     WC05*WC12wmax*duration[5] +
                     WC06*WC13wmax*duration[6],
         FDB_WkAvg = WC01*WC08wavg*duration[1]*BSA[1] + # weekly average bsa-adjusted time exposed
                       WC02*WC09wavg*duration[2]*BSA[2] +
                       WC03*WC10wavg*duration[3]*BSA[3] + 
                       WC04*WC11wavg*duration[4]*BSA[4] + 
                       WC05*WC12wavg*duration[5]*BSA[5] +
                       WC06*WC13wavg*duration[6]*BSA[6], 
         FDB_WkMin = WC01*WC08wmin*duration[1]*BSA[1] + # weekly average bsa-adjusted time exposed
                       WC02*WC09wmin*duration[2]*BSA[2] +
                       WC03*WC10wmin*duration[3]*BSA[3] + 
                       WC04*WC11wmin*duration[4]*BSA[4] + 
                       WC05*WC12wmin*duration[5]*BSA[5] +
                       WC06*WC13wmin*duration[6]*BSA[6],
         FDB_WkMax = WC01*WC08wmax*duration[1]*BSA[1] + # weekly average bsa-adjusted time exposed
                       WC02*WC09wmax*duration[2]*BSA[2] +
                       WC03*WC10wmax*duration[3]*BSA[3] + 
                       WC04*WC11wmax*duration[4]*BSA[4] + 
                       WC05*WC12wmax*duration[5]*BSA[5] +
                       WC06*WC13wmax*duration[6]*BSA[6],
         # activity-specific weekly frequencies 
         F_laundryWk = WC01*WC08wavg, 
         F_dishesWk = WC02*WC09wavg,
         F_collectWk = WC03*WC10wavg,
         F_irrigWk = WC04*WC11wavg,
         F_livestockWk = WC05*WC12wavg, 
         F_fishingWk = WC06*WC13wavg,
         # activity-specific weekly total duration 
         FD_laundryWk = WC01*WC08wavg*duration[1], 
         FD_dishesWk = WC02*WC09wavg*duration[2],
         FD_collectWk = WC03*WC10wavg*duration[3],
         FD_irrigWk = WC04*WC11wavg*duration[4],
         FD_livestockWk = WC05*WC12wavg*duration[5], 
         FD_fishingWk = WC06*WC13wavg*duration[6],
         # bsa-adjusted activity-specific weekly total duration 
         FDB_laundryWk = WC01*WC08wavg*duration[1]*BSA[1], 
         FDB_dishesWk = WC02*WC09wavg*duration[2]*BSA[2],
         FDB_collectWk = WC03*WC10wavg*duration[3]*BSA[3],
         FDB_irrigWk = WC04*WC11wavg*duration[4]*BSA[4],
         FDB_livestockWk = WC05*WC12wavg*duration[5]*BSA[5], 
         FDB_fishingWk = WC06*WC13wavg*duration[6]*BSA[6],
         # activity-specific monthly frequnecy
         F_laundryMo = WC01*WC08mavg,  
         F_dishesMo = WC02*WC09mavg,
         F_collectMo = WC03*WC10mavg,
         F_irrigMo = WC04*WC11mavg,
         F_livestockMo = WC05*WC12mavg, 
         F_fishingMo = WC06*WC13mavg, 
         # activity-specific monthly total duration
         FD_laundryMo = WC01*WC08mavg*duration[1],  
         FD_dishesMo = WC02*WC09mavg*duration[2],
         FD_collectMo = WC03*WC10mavg*duration[3],
         FD_irrigMo = WC04*WC11mavg*duration[4],
         FD_livestockMo = WC05*WC12mavg*duration[5], 
         FD_fishingMo = WC06*WC13mavg*duration[6], 
         # bsa-adjusted activity-specific monthly total duration 
         FDB_laundryMo = WC01*WC08mavg*duration[1]*BSA[1], 
         FDB_dishesMo = WC02*WC09mavg*duration[2]*BSA[2],
         FDB_collectMo = WC03*WC10mavg*duration[3]*BSA[3],
         FDB_irrigMo = WC04*WC11mavg*duration[4]*BSA[4],
         FDB_livestockMo = WC05*WC12mavg*duration[5]*BSA[5], 
         FDB_fishingMo = WC06*WC13mavg*duration[6]*BSA[6])

check <- kidsHH %>% 
  select(WC01:WC06, WC08mavg:WC13mavg, TE_MoAvg)

hist(kidsHH$FD_irrigWk)

# visualize distributions of exposure indices
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

################################
#
# exploratory analysis
#
################################

# correlation analysis
kidsHH <- kidsHH[!is.na(kidsHH),]
# overall
GGally::ggpairs(kidsHH, columns = c(136:137, 17, 77, 80, 81, 172:173, 220, 223, 226),
                ggplot2::aes(position = "dodge", alpha = 0.5),
                columnLabels = c("Sh med", "Sm med", "age", "wp visits/2wks", "wp visits/wk", 
                                 "field visits/wk", "distWP", "#activities", "FExp", "FDExp", "FDBExp"))
# by sex
GGally::ggpairs(kidsHH, columns = c(136:137, 17, 77, 80, 81, 172:173, 220, 223, 226),
                ggplot2::aes(col = factor(DEM01.x), position = "dodge", alpha = 0.5),
                columnLabels = c("Sh med", "Sm med", "age", "wp visits/2wks", "wp visits/wk", 
                                 "field visits/wk", "distWP", "#activities", "FExp", "FDExp", "FDBExp"))

# columns - Sh_median_18 (136), Sm_median_18 (137),
# DEM04 (17), WC07 (77), ENF02 (80), ENF03 (81), 
# distanceWatPt (172), nActivities (173), 
# F_WkAvg (220), FD_WkAvg (223), FDB_WkAvg (226)

# how many deaths? from what? 
table(indOld2018$DEMD3)
# Accident        Maladie Mort naturelle Mort Naturelle         Noyade 
#       2             59              1             15              1 

# where did others go? why? 
table(indOld2018$DEM21, indOld2018$DEM11)

# look at data on movement to field
table(kidsHH$MO06) # yes n = 236, no n = 847
table(kidsHH$MO06a) # nearby village n = 124, other n = 112
table(kidsHH$MO06ap1) # names of other villages
table(kidsHH$MO06ap2) # names of districts where other villages are located
table(kidsHH$MO06c) # highest frequency: men only (n = 135) followed by men/women (n = 59), 
# then women (n = 18), men/boys (n = 17), boys (n = 5) and everyone (n = 2)
table(kidsHH$MO06d) # highest frequency: several times a week but not every day (n = 134), 
# followed by every day (n = 69), a few times a week (n = 25) and once a week (n = 8)
 
unique(kidsHH$C_ConcessionId) # n = 601

# number of activities
table(kidsHH$nActivities)
#  0   1   2   3   4   5 
# 758 224 141  91   4   2 



################################
#
#  regression analysis
#
################################

### count outcomes
# Sh count predicted activity-specific total duration
shc.actFD <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 +
                         laundryFDwk + dishesFDwk + collectFDwk + 
                         irrigFDwk + livestockFDwk + fishingFDwk + 
                         (1|village.x/C_ConcessionId),
                       data = kidsHH, 
                       family = nbinom2(link = "log"))

summary(shc.actFD)
exp(confint(shc.actFD))

# Sm count predicted by activity-specific frequencies
smc.actFD <- glmmTMB(Sm_median_18 ~ factor(DEM01.x) + DEM04_18 +
                         laundryFDwk + dishesFDwk + collectFDwk + 
                         irrigFDwk + livestockFDwk + fishingFDwk + 
                         (1|village.x/C_ConcessionId),
                       data = kidsHH, 
                       family = nbinom2(link = "log"))

summary(smc.actFD)
exp(confint(shc.actFD))

# Sh count predicted activity-specific frequencies
shc.actfreq <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     laundryFwk + dishesFwk + collectFwk + 
                     irrigFwk + livestockFwk + fishingFwk + 
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))

summary(shc.actfreq)
exp(confint(shc.actfreq))

# Sm count predicted by activity-specific frequencies
smc.actfreq <- glmmTMB(Sm_median_18 ~ factor(DEM01.x) + DEM04_18 +
                         laundryFwk + dishesFwk + collectFwk + 
                         irrigFwk + livestockFwk + fishingFwk + 
                         (1|village.x/C_ConcessionId),
                       data = kidsHH, 
                       family = nbinom2(link = "log"))

summary(smc.actfreq)
exp(confint(smc.actfreq))

# Sh count predicted by presence of any exposure 
shc.any <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     anyActivity + 
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))
summary(shc.any)
exp(confint(shc.any))

# Sm count predicted by presence of any exposure 
smc.any <- glmmTMB(Sm_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     anyActivity + 
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))
summary(smc.any)
exp(confint(smc.any))


# Sh count predicted by activities
shc.act <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     WC01 + WC02 + WC03 + WC04 + WC05 + WC06 + 
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))
summary(shc.act)
exp(confint(shc.act))

# Sm count predicted by activities
smc.act <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     WC01 + WC02 + WC03 + WC04 + WC05 + WC06 +
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))
summary(smc.act)
exp(confint(smc.act))

# Sh count predicted by number of activities
shc.act <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     nActivities + 
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))
summary(shc.act)
exp(confint(shc.act))

# Sm count predicted by number of activities
smc.act <- glmmTMB(Sm_median_18 ~ factor(DEM01.x) + DEM04_18 +
                     nActivities + 
                     (1|village.x/C_ConcessionId),
                   data = kidsHH, 
                   family = nbinom2(link = "log"))
summary(smc.act)
exp(confint(smc.act))

# TE vs BSATE exposure indices
shc.bsate <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 + BSATE_WkAvg +
                       (1|village.x/C_ConcessionId),
                     data = kidsHH, family = nbinom2(link = "log"))
summary(shc.bsate)
exp(confint(shc.bsate))

shc.te <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 + TE_WkAvg + 
                    (1|village.x/C_ConcessionId),
                     data = kidsHH, family = nbinom2(link = "log"))
summary(shc.te)
exp(confint(shc.te))

# survey-derived water point visits
shc.WC07 <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 + WC07 +
                       (1|village.x/C_ConcessionId),
                     data = kidsHH, family = nbinom2(link = "log"))
summary(shc.WC07)
exp(confint(shc.WC07))

shc.ENF02 <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 + ENF02 +
                      (1|village.x/C_ConcessionId),
                    data = kidsHH, family = nbinom2(link = "log"))
summary(shc.ENF02)
exp(confint(shc.ENF02))

shc.ENF03 <- glmmTMB(Sh_median_18 ~ factor(DEM01.x) + DEM04_18 + ENF03 +
                       (1|village.x/C_ConcessionId),
                     data = kidsHH, family = nbinom2(link = "log"))
summary(shc.ENF03)
exp(confint(shc.ENF03))


### presence data
# Sh presence predicted by activities
shp.act <- glm(Sh_presence_18 ~ factor(DEM01.x) + DEM04 +
               WC01 + WC02 + WC03 + WC04 + WC05 + WC06,
             data = kidsHH, family = binomial)
summary(shp.act)
step(shp.act)
cbind(exp(coef(shp.act)), exp(confint(shp.act)))

shp.act.mix <- glmer(Sh_presence_18 ~ factor(DEM01.x) + DEM04 +
               WC01 + WC02 + WC03 + WC04 + WC05 + WC06 +
               (1|village/C_ConcessionId),
             data = kidsHH, family = binomial)
summary(shp.act.mix)


# Sm presence predicted by activities
smp.act <- glm(Sm_presence_18 ~ factor(DEM01.x) + DEM04 +
             WC01 + WC02 + WC03 + WC04 + WC05 + WC06,
           data = kidsHH, family = binomial)
summary(smp.act)

smp.mix <- glmer(Sm_presence_18 ~ factor(DEM01.x) + DEM04 +
                   WC01 + WC02 + WC03 + WC04 + WC05 + WC06 +
                   (1|village/C_ConcessionId),
                 data = kidsHH, family = binomial)
summary(smp.mix)

# coinf
coinf <- glm(coinf_18 ~ factor(DEM01.x) + DEM04 +
             WC01 + WC02 + WC03 + WC04 + WC05 + WC06,
           data = kidsHH, family = binomial)
summary(coinf)

# coinfection presence predicted by activities
coinf.mix <- glmer(coinf_18 ~ factor(DEM01.x) + DEM04 +
                   WC01 + WC02 + WC03 + WC04 + WC05 + WC06 +
                   (1|village/C_ConcessionId),
                 data = kidsHH, family = binomial)
summary(coinf.mix)

# Sh presence predicted by number of activities
shp.nact.mix <- glmer(Sh_presence_18 ~ factor(DEM01.x) + DEM04 +
                       nActivities + 
                       (1|village/C_ConcessionId),
                     data = kidsHH, family = binomial)
summary(shp.nact.mix)

# Sm presence predicted by number of activities
smp.nact.mix <- glmer(Sm_presence_18 ~ factor(DEM01.x) + DEM04 +
                        nActivities + 
                        (1|village/C_ConcessionId),
                      data = kidsHH, family = binomial)
summary(smp.nact.mix)

# Coinfection presence predicted by number of activities
coinf.nact.mix <- glmer(coinf_18 ~ factor(DEM01.x) + DEM04 +
                        nActivities + 
                        (1|village.y/C_ConcessionId),
                      data = kidsHH, family = binomial)
# warnings, max grad = 1.009, very large eigenvalue
summary(coinf.nact.mix)

# Sh presence predicted by nearest water point
shp.nwp.mix <- glmer(Sh_presence_18 ~ factor(DEM01.x) + DEM04 +
                 factor(nearestWatPt) + (1|village.y/C_ConcessionId),
               data = kidsHH, family = binomial)
# failure to converge in 10000 evaluations; unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue
summary(shp.nwp.mix)
# all significant associations between shp and individual water points are negative

# Sm presence predicted by nearest water point
smp.nwp.mix <- glmer(Sm_presence_18 ~ factor(DEM01.x) + DEM04 +
                       factor(nearestWatPt) + (1|village.y/C_ConcessionId),
                     data = kidsHH, family = binomial)
# failure to converge in 10000 evaluations; unable to evaluate scaled gradient; degenerate Hessian with 1 negative eigenvalue
summary(smp.nwp.mix)
# water points associated with significant increase in Sh presence
# Diokhor 2, Gankette, Malla, Merina 2, Ndiawdoune

# infection presence predicted by number of water point visits
# Sh
# in last two weeks
shp.visits1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + DEM04 + WC07 + 
                       (1|village.y/C_ConcessionId),
                     data = kidsHH, family = binomial)
summary(shp.visits1)

# in last week
shp.visits2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + DEM04 + 
                       ENF02 +
                       (1|village.y/C_ConcessionId),
                     data = kidsHH, family = binomial)
summary(shp.visits2)

# Sm
# in last two weeks
smp.visits1 <- glmer(Sm_presence_18 ~ factor(DEM01.x) + DEM04 + WC07 + 
                       (1|village.y/C_ConcessionId),
                     data = kidsHH, family = binomial)
summary(smp.visits1)

# in last week
smp.visits2 <- glmer(Sm_presence_18 ~ factor(DEM01.x) + DEM04 + 
                       ENF02 +
                       (1|village.y/C_ConcessionId),
                     data = kidsHH, family = binomial)
summary(smp.visits2)

# presence predicted by household frequency of water contact 
# Sh
shp.freq <- glmer(Sh_presence_18 ~ factor(DEM01.x) + DEM04 + 
                       factor(WC08) + factor(WC09) + factor(WC10) +
                       factor(WC11) + factor(WC12) + factor(WC12) + 
                       factor(WC13) + factor(WC14) +
                       (1|village.y/C_ConcessionId),
                     data = kidsHH, family = binomial)
# warning max grad = 0.00123672
summary(shp.freq)

# Sm
smp.freq <- glmer(Sm_presence_18 ~ factor(DEM01.x) + DEM04 + 
                    factor(WC08) + factor(WC09) + factor(WC10) +
                    factor(WC11) + factor(WC12) + factor(WC12) + 
                    factor(WC13) + factor(WC14) +
                    (1|village.y/C_ConcessionId),
                  data = kidsHH, family = binomial)
# warnings: failure to converge in 10000 evals, max grad = 0.0266952, large eigenvalue ratio
summary(smp.freq)

# visits and number of activities
table(kidsHH$ENF02)
mean(kidsHH$ENF02, na.rm = T) # 1.26
sd(kidsHH$ENF02, na.rm = T) # 1.29
table(kidsHH$nActivities)
mod <- glmer(ENF02 ~ nActivities + (1|village.y/C_ConcessionId),
      data = kidsHH, family = poisson(link = "log"))

IRR <- exp(fixef(mod)); IRR # = 1.15
# 15% increase in number of visits per week with every additional activity
CI <- exp(confint(mod)); CI # 1.08, 1.23
IRRCI <- cbind(OR, CI); IRRCI 

          