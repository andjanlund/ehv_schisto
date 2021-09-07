###################################################
#
#             regression analysis of
#        exposure - hazard - vulnerability
#       (excluding vegetation removal sites)
#               Code by: Andrea Lund
#         Last updated: September 7, 2021

# install packages 
library(tidyverse) 
library(lme4) # for mixed effects logistic regression
library(TMB)
library(glmmTMB) # for mixed effects negative binomial regression
library(AICcmodavg) # for model averaging of candidate set
library(car)
library(performance) # for binned residuals of binomial models
library(scales) # for plotting with logarithmic scales
library(GGally)
library(here)

# source script that calculates exposure indices
source(here("code/constructing exposure indices.R"))

# create scaled versions of all exposure indices
kidsHH <- kidsHH %>% 
  mutate(Sh_medianR_18 = round(Sh_median_18, 0),
         F_WkAvgC = scale(F_WkAvg),
         FD_WkAvgC = scale(FD_WkAvg),
         FDB_WkAvgC = scale(FDB_WkAvg),
         habitatArea_peakC = scale(habitatArea_peak),
         habitatAreaPeak_dC = scale(habitatAreaPeak_d),
         habitatArea_yrC = scale(habitatArea_yr),
         habitatAreaYr_dC = scale(habitatAreaYr_d),
         habitatAreaV_peakC = scale(habitatAreaV_peak),
         habitatAreaV_yrC = scale(habitatAreaV_yr)) %>% 
  select(-C_ConcessionId.y) %>% 
  rename("C_ConcessionId" = "C_ConcessionId.x")

# remove villages that underwent vegetation removal in 2017: 
#                           Guidick (n = 74), Lampsar (n = 71) and Mbane (n = 64)
kidsHH <- kidsHH %>% 
  filter(village %in% c("DF", "DT", "FS", "GG", "MA", "MB", "MD", "MG", "MO", "MT", "NE", "NM", "ST"))
  # n = 821 observations

# only 1030 observations for ENF02, so subset entire data set to make data set same size for all models
kidsHH <- kidsHH[!is.na(kidsHH$ENF02) & !is.na(kidsHH$F_WkAvg) & !is.na(kidsHH$Sh_medianR_18),]

########################################
#
#       S. HAEMATOBIUM PRESENCE
#
######## CONTACT-ONLY MODELS
# exposure = raw frequency
pexp1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(ENF02) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp1) #  849.7019
vif(pexp1)

# exposure = frequency only
pexp2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 F_WkAvgC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp2) # 851.1602
vif(pexp2)

# exposure = frequency + duration
pexp3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 FD_WkAvgC +
                 LakeYN + scale(Pop) + 
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp3) # 851.2009

# exposure = bsa-adjusted frequency + duration
pexp4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 FDB_WkAvgC +
                 LakeYN + scale(Pop) + 
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp4) # 851.1798

########     HAZARD-ONLY MODELS
# habitat area during transmission peak
phaz1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   habitatArea_peakC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))

AIC(phaz1) # 849.8695

# habitat area during entire year
phaz2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatArea_yrC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz2) # 848.3805

# distance-weighted habitat area during transmission peak
phaz3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaPeak_dC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz3) # 845.0784

# distance-weighted habitat area during entire year
phaz4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaYr_dC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz4) # 841.9315

# village habitat area during transmission peak
phaz5 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaV_peakC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz5) # 851.0572

# village habitat area over entire year
phaz6 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaV_yrC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz6) # 849.7349

######## VULNERABILITY-ONLY MODELS
# vulnerability = surface water dependence
pvul1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(surface) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

# max grad = 0.00211662
AIC(pvul1) # 843.296

# vulnerability = binary surface water dependence
pvul2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 surfaceYN +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pvul2) # 846.9709

# vulnerability = binary private vs shared/no sanitation
pvul3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 privateSan +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pvul3) # 847.404

# vulnerability = categorical sanitation (private vs shared vs no)
pvul4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(sanitation) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))


AIC(pvul4) # 848.2637

# vulnerability = asset index
pvul5 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(wealthQuintile1) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pvul5) # 852.7738

########     CONTACT-HAZARD MODELS
# simple exposure + distance-adjusted habitat area over year
peh1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(ENF02) + habitatAreaYr_dC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))
# max|grad| = 0.00139315
AIC(peh1) # 842.9508
vif(peh1)

# exposure frequency + distance-adjusted habitat area over year
peh2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +  
                   F_WkAvgC + habitatAreaYr_dC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))

AIC(peh2) #  843.5121
vif(peh2)

# exposure duration + distance-adjusted habitat area over year
peh3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +  
                  FD_WkAvgC + habitatAreaYr_dC +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = binomial(link = "logit"))

AIC(peh3) # 843.6321
vif(peh3)

# village  habitat area over entire year
peh4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +  
                  FDB_WkAvgC + habitatAreaYr_dC +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = binomial(link = "logit"))

AIC(peh4) #  843.8457
vif(peh4)


########     CONTACT-VULNERABILITY  MODELS
pev1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) + 
                  factor(ENF02) + factor(surface) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0037424
AIC(pev1) # 845.2792
vif(pev1)

pev2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                F_WkAvgC + factor(surface) +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
# max gad = 0.00290092
AIC(pev2) # 844.6025
vif(pev2)

pev3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                FD_WkAvgC + factor(surface) +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00119441
AIC(pev3) # 844.7865
vif(pev3)

pev4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                FDB_WkAvgC + factor(surface) +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
AIC(pev4) # 845.162
vif(pev4)

########     HAZARD-VULNERABILITY  MODELS
phv1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaYr_dC + factor(surface) +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
AIC(phv1) # 837.3656
vif(phv1)


########    (GLOBAL) CONTACT-HAZARD-VULNERABILITY MODELS 

# dependence on surface water
pehv1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) + 
                     factor(ENF02) + habitatAreaYr_dC + factor(surface) +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00429655
summary(pehv1)
AIC(pehv1) # 839.2104
vif(pehv1) # no concerning values here

# binary surface water dependence
pehv2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) + 
                     F_WkAvgC + habitatAreaYr_dC + factor(surface) +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00114572
AIC(pehv2) # 838.3239
vif(pehv2)


pehv3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) + 
                     FD_WkAvgC + habitatAreaYr_dC + factor(surface) +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
AIC(pehv3) # 838.574
vif(pehv3)

pehv4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   FDB_WkAvgC + habitatAreaYr_dC + factor(surface) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
AIC(pehv4) # 838.8878
vif(pehv4)


##################################################################
#
#           compute model averages for presence models

# create a list of all presence models
presenceModels <- list(pexp1, pexp2, pexp3, pexp4,
                       phaz1, phaz2, phaz3, phaz4, phaz5, phaz6,
                       pvul1, pvul2, pvul3, pvul4, pvul5,
                       peh1, peh2, peh3, peh4,
                       pev1, pev2, pev3, pev4,
                       phv1,
                       pehv1, pehv2, pehv3, pehv4)

modNamesP <- c('pexp1', 'pexp2', 'pexp3', 'pexp4',
                 'phaz1', 'phaz2', 'phaz3', 'phaz4', 'phaz5', 'phaz6',
                 'pvul1', 'pvul2', 'pvul3', 'pvul4', 'pvul5',
                 'peh1', 'peh2', 'peh3', 'peh4',
                 'pev1', 'pev2', 'pev3', 'pev4',
                 'phv1',
                 'pehv1', 'pehv2', 'pehv3', 'pehv4')

# construct model selection table 
presenceAIC <- aictab(presenceModels, modnames = modNamesP, second.ord = F, sort = F) # K is 3 greater than what I have in the spreadsheet - is this random effects?
write.csv(presenceAIC, here("output/presenceAIC28_factor_vegRem.csv"))
#presenceBIC <- bictab(presenceModels, modnames = modNamesP, second.ord = F, sort = F)
#write.csv(presenceBIC, "./Figures & Tables/presenceBIC22_factor.csv")
# determine confidence set
confset(presenceModels, modNamesP, method = 'raw', second.ord = F, level = 0.95) # 8 models
confSetP28 <- confset(presenceModels, modNamesP, method = 'ordinal', second.ord = F) 
  # substantial weight - 5 models, some weight - 5, little weight - 5 models, no weight - 7 models
write.csv(confSetP28$substantial, here("output/confsetPresence28_factor_vegRem.csv"))
confset(presenceModels, modNamesP, method = 'ratio', second.ord = F, delta = 10) # 15 models with delta = 10 cutoff

# compute model-averaged estimates of individual parameters
fraw1 <- modavg(presenceModels, parm = "factor(ENF02)1", modnames = modNamesP, second.ord = F); fraw1
fraw2 <- modavg(presenceModels, parm = "factor(ENF02)2", modnames = modNamesP, second.ord = F); fraw2
fraw3 <- modavg(presenceModels, parm = "factor(ENF02)3", modnames = modNamesP, second.ord = F); fraw3
fraw4 <- modavg(presenceModels, parm = "factor(ENF02)4", modnames = modNamesP, second.ord = F); fraw4
fsum <- modavg(presenceModels, parm = "F_WkAvgC", modnames = modNamesP, second.ord = F); fsum
fd <- modavg(presenceModels, parm = "FD_WkAvgC", modnames = modNamesP, second.ord = F); fd
fdb <- modavg(presenceModels, parm = "FDB_WkAvgC", modnames = modNamesP, second.ord = F); fdb
areaPeak_d <- modavg(presenceModels, parm = "habitatAreaPeak_dC", modnames = modNamesP, second.ord = F); areaPeak_d
areaYear_d <- modavg(presenceModels, parm = "habitatAreaYr_dC", modnames = modNamesP, second.ord = F); areaYear_d
areaPeak <- modavg(presenceModels, parm = "habitatArea_peakC", modnames = modNamesP, second.ord = F); areaPeak
areaYear <- modavg(presenceModels, parm = "habitatArea_yrC", modnames = modNamesP, second.ord = F); areaYear
areaPeakV <- modavg(presenceModels, parm = "habitatAreaV_peakC", modnames = modNamesP, second.ord = F); areaPeakV#areaYearV <- modavg(presenceModels, parm = "habitatAreaV_yrC", modnames = modNamesP, second.ord = F); areaYearV
surface1 <- modavg(presenceModels, parm = "factor(surface)1", modnames = modNamesP, second.ord = F); surface1
surface2 <- modavg(presenceModels, parm = "factor(surface)2", modnames = modNamesP, second.ord = F); surface2
surfaceYN <- modavg(presenceModels, parm = "surfaceYN", modnames = modNamesP, second.ord = F); surfaceYN
privateSan <- modavg(presenceModels, parm = "privateSan", modnames = modNamesP, second.ord = F); privateSan
sanitation1 <- modavg(presenceModels, parm = "factor(sanitation)1", modnames = modNamesP, second.ord = F); sanitation1
sanitation2 <- modavg(presenceModels, parm = "factor(sanitation)2", modnames = modNamesP, second.ord = F); sanitation2
sex <- modavg(presenceModels, parm = "factor(DEM01.x)2", modnames = modNamesP, second.ord = F); sex
age <- modavg(presenceModels, parm = "scale(DEM04)", modnames = modNamesP, second.ord = F); age
vPopn <- modavg(presenceModels, parm = "scale(Pop)", modnames = modNamesP, second.ord = F); vPopn
vLocn <- modavg(presenceModels, parm = "LakeYN", modnames = modNamesP, second.ord = F); vLocn

kidsHH$surfaceF <- as.factor(kidsHH$surface)


########################################
#
#       S. HAEMATOBIUM INTENSITY
#
######## CONTACT-ONLY MODELS
# exposure = raw frequency 
cexp1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(ENF02) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cexp1) # 5494.29


# exposure = frequency only
cexp2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   F_WkAvgC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cexp2) # 5511.231

# exposure = frequency + duration
cexp3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                   LakeYN + scale(Pop) + 
                   FD_WkAvgC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))

AIC(cexp3) # 5510.809

# exposure = bsa-adjusted frequency + duration
cexp4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   FDB_WkAvgC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cexp4) # 5510.88

######## HAZARD-ONLY MODELS
chaz1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatArea_peakC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz1) # 5508.681

chaz2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatArea_yrC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz2) # 5509.936

chaz3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaPeak_dC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz3) # 5509.578

chaz4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaYr_dC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz4) # 5508.982

chaz5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaV_peakC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz5) # 5511.311

chaz6 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaV_yrC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz6) # 5510.508

######## VULNERABILITY-ONLY MODELS
cvul1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(surface) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul1) # 5508.892

cvul2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   surfaceYN + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul2) # 5506.95

cvul3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   privateSan + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul3) # 5506.796

cvul4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(sanitation) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul4) # 5508.733

cvul5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(wealthQuintile1) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul5) # 5513.193

########     EXPOSURE-HAZARD MODELS
ceh1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_peakC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))
AIC(ceh1) # 5494.41

ceh2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_yrC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh2) # 5495.303

ceh3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + habitatAreaPeak_dC +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh3) # 5494.859

ceh4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + habitatAreaYr_dC +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh4) #  5494.455

ceh5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  factor(ENF02) + habitatAreaV_yrC +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh5) #  5495.498

########     CONTACT-VULNERABILITY MODELS
cev1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + surfaceYN +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))
AIC(cev1) # 5495.212

cev2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  factor(ENF02) + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(cev2) # 5494.491

cev3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  factor(ENF02) + factor(sanitation) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(cev3) # 5496.231

########     HAZARD-VULNERABILITY MODELS
chv1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_peakC + surfaceYN +
                 (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv1) # 5505.471

chv2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_yrC + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv2) # 5506.894

chv3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaYr_dC + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv3) # 5507.504

chv4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaPeak_dC + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv4) # 5507.947

chv5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_peakC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv5) # 5506.675

chv6 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_yrC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv6) # 5507.659

chv7 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaYr_dC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv7) # 5506.078

chv8 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaPeak_dC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv8) # 5506.648

chv9 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_peakC + factor(sanitation) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv9) # 5508.584

chv10 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_yrC + factor(sanitation) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv10) # 5509.576

chv11 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaYr_dC + factor(sanitation) +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv11) # 5508.008

chv12 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaPeak_dC + factor(sanitation) +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv12) # 5508.576


########    (GLOBAL) MODELS OF EXPOSURE-HAZARD-VULNERABILITY

# surface water dependence
cehv1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                       LakeYN + scale(Pop) +
                       factor(ENF02) + habitatArea_peakC + surfaceYN +
                       (1|village/C_ConcessionId), 
                     data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv1) # 5494.975
check_collinearity(cehv1)

cehv2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                       LakeYN + scale(Pop) +
                       factor(ENF02) + habitatArea_yrC + surfaceYN +
                       (1|village/C_ConcessionId), 
                     data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv2) # 5495.948
check_collinearity(cehv2)


cehv3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                      LakeYN + scale(Pop) +
                      factor(ENF02) + habitatAreaYr_dC + surfaceYN +
                      (1|village/C_ConcessionId), 
                    data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv3) # 5495.81
check_collinearity(cehv3)

cehv4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatAreaPeak_dC + surfaceYN +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv4) # 5496.161
check_collinearity(cehv4)

cehv5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_peakC + privateSan +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv5) # 5494.949
check_collinearity(cehv5)

cehv6 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_yrC + privateSan +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv6) # 5495.702
check_collinearity(cehv6)

cehv7 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   factor(ENF02) + habitatAreaYr_dC + privateSan +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv7) # 5494.456
check_collinearity(cehv7)

cehv8 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   factor(ENF02) + habitatAreaPeak_dC + privateSan +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv8) # 5494.864
check_collinearity(cehv8)

cehv9 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   factor(ENF02) + habitatArea_peakC + factor(sanitation) +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv9) # 5496.641
check_collinearity(cehv9)

cehv10 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   factor(ENF02) + habitatArea_yrC + factor(sanitation) +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv10) # 5497.408
check_collinearity(cehv10)

cehv11 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + habitatAreaYr_dC + factor(sanitation) +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv11) # 5496.187
check_collinearity(cehv11)

cehv12 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + habitatAreaPeak_dC + factor(sanitation) +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv12) # 5496.589
check_collinearity(cehv12)

##################################################################
#
#           compute model averages for intensity models

# create a list of all presence models
intensityModels <- list(cexp1, cexp2, cexp3, cexp4,
                       chaz1, chaz2, chaz3, chaz4, chaz5, chaz6,
                       cvul1, cvul2, cvul3, cvul4, cvul5,
                       ceh1, ceh2, ceh3, ceh4, ceh5,
                       cev1, cev2, cev3,
                       chv1, chv2, chv3, chv4, chv5, chv6,
                       chv7, chv8, chv9, chv10, chv11, chv12,
                       cehv1, cehv2, cehv3, cehv4, cehv5, cehv6,
                       cehv7, cehv8, cehv9, cehv10, cehv11, cehv12)

modNamesI <- c('cexp1', 'cexp2', 'cexp3', 'cexp4',
              'chaz1', 'chaz2', 'chaz3', 'chaz4', 'chaz5', 'chaz6',
              'cvul1', 'cvul2', 'cvul3', 'cvul4', 'cvul5',
              'ceh1', 'ceh2', 'ceh3', 'ceh4', 'ceh5',
              'cev1', 'cev2', 'cev3',
              'chv1', 'chv2', 'chv3', 'chv4', 'chv5', 'chv6',
              'chv7', 'chv8', 'chv9', 'chv10', 'chv11', 'chv12',
              'cehv1', 'cehv2', 'cehv3', 'cehv4', 'cehv5', 'cehv6',
              'cehv7', 'cehv8', 'cehv9', 'cehv10', 'cehv11', 'cehv12')

# construct model selection table 
intensityAIC <- aictab(intensityModels, modnames = modNamesI, second.ord = F, sort = F) # K is 3 greater than what I have in the spreadsheet - is this random effects?
write.csv(intensityAIC, here("output/intensityAIC47_vegRem.csv"))
#intensityBIC <- bictab(intensityModels, modnames = modNamesI, second.ord = F, sort = F)
#write.csv(intensityBIC, "./Figures & Tables/intensityBIC40_factor.csv")
# determine confidence set
confset(intensityModels, modNamesI, method = 'raw', second.ord = F, level = 0.95) # 13 models
intConfSet <- confset(intensityModels, modNamesI, method = 'ordinal', second.ord = F) 
write.csv(intConfSet$substantial, here("output/confsetIntensity47_vegRem.csv"))
# substantial weight - 11 models, some weight - 3, little weight - 6 models, no weight - 20 models
confset(intensityModels, modNamesI, method = 'ratio', second.ord = F, delta = 10) # 20 models with delta = 10 cutoff

# diagnostics on models with substantial weight

# compute model-averaged estimates of individual parameters
fraw1 <- modavg(intensityModels, parm = "factor(ENF02)1", modnames = modNamesI, second.ord = F); fraw1
fraw2 <- modavg(intensityModels, parm = "factor(ENF02)2", modnames = modNamesI, second.ord = F); fraw2
fraw3 <- modavg(intensityModels, parm = "factor(ENF02)3", modnames = modNamesI, second.ord = F); fraw3
fraw4 <- modavg(intensityModels, parm = "factor(ENF02)4", modnames = modNamesI, second.ord = F); fraw4
fsum <- modavg(intensityModels, parm = "F_WkAvgC", modnames = modNamesI, second.ord = F); fsum
fd <- modavg(intensityModels, parm = "FD_WkAvgC", modnames = modNamesI, second.ord = F); fd
fdb <- modavg(intensityModels, parm = "FDB_WkAvgC", modnames = modNamesI, second.ord = F); fdb
areaPeak <- modavg(intensityModels, parm = "habitatArea_peakC", modnames = modNamesI, second.ord = F); areaPeak
areaYear <- modavg(intensityModels, parm = "habitatArea_yrC", modnames = modNamesI, second.ord = F); areaYear
areaPeak_d <- modavg(intensityModels, parm = "habitatAreaPeak_dC", modnames = modNamesI, second.ord = F); areaPeak_d
areaYear_d <- modavg(intensityModels, parm = "habitatAreaYr_dC", modnames = modNamesI, second.ord = F); areaYear_d
areaYearV <- modavg(intensityModels, parm = "habitatAreaV_yrC", modnames = modNamesI, second.ord = F); areaYearV
areaPeakV <- modavg(intensityModels, parm = "habitatAreaV_peakC", modnames = modNamesI, second.ord = F); areaPeakV
surface1 <- modavg(intensityModels, parm = "factor(surface)1", modnames = modNamesI, second.ord = F); surface1
surface2 <- modavg(intensityModels, parm = "factor(surface)2", modnames = modNamesI, second.ord = F); surface2
surfaceYN <- modavg(intensityModels, parm = "surfaceYN", modnames = modNamesI, second.ord = F); surfaceYN
privateSan <- modavg(intensityModels, parm = "privateSan", modnames = modNamesI, second.ord = F); privateSan
sanitation1 <- modavg(intensityModels, parm = "factor(sanitation)1", modnames = modNamesI, second.ord = F); sanitation1
sanitation2 <- modavg(intensityModels, parm = "factor(sanitation)2", modnames = modNamesI, second.ord = F); sanitation2
sex <- modavg(intensityModels, parm = "factor(DEM01.x)2", modnames = modNamesI, second.ord = F); sex
age <- modavg(intensityModels, parm = "scale(DEM04)", modnames = modNamesI, second.ord = F); age
vPopn <- modavg(intensityModels, parm = "scale(Pop)", modnames = modNamesI, second.ord = F); vPopn
vLocn <- modavg(intensityModels, parm = "LakeYN", modnames = modNamesI, second.ord = F); vLocn


##################################################################
#
#           visualize effect sizes and precision for 
#           E, H and V components of risk

# import model-averaged estimates
allEst <- read.csv(here("output/allEstimates_vegRem.csv"))
allEst$Type <- factor(allEst$Type, 
                      levels = c("Vulnerability", "Hazard", "Exposure", "Demographics"))
allEst$est <- factor(allEst$est, 
                      levels = c("Presence", "Intensity"))
allEst <- filter(allEst, !is.na(nmod))

# generate forest plot for all Sh presence estimates
avgPlot <- ggplot(allEst,
       aes(x = Variable, y = OR, ymin = orLL, ymax = orUL)) +
  geom_pointrange(aes(col = Type, pch = Type),
                  position = position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = orLL, ymax = orUL, col = Type, pch = Type), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab("Variable") +
  ylab("Odds Ratio (Presence) / Rate Ratio (Intensity)") +
    scale_y_continuous(trans = 'log2', limits = c(0.03, 46))+ 
  #                   breaks = trans_breaks('log2', function(x) 2^x)) +
  facet_grid(est ~ .) +
  coord_flip() +
  #scale_linetype_manual(values = c("longdash", "dotdash", "solid")) +
  scale_colour_manual(values = c("steelblue4", "turquoise4", "black", "darkorange2")) +
  scale_shape_manual(values = c(15:18)) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("Figure 2.tiff", plot = avgPlot, device = "tiff",
       path = here("output"), dpi = 300, 
       width = 7, height = 5.4, units = "in")

##################################################################
#
#             export aggregated data for Dryad
#

kidsHH_agg <- kidsHH %>%   # select variables used in primary, secondary and sensitivity analyses
  select(village, LakeYN, Pop, OD10p, concession,
         DEM01.x, DEM04,
         Sh_presence_18, Sh_mean_18, Sh_medianR_18,
         ENF02, F_WkAvg, FD_WkAvg, FDB_WkAvg,
         habitatArea_peak, habitatAreaPeak_d, 
         habitatArea_yr, habitatAreaYr_d, 
         habitatAreaV_peak, habitatAreaV_yr,
         surface, surfaceYN, 
         sanitationYN, sanitation, privateSan,
         wealthQuintile1) %>%
  mutate(exp_cat0 = ifelse(ENF02 == 0, 1, 0),
         exp_cat1 = ifelse(ENF02 == 1, 1, 0),
         exp_cat2 = ifelse(ENF02 == 2, 1, 0),
         exp_cat3 = ifelse(ENF02 == 3, 1, 0),
         exp_cat4 = ifelse(ENF02 == 4, 1, 0),
         surface0 = ifelse(surface == 0, 1, 0),
         surface1 = ifelse(surface == 1, 1, 0),
         surface2 = ifelse(surface == 2, 1, 0),
         sanitation0 = ifelse(sanitation == 0, 1, 0),
         sanitation1 = ifelse(sanitation == 1, 1, 0),
         sanitation2 = ifelse(sanitation == 2, 1, 0),
         ses1 = ifelse(wealthQuintile1 == 1, 1, 0),
         ses2 = ifelse(wealthQuintile1 == 2, 1, 0), 
         ses3 = ifelse(wealthQuintile1 == 3, 1, 0),
         ses4 = ifelse(wealthQuintile1 == 4, 1, 0),
         ses5 = ifelse(wealthQuintile1 == 5, 1, 0)) %>% 
  group_by(village, DEM01.x) %>% 
  # calculate means of numeric variables, sum binary variables
  summarise(nInd = n_distinct(OD10p),
            nHH = n_distinct(concession),
            age_mean = mean(DEM04, na.rm = T),
            village_pop = mean(Pop, na.rm = T),
            Sh_cases = sum(Sh_presence_18, na.rm = T),
            Sh_mean = mean(Sh_mean_18, na.rm = T),
            Sh_median = median(Sh_medianR_18, na.rm = T),
            nExp_0 = sum(exp_cat0, na.rm = T),
            nExp_1 = sum(exp_cat1, na.rm = T),
            nExp_2 = sum(exp_cat2, na.rm = T),
            nExp_3 = sum(exp_cat3, na.rm = T),
            nExp_4 = sum(exp_cat4, na.rm = T),
            F_mean = mean(F_WkAvg, na.rm = T), 
            FD_mean = mean(FD_WkAvg, na.rm = T), 
            FDB_mean = mean(FDB_WkAvg, na.rm = T),
            habitatPeak_mean = mean(habitatArea_peak, na.rm = T),
            habitatPeakD_mean = mean(habitatAreaPeak_d, na.rm = T),
            habitatYear_mean = mean(habitatArea_yr, na.rm = T),
            habitatYearD_mean = mean(habitatAreaYr_d, na.rm = T),
            habitatPeakV_mean = mean(habitatAreaV_peak, na.rm = T),
            habitatYearV_mean = mean(habitatAreaV_yr, na.rm = T),
            nSurface_0 = sum(surface0, na.rm = T),
            nSurface_1 = sum(surface1, na.rm = T),
            nSurface_2 = sum(surface2, na.rm = T),
            nSurface_Any = sum(surfaceYN, na.rm = T),
            nSanit_0 = sum(sanitation0, na.rm = T),
            nSanit_1 = sum(sanitation1, na.rm = T),
            nSanit_2 = sum(sanitation2, na.rm = T),
            nSanit_Any = sum(sanitationYN, na.rm = T),
            nPrivateSan = sum(privateSan, na.rm = T),
            ses_mean = mean(wealthQuintile1, na.rm = T),
            nSes_1 = sum(ses1, na.rm = T),
            nSes_2 = sum(ses2, na.rm = T),
            nSes_3 = sum(ses3, na.rm = T),
            nSes_4 = sum(ses4, na.rm = T),
            nSes_5 = sum(ses5, na.rm = T)) %>% 
  rename("sex" = "DEM01.x")

write.csv(kidsHH_agg, file = "data/ehvData_agg.csv")


