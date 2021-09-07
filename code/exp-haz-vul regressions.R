###################################################
#
#             regression analysis of
#        exposure - hazard - vulnerability
#               Code by: Andrea Lund
#         Last updated: April 23, 2020

# install packages 
library(tidyverse) 
library(lme4) # for mixed effects logistic regression 
library(glmmTMB) # for mixed effects negative binomial regression
library(AICcmodavg) # for model averaging of candidate set
library(performance) # for binned residuals of binomial models
library(ggplot2)
library(scales) # for plotting with logarithmic scales

# set working directory 
setwd("C:/Users/andre/Box Sync/Papers/Ch3Prawn-CEA")

# source script that calculates exposure indices
source("./Code/constructing exposure indices.R")

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

# only 1030 observations for ENF02, so subset entire data set to make data set same size for all models
kidsHH <- kidsHH[!is.na(kidsHH$ENF02) & !is.na(kidsHH$F_WkAvg) & !is.na(kidsHH$Sh_medianR_18),]

##############   FUNCTIONS  ##############  
## calculate confidence intervals for fitted models
# https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
ci_calc <- function(mod) {
  se <- sqrt(diag(vcov(mod)))
  estCI <- cbind(OR = fixef(mod), 
                 LL = fixef(mod)-1.96*se, 
                 UL = fixef(mod)+1.96*se)
  orCI <- exp(estCI)
  return(orCI)
}

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

AIC(pexp1) # 1058.997
ci_calc(pexp1)

# exposure = frequency only
pexp2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 F_WkAvgC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp2) # 1061.193

# exposure = frequency + duration
pexp3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 FD_WkAvgC +
                 LakeYN + scale(Pop) + 
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp3) # 1061.233

# exposure = bsa-adjusted frequency + duration
pexp4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 FDB_WkAvgC +
                 LakeYN + scale(Pop) + 
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pexp4) # 1061.386

########     HAZARD-ONLY MODELS
# habitat area during transmission peak
phaz1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   habitatArea_peakC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))

AIC(phaz1) # 1060.099

# habitat area during entire year
phaz2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatArea_yrC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz2) # 1059.592

# distance-weighted habitat area during transmission peak
phaz3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaPeak_dC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz3) # 1055.312

# distance-weighted habitat area during entire year
phaz4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaYr_dC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0011511

AIC(phaz4) # 1052.944

# village habitat area during transmission peak
phaz5 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaV_peakC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz5) # 1060.584

# village habitat area over entire year
phaz6 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 habitatAreaV_yrC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(phaz6) # 1058.493

######## VULNERABILITY-ONLY MODELS
# vulnerability = surface water dependence
pvul1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(surface) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

# max grad = 0.00211662
AIC(pvul1) # 1055.452

# vulnerability = binary surface water dependence
pvul2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 surfaceYN +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pvul2) # 1057.523

# vulnerability = binary private vs shared/no sanitation
pvul3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 privateSan +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pvul3) # 1057.383

# vulnerability = categorical sanitation (private vs shared vs no)
pvul4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(sanitation) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0236004

AIC(pvul4) # 1058.963

# vulnerability = asset index
pvul5 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(wealthQuintile1) +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))

AIC(pvul5) # 1063.277

########     CONTACT-HAZARD MODELS
# habitat area during transmission peak
peh1.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                 LakeYN + scale(Pop) + 
                 factor(ENF02) + habitatArea_peakC +
                 (1|village/C_ConcessionId), 
               data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0120003

AIC(peh1.1) # 1059.233

# total habitat area over the year
peh2.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +  
                   factor(ENF02) + habitatArea_yrC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.018522
AIC(peh2.1) # 1058.783

# village  habitat area at transmission peak
peh3.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +  
                  factor(ENF02) + habitatAreaV_peakC +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = binomial(link = "logit"))

AIC(peh3.1) # 1060.18

# village  habitat area over entire year
peh4.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +  
                  factor(ENF02) + habitatAreaV_yrC +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = binomial(link = "logit"))

AIC(peh4.1) # 1058.265

# distance-weighted habitat area at transmission peak
pehd1.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +  
                   factor(ENF02) + habitatAreaPeak_dC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00437232
AIC(pehd1.1) # 1055.593
ci_calc(pehd1.1)

# distance-weighted habitat area over entire year 
pehd2.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +                     
                    factor(ENF02) + habitatAreaYr_dC +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00240995
AIC(pehd2.1) # 1053.311
ci_calc(pehd2.1)

########     CONTACT-VULNERABILITY  MODELS
pev1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) + 
                  factor(ENF02) + factor(surface) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00346483
AIC(pev1) # 1056.22

pev2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                factor(ENF02) + surfaceYN +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
AIC(pev2) # 1058.974

pev3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                factor(ENF02) + privateSan +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
AIC(pev3) # 1058.279

########     HAZARD-VULNERABILITY  MODELS
phv1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaPeak_dC + factor(surface) +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
AIC(phv1) # 1051.982

phv2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaYr_dC + factor(surface) +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.002079
AIC(phv2) # 1050.23

phv3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaPeak_dC + surfaceYN +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))

AIC(phv3) # 1054.754

phv4 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaYr_dC + surfaceYN +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00121944
AIC(phv4) # 1052.525

phv5 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaPeak_dC + privateSan +
                (1|village/C_ConcessionId), 
              data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00160567
AIC(phv5) # 1052.789

phv6 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                LakeYN + scale(Pop) + 
                habitatAreaYr_dC + privateSan +
                (1|village/C_ConcessionId),
              data = kidsHH, family = binomial(link = "logit"))

AIC(phv6) # 1050.454


########    (GLOBAL) CONTACT-HAZARD-VULNERABILITY MODELS 

# dependence on surface water
pehv1.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) + 
                     factor(ENF02) + habitatAreaYr_dC + factor(surface) +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00859006
summary(pehv1.1)
AIC(pehv1.1) # 1051.261
ci_calc(pehv1.1)
vif(pehv1.1) # no concerning values here

# binary surface water dependence
pehv2.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) + 
                     factor(ENF02) + habitatAreaYr_dC + surfaceYN +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00246493
AIC(pehv2.1) # 1054.195
ci_calc(pehv2.1)

pehv2.2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) + 
                     factor(ENF02)*surfaceYN + habitatAreaYr_dC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00641314
AIC(pehv2.2) # 1045.879

anova(pehv2.1, pehv2.2, test = "Chisq") # Chisq = 61.021, df = 3, p = 0.0135

pehv2.3 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   factor(ENF02)*surfaceYN*habitatAreaYr_dC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0175674
AIC(pehv2.3) # 1058.649

# private sanitation (y/n)
pehv3.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   factor(ENF02) + habitatAreaYr_dC + privateSan +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00213937
AIC(pehv3.1) # 1052.092
ci_calc(pehv3.1)

pehv3.2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   factor(ENF02)*privateSan + habitatAreaYr_dC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00249497
AIC(pehv3.2) # 1059.13

pehv4.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   factor(ENF02) + habitatAreaPeak_dC + factor(surface) +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0017862
AIC(pehv4.1) # 1052.938
ci_calc(pehv4.1)

pehv4.2 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   factor(ENF02)*factor(surface) + habitatAreaPeak_dC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00560929
AIC(pehv4.2) # 1051.52

pehv5.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   factor(ENF02) + habitatAreaPeak_dC + surfaceYN +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.0010271
AIC(pehv5.1) # 1056.394

pehv6.1 <- glmer(Sh_presence_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   ENF02 + habitatAreaPeak_dC + privateSan +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = binomial(link = "logit"))
# max grad = 0.00183034
AIC(pehv6.1) # 1051.81

##################################################################
#
#           compute model averages for presence models

# create a list of all presence models
presenceModels <- list(pexp1, pexp2, pexp3, pexp4,
                       phaz1, phaz2, phaz3, phaz4, phaz5, phaz6,
                       pvul1, pvul2, pvul3, pvul4, pvul5,
                       pehd2.1,
                       pev1, pev3,
                       phv2, phv6,
                       pehv1.1, pehv3.1)

modNamesP <- c('pexp1', 'pexp2', 'pexp3', 'pexp4',
                 'phaz1', 'phaz2', 'phaz3', 'phaz4', 'phaz5', 'phaz6',
                 'pvul1', 'pvul2', 'pvul3', 'pvul4', 'pvul5',
                 'pehd2.1',
                 'pev1', 'pev3',
                 'phv2', 'phv6',
                 'pehv1.1', 'pehv3.1')

# construct model selection table 
presenceAIC <- aictab(presenceModels, modnames = modNamesP, second.ord = F, sort = F) # K is 3 greater than what I have in the spreadsheet - is this random effects?
write.csv(presenceAIC, "./Figures & Tables/presenceAIC22_factor.csv")
presenceBIC <- bictab(presenceModels, modnames = modNamesP, second.ord = F, sort = F)
write.csv(presenceBIC, "./Figures & Tables/presenceBIC22_factor.csv")
# determine confidence set
confset(presenceModels, modNamesP, method = 'raw', second.ord = F, level = 0.95) # 8 models
confSetP22 <- confset(presenceModels, modNamesP, method = 'ordinal', second.ord = F) 
  # substantial weight - 5 models, some weight - 5, little weight - 5 models, no weight - 7 models
write.csv(confSetP22$substantial, "./Figures & Tables/confsetPresence22_factor.csv")
confset(presenceModels, modNamesP, method = 'ratio', second.ord = F, delta = 10) # 15 models with delta = 10 cutoff

# diagnostics on models with substantial weight
binned_residuals(pehv3.1) # 66%
binned_residuals(pehv1.1) # 66%
binned_residuals(phv2) # 66%
binned_residuals(phv6) # 75
#binned_residuals(pehd2.1) # 72%
#binned_residuals(pehv4.1) # 81%
#binned_residuals(pehv2.1) # 75%
#binned_residuals(pehv4.2) # 75%

# compute model-averaged estimates of individual parameters
fraw1 <- modavg(presenceModels, parm = "factor(ENF02)1", modnames = modNamesP, second.ord = F); fraw1
fraw2 <- modavg(presenceModels, parm = "factor(ENF02)2", modnames = modNamesP, second.ord = F); fraw2
fraw3 <- modavg(presenceModels, parm = "factor(ENF02)3", modnames = modNamesP, second.ord = F); fraw3
fraw4 <- modavg(presenceModels, parm = "factor(ENF02)4", modnames = modNamesP, second.ord = F); fraw4
#fsum <- modavg(presenceModels, parm = "F_WkAvg", modnames = modNamesP, second.ord = F); fsum
#fd <- modavg(presenceModels, parm = "FD_WkAvgC", modnames = modNamesP, second.ord = F); fd
#fdb <- modavg(presenceModels, parm = "FDB_WkAvgC", modnames = modNamesP, second.ord = F); fdb
areaPeak_d <- modavg(presenceModels, parm = "habitatAreaPeak_dC", modnames = modNamesP, second.ord = F); areaPeak_d
areaYear_d <- modavg(presenceModels, parm = "habitatAreaYr_dC", modnames = modNamesP, second.ord = F); areaYear_d
areaPeak <- modavg(presenceModels, parm = "habitatArea_peakC", modnames = modNamesP, second.ord = F); areaPeak
areaYear <- modavg(presenceModels, parm = "habitatArea_yrC", modnames = modNamesP, second.ord = F); areaYear
areaPeakV <- modavg(presenceModels, parm = "habitatAreaV_peakC", modnames = modNamesP, second.ord = F); areaPeakV
areaYearV <- modavg(presenceModels, parm = "habitatAreaV_yrC", modnames = modNamesP, second.ord = F); areaYearV
surface1 <- modavg(presenceModels, parm = "factor(surface)1", modnames = modNamesP, second.ord = F); surface1
surface2 <- modavg(presenceModels, parm = "factor(surface)2", modnames = modNamesP, second.ord = F); surface2
surfaceYN <- modavg(presenceModels, parm = "surfaceYN", modnames = modNamesP, second.ord = F); surfaceYN
privateSan <- modavg(presenceModels, parm = "privateSan", modnames = modNamesP, second.ord = F); privateSan
sex <- modavg(presenceModels, parm = "factor(DEM01.x)2", modnames = modNamesP, second.ord = F); sex
age <- modavg(presenceModels, parm = "scale(DEM04)", modnames = modNamesP, second.ord = F); age
vPopn <- modavg(presenceModels, parm = "scale(Pop)", modnames = modNamesP, second.ord = F); vPopn
vLocn <- modavg(presenceModels, parm = "LakeYN", modnames = modNamesP, second.ord = F); vLocn

kidsHH$surfaceF <- as.factor(kidsHH$surface)

pehv4.1_int <- glmer(Sh_presence_18 ~ DEM01.x + DEM04 + 
                     LakeYN + scale(Pop) + 
                     ENF02*surfaceF + habitatAreaPeak_dC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = binomial(link = "logit"))
summary(pehv4.1_int)
anova(pehv4.1_int, pehv4.1, test = "Chisq") # Chisq = 7.8138, df = 2, p = 0.0201

# calculate confidence intervals for surface strata
# ENF02 [6] surfaceF1 [7] surfaceF2 [8] ENF02:surfaceF1 [10] ENF02:surfaceF2 [11]
# 1 vs 0
coefs <- fixef(pehv4.1_int)
se <- sqrt(diag(vcov(pehv4.1_int)))
eBeta <- coefs[6]; eSE <- se[6]
vBeta <- coefs[7]; vSE <- se[7]
evBeta <- coefs[10]; evSE <- se[10]
est1v0 <- vBeta + evBeta
ll_1v0 <- (vBeta - vSE) + (evBeta - evSE)
ul_1v0 <- (vBeta + vSE) + (evBeta + evSE)
estCI_1v0 <- cbind(est1v0, ll_1v0, ul_1v0)
ORCI_1v0 <- exp(estCI_1v0); ORCI_1v0

# 2 vs 0
eBeta <- coefs[6]; eSE <- se[6]
vBeta <- coefs[8]; vSE <- se[8]
evBeta <- coefs[11]; evSE <- se[11]
est2v0 <- vBeta + evBeta
ll_2v0 <- (vBeta - vSE) + (evBeta - evSE)
ul_2v0 <- (vBeta + vSE) + (evBeta + evSE)
estCI_2v0 <- cbind(est2v0, ll_2v0, ul_2v0)
ORCI_2v0 <- exp(estCI_2v0); ORCI_2v0

stratOR <- rbind(ORCI_1v0, ORCI_2v0); stratOR


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
AIC(cexp1) # 7112.272


# exposure = frequency only
cexp2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   F_WkAvgC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cexp2) # 7130.451

# exposure = frequency + duration
cexp3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                   LakeYN + scale(Pop) + 
                   FD_WkAvgC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))

AIC(cexp3) # 7129.78

# exposure = bsa-adjusted frequency + duration
cexp4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) + 
                   FDB_WkAvgC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cexp4) # 7129.356

######## HAZARD-ONLY MODELS
chaz1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatArea_peakC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz1) # 7126.897

chaz2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatArea_yrC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz2) # 7128.062

chaz3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaPeak_dC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz3) # 7129.43

chaz4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaYr_dC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz4) # 7128.901

chaz5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaV_peakC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz5) # 7130.122

chaz6 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   habitatAreaV_yrC + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chaz6) # 7128.812

######## VULNERABILITY-ONLY MODELS
cvul1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(surface) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul1) # 7129.315

cvul2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   surfaceYN + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul2) # 7127.337

cvul3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   privateSan + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul3) # 7129.115

cvul4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(sanitation) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul4) # 7131.096

cvul5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +                    
                   factor(wealthQuintile1) + 
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(cvul5) # 7134.819

########     CONTACT-HAZARD MODELS
# habitat area during transmission peak
ceh1.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_peakC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))
AIC(ceh1.1) # 7111.857
exp(confint(ceh1.1))

# habitat area during entire year
ceh2.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_yrC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh2.1) # 7112.662
exp(confint(ceh2.1))

# habitat area during transmission peak
ceh3.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + habitatAreaPeak_dC +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh3.1) # 7113.293
exp(confint(ceh3.1))

# village habitat area during entire year
ceh4.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + habitatAreaV_yrC +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))

AIC(ceh4.1) #  7112.556
exp(confint(ceh4.1))


# distance-weighted habitat area during transmission peak
cehd1.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatAreaPeak_d +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))
AIC(cehd1.1) # 7113.293
exp(confint(cehd1.1))

# distance-weighted habitat area over entire year
cehd2.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                      LakeYN + scale(Pop) +
                      factor(ENF02) + habitatAreaYr_d +
                      (1|village/C_ConcessionId), 
                    data = kidsHH, family = nbinom2(link = "log"))
AIC(cehd2.1) # 7112.815
exp(confint(cehd2.1))

########     CONTACT-VULNERABILITY MODELS
cev1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                    LakeYN + scale(Pop) +
                    factor(ENF02) + factor(surface) +
                    (1|village/C_ConcessionId), 
                  data = kidsHH, family = nbinom2(link = "log"))
AIC(cev1) # 7115.648

cev2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  factor(ENF02) + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(cev2) # 7118.54

cev3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  factor(ENF02) + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(cev3) # 7113.97

cev4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  factor(ENF02) + factor(sanitation) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(cev4) # 7115.969

########     HAZARD-VULNERABILITY MODELS
chv1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_peakC + factor(surface) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv1) # 7127.258

chv2 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_yrC + factor(surface) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv2) # 7128.498

chv3 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaV_yrC + factor(surface) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv3) # 7129.745

chv4 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_peakC + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv4) # 7125.271

chv5 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_yrC + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv5) # 7126.538

chv6 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaV_yrC + surfaceYN +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv6) # 7127.749

chv7 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_peakC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv7) # 7127.679

chv8 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatArea_yrC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv8) # 7128.646

chv9 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaV_yrC + privateSan +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv9) # 7129.303

chv10 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                  LakeYN + scale(Pop) +
                  habitatAreaPeak_dC + factor(surface) +
                  (1|village/C_ConcessionId), 
                data = kidsHH, family = nbinom2(link = "log"))
AIC(chv10) # 7130.581

chv11 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaYr_dC + factor(surface) +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv11) # 7130.132

chv12 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaPeak_dC + surfaceYN +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv12) # 7128.619

chv13 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaYr_dC + surfaceYN +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv13) # 7128.157

chv14 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaPeak_dC + privateSan +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv14) # 7129.776

chv15 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                   LakeYN + scale(Pop) +
                   habitatAreaYr_dC + privateSan +
                   (1|village/C_ConcessionId), 
                 data = kidsHH, family = nbinom2(link = "log"))
AIC(chv15) # 7129.273

########     CONTACT-HAZARD MODELS STRATIFIED BY VULNERABILITY

# surface water dependence
cehv1.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                       LakeYN + scale(Pop) +
                       factor(ENF02) + habitatArea_peakC + factor(surface) +
                       (1|village/C_ConcessionId), 
                     data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv1.1) # 7115.084
exp(confint(cehv1.1))


cehv2.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                       LakeYN + scale(Pop) +
                       factor(ENF02) + habitatArea_yrC + factor(surface) +
                       (1|village/C_ConcessionId), 
                     data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv2.1) # 7115.925


cehv3.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                      LakeYN + scale(Pop) +
                      factor(ENF02) + habitatArea_peakC + surfaceYN +
                      (1|village/C_ConcessionId), 
                    data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv3.1) # 7113.084

cehv4.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_yrC + surfaceYN +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv4.1) # 7118.852
exp(confint(cehv4.1))

cehv5.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) +
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_peakC + privateSan +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv5.1) # 7120.011


cehv6.1 <- glmmTMB(Sh_medianR_18 ~ factor(DEM01.x) + scale(DEM04) + 
                     LakeYN + scale(Pop) +
                     factor(ENF02) + habitatArea_yrC + privateSan +
                     (1|village/C_ConcessionId), 
                   data = kidsHH, family = nbinom2(link = "log"))

AIC(cehv6.1) # 7118.852

##################################################################
#
#           compute model averages for presence models

# create a list of all presence models
intensityModels <- list(cexp1, cexp2, cexp3, cexp4,
                       chaz1, chaz2, chaz3, chaz4, chaz5, chaz6,
                       cvul1, cvul2, cvul3, cvul4, cvul5,
                       ceh1.1, ceh2.1, cehd2.1, ceh4.1,
                       cev1, cev2, cev3,
                       chv1, chv2, chv3, chv4, chv5, chv6,
                       chv7, chv8, chv9, chv11, chv13, chv15,
                       cehv1.1, cehv2.1, cehv3.1,
                       cehv4.1, cehv5.1, cehv6.1)

modNamesI <- c('cexp1', 'cexp2', 'cexp3', 'cexp4',
              'chaz1', 'chaz2', 'chaz3', 'chaz4', 'chaz5', 'chaz6',
              'cvul1', 'cvul2', 'cvul3', 'cvul4', 'cvul5',
              'ceh1.1', 'ceh2.1', 'cehd2.1', 'ceh4.1',
              'cev1', 'cev2', 'cev3',
              'chv1', 'chv2', 'chv3', 'chv4', 'chv5', 'chv6',
              'chv7', 'chv8', 'chv9', 'chv11', 'chv13', 'chv15',
              'cehv1.1', 'cehv2.1', 'cehv3.1',
              'cehv4.1', 'cehv5.1', 'cehv6.1')

# construct model selection table 
intensityAIC <- aictab(intensityModels, modnames = modNamesI, second.ord = F, sort = F) # K is 3 greater than what I have in the spreadsheet - is this random effects?
write.csv(intensityAIC, "./Figures & Tables/intensityAIC40_factor.csv")
intensityBIC <- bictab(intensityModels, modnames = modNamesI, second.ord = F, sort = F)
write.csv(intensityBIC, "./Figures & Tables/intensityBIC40_factor.csv")
# determine confidence set
confset(intensityModels, modNamesI, method = 'raw', second.ord = F, level = 0.95) # 13 models
intConfSet <- confset(intensityModels, modNamesI, method = 'ordinal', second.ord = F) 
write.csv(intConfSet$substantial, "./Figures & Tables/confsetIntensity40_factor.csv")
# substantial weight - 11 models, some weight - 3, little weight - 6 models, no weight - 20 models
confset(intensityModels, modNamesI, method = 'ratio', second.ord = F, delta = 10) # 20 models with delta = 10 cutoff

# diagnostics on models with substantial weight

# compute model-averaged estimates of individual parameters
fraw1 <- modavg(intensityModels, parm = "factor(ENF02)1", modnames = modNamesI, second.ord = F); fraw1
fraw2 <- modavg(intensityModels, parm = "factor(ENF02)2", modnames = modNamesI, second.ord = F); fraw2
fraw3 <- modavg(intensityModels, parm = "factor(ENF02)3", modnames = modNamesI, second.ord = F); fraw3
fraw4 <- modavg(intensityModels, parm = "factor(ENF02)4", modnames = modNamesI, second.ord = F); fraw4
areaPeak <- modavg(intensityModels, parm = "habitatArea_peakC", modnames = modNamesI, second.ord = F); areaPeak
areaYear <- modavg(intensityModels, parm = "habitatArea_yrC", modnames = modNamesI, second.ord = F); areaYear
areaPeak_d <- modavg(intensityModels, parm = "habitatAreaPeak_dC", modnames = modNamesI, second.ord = F); areaPeak_d
areaYear_d <- modavg(intensityModels, parm = "habitatAreaYr_dC", modnames = modNamesI, second.ord = F); areaYear_d
areaYearV <- modavg(intensityModels, parm = "habitatAreaV_yrC", modnames = modNamesI, second.ord = F); areaYearV
surface1 <- modavg(intensityModels, parm = "factor(surface)1", modnames = modNamesI, second.ord = F); surface1
surface2 <- modavg(intensityModels, parm = "factor(surface)2", modnames = modNamesI, second.ord = F); surface2
surfaceYN <- modavg(intensityModels, parm = "surfaceYN", modnames = modNamesI, second.ord = F); surfaceYN
privateSan <- modavg(intensityModels, parm = "privateSan", modnames = modNamesI, second.ord = F); privateSan
#sanitation1 <- modavg(intensityModels, parm = "factor(sanitation)1", modnames = modNamesI, second.ord = F); sanitation1
#sanitation2 <- modavg(intensityModels, parm = "factor(sanitation)2", modnames = modNamesI, second.ord = F); sanitation2
sex <- modavg(intensityModels, parm = "factor(DEM01.x)2", modnames = modNamesI, second.ord = F); sex
age <- modavg(intensityModels, parm = "scale(DEM04)", modnames = modNamesI, second.ord = F); age
vPopn <- modavg(intensityModels, parm = "scale(Pop)", modnames = modNamesI, second.ord = F); vPopn
vLocn <- modavg(intensityModels, parm = "LakeYN", modnames = modNamesI, second.ord = F); vLocn



# evidence ratios
evidence(intensityAIC)

##################################################################
#
#           visualize effect sizes and precision for 
#           E, H and V components of risk

# import model-averaged estimates
allEst <- read.csv("./Figures & Tables/allEstimates_factor_allVars.csv")
allEst$Type <- factor(allEst$Type, 
                      levels = c("Vulnerability", "Hazard", "Exposure", "Demographics"))
allEst$est <- factor(allEst$est, 
                      levels = c("Presence", "Intensity"))

library(scales)

# generate forest plot for all Sh presence estimates
avgPlot <- ggplot(allEst,
       aes(x = Variable, y = OR, ymin = orLL, ymax = orUL)) +
  geom_pointrange(aes(col = Type, pch = Type),
                  position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = orLL, ymax = orUL, col = Type, pch = Type), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab("Variable") +
  ylab("Odds Ratio (Presence) / Rate Ratio (Intensity)") +
  scale_y_continuous(trans = 'log2', limits = c(0.25, 18), 
                     breaks = trans_breaks('log2', function(x) 2^x)) +
  facet_grid(est ~ .) +
  coord_flip() +
  #scale_linetype_manual(values = c("longdash", "dotdash", "solid")) +
  scale_colour_manual(values = c("steelblue4", "turquoise4", "black", "darkorange2")) +
  scale_shape_manual(values = c(15:18)) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("Figure 2.tiff", plot = avgPlot, device = "tiff",
       path = "./Figures & Tables", dpi = 300)

# subset data into vulnerable and not vulnerable strata 
# based on binary surface water dependence variable 

# more vulnerable
haz3d_V <- glmer(Sh_medianR_18 ~ DEM01.x + DEM04 + 
                   ENF02 + pctHabitat_dC +
                   (1|village/C_ConcessionId), 
                 data = kidsHH_V, family = nbinom2(link = "log"))
summary(haz3d_V)
AIC(haz3d_V) # Inf

# less vulnerable
haz3d_noV <- glmer(Sh_medianR_18 ~ DEM01.x + DEM04 + 
                     ENF02 + pctHabitat_dC +
                     (1|village/C_ConcessionId), 
                   data = kidsHH_noV, family = nbinom2(link = "log"))
summary(haz3d_noV)
AIC(haz3d_noV) # Inf