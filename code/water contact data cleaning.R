#############################################################
#
#           cleaning both rounds of survey data
#           + all infection data
#           + household location data
#           + ecological data
#           for
#           analysis of water contact behavior
#                   Dissertation chapter 3
#                   Code by: Andrea Lund
#                   Created: 14 June 2019
#                Last updated: 7 August 2020
#
#############################################################

# install packages
library(readxl)
library(tidyverse)
library(lubridate)
library(here)


#########################################
#
#       import data
#

# import survey data
indOld2018 <- read_excel(here("data/S11.xlsx"))
indOld2016 <- read_excel(here("data/S1-resolved.xls"))
indNew2018 <- read_excel(here("data/S12.xlsx"))
enfants2018 <- read_excel(here("data/S14.xlsx"))
watercontact <- read_excel(here("data/S51.xlsx")) # WCO1-WC07
watercontact2 <- read_excel(here("data/S51_S5ActiviteEau.xlsx"))
household <- read_excel(here("data/household.xlsx"))

# import infection data
infect18 <- read_excel(here("data/infect18.xlsx"))
allInfect <- read_csv(here("data/allInfection.csv"))

# import household location data 
hhloc <- read_csv(here("data/householdGPS_2016_wDistance.csv"))
hhloc <- rename(hhloc, "concession" = "C_concessi", "nearestWatPt" = "Distance matrix_TargetID", "distanceWatPt" = "Distance matrix_Distance")

# import 2016 wealth indices 
wealth2016 <- read_csv(here("data/wealthData2016.csv"))

wealth2016 <- wealth2016 %>% 
  rename("concession" = "C_concession")

# import distance between household and nearest water point for snail sampling sites only
hhloc_snail <- read_csv(here("data/snailHH18_distanceMatrix.csv"))

# import water point size and snail metrics from ecological data
enviro <- read_csv(here("data/SiteData_2years.csv"))

#########################################
#
#       process household data
#

#### survey data #### 
# use questions on reported water sources to calculate surface water dependence
for (i in 1:nrow(household)) {
  household$surface[i] = ifelse((household$HH22[i] == 5 && household$HH23[i] == 5), 2, 
                              ifelse((household$HH22[i] == 5 && household$HH23[i] %in% c(1:4,6)), 1,
                                     ifelse((household$HH22[i] %in% c(1:4,6) && household$HH23[i] == 5), 1,
                                            ifelse((household$HH22[i] %in% c(1:4,6) && household$HH23[i] %in% c(1:4,6)), 0))))
}

# calculate binary variable for dependence on surface water 
household$surfaceYN <- ifelse(household$surface > 0, 1, 0)

# calculate access to private, shared or no sanitation infrastructure
household$sanitationYN <- ifelse(is.na(household$HH25), NA,
                               ifelse(household$HH25 %in% c(1:3), 1, 0))

for (i in 1:nrow(household)) {
  household$sanitation[i] = ifelse(is.na(household$sanitationYN[i]), NA,
                                 ifelse((household$sanitationYN[i] == 1 && household$HH26[i] == 0), 2, 
                                        ifelse((household$sanitationYN[i] == 1 && household$HH26[i] == 1), 1,
                                               ifelse((household$sanitationYN[i] == 0), 0, -9))))
}

household$privateSan <- ifelse(household$sanitation == 2, 1, 0)

#### location data #### 
# investigate repeated cases of CC/GK/10
# two rows for CC/GK/10 with the same lat/lon coordinates and concession ID, apparently imported from different text files
hhloc <- hhloc %>% 
  distinct(concession, .keep_all = T) %>%  # remove row where C_concession == "CC/GK/10" & file = "Device9-ME-SS.txt"
  filter(concession != "CC/ME/56") # remove concession whose GPS coordinates map to another village

# combine location data
hhloc_all <- dplyr::full_join(hhloc, hhloc_snail, by = "concession") %>% 
              select(concession, lat, lon, nearestWatPt, distanceWatPt, nearestSnailSite, distanceSnail)

#########################################
#
#       process ecological data
#

snail_peak <- enviro %>% 
  filter(FM == 5) %>% 
  select(Site, FM, Date, lat, long, AreaPrawn, PercOther,
         BulinusAvgDensSite, ShSnailPrevSite, 
         BulinusTotalSite, ShSnailTotalSite, ShSnailDensSite,
         LakeYN, Pop) %>% 
  mutate(habitatArea_peak = AreaPrawn*PercOther) %>% 
  rename("nearestSnailSite" = "Site", "lat_site" = "lat", "lon_site" = "long",
         "siteArea_peak" = "AreaPrawn", "pctHabitat_peak" = "PercOther",
         "snailDens_peak" = "BulinusAvgDensSite", "snailPrev_peak" = "ShSnailPrevSite", 
         "snailTot_peak" = "BulinusTotalSite", "infSnail_peak" = "ShSnailTotalSite",
         "infSnailDens_peak" = "ShSnailDensSite")

snail_yr <- enviro %>% 
  filter(FM %in% c(4:6)) %>% 
  select(Site, FM, Date, AreaPrawn, PercOther,
         BulinusAvgDensSite, ShSnailPrevSite, VegMassSite, 
         BulinusTotalSite, ShSnailTotalSite, ShSnailDensSite,
         LakeYN, Pop) %>% 
  mutate(habitatArea = AreaPrawn*PercOther) %>% 
  group_by(Site) %>% 
  summarize(habitatArea_yr = sum(habitatArea),
            siteArea_yr = sum(AreaPrawn)) %>%
  rename("nearestSnailSite" = "Site")  

snailV_yr <- enviro %>% 
  filter(FM %in% c(4:6)) %>% 
  select(Site, FM, Date, AreaPrawn, PercOther, LakeYN, Village, VillageCode) %>% 
  mutate(habitatArea = AreaPrawn*PercOther) %>% 
  group_by(VillageCode) %>% 
  summarize(habitatAreaV_yr = sum(habitatArea))

snailV_peak <- enviro %>% 
  filter(FM == 5) %>% 
  select(Site, FM, Date, AreaPrawn, PercOther, LakeYN, Village, VillageCode) %>% 
  mutate(habitatArea = AreaPrawn*PercOther) %>% 
  group_by(VillageCode) %>% 
  summarize(habitatAreaV_peak = sum(habitatArea))

snail_vars <- full_join(snail_peak, snail_yr, by = "nearestSnailSite")

snailV_vars <- full_join(snailV_peak, snailV_yr, by = "VillageCode")

#########################################
#
#       process individual-level data
#

# subset 2018 data on 2016 household members that are still present, with relevant variables
oldMemPresent <- indOld2018 %>% 
  filter(DEM18 == 1) %>% 
  select(concession, Rang, DEM19, DEM20, OD01, OD01p, OD02, OD03, OD03p, OD04, OD04p) %>% 
  rename_at(vars(-"concession", -"Rang"), function(x) paste0(x, "_18"))
#summary(oldMemPresent) # n = 8172

# subset 2016 household members that have left, with relevant varialbes
oldMemGone <- indOld2018 %>% 
  filter(DEM18 == 0) %>% 
  select(concession, VillageId, Rang, DEM18, DEM21, DEM21p, DEM11, DEM11p, DEMD3)
#summary(oldMemGone) # n = 2115

# subset 2016 data on 2016 observations only
indOld2016_16 <- indOld2016 %>% 
  mutate(yr = year(ymd(date1))) %>% 
  filter(yr == 2016) %>% # n = 10987
  select(Libelle, Initial, C_ConcessionId, concession, date1, lang, 
         enqueteur, remarks, Rang, DEM01:DEMD3) %>% 
  mutate(DEM04_18 = DEM04 + 2)

# select variables for new 2018 members
indNew2018_ <- indNew2018 %>% 
  select(-ConcessionId, -VillageId, -S12Id, 
         -WC02, -WC03, -WC04, -WC04, -WC05, 
         -WC06, -WC07, -DEM07p) %>% 
  rename("Rang" = "rang") # n = 1336

#########################################
#
#           combine data sets
#
#

#### FOR CHAPTER 3 ANALYSIS ####
# combine all individuals into single data frame
oldNew <- plyr::rbind.fill(indOld2016_16, indNew2018_)

watercontact <- select(watercontact, -EnqueteId)

allInd <- full_join(oldNew, oldMemPresent, by = c("concession", "Rang")) %>% 
  full_join(., watercontact, by = c("concession", "Rang")) %>% 
  filter(!is.na(WC01))

# 2016 data from parasito kids
enfants2016 <- indOld2016_16 %>%
  filter(!is.na(OD10p)) %>% 
  select(Libelle, Initial, concession, date1, lang,
         enqueteur, remarks, Rang, DEM01:DEM11b, DEM04_18, OD10p) %>% 
  rename_at(vars("date1":"remarks"), function(x) paste0(x, "_16"))

# 2018 data from parasito kids
enfants2018 <- enfants2018 %>% 
  select(concession, Rang, OD10p, ENF00:OD08p) %>% 
  rename_at(vars("OD05", "OD06", "OD07", "OD08"), function(x) paste0(x, "_18"))

kidsOnly18 <- filter(allInd, !is.na(OD10p)) # n = 1216

# combine data on parasito kids at both time points
allKids <- full_join(kidsOnly18, enfants2018, by = c("OD10p", "Rang", "concession")) %>% 
  full_join(., allInfect, by = "OD10p") 

# extract household data to combine with data on parasito kids
hhMerge <- select(household, ConcessionId, MO06:MO06d, MO06ap1:MO06ap2, WC08:WC14, HH22:HH23, HH25:HH26)

# merge allKids and hhMerge by concession
kidsHH <- right_join(allKids, hhMerge, by = "ConcessionId") %>% 
  rename("VillageCode" = "Initial") %>% 
  right_join(., wealth2016, by = "concession") %>% 
  right_join(., hhloc, by = "concession") %>% 
  right_join(., hhloc_snail, by = "concession") %>% 
  right_join(., snail_vars, by = "nearestSnailSite") %>% 
  right_join(., snailV_vars, by = "VillageCode") %>% 
  mutate(nActivities = WC01 + WC02 + WC03 + WC04 + WC05 + WC06,
         anyActivity = ifelse(nActivities > 0, 1, 0),
         siteAreaPeak_d = siteArea_peak/distanceSnail,
         siteAreaYr_d = siteArea_yr/distanceSnail,
         pctHabitatPeak_d = pctHabitat_peak/distanceSnail, 
         habitatAreaPeak_d = habitatArea_peak/distanceSnail,
         habitatAreaYr_d = habitatArea_yr/distanceSnail)


# check for missing village data
#missingVillage <- kidsHH[is.na(kidsHH$village.x),] # n = 97
#missingVillage$OD10p # ~ 10 of those missing have IDs
kidsHH$village <- substr(kidsHH$OD10p, 1, 2) # redefine village variable
missingVillage <- kidsHH[is.na(kidsHH$village),] # n = 87, households with GPS data but nothing else
kidsHH <- kidsHH[!is.na(kidsHH$village),] # n = 1264

write.csv(kidsHH, "data/kidsHHeco.csv")



