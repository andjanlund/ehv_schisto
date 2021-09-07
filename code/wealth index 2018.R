###########################################################################
#   generating wealth index from asset data from 2018 survey
#
#   use principal componenets analysis to reduce many asset variables  
#   to a single index of wealth for each household
#
#   Code by: Andrea Lund
#   Last updated: April 23, 2020
##########################################################################

# install packages
library(readxl)
library(tidyverse)
library(lubridate)
library(corrplot)

# set working directory
setwd(dir = "C:/Users/andre/Box Sync/Schisto (alund2@stanford.edu)/Data")

#### import data #####
# household level survey data
household <- read_excel("./Andrea's Field Data/2018 Survey Data/household.xlsx")

# livestock data from survey
livestock <- read_excel("./Andrea's Field Data/2018 Survey Data/S42.xlsx")

# establish concession key
concessionKey <- livestock %>% 
  group_by(ConcessionId) %>% 
  select(ConcessionId, concession)

#### summarize data ####
# unique concession values in different data sets
unique(factor(household$ConcessionId)) # 681
unique(factor(livestock$concession)) # 624
unique(factor(livestock$ConcessionId)) # 624

# review column names
colnames(household)
colnames(livestock)

# review unique values
unique(livestock$Pl)

#### extract and combine relevant data #####
### extract asset data
assets <- household %>% 
  select(ConcessionId, HH01:HH30) %>% 
  full_join(concessionKey, ., by = "ConcessionId") %>% 
  distinct()

## dichotomize categorical variables (for H22:H29, spread categories in columns)
# H22 - drinking water source
ddichot <- assets %>% 
  mutate(v = 1, drink = HH22) %>% 
  spread(drink, v, fill = 0) %>% 
  select(ConcessionId, c("1", "2", "4", "5")) %>% 
  rename("dpiped" = "1", "dforage" = "2", "dwell" = "4", "dsurface" = "5")

# H23 - laundry water source
ldichot <- assets %>%
  mutate(v = 1, linge = HH23) %>% 
  spread(linge, v, fill = 0) %>% 
  select(ConcessionId, c("1", "2", "3", "4", "5", "6")) %>% 
  rename("lpiped" = "1", "lforage" = "2", "lwellp" = "3", "lwellu" = "4", "lsurface" = "5", "lother" = "6")

# HH25 - toilet type
assets$toilet <- ifelse(assets$HH25 == 1, "ttoilet",
                        ifelse(assets$HH25 %in% c(2,3), "tlatrine", "tother"))
tdichot <- assets %>%
  mutate(v = 1, toilet = toilet) %>% 
  spread(toilet, v, fill = 0) %>% 
  select(ConcessionId, "tlatrine", "ttoilet", "tother")

# HH27 - floor material
fdichot <- assets %>% 
  mutate(v = 1, floor = HH27) %>% 
  spread(floor, v, fill = 0) %>% 
  select(ConcessionId, c("1", "2", "3", "4")) %>% 
  rename("fdirt" = "1", "fbanco" = "2", "fcement" = "3", "ftile" = "4")

# HH28 - roof material
rdichot <- assets %>% 
  mutate(v = 1, roof = HH28) %>% 
  spread(roof, v, fill = 0) %>% 
  select(ConcessionId, "1", "2", "3", "4", "5") %>% 
  rename("rstraw" = "1", "rwood" = "2", "rmetal" = "3", "recement" = "4", "rtile" = "5")

# HH29 - material of external walls
assets$wallmat <- ifelse(assets$HH29 == 1, "wcement",
                         ifelse(assets$HH29 == 3, "wstraw", "wother"))
wdichot <- assets %>%
  mutate(v = 1, wall = wallmat) %>% 
  spread(wall, v, fill = 0) %>% 
  select(ConcessionId, wcement, wstraw, wother)

# combine all dichotomized data by household ID
catvars <- full_join(concessionKey, ddichot, by = "ConcessionId") %>% 
           full_join(., ldichot, by = "ConcessionId") %>% 
           full_join(., tdichot, by = "ConcessionId") %>% 
           full_join(., fdichot, by = "ConcessionId") %>% 
           full_join(., rdichot, by = "ConcessionId") %>% 
           full_join(., wdichot, by = "ConcessionId") %>% 
           distinct()

## extract livestock data
cows <- livestock %>% 
  group_by(ConcessionId) %>%
  filter(Pl %in% c("Boeuf", "vache", "Vache")) %>% 
  select(concession, Pl05) %>% 
  rename("nCows" = "Pl05")

sheep <- livestock %>% 
  group_by(ConcessionId) %>% 
  filter(Pl == "Mouton") %>% 
  select(Pl05) %>% 
  rename("nSheep" = "Pl05")

donkey <- livestock %>% 
  group_by(ConcessionId) %>% 
  filter(Pl == "Ane") %>% 
  select(Pl05) %>% 
  rename("nDonkeys" = "Pl05")

goat <- livestock %>% 
  group_by(ConcessionId) %>% 
  filter(Pl %in% c("chevre", "Chèvre")) %>% 
  select(Pl05) %>% 
  rename("nGoats" = "Pl05")

# combine livestock data by household ID
animals <- full_join(concessionKey, cows, by = "ConcessionId", fill = 0) %>% 
           full_join(., sheep, by = "ConcessionId", fill = 0) %>% 
           full_join(., donkey, by = "ConcessionId", fill = 0) %>% 
           full_join(., goat, by = "ConcessionId", fill = 0) %>% 
           select(-concession.y) %>% 
           rename("concession" = "concession.x") %>% 
           distinct()

# recode missing values as 0        
animals[is.na(animals)] <- 0    

## merge assets, dichotomized categorical variables and livestock data
allAssets <- full_join(assets, catvars, by = c("concession", "ConcessionId")) %>% 
             full_join(., animals, by = c("concession", "ConcessionId")) %>% 
             select(-VillageId, -toilet, -wallmat, -HH22:-HH29) %>% 
             distinct()

# remove 57 observations with missing livestock data 
allAssets <- allAssets[!is.na(allAssets$nCows),]

# remove 1 observation with missing data for fan
allAssets <- allAssets[!is.na(allAssets$fan),]

# name columns
colnames(allAssets)[3:21] <- c("elec", "radio", "tv", "cable", "lave", "refri",
                              "stove", "fan", "computer", "internet", "mobile",
                              "satdish", "bici", "moto", "car", "cart",
                              "plow", "canoe", "pump")

allAssets <- as.numeric(allAssets[,3:50])

#### calculate descriptive statistics for all asset variables #### 
mean_sd <- summarise_if(allAssets, is.numeric, funs(mean, sd))
sds <- summarise_if(allAssets, is.numeric, funs(sd))


# how much does sd of asset variables vary? if it varies by orders of magnitude, need to scale in pca
summary(mean_sd$sd) # min = 0.03887, mean = 0.96744, median = 0.44450, max = 9.64724 - scaling seems necessary
# remove variables with very low variability (SD lt 0.2)
# ... dforage, dwell, fbanco, lave, lforage, lother, lwellp, lwellu, rwood, tother, wother

#### visualize asset data in corrplot ####
M <- cor(allAssets[,3:50])
# visualize correlation between asset variables
corrplot(M, method = "circle", type = "upper")
corrplot(M, method = "circle", type = "upper", order = "FPC") # order by first principal component

# run principal components analysis on all assets
pca <- prcomp(allAssets[,c(3:21, 46:50)], scale = T)
summary(pca) # 27.9% of variance accounted for in first component

names(pca) # center = mean of the variables; sdev = sd of the variables
pca$center # the centering and scaling used (same as means calculated above)
pca$sdev # standard deviations of the principal components

plot(pca) # scree plot

# process and export loading data for index creation
loadings <- as.data.frame(pca$rotation[,1]) # variable loadings for first principal component. 
## ^ these loadings ^ are very counterintuitive 
# (e.g. positive values for surface water use, negative values for tv/computer ownership)
# exactly the opposite of what you'd expect, and opposite of what we see in 2016 survey data
# counterintuitive loadings persisted regardless of subset of asset variables used
# (e.g. assets only vs assets + dichotomized vars vs assets + dichotomized variables + livestock)

loadings2 <- loadings %>%
  add_rownames() %>%
  rename(asset = rowname, fpc = `pca$rotation[, 1]`)
# join variable loadings to md_SD table
assetDescriptives <- full_join(mean_sd, loadings2, by = "asset") %>% 
  write_csv("asset descriptives.csv")

# extract the first column for loadings of the first principal component to calculate wealth scores for each household
biplot(pca)
boxplot(pca$x)
pca$x

# pca just durable assets, not water, hhsize or construction materials
pca_redux <- prcomp(fullSurvey[,64:83], scale = T)
summary(pca_redux) # 25.96% of variation explained by first principal component
plot(pca_redux)
pca_redux$rotation

# pca excluding assets with low variation
pca_hisd <- prcomp(fullSurvey[,c(51:54,56:74,77:78,82,84:85,87,89:91,93:97,99:101)], scale = T)
summary(pca_hisd) # 15.67% of variation explained by first principal component
plot(pca_hisd)
pca_hisd$rotation[,1]
