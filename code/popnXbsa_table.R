#############################################################
#
#             Population x BSA X Duration Table 
#       Schistosomiasis risk by water contact activity
#                   Dissertation chapter 3
#                   Code by: Andrea Lund
#                   Created: 14 June 2019
#                Last updated: 6 April 2020
#
#############################################################

# import household data set

# import allInd data set

# import allKids data set

summary(household$WC08) # n = 4 missing
summary(household$WC09) # n = 4 missing
summary(household$WC10)
summary(household$WC11)
summary(household$WC12)
summary(household$WC13)
summary(household$WC14)

householdWC <- household %>%
  filter(!is.na(household$WC08))

summary(householdWC$WC08)

# summarize household frequencies of visitation for all activities
hlaundry <- table(householdWC$WC08)
hdishes <- table(householdWC$WC09)
hwater <- table(householdWC$WC10)
hirrig <- table(householdWC$WC11)
hlivestock <- table(householdWC$WC12)
hfishing <- table(householdWC$WC13)
hbathing <- table(householdWC$WC14)

# distinguish canoe fishing from shore fishing
householdWC$canoe <- ifelse(householdWC$HH19 > 0, 1, 0)
table(householdWC$canoe, householdWC$HH19)
fishingType <- table(householdWC$canoe, householdWC$WC13); fishingType

# combine frequencies into single data frame
householdFreq <- rbind(hlaundry, hdishes, hwater, hirrig, hlivestock, hfishing, hbathing)
householdFreq <- as.data.frame(householdFreq)
colnames(householdFreq) <- c("None", "<1x/2wks", "1x/2wks", ">1x/2wks", "1x/day", ">1x/day"); householdFreq
write.csv(householdFreq, file = "householdFreq.csv")

# summarize proportion of population exposed via each activity
summary(allInd$WC01) # n = 0                        10928
summary(allInd$WC02) # n = 1 missing                10927
summary(allInd$WC03) # n = 2 missing                10926
summary(allInd$WC04) # n = 2 missing                10926 
summary(allInd$WC05) # n = 2 missing                10926
summary(allInd$WC06) # n = 1 missing                10927

ilaundry <- table(allInd$WC01)
idishes <- table(allInd$WC02)
iwater <- table(allInd$WC03)
iirrig <- table(allInd$WC04)
ilivestock <- table(allInd$WC05)
ifishing <- table(allInd$WC06)

# combine frequencies into single 
indFreq <- rbind(ilaundry, idishes, iwater, iirrig, ilivestock, ifishing)
indFreq <- as.data.frame(indFreq)
rownames(indFreq) <- c("laundry", "dishes", "water", "irrigation", "livestock", "fishing")
colnames(indFreq) <- c("E-", "E+"); indFreq

# male frequencies
males <- allInd %>% 
  filter(DEM01 == 1) # n = 4711

summary(males$WC01) # n = 5206
summary(males$WC02) # n = 5205 
summary(males$WC03) # n = 5205
summary(males$WC04) # n = 5204
summary(males$WC05) # n = 5205
summary(males$WC06) # n = 5206

mlaundry <- table(males$WC01)
mdishes <- table(males$WC02)
mwater <- table(males$WC03)
mirrig <- table(males$WC04)
mlivestock <- table(males$WC05)
mfishing <- table(males$WC06)

# combine frequencies into single 
mFreq <- rbind(mlaundry, mdishes, mwater, mirrig, mlivestock, mfishing)
mFreq <- as.data.frame(mFreq)
rownames(mFreq) <- c("laundry", "dishes", "water", "irrigation", "livestock", "fishing")
colnames(mFreq) <- c("E-", "E+"); mFreq

# female frequencies
females <- allInd %>% 
  filter(DEM01 == 2) # n = 5197

summary(females$WC01) # n = 5197
summary(females$WC02) # n = 5197  
summary(females$WC03) # n = 5196
summary(females$WC04) # n = 5197 
summary(females$WC05) # n = 5196
summary(females$WC06) # n = 5196

flaundry <- table(females$WC01)
fdishes <- table(females$WC02)
fwater <- table(females$WC03)
firrig <- table(females$WC04)
flivestock <- table(females$WC05)
ffishing <- table(females$WC06)

# combine frequencies into single data frame
fFreq <- rbind(flaundry, fdishes, fwater, firrig, flivestock, ffishing)
fFreq <- as.data.frame(fFreq)
rownames(fFreq) <- c("laundry", "dishes", "water", "irrigation", "livestock", "fishing")
colnames(fFreq) <- c("E-", "E+"); fFreq

# summarize proportion of population exposed via each activity
summary(allKids$WC01) # n = 398 missing              1217
summary(allKids$WC02) # n = 398                      1217
summary(allKids$WC03) # n = 399                      1216
summary(allKids$WC04) # n = 399                      1216
summary(allKids$WC05) # n = 398                      1217
summary(allKids$WC06) # n = 398                      1217

kidsLaundry <- table(allKids$WC01)
kidsDishes <- table(allKids$WC02)
kidsWater <- table(allKids$WC03)
kidsIrrig <- table(allKids$WC04)
kidsLivestock <- table(allKids$WC05)
kidsFishing <- table(allKids$WC06)

# combine frequencies into single 
kidsFreq <- rbind(kidsLaundry, kidsDishes, kidsWater, kidsIrrig, kidsLivestock, kidsFishing)
kidsFreq <- as.data.frame(kidsFreq)
rownames(kidsFreq) <- c("laundry", "dishes", "water", "irrigation", "livestock", "fishing")
colnames(kidsFreq) <- c("E-", "E+"); kidsFreq

# boy frequencies
boys <- allKids %>% 
  filter(DEM01.x == 1) # n = 630

summary(boys$WC01) # n = 630
summary(boys$WC02) # n = 630 
summary(boys$WC03) # n = 629
summary(boys$WC04) # n = 629
summary(boys$WC05) # n = 630
summary(boys$WC06) # n = 630

bLaundry <- table(boys$WC01)
bDishes <- table(boys$WC02)
bWater <- table(boys$WC03)
bIrrig <- table(boys$WC04)
bLivestock <- table(boys$WC05)
bFishing <- table(boys$WC06)

# combine frequencies into single 
bFreq <- rbind(bLaundry, bDishes, bWater, bIrrig, bLivestock, bFishing)
bFreq <- as.data.frame(bFreq)
rownames(bFreq) <- c("laundry", "dishes", "water", "irrigation", "livestock", "fishing")
colnames(bFreq) <- c("E-", "E+"); bFreq

# girl frequencies
girls <- allKids %>% 
  filter(DEM01.x == 2) # n = 587

summary(girls$WC01) # n = 587
summary(girls$WC02) # n = 587 
summary(girls$WC03) # n = 587
summary(girls$WC04) # n = 587
summary(girls$WC05) # n = 587
summary(girls$WC06) # n = 587

gLaundry <- table(girls$WC01)
gDishes <- table(girls$WC02)
gWater <- table(girls$WC03)
gIrrig <- table(girls$WC04)
gLivestock <- table(girls$WC05)
gFishing <- table(girls$WC06)

# combine frequencies into single 
gFreq <- rbind(gLaundry, gDishes, gWater, gIrrig, gLivestock, gFishing)
gFreq <- as.data.frame(gFreq)
rownames(gFreq) <- c("laundry", "dishes", "water", "irrigation", "livestock", "fishing")
colnames(gFreq) <- c("E-", "E+"); gFreq