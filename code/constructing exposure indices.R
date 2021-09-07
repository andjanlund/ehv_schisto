#############################################################
#
#                 Creation of exposure indices 
#       Schistosomiasis risk by water contact activity
#                   Dissertation chapter 3
#                   Code by: Andrea Lund
#                   Created: 14 June 2019
#                Last updated: 7 September 2021
#
#############################################################

library(tidyr)
library(GGally)  
library(here)

# import data
kidsHH <- read.csv("data/kidsHHeco.csv") # n = 1264

###################################################
#
#         create new frequency variables
#

# recode survey reponse categories of WC08-WC14:
# these variables represent household-level frequency of water contact
# but each category represents a range, so converted to avg, min and max 
# code resource: http://www.cookbook-r.com/Manipulating_data/Recoding_data/
oldvalues <- c(0:5) # numbers the represent categories in the survey 
wkavgvals <- c(0, 0, 1, 4, 7, 14) # converted to approx weekly frequency
wkminvals <- c(0, 0, 1, 3, 7, 10) # converted to minimum weekly frequency
wkmaxvals <- c(0, 0, 1, 5, 7, 21) # converted to maximum weekly frequency
moavgvals <- c(0, 2, 4, 14, 30, 90) # converted to approx monthly frequency
mominvals <- c(0, 1, 4, 8, 30, 60) # converted to min monthly frequency
momaxvals <- c(1,3, 4, 21, 90, 120) # converted to max monthly frequency

### weekly frequencies by activity
# average
kidsHH$WC08wavg <- wkavgvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09wavg <- wkavgvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10wavg <- wkavgvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11wavg <- wkavgvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12wavg <- wkavgvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13wavg <- wkavgvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC14wavg <- wkavgvals[match(kidsHH$WC14, oldvalues)]

# minimum
kidsHH$WC08wmin <- wkminvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09wmin <- wkminvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10wmin <- wkminvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11wmin <- wkminvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12wmin <- wkminvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13wmin <- wkminvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC14wmin <- wkminvals[match(kidsHH$WC14, oldvalues)]

# maximum
kidsHH$WC08wmax <- wkmaxvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09wmax <- wkmaxvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10wmax <- wkmaxvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11wmax <- wkmaxvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12wmax <- wkmaxvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13wmax <- wkmaxvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC14wmax <- wkmaxvals[match(kidsHH$WC14, oldvalues)]

### monthly frequencies by activity
# average
kidsHH$WC08mavg <- moavgvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09mavg <- moavgvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10mavg <- moavgvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11mavg <- moavgvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12mavg <- moavgvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13mavg <- moavgvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC13mavg <- moavgvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC14mavg <- moavgvals[match(kidsHH$WC14, oldvalues)]

# minimum
kidsHH$WC08mmin <- mominvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09mmin <- mominvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10mmin <- mominvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11mmin <- mominvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12mmin <- mominvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13mmin <- mominvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC14mmin <- mominvals[match(kidsHH$WC14, oldvalues)]

# maximum
kidsHH$WC08mmax <- momaxvals[match(kidsHH$WC08, oldvalues)]
kidsHH$WC09mmax <- momaxvals[match(kidsHH$WC09, oldvalues)]
kidsHH$WC10mmax <- momaxvals[match(kidsHH$WC10, oldvalues)]
kidsHH$WC11mmax <- momaxvals[match(kidsHH$WC11, oldvalues)]
kidsHH$WC12mmax <- momaxvals[match(kidsHH$WC12, oldvalues)]
kidsHH$WC13mmax <- momaxvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC13mmax <- momaxvals[match(kidsHH$WC13, oldvalues)]
kidsHH$WC14mmax <- momaxvals[match(kidsHH$WC14, oldvalues)]


###################################################
#
#         create duration and bsa variables
#
kidsHH$bathing <- 1 # assume all kids have bathed/swam in surface water in last two weeks
kidsBSA = 1.130 # sq meters of body surface area for children under 14, derived from Mosteller 1987 and cited in Sudat et al 2010
duration = c(13.3, 13.3, 2.0, 7.9, 9.3, 23.8, 11.5) # derived from Sow et al 2011
pctBSA = c(0.358, 0.318, 0.453, 0.225, 0.415, 0.476, 1.0) # derived from BSA interviews - ask Andrea for raw data, if needed
BSA <- kidsBSA*pctBSA # calculate estimated m2 exposed per activity for children

###################################################
#
#         calculate exposure indices
#

kidsHH <- kidsHH %>% 
  dplyr::mutate(F_MoAvg = WC01*WC08mavg +  # monthly average frequency
                   WC02*WC09mavg +
                   WC03*WC10mavg + 
                   WC04*WC11mavg + 
                   WC05*WC12mavg +
                   WC06*WC13mavg + 
                   bathing*WC14mavg,
         F_MoMin = WC01*WC08mmin + # monthly minimum frequency
                   WC02*WC09mmin +
                   WC03*WC10mmin + 
                   WC04*WC11mmin + 
                   WC05*WC12mmin +
                   WC06*WC13mmin +
                   bathing*WC14mmin,
         F_MoMax = WC01*WC08mmax + # monthly maximum frequency
                   WC02*WC09mmax +
                   WC03*WC10mmax + 
                   WC04*WC11mmax + 
                   WC05*WC12mmax +
                   WC06*WC13mmax +
                   bathing*WC14mmax,
         FD_MoAvg = WC01*WC08mavg*duration[1] +  # monthly average time exposed
                    WC02*WC09mavg*duration[2] +
                    WC03*WC10mavg*duration[3] + 
                    WC04*WC11mavg*duration[4] + 
                    WC05*WC12mavg*duration[5] +
                    WC06*WC13mavg*duration[6] +
                    bathing*WC14mavg*duration[7],
         FD_MoMin = WC01*WC08mmin*duration[1] + # monthly minimum time exposed
                    WC02*WC09mmin*duration[2] +
                    WC03*WC10mmin*duration[3] + 
                    WC04*WC11mmin*duration[4] + 
                    WC05*WC12mmin*duration[5] +
                    WC06*WC13mmin*duration[6] +
                    bathing*WC14mmin*duration[7],
         FD_MoMax = WC01*WC08mmax*duration[1] + # monthly maximum time exposed
                    WC02*WC09mmax*duration[2] +
                    WC03*WC10mmax*duration[3] + 
                    WC04*WC11mmax*duration[4] + 
                    WC05*WC12mmax*duration[5] +
                    WC06*WC13mmax*duration[6] +
                    bathing*WC14mmax*duration[7],
         FDB_MoAvg = WC01*WC08mavg*duration[1]*BSA[1] + # monthly average bsa-adjusted time exposed
                       WC02*WC09mavg*duration[2]*BSA[2] +
                       WC03*WC10mavg*duration[3]*BSA[3] + 
                       WC04*WC11mavg*duration[4]*BSA[4] + 
                       WC05*WC12mavg*duration[5]*BSA[5] +
                       WC06*WC13mavg*duration[6]*BSA[6] +
                       bathing*WC14mavg*duration[7]*BSA[7], 
         FDB_MoMin = WC01*WC08mmin*duration[1]*BSA[1] + # monthly average bsa-adjusted time exposed
                       WC02*WC09mmin*duration[2]*BSA[2] +
                       WC03*WC10mmin*duration[3]*BSA[3] + 
                       WC04*WC11mmin*duration[4]*BSA[4] + 
                       WC05*WC12mmin*duration[5]*BSA[5] +
                       WC06*WC13mmin*duration[6]*BSA[6] +
                       bathing*WC14mmin*duration[7]*BSA[7],
         FDB_MoMax = WC01*WC08mmax*duration[1]*BSA[1] + # monthly average bsa-adjusted time exposed
                       WC02*WC09mmax*duration[2]*BSA[2] +
                       WC03*WC10mmax*duration[3]*BSA[3] + 
                       WC04*WC11mmax*duration[4]*BSA[4] + 
                       WC05*WC12mmax*duration[5]*BSA[5] +
                       WC06*WC13mmax*duration[6]*BSA[6] +
                       bathing*WC14mmax*duration[7]*BSA[7],
         F_WkAvg = WC01*WC08wavg +  # weekly average frequency
                   WC02*WC09wavg +
                   WC03*WC10wavg + 
                   WC04*WC11wavg + 
                   WC05*WC12wavg +
                   WC06*WC13wavg +
                   bathing*WC14wavg,
         F_WkMin = WC01*WC08wmin + # weekly minimum frequency
                   WC02*WC09wmin +
                   WC03*WC10wmin + 
                   WC04*WC11wmin + 
                   WC05*WC12wmin +
                   WC06*WC13wmin +
                   bathing*WC14wmin,
         F_WkMax = WC01*WC08wmax + # weekly maximum frequency
                   WC02*WC09wmax +
                   WC03*WC10wmax + 
                   WC04*WC11wmax + 
                   WC05*WC12wmax +
                   WC06*WC13wmax +
                   bathing*WC14wmax,
         FD_WkAvg = WC01*WC08wavg*duration[1] +  # weekly average time exposed
                    WC02*WC09wavg*duration[2] +
                    WC03*WC10wavg*duration[3] + 
                    WC04*WC11wavg*duration[4] + 
                    WC05*WC12wavg*duration[5] +
                    WC06*WC13wavg*duration[6] +
                    bathing*WC14wavg*duration[7],
         FD_WkMin = WC01*WC08wmin*duration[1] + # weekly minimum time exposed
                    WC02*WC09wmin*duration[2] +
                    WC03*WC10wmin*duration[3] + 
                    WC04*WC11wmin*duration[4] + 
                    WC05*WC12wmin*duration[5] +
                    WC06*WC13wmin*duration[6] +
                    bathing*WC14wmin*duration[7],
         FD_WkMax = WC01*WC08wmax*duration[1] + # weekly maximum time exposed
                     WC02*WC09wmax*duration[2] +
                     WC03*WC10wmax*duration[3] + 
                     WC04*WC11wmax*duration[4] + 
                     WC05*WC12wmax*duration[5] +
                     WC06*WC13wmax*duration[6] +
                     bathing*WC14wmax*duration[7],
         FDB_WkAvg = WC01*WC08wavg*duration[1]*BSA[1] + # weekly average bsa-adjusted time exposed
                       WC02*WC09wavg*duration[2]*BSA[2] +
                       WC03*WC10wavg*duration[3]*BSA[3] + 
                       WC04*WC11wavg*duration[4]*BSA[4] + 
                       WC05*WC12wavg*duration[5]*BSA[5] +
                       WC06*WC13wavg*duration[6]*BSA[6] +
                       bathing*WC14wavg*duration[7]*BSA[7], 
         FDB_WkMin = WC01*WC08wmin*duration[1]*BSA[1] + # weekly average bsa-adjusted time exposed
                       WC02*WC09wmin*duration[2]*BSA[2] +
                       WC03*WC10wmin*duration[3]*BSA[3] + 
                       WC04*WC11wmin*duration[4]*BSA[4] + 
                       WC05*WC12wmin*duration[5]*BSA[5] +
                       WC06*WC13wmin*duration[6]*BSA[6] +
                       bathing*WC14wmin*duration[7]*BSA[7],
         FDB_WkMax = WC01*WC08wmax*duration[1]*BSA[1] + # weekly average bsa-adjusted time exposed
                       WC02*WC09wmax*duration[2]*BSA[2] +
                       WC03*WC10wmax*duration[3]*BSA[3] + 
                       WC04*WC11wmax*duration[4]*BSA[4] + 
                       WC05*WC12wmax*duration[5]*BSA[5] +
                       WC06*WC13wmax*duration[6]*BSA[6] +
                       bathing*WC14wmax*duration[7]*BSA[7],
         # activity-specific weekly frequencies 
         F_laundryWk = WC01*WC08wavg, 
         F_dishesWk = WC02*WC09wavg,
         F_collectWk = WC03*WC10wavg,
         F_irrigWk = WC04*WC11wavg,
         F_livestockWk = WC05*WC12wavg, 
         F_fishingWk = WC06*WC13wavg,
         F_bathingWk = bathing*WC14wavg,
         # activity-specific minimum frequencies 
         F_laundryWk_min = WC01*WC08wmin, 
         F_dishesWk_min = WC02*WC09wmin,
         F_collectWk_min = WC03*WC10wmin,
         F_irrigWk_min = WC04*WC11wmin,
         F_livestockWk_min = WC05*WC12wmin, 
         F_fishingWk_min = WC06*WC13wmin,
         F_bathingWk_min = bathing*WC14wmin,
         # activity-specific maximum frequencies 
         F_laundryWk_max = WC01*WC08wmax, 
         F_dishesWk_max = WC02*WC09wmax,
         F_collectWk_max = WC03*WC10wmax,
         F_irrigWk_max = WC04*WC11wmax,
         F_livestockWk_max = WC05*WC12wmax, 
         F_fishingWk_max = WC06*WC13wmax,
         F_bathingWk_max = bathing*WC14wmax,
         # activity-specific weekly total duration 
         FD_laundryWk = WC01*WC08wavg*duration[1], 
         FD_dishesWk = WC02*WC09wavg*duration[2],
         FD_collectWk = WC03*WC10wavg*duration[3],
         FD_irrigWk = WC04*WC11wavg*duration[4],
         FD_livestockWk = WC05*WC12wavg*duration[5], 
         FD_fishingWk = WC06*WC13wavg*duration[6],
         FD_bathingWk = bathing*WC14wavg*duration[7],
         # bsa-adjusted activity-specific weekly total duration 
         FDB_laundryWk = WC01*WC08wavg*duration[1]*BSA[1], 
         FDB_dishesWk = WC02*WC09wavg*duration[2]*BSA[2],
         FDB_collectWk = WC03*WC10wavg*duration[3]*BSA[3],
         FDB_irrigWk = WC04*WC11wavg*duration[4]*BSA[4],
         FDB_livestockWk = WC05*WC12wavg*duration[5]*BSA[5], 
         FDB_fishingWk = WC06*WC13wavg*duration[6]*BSA[6],
         FDB_bathingWk = bathing*WC14wavg*duration[7]*BSA[7],
         # activity-specific monthly frequnecy
         F_laundryMo = WC01*WC08mavg,  
         F_dishesMo = WC02*WC09mavg,
         F_collectMo = WC03*WC10mavg,
         F_irrigMo = WC04*WC11mavg,
         F_livestockMo = WC05*WC12mavg, 
         F_fishingMo = WC06*WC13mavg, 
         F_bathingMo = bathing*WC14mavg,
         # activity-specific monthly total duration
         FD_laundryMo = WC01*WC08mavg*duration[1],  
         FD_dishesMo = WC02*WC09mavg*duration[2],
         FD_collectMo = WC03*WC10mavg*duration[3],
         FD_irrigMo = WC04*WC11mavg*duration[4],
         FD_livestockMo = WC05*WC12mavg*duration[5], 
         FD_fishingMo = WC06*WC13mavg*duration[6], 
         FD_bathingMo = bathing*WC14mavg*duration[7],
         # bsa-adjusted activity-specific monthly total duration 
         FDB_laundryMo = WC01*WC08mavg*duration[1]*BSA[1], 
         FDB_dishesMo = WC02*WC09mavg*duration[2]*BSA[2],
         FDB_collectMo = WC03*WC10mavg*duration[3]*BSA[3],
         FDB_irrigMo = WC04*WC11mavg*duration[4]*BSA[4],
         FDB_livestockMo = WC05*WC12mavg*duration[5]*BSA[5], 
         FDB_fishingMo = WC06*WC13mavg*duration[6]*BSA[6],
         FDB_bathingMo = bathing*WC14mavg*duration[7]*BSA[7])

