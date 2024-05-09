
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Subset data based on demographics
#
# Last Update: 09 Jun 2023
#
################################################################################


packages <- c("dplyr", "magrittr", "readr", "lubridate")
for (package in packages) {
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

dat <- read_csv("./screening/data/screen_overall.csv")


# Race --------------------------------------------------------------------

white <- filter(dat, race==1)
write_csv(white, paste0("./screening/data/screen_white.csv"))
rm(white)

black <- filter(dat, race==2)
write_csv(black, paste0("./screening/data/screen_black.csv"))
rm(black)

hisp <- filter(dat, race==3)
write_csv(hisp, paste0("./screening/data/screen_hisp.csv"))
rm(hisp)


# Sex ---------------------------------------------------------------------

male <- filter(dat, sex=="M")
write_csv(male, paste0("./screening/data/screen_male.csv"))
rm(male)

female <- filter(dat, sex=="F")
write_csv(female, paste0("./screening/data/screen_female.csv"))
rm(female)
