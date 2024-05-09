
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Wrangle results into tables
#
# Last Update: 28 Aug 2023
#
################################################################################


library("tidyverse")

analysis <- "main" # {main, post2010, related6mo, colonoscopy}
weights <- "10yr" # {2yr, 4yr, 10yr, hivwt}


# Read in results ---------------------------------------------------------

if (analysis=="main") {
  res <- read_csv(paste0("../results/screen_res_", weights, ".csv"))
} else {
  res <- read_csv(paste0("../results/screen_res_", analysis, "-", weights, ".csv")) 
}


# Summarize weights -------------------------------------------------------

res %>% 
  summarize(mean_avg_wt = mean(avg_wt, na.rm=T),
            min_avg_wt = min(avg_wt, na.rm=T),
            max_avg_wt = max(avg_wt, na.rm=T),
            mean_min_wt = mean(min_wt, na.rm=T),
            min_min_wt = min(min_wt, na.rm=T),
            max_min_wt = max(min_wt, na.rm=T),
            mean_max_wt = mean(max_wt, na.rm=T),
            min_max_wt = min(max_wt, na.rm=T),
            max_max_wt = max(max_wt, na.rm=T))


# Functions to process results --------------------------------------------

risk_age <- function(dat, ages) {
  
  by_age <- function(t) {
    temp <- filter(dat, floor(time)==t) %>% 
      mutate(hiv_risk = round(risk_colon1*100, 1),
             hiv_CI = paste0("(", format(round(risk_colon1_ll*100, 1), nsmall=1), ", ",
                             format(round(risk_colon1_ul*100, 1), nsmall=1),")"),

             nohiv_risk = round(risk_colon0*100, 1),
             nohiv_CI = paste0("(", format(round(risk_colon0_ll*100, 1), nsmall=1), ", ",
                               format(round(risk_colon0_ul*100, 1), nsmall=1),")"),

             rd = round(rd_colon*100, 1),
             rd_CI = paste0("(", format(round(rd_colon_ll*100, 1), nsmall=1), ", ",
                            format(round(rd_colon_ul*100, 1), nsmall=1),")"),
             rr = round(rr_colon, 2),
             rr_CI = paste0("(", format(round(rr_colon_ll, 2), nsmall=2), ", ",
                            format(round(rr_colon_ul, 2), nsmall=2),")"),
             age = t) %>% 
      filter(!duplicated(age, fromLast=T)) %>% 
      select(age, hiv_risk, hiv_CI, nohiv_risk, nohiv_CI, rd, rd_CI, rr, rr_CI)
  }
  new.dat <- bind_rows(lapply(ages, function(x){by_age(x)}))
  
  return(new.dat)
}

death_age <- function(dat, ages) {
  
  by_age <- function(t) {
    temp <- filter(dat, floor(time)==t) %>% 
      mutate(hiv_risk = round(risk_death1*100, 1),
             hiv_CI = paste0("(", format(round(risk_death1_ll*100, 1), nsmall=1), ", ",
                             format(round(risk_death1_ul*100, 1), nsmall=1),")"),
             
             nohiv_risk = round(risk_death0*100, 1),
             nohiv_CI = paste0("(", format(round(risk_death0_ll*100, 1), nsmall=1), ", ",
                               format(round(risk_death0_ul*100, 1), nsmall=1),")"),
             
             rd = round(rd_death*100, 1),
             rd_CI = paste0("(", format(round(rd_death_ll*100, 1), nsmall=1), ", ",
                            format(round(rd_death_ul*100, 1), nsmall=1),")"),
             rr = round(rr_death, 1),
             rr_CI = paste0("(", format(round(rr_death_ll, 1), nsmall=1), ", ",
                            format(round(rr_death_ul, 1), nsmall=1),")"),
             age = t) %>% 
      filter(!duplicated(age, fromLast=T)) %>% 
      select(age, hiv_risk, hiv_CI, nohiv_risk, nohiv_CI, rd, rd_CI, rr, rr_CI)
  }
  new.dat <- bind_rows(lapply(ages, function(x){by_age(x)}))
  
  return(new.dat)
}


# Weighted analysis -------------------------------------------------------

# Colon cancer
ages <- risk_age(res, c(55, 60, 65))
  
if (analysis=="main") {
  name <- paste0("../tables/table_colon_", weights, ".csv")
} else {
  name <- paste0("../tables/table_colon_", analysis, "-", weights, ".csv")
}
write_csv(ages, name)

# Death
ages <- death_age(res, c(55, 60, 65))  

if (analysis=="main") {
  name <- paste0("../tables/table_death_", weights, ".csv")
} else {
  name <- paste0("../tables/table_death_", analysis, "-", weights, ".csv")
}
write_csv(ages, name)
  

# Crude analysis ----------------------------------------------------------

# Colon cancer  
for (age in c(55, 60, 65)) {

  res <- read_csv("../results/screen_res_crude.csv") %>% 
    filter(state==1)
  
  hiv <- res %>% 
    filter(strata=="hiv_base=1",
           floor(time)==age) %>% 
    mutate(age=floor(time)) %>% 
    filter(!duplicated(age, fromLast=T))
  
  hiv0 <- res %>% 
    filter(strata=="hiv_base=0",
           floor(time)==age) %>% 
    mutate(age=floor(time)) %>% 
    filter(!duplicated(age, fromLast=T))
  
  if (age==55) {
    crude.res <- bind_rows(hiv, hiv0)
  } else {
    crude.res <- bind_rows(crude.res, hiv, hiv0)
  }
  
}

write_csv(crude.res, "../tables/table_colon_crude.csv")


# Death
for (age in c(55, 60, 65)) {
  
  res <- read_csv("../results/screen_res_crude.csv") %>% 
    filter(state==2)
  
  hiv <- res %>% 
    filter(strata=="hiv_base=1",
           floor(time)==age) %>% 
    mutate(age=floor(time)) %>% 
    filter(!duplicated(age, fromLast=T))
  
  hiv0 <- res %>% 
    filter(strata=="hiv_base=0",
           floor(time)==age) %>% 
    mutate(age=floor(time)) %>% 
    filter(!duplicated(age, fromLast=T))
  
  if (age==55) {
    crude.res <- bind_rows(hiv, hiv0)
  } else {
    crude.res <- bind_rows(crude.res, hiv, hiv0)
  }
  
}

write_csv(crude.res, "../tables/table_death_crude.csv")

