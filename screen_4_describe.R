
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Describe data 
#
# Last Update: 26 Apr 2023
#
################################################################################


packages <- c("dplyr", "tidyr", "magrittr", "readr", "lubridate", "parallel", "gmodels")
for (package in packages){
  library(package, character.only=T)
} 

# Program settings
analysis <- "main"


# Read in data ------------------------------------------------------------

if (analysis=="main") {
  dat <- read_csv("./screening/data/screen_overall.csv", 
                  col_types=c("cDDDfdfffffDdddDdDdfdfDdf"))
} else if (analysis=="post2010") {
  dat <- read_csv("./screening/data/screen_post2010.csv", 
                  col_types=c("cDDDfdfffffDdddDdDdfdfDdf"))
}


# Sample data -------------------------------------------------------------

p <- 0.25 # Proportion of unexposed to sample

# Sample at the person level
base <- dat %>% 
  filter(!duplicated(bene_id)) %>% 
  select(bene_id, hiv_base)

set.seed(123)
samp1 <- filter(base, hiv_base==1)
samp0 <- filter(base, hiv_base==0)[rbinom(nrow(filter(base, hiv_base==0)), 1, p)==1, ] 
include <- bind_rows(samp1, samp0)
samp <- filter(dat, bene_id %in% include$bene_id) %>% 
  mutate(wt_s = ifelse(hiv_base==1, 1, 1/p))


# Manage data -------------------------------------------------------------

# Run in parallel to speed up process
nsplits <- 1000
s <- sample(rep(1:nsplits, ceiling(nrow(samp)/nsplits))[1:nrow(samp)])
ids <- unique(samp$bene_id)

set.up <- function(split) {
  include <- ids[s==split]
  new.dat <- filter(samp, bene_id %in% include)
  
  new.dat <- new.dat %>% 
    filter(!is.na(race) & !is.na(sex)) %>% 
    group_by(bene_id) %>% 
    mutate(
      # Set up confounders
      n_cond = factor(ifelse(n_cond==0, 0,
                             ifelse(n_cond==1, 1, 2))),
      n_cond_1 = lag(n_cond),
      n_cond_1 = factor(ifelse(is.na(n_cond_1), n_cond, n_cond_1)),
      n_cond_2 = lag(n_cond, n=2),
      n_cond_2 = factor(ifelse(is.na(n_cond_2), n_cond_1, n_cond_2)),
      enroll_period = factor(case_when(
        yr_base %in% 2001:2005 ~ 1,
        yr_base %in% 2006:2010 ~ 2,
        T ~ 3)),
      
      # Screening protocol
      screen = case_when(
        # Censor if it's been more than X years since last endoscopy or baseline
        years_between_endo>10 ~ 0, 
        # Censor if it's been less than X years between endoscopies
        endo==1 & prior_endo==1 & years_between_endo<10 ~ 0, 
        T ~ 1),
      screen_1 = lag(screen, default=1),
      screen_2 = lag(screen, n=2, default=1),
      cum_screen = cumprod(screen),
      
      # Set up outcome
      event = ifelse(int!=event_int, 0, 
                     ifelse(first_event=="censor", 0,
                            ifelse(first_event=="colon", 1, 2))),
      event = factor(event)) %>% 
    ungroup()
  
  return(new.dat)
}

cores <- detectCores()
all.dat <- mclapply(1:nsplits, function(x){set.up(split=x)}, mc.cores=cores/4, mc.set.seed=F)
dat <- bind_rows(all.dat)

base <- filter(dat, !duplicated(bene_id))


# Describe data -----------------------------------------------------------

# Descriptive statistics
  # Overall
  summary(base$age_start)
  CrossTable(x=base$sex)
  CrossTable(x=base$race)
  CrossTable(x=base$enroll_period)
  CrossTable(x=base$first_elig_state)
  CrossTable(base$n_cond)

  # By HIV status
  CrossTable(x=base$hiv_base)
  summary(base$age_start[base$hiv_base==1])
  summary(base$age_start[base$hiv_base==0])
  CrossTable(x=base$sex, y=base$hiv_base, prop.r=F, prop.t=F, prop.chisq=F)
  CrossTable(x=base$race, y=base$hiv_base, prop.r=F, prop.t=F, prop.chisq=F)
  CrossTable(x=base$enroll_period, y=base$hiv_base, prop.r=F, prop.t=F, prop.chisq=F)
  CrossTable(x=base$first_elig_state, y=base$hiv_base, prop.r=F, prop.t=F, prop.chisq=F)
  CrossTable(x=base$n_cond, y=base$hiv_base, prop.r=F, prop.t=F, prop.chisq=F)
  
# Describe follow-up
base %>% 
  mutate(length = (baseline %--% end_date)/months(1)) %>% 
  summarize(tot_length = sum(length),
            med_length = median(length))
CrossTable(x=base$first_event)

  
# When did people stop following protocol?
protocol <- dat %>% 
  group_by(bene_id) %>% 
  mutate(age_int = floor(age_start),
         deviate = as.numeric(screen==0),
         cum_dev = cumsum(cumsum(deviate))) %>% 
  filter(cum_dev<=1) %>% 
  ungroup()
protocol_hiv <- filter(protocol, hiv_base==1)
protocol_nohiv <- filter(protocol, hiv_base==0)

CrossTable(x=protocol_hiv$age_int, y=protocol_hiv$deviate)
CrossTable(x=protocol_nohiv$age_int, y=protocol_nohiv$deviate)
