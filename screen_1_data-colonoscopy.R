
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Manage raw data 
#
# Last Update: 10 Jul 2023
#
################################################################################


packages <- c("dplyr", "tidyr", "magrittr", "readr", "lubridate", "parallel", "zoo")
for (package in packages){
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

# Individual data
dat <- read_csv(file="../../data/MedicanDataRequest/jackie/20230622_jackie_crc_v4/crc_file1.csv", 
              col_types=paste0("cfffffcc", paste(rep("D", 23), collapse=""))) %>% 
  rename(start_date = first_elig_period_start_date_65,
         end_date = first_elig_end_date_65)

# Endoscopy data
proc <- read_csv(file="../../data/MedicanDataRequest/jackie/20230622_jackie_crc_v4/crc_file2.csv",
                 col_types=c("cfD")) %>% 
  filter(proc_type==1) # Only colonoscopies


# Sample split for efficiency ---------------------------------------------

set.seed(123)
nsplits <- 1000
s <- sample(rep(1:nsplits, ceiling(nrow(dat)/nsplits))[1:nrow(dat)])
ids <- unique(dat$bene_id)

set.up <- function(split) {
  include <- ids[s==split]
  new.dat <- filter(dat, bene_id %in% include)
  

# Manage individual data --------------------------------------------------

dat2 <- new.dat %>% 
  mutate(baseline = start_date + months(6),
         
         # Recode leap years birthdays
         dob = ifelse(month(dob)==2 & day(dob)==29, dob - days(1), dob),
         dob = as.Date(dob, origin=ymd("1970-01-01")),
         
         # Find 50th birthday and move baseline
         date50 = dob + years(50),
         baseline = ifelse(date50>baseline, date50, baseline),
         baseline = as.Date(baseline, origin=ymd("1970-01-01")),
         
         # Censor on Sept 30, 2015 (UNTIL ICD-10 FULLY IMPLEMENTED)
	       end_date = ifelse(end_date>ymd("2015-09-30"), ymd("2015-09-30"), end_date),
         end_date = as.Date(end_date, origin=ymd("1970-01-01")),

         # HIV status at baseline
         hiv_base = ifelse(is.na(hiv_cms_date_type2), 0, 
                           as.numeric(hiv_cms_date_type2<=baseline))) %>% 
  # Remove ineligible beneficiaries:
    # Died, had cancer, or dropped out before baseline
  filter(is.na(dod) | dod>baseline,
         is.na(any_cx_data_type1) | any_cx_data_type1>baseline,
         end_date>baseline)

# Determine first event
first <- as_date(apply(dat2[, c("end_date", "colon_cx_data_type2", "any_cx_data_type2", "dod")],
                         1, min, na.rm=T))
dat2 <- bind_cols(dat2, "event_date"=first)
event_type <- rep(NA, nrow(dat2))
event_type[dat2$event_date == dat2$end_date] <- "censor"
event_type[dat2$event_date == dat2$any_cx_data_type2] <- "censor" # Censor at other cancers
event_type[dat2$event_date == dat2$dod] <- "dead"
event_type[dat2$event_date == dat2$colon_cx_data_type2] <- "colon"
dat2 <- bind_cols(dat2, "first_event"=event_type)

# Calculate times  
dat3 <- dat2 %>% 
  mutate(age_base = floor((dob %--% baseline)/years(1)),
         yr_base = year(baseline),
         
         # Time to event
         event_time = (baseline %--% event_date)/months(3),
         event_int = ceiling(event_time),
         # If event_int is an integer, then event date occurs in next interval
         event_int = ifelse(event_time%%1==0, event_int+1, event_int),
         
         # If they had anal cancer w/in 60 days of colon cancer, censor
         days_anal = (colon_cx_data_type2 %--% anal_cx_data_type2)/days(1),
         anal = as.numeric(days_anal>=0 & days_anal<=60 & !is.na(days_anal)),
         first_event = ifelse(first_event=="colon" & anal==1, "censor", first_event),

         # Charlson comorbidites
         t_mi = ceiling((baseline %--% myocard_infarct_date)/months(3)),
         t_chf = ceiling((baseline %--% cong_heart_fail_date)/months(3)),
         t_pvas = ceiling((baseline %--% periph_vas_dis_date)/months(3)),
         t_cvas = ceiling((baseline %--% cereb_vas_dis_date)/months(3)), 
         t_dem = ceiling((baseline %--% dementia_date)/months(3)),
         t_pulm = ceiling((baseline %--% chron_pulm_dis_date)/months(3)),
         t_rheum = ceiling((baseline %--% rheumatic_dis_date)/months(3)),
         t_ulc =  ceiling((baseline %--% peptic_ulcer_date)/months(3)),
         t_miliv = ceiling((baseline %--% mild_liv_dis_date)/months(3)),
         t_modliv = ceiling((baseline %--% mod_sevr_liver_dis_date)/months(3)),
         t_diabwo = ceiling((baseline %--% diabt_wo_chron_complict_date)/months(3)),
         t_diabcomp = ceiling((baseline %--% diabt_chron_complict_date)/months(3)),
         t_pleg = ceiling((baseline %--% hemiplegia_paraplegia_date)/months(3)),
         t_renal = ceiling((baseline %--% renal_dis_date)/months(3))) %>% 
  select(bene_id, dob, baseline, end_date, hiv_base, age_base, yr_base, first_elig_state, 
         race, sex, first_event, event_date, event_time, event_int, 
         t_mi, t_chf, t_pvas, t_cvas, t_dem, t_pulm, t_rheum, t_ulc, t_miliv, 
         t_modliv, t_diabwo, t_diabcomp, t_pleg, t_renal) %>% 
  filter(event_time>=0) 
  
eligible.ids <- unique(dat3$bene_id)
  

# Manage endoscopies ------------------------------------------------------

# Grab endoscopies for eligible beneficiaries
proc2 <- proc %>% 
  filter(bene_id %in% eligible.ids) %>% 
  arrange(bene_id, proc_date) %>% 
  # Grab enrollment information
  left_join(select(dat2, bene_id, baseline, end_date, date50), by="bene_id") %>% 
  # Only keep endoscopies observed before or during follow-up
  filter(proc_date<=end_date)

# Determine who had endoscopy before age 50 (n=141,323 of full sample)
exclude <- proc2 %>% 
  filter(proc_date < date50,
         !duplicated(bene_id))
exclude.ids <- exclude$bene_id

# Exclude beneficiaries with prior endo and deal with related endoscopies
proc3 <- proc2 %>% 
  filter(!(bene_id %in% exclude.ids)) %>% 
  group_by(bene_id) %>% 
  # If endoscopies are within 1 year, they are "related"
    # Use date of first procedure
  mutate(unique = as.numeric(!duplicated(bene_id)),
         t_between = (lag(proc_date) %--% proc_date)/years(1),
         t_between = ifelse(is.na(lag(proc_date)), 0, t_between),
         cum_t_between = cumsum(t_between)) %>% 
  filter(unique==1 | cum_t_between>1) 

# Repeat process until only unrelated procedures
repeat {
  proc3 <- proc3 %>%   
    mutate(unique = ifelse(lag(unique)==1 | is.na(lag(unique)), 1, 0),
           t_between = ifelse(unique==1, 0, t_between),
           cum_t_between = cumsum(t_between)) %>% 
    filter(unique==1 | cum_t_between>1)
  temp <- filter(proc3, unique==0 & t_between<1)
  n <- nrow(temp)
  if (n==0) {break}
}  
  
# Determine time to endoscopy
proc4 <- proc3 %>% 
  mutate(endo_time = (baseline %--% proc_date)/months(3),
         int = ceiling(endo_time), 
         # If int is an integer, then endoscopy date occurs in next interval
         int = ifelse(endo_time%%1==0, int+1, int),
         t_between = (lag(proc_date) %--% proc_date)/years(1),
         flag = 1) %>% 
  select(-c(unique, cum_t_between, baseline, end_date, date50)) 


# Make data long ----------------------------------------------------------

# Create 1 record for each 3-month interval of follow-up
dat4 <- bind_rows(lapply(dat3, rep, dat3$event_int)) %>%
  group_by(bene_id) %>% 
  mutate(# Set up time
         n = 1,
         int = cumsum(n),
         drop = 0, # Has not yet dropped
         
         # Start date of 3-month interval
         int_start = baseline + months((int-1)*3),
         int_start = ifelse(is.na(int_start), 
                            ymd(substr(baseline + dmonths((int-1)*3), 1, 10)),
                            int_start),
         int_start = as.Date(int_start, origin=ymd("1970-01-01")),
         age_start = (dob %--% int_start)/years(1),
         
         # End date of 3-month interval
         int_end = lead(int_start) - days(1), # For everything except last record
         int_end = ifelse(is.na(int_end), event_date, int_end),
         int_end = as.Date(int_end, origin=ymd("1970-01-01")),
         age_end = lead(age_start),
         age_end = ifelse(is.na(age_end), (dob %--% int_end)/years(1), age_end),
         
         # Deal with rare cases where age_start and age_end are the same day
         event_int = ifelse(event_int==int & age_start>=age_end, event_int - 1, event_int),
         
         # Flag and count comorbidities
         across(c(t_mi, t_chf, t_pvas, t_cvas, t_dem, 
                  t_pulm, t_rheum, t_ulc, t_miliv, t_modliv,
                  t_diabwo, t_diabcomp, t_pleg, t_renal), 
                ~ ifelse(is.na(.x), 0, ifelse(int>=.x, 1, 0))),
         n_cond = t_mi + t_chf + t_pvas + t_cvas + t_dem + 
           t_pulm + t_rheum + t_ulc + t_miliv + t_modliv +
           t_diabwo + t_diabcomp + t_pleg + t_renal) %>%
  select(-c(n, t_mi, t_chf, t_pvas, t_cvas, t_dem, 
            t_pulm, t_rheum, t_ulc, t_miliv, t_modliv, t_diabwo, 
            t_diabcomp, t_pleg, t_renal)) %>% 
  filter(int<=event_int)

# For those censored, build artificial extra day of follow-up
censor <- dat4 %>% 
  filter(int==event_int & first_event=="censor") %>% # Censoring events
  filter(int_end<(dob + years(65) - months(1)) & int_end<ymd("2015-09-30")) %>% # But not admin censoring
  mutate(int = int + 1,
         int_start = int_end,
         age_start = age_end,
         int_end = int_end + days(1),
         age_end = (dob %--% int_end)/years(1),
         drop = 1)

dat5 <- bind_rows(dat4, censor) %>% 
  arrange(bene_id, int)

# Merge data sets
dat6 <- left_join(dat5, proc4, by=c("bene_id", "int")) %>% 
  filter(!(bene_id %in% exclude.ids)) %>%
  # Mark if they had an endoscopy 
  mutate(endo = ifelse(is.na(flag), 0, 1)) %>% 
  select(-flag)

date.fill <- dat6 %>% 
  select(bene_id, int, proc_date) %>% 
  group_by(bene_id) %>%
  mutate_at(vars(-group_cols()), list( ~ na.locf(., na.rm = FALSE))) %>% 
  rename(date_last_endo = proc_date)

dat7 <- left_join(dat6, date.fill, by=c("bene_id", "int")) %>% 
  mutate(date_last_endo = ifelse(endo==1, lag(date_last_endo), date_last_endo),
         date_last_endo = as.Date(date_last_endo, origin=ymd("1970-01-01")),
         
         # Years since last endoscopy
         years_between_endo = case_when(
                                is.na(date_last_endo) & endo==0 ~ (baseline %--% int_end)/years(1),
                                is.na(date_last_endo) & endo==1 ~ (baseline %--% proc_date)/years(1),
                                endo==0 ~ (date_last_endo %--% int_end)/years(1),
                                endo==1 ~ (date_last_endo %--% proc_date)/years(1)),
           
         # Indicator of whether beneficiary had prior endoscopy
         prior_endo = as.numeric(!is.na(date_last_endo))) %>% 
  select(-c(t_between, proc_date))

return(dat7)

}

all.dat <- mclapply(1:nsplits, function(x){set.up(split=x)}, mc.cores=10, mc.set.seed=F)
all.dat <- do.call(rbind, all.dat)


# Output data -------------------------------------------------------------

write_csv(all.dat, "./screening/data/screen_colonoscopy.csv")





