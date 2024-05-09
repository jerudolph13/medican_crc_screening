
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Run analysis 
#
# Last Update: 01 Sep 2022
#
################################################################################

packages <- c("dplyr", "tidyr", "magrittr", "readr", "broom", "splines", 
              "survival", "parallel", "zoo", "tictoc")
for (package in packages){
  library(package, character.only=T)
} 

# Number of bootstrap resamples
nboot <- 0 #500


# Read in data ------------------------------------------------------------

dat <- read_csv("./screening/data/screen_exp-sample.csv")


# Manage data -------------------------------------------------------------

dat <- dat %>% 
  filter(!is.na(age) & !is.na(race) & !is.na(sex)) %>% 
  group_by(bene_id) %>% 
  mutate(# Set up confounders
         race = factor(race),
         sex = factor(sex),
         state = factor(first_elig_state),
         n_cond = factor(ifelse(n_cond==0, 0,
                         ifelse(n_cond==1, 1, 2))),
         n_cond_1 = lag(n_cond),
         n_cond_1 = factor(ifelse(is.na(n_cond_1), n_cond, n_cond_1)),
         n_cond_2 = lag(n_cond, n=2),
         n_cond_2 = factor(ifelse(is.na(n_cond_2), n_cond_1, n_cond_2)),
         age_f = factor(age),
         enroll_period = factor(case_when(
           yr_base %in% 2001:2005 ~ 1,
           yr_base %in% 2006:2010 ~ 2,
           T ~ 3)),
         
         # Assign screening protocol:
          # No endoscopies under 50
          # Required to have 1 endoscopy before 60
          # Can get 2nd endoscopy after 60 but no more
         cum_endo = cumsum(endo),
         screen = case_when(
           age<50 & cum_endo>0 ~ 0,
           age>=50 & age<60 & cum_endo>1 ~ 0,
           age==60 & cum_endo==0 ~ 0,
           age>60 & cum_endo>2 ~ 0,
           T ~ 1),
         screen_1 = lag(screen, default=1),
         screen_2 = lag(screen, n=2, default=1),
         cum_screen = cumprod(screen),
         
         # Set up outcome
         event = factor(event)) %>% 
  ungroup()

base <- filter(dat, !duplicated(bene_id))


# Bootstrap ---------------------------------------------------------------

boot_rep <- function(r) {
  
  tic(paste0("bootstrap ", r))
  set.seed(r+1)
  samp <- table(base[sample(1:nrow(base), nrow(base), replace=T), (names(base)=="bene_id")])
 
  # The step below pulls in the simulated data for boot=0; 
  # otherwise grabs all records for the resampled observations
  boot <- NULL
  if(r==0){
    boot <- dat %>% 
      rename(bid = bene_id)
  } else{
    for(zzz in 1:max(samp)){ 
      cc <- dat[dat$bene_id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid <- paste0(cc$bene_id, zzz)
      boot <- rbind(boot, cc)
    }
    boot <- select(boot, -bene_id)
  }
  
  boot_base <- boot %>% 
    mutate(first = !duplicated(bid)) %>% 
    filter(first==1)
  
  
# IP-weights --------------------------------------------------------------

  # For HIV status
  mod_hiv <- glm(as.factor(hiv_base) ~ bs(age, df=3) + race + sex + state + enroll_period, 
                 family=binomial(link="logit"), data=boot_base)$fitted.values
  den_hiv <- boot_base$hiv_base*mod_hiv + (1 - boot_base$hiv_base)*(1 - mod_hiv)
  
  mod_hiv <- glm(as.factor(hiv_base) ~ 1,
                 family=binomial(link="logit"), data=boot_base)$fitted.values
  num_hiv <- boot_base$hiv_base*mod_hiv + (1 - boot_base$hiv_base)*(1 - mod_hiv)
  
  boot_base$wt_hiv <- num_hiv/den_hiv
  boot_base <- select(boot_base, bid, wt_hiv)

  rm(mod_hiv)
  rm(den_hiv)
  rm(num_hiv)

  # For screening
  mod_screen <- glm(as.factor(screen) ~ age_f + screen_1 + screen_2 + race + sex + 
                      state + n_cond + n_cond_1 + enroll_period,
                     family=binomial(link="logit"), data=boot)$fitted.values
  den_screen <- boot$screen*mod_screen + (1 - boot$screen)*(1 - mod_screen)
  
  mod_screen <- glm(as.factor(screen) ~ age_f + screen_1 + screen_2,
                     family=binomial(link="logit"), data=boot)$fitted.values
  num_screen <- boot$screen*mod_screen + (1 - boot$screen)*(1 - mod_screen)
  
  boot$wt_screen <- num_screen/den_screen

  rm(mod_screen)
  rm(den_screen)
  rm(num_screen)
  
  # Combine
  boot <- boot %>% 
    left_join(boot_base, by="bid") %>% 
    group_by(bid) %>% 
    mutate(cum_wt_screen = cumprod(wt_screen),
           wt_final = cum_screen*wt_hiv*cum_wt_screen) %>% 
    select(bid, age, hiv_base, event, cum_screen, wt_final)
  
  # Get mean, min, and max of weights
  wt.summ <- boot %>% 
    filter(cum_screen==1) %>% 
    group_by(age) %>% 
    summarize(avg_wt = mean(wt_final),
              min_wt = min(wt_final),
              max_wt = max(wt_final)) %>% 
    rename(time = age)
  
  # Crude cumulative incidence
  fit <- summary(survfit(Surv(age-1, age, event) ~ hiv_base, data=boot, id=bid,
                         conf.type="none", se.fit=F))

  # Weighted cumulative incidence
  # fit <- summary(survfit(Surv(age, age+1, event) ~ hiv_base, data=boot, id=bid,
  #                     weights=wt_final, conf.type="none", se.fit=F))
  
  time0 = fit$time[fit$strata=="hiv_base=0"]
  res0 <- tibble(time = c(18, time0),
                 risk_colon0 = c(0, fit$pstate[1:length(time0), 2]),
                 risk_other0 = c(0, fit$pstate[1:length(time0), 3]),
                 risk_death0 = c(0, fit$pstate[1:length(time0), 4]))
  res1 <- tibble(time = c(18, fit$time[fit$strata=="hiv_base=1"]),
                 risk_colon1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 2]),
                 risk_other1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 3]),
                 risk_death1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 4]))
  
  res<- full_join(res0, res1, by="time") %>% 
    na.locf() %>% 
    mutate(rd_colon = risk_colon1 - risk_colon0,
           rd_other = risk_other1 - risk_other0,
           rd_death = risk_death1 - risk_death0,
           rep = r) %>% 
    left_join(wt.summ, by="time")
  
  toc()
  return(res)

}

start_time <- Sys.time()

#all.boot <- lapply(0:0, function(tt) {boot_rep(tt)})
all.boot <- mclapply(0:nboot, function(tt) {boot_rep(tt)}, mc.cores=10, mc.set.seed=F) 
all.boot <- do.call(rbind, all.boot)

end_time <- Sys.time()
print(end_time - start_time)


# Aggregate results -------------------------------------------------------

# For point estimate, pull out results where boot=0
est <- filter(all.boot, rep==0)

# Summarize over bootstraps
boot.summ <- all.boot %>% 
  group_by(time) %>% 
  summarize(se_colon0 = sd(risk_colon0, na.rm=T),
            se_colon1 = sd(risk_colon1, na.rm=T),
            se_colon_rd = sd(rd_colon, na.rm=T),
            se_other0 = sd(risk_other0, na.rm=T),
            se_other1 = sd(risk_other1, na.rm=T),
            se_other_rd = sd(rd_other, na.rm=T),
            se_death0 = sd(risk_death0, na.rm=T),
            se_death1 = sd(risk_death1, na.rm=T),
            se_death_rd = sd(rd_death, na.rm=T))

# Merge back to point estimates
est <- left_join(est, boot.summ, by="time") %>% 
  select(-rep)
write_csv(est, "./screening/results/screen_res_exp-crude.csv")




