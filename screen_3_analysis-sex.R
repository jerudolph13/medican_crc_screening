
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Run sex-stratified analysis 
#
# Last Update: 10 Jul 2023
#
################################################################################


packages <- c("dplyr", "tidyr", "magrittr", "readr", "lubridate", "broom", "splines", 
              "survival", "parallel", "zoo", "tictoc")
for (package in packages){
  library(package, character.only=T)
} 

# Program settings
args <- commandArgs(trailingOnly=T)
nboot <- as.numeric(args[1])
def <- as.numeric(args[2])
sex <- as.character(args[3])


# Read in data ------------------------------------------------------------

dat <- read_csv(paste0("./screening/data/screen_", sex, ".csv"), 
                col_types=c("cDDDfdfffffDdddfDdDdffdfDdf"))


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
        years_between_endo>(def+0.25) ~ 0, 
        # Censor if it's been less than X years between endoscopies
        endo==1 & prior_endo==1 & years_between_endo<(def-0.25) ~ 0, 
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


# Bootstrap ---------------------------------------------------------------

boot_rep <- function(r, analysis) {
  
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
  
  boot.base <- boot %>% 
    filter(!duplicated(bid))
  
  
# IP-weights --------------------------------------------------------------

  # For HIV status
  mod_hiv <- glm(as.factor(hiv_base) ~ bs(age_start, df=3) + race + first_elig_state + enroll_period, 
                 family=binomial(link="logit"), data=boot.base)$fitted.values
  den_hiv <- ifelse(boot.base$hiv_base==1, mod_hiv, 1 - mod_hiv)
  
  mod_hiv <- glm(as.factor(hiv_base) ~ 1,
                 family=binomial(link="logit"), data=boot.base)$fitted.values
  num_hiv <- ifelse(boot.base$hiv_base==1, mod_hiv, 1 - mod_hiv)
  
  boot.base$wt_hiv <- num_hiv/den_hiv
  boot.base <- select(boot.base, bid, wt_hiv)

  rm(mod_hiv)
  rm(den_hiv)
  rm(num_hiv)

  # For screening
  boot$wt_screen <- rep(NA, nrow(boot))
    # Among those with HIV
    boot.hiv <- filter(boot, hiv_base==1)
    mod_screen <- glm(as.factor(screen) ~ bs(age_start, df=3) + screen_1 + screen_2 +  
                      race + first_elig_state + n_cond + n_cond_1 + enroll_period,
                     family=binomial(link="logit"), data=boot.hiv)$fitted.values
    den_screen <- ifelse(boot.hiv$screen==1, mod_screen, 1 - mod_screen)
  
    mod_screen <- glm(as.factor(screen) ~ bs(age_start, df=3) + screen_1 + screen_2,
                     family=binomial(link="logit"), data=boot.hiv)$fitted.values
    num_screen <- ifelse(boot.hiv$screen==1, mod_screen, 1 - mod_screen)
  
    boot$wt_screen[boot$hiv_base==1] <- num_screen/den_screen
    
    # Among those without HIV
    boot.hiv <- filter(boot, hiv_base==0)
    mod_screen <- glm(as.factor(screen) ~ bs(age_start, df=3) + screen_1 + screen_2 +  
                        race + first_elig_state + n_cond + n_cond_1 + enroll_period,
                      family=binomial(link="logit"), data=boot.hiv)$fitted.values
    den_screen <- ifelse(boot.hiv$screen==1, mod_screen, 1 - mod_screen)
    
    mod_screen <- glm(as.factor(screen) ~ bs(age_start, df=3) + screen_1 + screen_2,
                      family=binomial(link="logit"), data=boot.hiv)$fitted.values
    num_screen <- ifelse(boot.hiv$screen==1, mod_screen, 1 - mod_screen)
    
    boot$wt_screen[boot$hiv_base==0] <- num_screen/den_screen
  
  rm(boot.hiv)
  rm(mod_screen)
  rm(den_screen)
  rm(num_screen)
  
  # For drop out
  boot$wt_drop <- rep(NA, nrow(boot))
  
  # Among those with HIV
  boot.hiv <- filter(boot, hiv_base==1)
  mod_drop <- glm(I(drop==0) ~ bs(age_start, df=3) + screen + screen_1 + screen_2 +  
                      race + first_elig_state + n_cond + n_cond_1 + enroll_period,
                    family=binomial(link="logit"), data=boot.hiv)$fitted.values
  den_drop <- ifelse(boot.hiv$drop==0, mod_drop, 1 - mod_drop)
  
  mod_drop <- glm(I(drop==0) ~ bs(age_start, df=3) + screen + screen_1 + screen_2,
                    family=binomial(link="logit"), data=boot.hiv)$fitted.values
  num_drop <- ifelse(boot.hiv$drop==0, mod_drop, 0)
  
  boot$wt_drop[boot$hiv_base==1] <- num_drop/den_drop
  
  # Among those without HIV
  boot.hiv <- filter(boot, hiv_base==0)
  mod_drop <- glm(I(drop==0) ~ bs(age_start, df=3) + screen + screen_1 + screen_2 +  
                    race + first_elig_state + n_cond + n_cond_1 + enroll_period,
                  family=binomial(link="logit"), data=boot.hiv)$fitted.values
  den_drop <- ifelse(boot.hiv$drop==0, mod_drop, 1 - mod_drop)
  
  mod_drop <- glm(I(drop==0) ~ bs(age_start, df=3) + screen + screen_1 + screen_2,
                  family=binomial(link="logit"), data=boot.hiv)$fitted.values
  num_drop <- ifelse(boot.hiv$drop==0, mod_drop, 0)
  
  boot$wt_drop[boot$hiv_base==0] <- num_drop/den_drop
  
  # Combine
  boot <- boot %>% 
    left_join(boot.base, by="bid") %>% 
    group_by(bid) %>% 
    mutate(cum_wt_screen = cumprod(wt_screen),
           cum_wt_drop = cumprod(wt_drop),
           wt_final = cum_screen*wt_hiv*cum_wt_screen*cum_wt_drop,
           wt_final_hiv = wt_hiv*cum_wt_drop) %>% 
    ungroup()
  
  # Get mean, min, and max of weights
    # We will need to truncate weights
  wt.summ <- boot %>% 
    filter(cum_screen==1) %>% 
    ungroup() %>% 
    group_by(age_start) %>% 
    summarize(n = n(),
              avg_wt = mean(wt_final),
              min_wt = min(wt_final),
              p5_wt = quantile(wt_final, p=0.05),
              p95_wt = quantile(wt_final, p=0.95),
              max_wt = max(wt_final))
  
  # If necessary, truncate weights at 5th, 95th percentile
  boot <- boot %>% 
    left_join(wt.summ, by="age_start") %>% 
    mutate(wt_final = ifelse(wt_final>100, p95_wt, wt_final))
  
  # Estimate cumulative incidence functions
  if (analysis=="screen weight") {

    fit <- summary(survfit(Surv(age_start, age_end, event) ~ hiv_base, data=boot, id=bid,
                           weights=wt_final, conf.type="none", se.fit=F, timefix=F))
    
    time0 = fit$time[fit$strata=="hiv_base=0"]
    res0 <- tibble(time = c(50, time0),
                   risk_colon0 = c(0, fit$pstate[1:length(time0), 3]),
                   risk_death0 = c(0, fit$pstate[1:length(time0), 2]))
    res1 <- tibble(time = c(50, fit$time[fit$strata=="hiv_base=1"]),
                   risk_colon1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 3]),
                   risk_death1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 2]))
    
    wt.summ <- wt.summ %>% 
      rename(time = age_start)
    
    res<- left_join(res0, res1, by="time") %>% 
      na.locf() %>% 
      mutate(rd_colon = risk_colon1 - risk_colon0,
             rr_colon = risk_colon1/risk_colon0,
             rd_death = risk_death1 - risk_death0,
             rr_death = risk_death1/risk_death0,
             rep = r) %>% 
      left_join(wt.summ, by="time")
    
  } else if (analysis=="hiv weight") {

    fit <- summary(survfit(Surv(age_start, age_end, event) ~ hiv_base, data=boot, id=bid,
                           weights=wt_final_hiv, conf.type="none", se.fit=F, timefix=F))
    
    time0 = fit$time[fit$strata=="hiv_base=0"]
    res0 <- tibble(time = c(50, time0),
                   risk_colon0 = c(0, fit$pstate[1:length(time0), 3]),
                   risk_death0 = c(0, fit$pstate[1:length(time0), 2]))
    res1 <- tibble(time = c(50, fit$time[fit$strata=="hiv_base=1"]),
                   risk_colon1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 3]),
                   risk_death1 = c(0, fit$pstate[(length(time0) + 1):length(fit$time), 2]))
    
    res<- left_join(res0, res1, by="time") %>% 
      na.locf() %>% 
      mutate(rd_colon = risk_colon1 - risk_colon0,
             rr_colon = risk_colon1/risk_colon0,
             rd_death = risk_death1 - risk_death0,
             rr_death = risk_death1/risk_death0,
             rep = r) 
    
  } else if (analysis=="crude") {
    
    fit <- survfit(Surv(age_start, age_end, event) ~ hiv_base, data=boot, id=bid,
                           conf.type="log-log", timefix=F)
    res <- tidy(fit) %>% 
      filter(state!="(s0)") %>% 
      select(time, estimate, state, strata, conf.high, conf.low)
    
  }

  
  toc()
  return(res)

}

start_time <- Sys.time()

# all.boot <- boot_rep(0, analysis="crude")
# all.boot <- mclapply(0:nboot, function(tt) {boot_rep(tt, analysis="hiv weight")}, mc.cores=10, mc.set.seed=F)
all.boot <- mclapply(0:nboot, function(tt) {boot_rep(tt, analysis="screen weight")}, mc.cores=10, mc.set.seed=F)
all.boot <- do.call(rbind, all.boot)

end_time <- Sys.time()
print(end_time - start_time)


# Aggregate results -------------------------------------------------------


# For point estimate, pull out results where boot=0
est <- filter(all.boot, rep==0)

# Summarize over bootstraps
boot.summ <- all.boot %>%
  group_by(time) %>%
  summarize(risk_colon0_ll = quantile(risk_colon0, p=0.025, na.rm=T),
            risk_colon0_ul = quantile(risk_colon0, p=0.975, na.rm=T),
            risk_colon1_ll = quantile(risk_colon1, p=0.025, na.rm=T),
            risk_colon1_ul = quantile(risk_colon1, p=0.975, na.rm=T),

            rd_colon_ll = quantile(rd_colon, p=0.025, na.rm=T),
            rd_colon_ul = quantile(rd_colon, p=0.975, na.rm=T),

            rr_colon_ll = quantile(rr_colon, p=0.025, na.rm=T),
            rr_colon_ul = quantile(rr_colon, p=0.975, na.rm=T),

            risk_death0_ll = quantile(risk_death0, p=0.025, na.rm=T),
            risk_death0_ul = quantile(risk_death0, p=0.975, na.rm=T),

            risk_death1_ll = quantile(risk_death1, p=0.025, na.rm=T),
            risk_death1_ul = quantile(risk_death1, p=0.975, na.rm=T),

            rd_death_ll = quantile(rd_death, p=0.025, na.rm=T),
            rd_death_ul = quantile(rd_death, p=0.975, na.rm=T),

            rr_death_ll = quantile(rr_death, p=0.025, na.rm=T),
            rr_death_ul = quantile(rr_death, p=0.975, na.rm=T))

# Merge back to point estimates
est <- left_join(est, boot.summ, by="time") %>%
  select(-rep)


# Output results ----------------------------------------------------------

if (nboot==0) {
  write_csv(est, paste0("./screening/results/screen_res_", def, "yr-ptest.csv"))
} else {
  write_csv(est, paste0("./screening/results/screen_res_", def, "yr_", sex, ".csv"))
}

# write_csv(est, "./screening/results/screen_res_hivwt.csv")
# write_csv(all.boot, "./screening/results/screen_res_crude.csv")



