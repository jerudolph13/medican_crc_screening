
################################################################################
#
# Project: MediCan CRC Screening & Incidence
#
# Purpose: Visualize results 
#
# Last Update: 28 Aug 2023
#
################################################################################

library("tidyverse")

analysis <- "main" # {main, post2010, related6mo, colonoscopy}
weights <- "2yr" # {2yr, 4yr, 10yr, hivwt}


# Read in results ---------------------------------------------------------

if (analysis=="main") {
  res <- read_csv(paste0("../results/screen_res_", weights, ".csv"))
} else {
  res <- read_csv(paste0("../results/screen_res_", analysis, "-", weights, ".csv")) 
}

crude <- read_csv("../results/screen_res_crude.csv")

res <- res %>% 
  mutate(risk_colon1_ll = ifelse(risk_colon1_ll<0, 0, risk_colon1_ll),
         risk_colon0_ll = ifelse(risk_colon0_ll<0, 0, risk_colon0_ll))


# Weighted results --------------------------------------------------------

# Just colon
if (analysis=="main") {
  name <- paste0("../results/fig_colon_", weights, ".csv")
} else {
  name <- paste0("../results/fig_colon_", analysis, "-", weights, ".csv")
}

jpeg(name, width=5, height=5, units="in", res=300)
ggplot(data=res, aes(x=time)) +
  labs(x="\nAge (years)", y="Risk of colon cancer\n", color="", fill="") +
  theme_classic() +
  theme(legend.position=c(0.25, 0.85),
        legend.text=element_text(color="black", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title=element_text(color="black", size=14)) +
  geom_ribbon(aes(ymin=risk_colon1_ll,
                  ymax=risk_colon1_ul,
                  fill="With HIV"), alpha=0.3) +
  geom_ribbon(aes(ymin=risk_colon0_ll,
                  ymax=risk_colon0_ul,
                  fill="Without HIV"), alpha=0.3) +
  geom_step(aes(y=risk_colon0, color="Without HIV", group=1)) +
  geom_step(aes(y=risk_colon1, color="With HIV", group=1)) +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(limits=c(0, 0.04))
dev.off()

# With death
if (analysis=="main") {
  name <- paste0("../results/fig_death_", weights, ".csv")
} else {
  name <- paste0("../results/fig_death_", analysis, "-", weights, ".csv")
}

jpeg(name, width=5, height=5, units="in", res=300)
ggplot(data=res, aes(x=time)) +
  labs(x="\nAge (years)", y="Risk\n", color="", fill="", linetype="") +
  theme_classic() +
  theme(legend.position=c(0.25, 0.85),
        legend.text=element_text(color="black", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title=element_text(color="black", size=14)) +
  geom_ribbon(aes(ymin=risk_colon1_ll,
                  ymax=risk_colon1_ul,
                  fill="With HIV"), alpha=0.3) +
  geom_ribbon(aes(ymin=risk_colon0_ll,
                  ymax=risk_colon0_ul,
                  fill="Without HIV"), alpha=0.3) +
  geom_ribbon(aes(ymin=risk_death1_ll,
                  ymax=risk_death1_ul,
                  fill="With HIV"), alpha=0.3) +
  geom_ribbon(aes(ymin=risk_death0_ll,
                  ymax=risk_death0_ul,
                  fill="Without HIV"), alpha=0.3) +
  geom_step(aes(y=risk_colon0, color="Without HIV", linetype="Colon cancer", group=1)) +
  geom_step(aes(y=risk_colon1, color="With HIV", linetype="Colon cancer", group=1)) +
  geom_step(aes(y=risk_death0, color="Without HIV", linetype="Death", group=1)) +
  geom_step(aes(y=risk_death1, color="With HIV", linetype="Death", group=1)) +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  scale_linetype_manual(values=c(1, 3)) +
  scale_y_continuous(limits=c(0, 0.6))
dev.off()


# Crude results -----------------------------------------------------------

colon <- crude %>% 
  filter(state==1) %>% 
  mutate(hiv = ifelse(strata=="hiv_base=0", "Without HIV", "With HIV"))

# Just colon cancer
jpeg("../figures/fig_colon_crude.jpeg", width=5, height=5, units="in", res=300)
ggplot(data=colon, aes(x=time, y=estimate, color=as.factor(hiv), fill=as.factor(hiv))) +
  labs(x="\nAge (years)", y="Risk of colon cancer\n", color="", fill="") +
  theme_classic() +
  theme(legend.position=c(0.25, 0.85),
        legend.text=element_text(color="black", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title=element_text(color="black", size=14)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), color=NA, alpha=0.3) +
  geom_step() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(limits=c(0, 0.03))
dev.off()

# With death
death <- crude %>% 
  filter(state %in% c(1,2)) %>% 
  mutate(hiv = ifelse(strata=="hiv_base=0", "Without HIV", "With HIV"),
         outcome = ifelse(state==1, "Colon cancer", "Death"))

jpeg(paste0("../figures/fig_death_crude.jpeg"), width=5, height=5, units="in", res=300)
ggplot(data=death, aes(x=time, y=estimate, color=as.factor(hiv), fill=as.factor(hiv), 
                     linetype=as.factor(outcome))) +
  labs(x="\nAge (years)", y="Risk of colon cancer\n", color="", fill="", linetype="") +
  theme_classic() +
  theme(legend.position=c(0.25, 0.85),
        legend.text=element_text(color="black", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title=element_text(color="black", size=14)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), color=NA, alpha=0.3) +
  geom_step() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(limits=c(0, 0.4))
dev.off()
