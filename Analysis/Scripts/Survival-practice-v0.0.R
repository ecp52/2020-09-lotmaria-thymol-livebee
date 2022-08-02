## 2020-10-16 ##
## Review survival analysis ##
## Two tutorials ##
#One for R with Kassambara
#And a Data Camp tutorial
#https://www.datacamp.com/community/tutorials/survival-analysis-R

#Background:
#We are evaluating p(survival) to time (t)
#and testing for differences in p(survival) across groups
#In cox proportional hazards model

#Data analysis: Orders the data by longest-surviving
#log-rank test compares survivorship of two different groups

#assumptions of proportional hazards:
#relative p(death) in the two groups are constant over time
# (i.e., curves don't cross)

library(survival)
data(ovarian)
library(dplyr)
library(survminer)

#install.packages("ggeffects")
library(ggeffects)

glimpse(ovarian)
#futime and fustat (outcome) are recorded
#rx is the treatment
#ecog.ps is performance
#resid.ds is disease regression

#relabel variables
# Dichotomize age and change data labels
ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "yes"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("good", "bad"))


#dichotomize age
hist(ovarian$age)
ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "old", "young"))
ovarian$age_group <- factor(ovarian$age_group)

surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object 
#+ means alive at last observation
plot(surv_object)

fit1 <- survfit(surv_object ~ rx, data = ovarian)
summary(fit1)
plot(fit1)
str(fit1)

library(ggplot2)
survminer::ggsurvplot(fit1, data = ovarian, pval = TRUE)
#this is p-value for more conservative log-rank test

fit2 <- survfit(surv_object ~ resid.ds, data = ovarian)
ggsurvplot(fit2, data = ovarian, pval = TRUE)
#vertical lines for censored observations

#full model and plots of hazard ratios
fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)
summary(fit.coxph)
ggforest(fit.coxph, data = ovarian)
#that's informative:
#rx B has lower hazard
#persistent disease --> higher hazard
#young age <50 years --> lower hazard
#proportional hazards test is more powerful and controls for covariates
#log-rank test is more conservative

#Check the assumptions:
ftest <- cox.zph(fit.coxph)
ftest
#check for homogeneity
#if proportional hazards assumption holds, then 
#beta(t) function should be horizontal line
survminer::ggcoxzph(ftest)
#p-value tests null hypothesis of slope = 0

## Work with output
car::Anova(fit.coxph)

library(emmeans)
emm.ph <- emmeans(fit.coxph, ~age_group + rx, 
                  type = "response")
emm.df <- as_tibble(summary(emm.ph))
ggplot(emm.df, aes(x = age_group, y = hazard, shape = rx)) +
  geom_pointrange(aes(ymin = hazard - asymp.LCL, ymax = hazard + asymp.LCL ),
                  position = position_dodge2(width = 0.7)) +
  ylab("Relative risk")


## Try changing default contrasts 
#Change contrasts
contrasts <- c("contr.sum", "contr.poly")

fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)
summary(fit.coxph)

fit.int <- coxph(surv_object ~ rx * age_group, 
                 data = ovarian) ##warning

ggforest(fit.coxph, data = ovarian)

emm.ph <- emmeans(fit.coxph, ~age_group + rx, 
                  type = "response")
emm.df <- as_tibble(summary(emm.ph))
ggplot(emm.df, aes(x = age_group, y = hazard, shape = rx)) +
  geom_pointrange(aes(ymin = hazard - asymp.LCL, ymax = hazard + asymp.LCL ),
                  position = position_dodge2(width = 0.7)) +
  ylab("Relative risk")

emm.contrast <- emmeans(fit.int, ~rx|age_group, 
                        type = "response")
contrast(emm.contrast)

#This can also be directly ported to data frame
#install.packages("ggeffects")
# library(ggeffects)
?ggemmeans
ggemm <- ggemmeans(fit.int, ~rx|age_group, 
          type = "fixed")
plot(ggemm)

#Manually compute hazard rate as inverse logit
emm.contrast <- 
  summary(emmeans(fit.int, ~rx|age_group)) %>%
  mutate(Rate = boot::inv.logit(emmean),
         Lower = boot::inv.logit(emmean - SE),
         Upper = boot::inv.logit(emmean + SE))
ggplot(emm.contrast,
       aes(x = age_group, y = Rate, color = rx)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge2(width = 0.7))

###################################################################
#More plotting tips for ggsurvplot and survminer
#https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf
###############################################################
?ggsurvplot
ggsurvplot(fit2)

?surv_summary
summary_alt <- surv_summary(fit2) #data frame rather than list
str(summary_alt)

#Survminer plot
str(summary_alt)
ggsurvplot(summary_alt, data = ovarian, 
           linetype = "strata",
           conf.int = TRUE,
           pval = TRUE, fun = "pct",
           legend.title = "Disease status",
           legend.labs = c("Remission", "Active"),
           surv.median.line = "h")
?ggsurvplot
#Manual plot
ggplot(summary_alt, aes(x = time, y = surv, 
                        color = resid.ds,
                        linetype = resid.ds))+
  geom_line() +
  geom_ribbon(aes(ymin = surv - std.err, ymax = surv + std.err),
              alpha = 0.2) +
  ylab("Survival probability") +
  theme_bw()


#################################################
#### Kassambara tutorial with ggsurvplot #####
######################################################
data("lung")
glimpse(lung)

#Survival curves:
(fit <- survfit(Surv(time, status) ~ sex, data = lung))

summary(fit)
summary(fit)$table #max, min, median survival times

#Coax into a data frame
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
                )
head(d)

#or
(d <- surv_summary(fit))

d.summ <- attr(d, "table") 
str(d.summ)
Sex = tibble(Sex = row.names(d.summ))
d.summ.df <-
  bind_cols(Sex, d.summ)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           legend.title = "Sex",
           legend.labs = c("Male", "Female"))

## More than two curves
require("survival")
fit2 <- survfit( Surv(time, status) ~ sex + rx + adhere,
                 data = colon )
# Plot survival curves by sex and facet by rx and adhere
ggsurv <- ggsurvplot(fit2, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw(),
                     legend.title = "Adherence",
                     legend.labs = c("O", "L", "F/u"),
                     facet.by = c("sex","adhere"),
                     short.panel.labs = TRUE,
                     panel.labs = list(sex = c("Male", "Female")))
ggsurv
?ggsurvplot

##Tweak theme using $plot as base plot
ggsurv$plot +
  theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(rx ~ adhere)

#Extended practice 
#https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
