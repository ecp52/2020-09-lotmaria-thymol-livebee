############################################################
### 2020.09 Lotmaria thymol -- Analysis of consumption  ###
###########################################################

## Do diet treatments affect consumption?
## Update 2022-04: change comparisons to pairwise vs 0-P group

### Main data ###
setwd("C:/Users/Evan/OneDrive - University of California, Riverside/Lotmaria/2020-09-lotmaria-thymol-livebee/Analysis")

## Using shared folder:
setwd("C:/Users/Evan.Palmer-Young/OneDrive/2020-09-lotmaria-thymol-livebee/Analysis")
setwd("C:/Users/Evan/OneDrive/2020-09-lotmaria-thymol-livebee/Analysis")

## Change log:
## 2022.04 **to do**
#Change working directory and 
 ##update comparisons to use 0-P (no compound + parasite) as ref group
## Update emmeans using regrid() to ratios
 ##https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html
  ## piglog.emm.s <- regrid(emmeans(pigroot.lm, "source"), transform = "log")
##Export table with Relative consumption\SE\T\P
## Export raw data for supplement


#Packages
library(readxl)
library(plyr); library(dplyr)
library(tidyverse)
library(directlabels)
library(cowplot)
library(lme4)
library(emmeans)
library(multcomp)

###########################      INPUT DATA      ###################################
####################################################################################

D <- read_excel('../2020-09-lotmaria-thymol-livebee_data_v0.xlsx',
                sheet = "Consumption")

str(D)

#Evaporation controls--###
 ##control for day-to-day humidity changes
Df.e <- read_excel('../2020-09-lotmaria-thymol-livebee_data_v0.xlsx',
                  sheet = "Consumption_controls")
str(Df.e)

Df.e.s <-
  Df.e %>%
  mutate(Start = lubridate::ymd("2020-09-22"),
         Date1 = lubridate::ymd(Date),
         Day = Date1 -Start) %>%
  group_by(Replicate) %>%
  mutate(Mf = lead(Start_mass),
         Dm = Start_mass - Mf) %>%
  ungroup() %>%
  group_by(Day) %>%
  dplyr::summarize(Dm.evap = mean(Dm)) %>%
  filter(!is.na(Dm.evap)) %>%
  mutate(Day = as.numeric(Day))

str(Df.e.s) ## Merge back to data with by = "Day"

ggplot(Df.e.s, aes(x = Day, y = Dm.evap)) +
  geom_point()## check for anomalies

##Merge evaporation data back to consumption data
D.net <- 
  D %>%
  left_join(Df.e.s, by = "Day")

### Practice ###
#We want to compute change in mass, dM, as difference from mass of previous day
## Except for one instance, mass of start equals mass of previous day
##or mass of end equals mass of next day

#Dplyr lag function
#https://dplyr.tidyverse.org/reference/lead-lag.html

lead(1:5)

Df.test <-
  tibble(Name = rep(c("Tom", "Mary"), 3),
         Score = letters[1:6]) %>%
  group_by(Name) %>%
  mutate(next.score = lead(Score))
Df.test
## This works fine on grouped data?

glimpse(D)

D.real <-
  D.net %>%
  group_by(Cup_ID, Day) %>%
  mutate(M0 = ifelse(is.na(M.new), M, M.new)) %>%
  arrange(Day, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Mf = dplyr::lead(M, order_by = Cup_ID, n = 1, default = NA),
         Nf = dplyr::lead(N, order_by = Cup_ID),
         Dm = M0 - Mf - Dm.evap, ##corrects for day-specific evaporation
         N.mean = (N + Nf)/2,
         Dm.norm = Dm/N.mean,
         ) 
## Caution-- no results without 'order_by', wrong results with 'order_by = Day'!

Df.c <-
  D.real %>% 
  dplyr::select(Parasite:Day, M0:Dm.norm)%>%
  filter(!(is.na(Dm.norm)|N.mean<5)) %>%
  mutate(Treatment = ifelse(Compound == "X", 
                            paste(Conc, Parasite, sep = "-"),
                            paste(Compound, Conc, sep = "-")),
         Cup = paste(Treatment, Rep, sep = "_")) 
Df.c$Treatment
unique(Df.c$Conc)

Highs <- unique(Df.c$Treatment[str_detect(Df.c$Treatment, "H")])
Lows <- unique(Df.c$Treatment[str_detect(Df.c$Treatment, "L")])

Treatment.order = c("0-S", "0-P", Lows, Highs)

unique(Df.c$Treatment)

Df.c.ord <-
  Df.c %>%
  mutate(Treatment.o = factor(Treatment, levels = Treatment.order),
         Conc.o = factor(Conc, levels = c("0", "L", "H")),
         Parasite.o = factor(Parasite, levels = c("S", "P")),
         Compound.o = factor(Compound, levels = c("X", "C", "E", "N", "T"))
  )
levels(Df.c.ord$Compound.o)[1:5]

Cons.summ <-
  Df.c.ord %>%
  group_by(Parasite.o, Compound.o, Conc.o, Treatment.o) %>%
  summarise(Mean = mean(Dm.norm),
            SD = sd(Dm.norm),
            N = length(Dm.norm),
            SE = SD/sqrt(N)) 

Ylab <- expression(Consumption~(mg~bee^-1~d^-1))

ggplot(Cons.summ, aes(x = Treatment.o, y = Mean, 
                      color = Parasite.o, shape = Conc.o)) +
  geom_pointrange(aes(ymin = Mean - 2*SE, ymax = Mean + 2*SE, 
                      color = Parasite.o, shape = Conc.o))+
  scale_color_discrete(name = "Inoculation\ntreatment",
                       labels = c("Sham", "Parasite"),
                       direction = -1)+
  scale_shape_discrete(name = "Chemical\nconcentration", 
                       labels = c("0", "Low", "High"))+
  
  #ylim(0, 1.05*max(Cons.summ$Mean + Cons.summ$SE))+
  ylab(Ylab) +
  xlab("Treatment") +
  geom_point(data = Df.c.ord %>% filter(Dm.norm < 0.05),
             aes(y = Dm.norm),
             position = position_jitter(width = 0.3), 
             size = 1, alpha = 1/3) +
  theme_bw()

getwd()
ggsave("Consumption_v1.1_raw_no_out.pdf", height = 3, width = 5)

#Parasite inoculation looks to have reduced consumption
#Among parasite-inoculated, Best consumption of eugenol low 
#Carvacrol, eugenol, cinnamaldehyde decreased consump at hi concs
#Thymol reduced consumption at low but not high conc?

#Regroup - cluster by treatment
Compounds <- c("None", "Carvacrol", "Eugenol", "Cinnamaldehyde", "Thymol")

Xlab <- Compounds

ggplot(Cons.summ, aes(x = Compound.o, y = Mean, 
                      color = Parasite.o, 
                      shape = Compound.o,
                      size = Conc.o)) +
  geom_pointrange(aes(ymin = Mean - 2*SE, ymax = Mean + 2*SE),
                  position = position_dodge(width = 0.75))+
  scale_color_discrete(name = "Inoculation",
                       labels = c("Sham", "Parasite"),
                       direction = -1)+
  scale_shape_discrete(name = "Compound", 
                       labels = Compounds)+
  scale_size_manual(name = "Concentration", 
                      labels = c("0", "Low", "High"),
                      values = (c(0.3, 0.8, 1.3)^2))+
  #ylim(0, 1.05*max(Cons.summ$Mean + Cons.summ$SE))+
  scale_x_discrete(labels = Compounds) +
  ylab(Ylab) +
  xlab("Compound") +
  geom_point(data = Df.c.ord %>% filter(Dm.norm < 0.05),
             aes(y = Dm.norm, 
                 color = Parasite.o, 
                 shape = Compound.o,
                 size = Conc.o),
             position = position_jitterdodge(jitter.width = 0.3,
                                             dodge.width = 0.75),
             alpha = 1/3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5,
                                   vjust = 0.5, color = "black"))

##Third try -- don't need shape for compound now that this is on x-axis
#Regroup - cluster by treatment
Compounds <- c("None", "Carvacrol", "Eugenol", "Cinnamaldehyde", "Thymol")
Xlab <- Compounds

ggplot(Cons.summ, aes(x = Compound.o, y = Mean, 
                      color = Parasite.o, 
                      shape = Conc.o)) +
  geom_pointrange(aes(ymin = Mean - 2*SE, ymax = Mean + 2*SE),
                  position = position_dodge(width = 0.75))+
  scale_color_discrete(name = "Inoculation",
                       labels = c("Sham", "Parasite"),
                       direction = -1)+
  scale_shape_discrete(name = "Concentration", 
                    labels = c("0", "Low", "High"))+
  #ylim(0, 1.05*max(Cons.summ$Mean + Cons.summ$SE))+
  scale_x_discrete(labels = Compounds) +
  ylab(Ylab) +
  xlab("Compound") +
  geom_point(data = Df.c.ord %>% filter(Dm.norm < 0.05),
             aes(y = Dm.norm, 
                 color = Parasite.o, 
                 shape = Conc.o),
             position = position_jitterdodge(jitter.width = 0.5,
                                             dodge.width = 0.75),
             alpha = 1/5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1,
                                   vjust = 1, color = "black"))

ggsave("Consumption_v1.2_bycompound.pdf", height = 3, width = 5)

####################################
############ MODELS #############
##########################################

## Model should assess variation by treatment w/ pairwise comparisons
str(Df.c.ord)

Df.c.noout <- 
  Df.c.ord %>% 
  filter(Dm.norm < 0.05)
M.c <- lmer(Dm.norm ~ Treatment + Day + (1|Cup), data = Df.c.noout) 
summary(M.c)
car::Anova(M.c)  
## Treatment effect is stronger after controlling for day-to-day changes in evaporation

library(emmeans)
M.c.emm <- emmeans(M.c, ~Day, at = list(Day = unique(Df.c$Day)))
M.c.emm
str(M.c.emm)
emm.plot <- as.data.frame(M.c.emm)
ggplot(emm.plot, aes(x = Day, y = emmean)) +
  geom_point() ##Better correct for evaporation if humidity decreasing
##Interesting trend of increased consumption over time 
## Perhaps reflects selective survival of bees that found feeder?

##Plot original data
ggplot(Df.c, aes(x = Day, y = Dm.norm)) +
         geom_point() +
         stat_smooth()
##Nonlinear increase of consumption over time...
##Consumption increases over first 4 d, then stable

#Try squared term
M.c.sq <- lmer(Dm.norm ~ Treatment.o + Day + I(Day^2) + (1|Cup), data = Df.c.noout) 
summary(M.c.sq)
car::Anova(M.c.sq) ## Significant

M.c.sq.emm <- emmeans(M.c.sq, ~Day, at = list(Day = unique(Df.c$Day)))
M.c.sq.emm
str(M.c.sq.emm)
emm.plot <- as_tibble(M.c.sq.emm)
ggplot(emm.plot, aes(x = Day, y = emmean)) +
  geom_point()

#Plot emmeans by treatment
#Get pairwise contrasts
?emmeans
?contrast
(M.c.sq.emm <- emmeans(M.c.sq, ~Treatment.o))
(Contrasts <- 
  contrast(M.c.sq.emm, method = "pairwise") %>%
  as_tibble() %>%
  arrange(p.value)) ##E-L & 0-S > N-H; 

write_excel_csv(Contrasts, "Consumption_contrasts_v0.1_maxpoint05.csv")

emmeans::CLD(M.c.sq.emm, details = FALSE, sort = TRUE, by = "Treatment.o",
    alpha = 0.05, Letters = letters,
    reversed = FALSE) ##use is deprecated

?cld
CLD.mc <- multcomp::cld(M.c.sq.emm) %>% ## still works
  mutate(group = as.numeric(.group),
         #letter.lab = letters[group],
         letter.lab = mapvalues(group, 
                                from = c("1", "12", "2"),
                                to = c("a", "ab", "b"))) 
CLD.mc

## Should plot the model results and CI's to account for 
 ##controlling for time and
 ##accounting for psuedo-replication/non-independence within cups
## and add text layer for sample size
glimpse(Df.c.noout)
(CLD.n <-
  Df.c.noout %>%
  group_by(Treatment.o, Parasite.o, Compound.o, Conc.o) %>%
  dplyr::summarize(N = length(!is.na(Dm.norm))) %>%
    left_join(CLD.mc, by = "Treatment.o")
)

## Replot ##
Compounds <- c("None", "Carvacrol", "Eugenol", "Cinnamaldehyde", "Thymol")
Xlab <- Compounds
pd <- position_jitterdodge(jitter.width = 0.5,
                           dodge.width = 0.75)
ggplot(CLD.n, aes(x = Compound.o, y = emmean*1000, 
                      color = Parasite.o, 
                      shape = Conc.o)) +
  geom_pointrange(aes(ymin = lower.CL*1000, ymax = upper.CL*1000),
                  position = pd)+
  scale_color_discrete(name = "Inoculation",
                       labels = c("Sham", "Parasite"),
                       direction = -1)+
  scale_shape_discrete(name = "Concentration", 
                       labels = c("0", "Low", "High"))+
  #ylim(0, 1.05*max(Cons.summ$Mean + Cons.summ$SE))+
  scale_x_discrete(labels = Compounds) +
  ylab(Ylab) +
  xlab("Compound") +
  geom_point(data = Df.c.ord %>% filter(Dm.norm < 0.05),
             aes(y = Dm.norm*1000, 
                 color = Parasite.o, 
                 shape = Conc.o),
             position = pd,
             alpha = 1/5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1,
                                   vjust = 1, color = "black")) +
  geom_text(aes(y = upper.CL*1000, label = letter.lab), 
            vjust = -0.5, #color = "black",
            position = pd)

ggsave("Consumption_v1.3_emmean.pdf", height = 3, width = 5)

## Review model results ##
car::Anova(M.c.sq)
summary(M.c.sq)

summary(M.c.sq.emm) %>% arrange(desc(emmean))

Contrasts %>% arrange(p.value)

## Overall summary ##
Df.c.noout %>%
  ungroup() %>%
  dplyr::summarise(N = length(!is.na(Dm.norm)),
                   Mean = mean(Dm.norm, na.rm = TRUE),
                   SD = sd(Dm.norm, na.rm = TRUE),
                   SE = SD/sqrt(N))

#############################################
####           Write data           ##########
###############################################

str(Df.c.ord)
Df.exp <-
  Df.c.ord %>%
  dplyr::select(Parasite:Dm.norm) %>%
  write_excel_csv(path = "Consumption-data-maxpoint05.csv")

Df.exp <-
  Df.c %>%
  dplyr::select(Parasite:Dm.norm) %>%
  write_excel_csv(path = "Consumption-data-all.csv")

## More informative data output
#Optional: rename treatments #
(Originals = unique(Df.c.ord$Compound))
Full_name = c( "Carvacrol", "Eugenol", "Cinnamaldehyde", "Thymol", "None")
Treatment.short <- c("P", "S")
Treatment.long <- c("Parasite", "Sham")
Concentration.short <- c("L", "H", "0")
Concentration.long <- c("Low", "High", "0")

Df.exp.named <- 
  Df.c %>%
  mutate(Compound.name = mapvalues(Compound,
                                   from = Originals,
                                   to = Full_name),
         Parasite.treatment = mapvalues(Parasite,
                                        from = Treatment.short,
                                        to = Treatment.long),
         Concentration = mapvalues(Conc,
                                   from = Concentration.short,
                                   to = Concentration.long))
Df.exp.named %>%
  dplyr::select(Parasite.treatment,
                Compound.name,
                Concentration,
                Rep:Mf,
                N.mean:Dm.norm) %>%
  write_excel_csv("Consumption_data_named.csv")

