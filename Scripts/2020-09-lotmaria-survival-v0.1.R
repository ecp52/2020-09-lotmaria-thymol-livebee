############################################################
### 2020.09 Lotmaria thymol -- Analysis of mortality  ###
###########################################################

## Does parasite inoculation affect survival?
## Do diet treatments affect survival?

#v0.1-- 
#add survfit curves
#remove inverse logit back transformation

#Packages
library(readxl)
library(lubridate)
library(plyr); library(dplyr)
library(tidyverse)
library(coxme)
library(emmeans)
library(multcomp)
#library(devtools)
#install_github("kassambara/survminer") #get developmental version
library(survminer)
library(`survminer`)
library(directlabels)
library(cowplot)

###########################      INPUT DATA      ###################################
####################################################################################

#epy:
setwd("C:/Users/Evan/OneDrive - University of California, Riverside/Lotmaria/2020-09-lotmaria-thymol-livebee/Analysis")
## new shared folder
setwd("C:/Users/Evan.Palmer-Young/OneDrive/2020-09-lotmaria-thymol-livebee/Analysis")

D <- read_excel('../2020-09-lotmaria-thymol-livebee_data_v0.xlsx',
                sheet = "Inocs")
str(D$End_date)

Df <-
  D %>%
  dplyr::select(Parasite:Number, Dead_before, End_date) %>%
  mutate(End_date = ymd(End_date),
         Inoc_date = ymd(as.character("2020-09-22")))

Start.15 <- 
  Df %>%
  filter(Compound == "C" & Conc == "L" & Rep == 3 & Number == 15)

Df <-
  Df %>%
  filter(!Number == 15) %>%
  bind_rows(Start.15)%>%
  mutate(End.time = (End_date - Inoc_date)) 
str(Df)

## Relevel treatments ##
Df.c <-
  Df %>% 
  dplyr::select(Parasite:Number, Dead_before, End.time)%>%
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

Df.impute <- 
  Df.c.ord %>%
  mutate(End_time = ifelse(is.na(End.time), 15, End.time),
         Dead = ifelse(is.na(Dead_before), 0, Dead_before))


###################################
####    MODELS   ####################
##################################
contrasts <- c("contr.sum", "contr.poly")
fit <- coxme(Surv(End_time, Dead) ~ Treatment.o + (1|Rep), 
              data= Df.impute)
summary(fit)
broom::tidy(fit)
zph <- cox.zph(fit, transform = "identity")
zph
car::Anova(fit)
# Df  Chisq Pr(>Chisq)    
# Treatment.o  74.334  2.139e-12 ***

#Manually compute hazard rate as inverse logit?
#no, R works on relative ratios
#The Cox proportional hazards model does not estimate the baseline hazard
#https://stats.stackexchange.com/questions/367078/estimate-the-cumulative-baseline-hazard-function-from-a-cox-model-coefficients
(emm.contrast <-emmeans(fit, ~Treatment.o)) 
CLD(emm.contrast)

(emm.contrast <-emmeans(fit, ~Treatment.o, type = "response")) 

contrast(emm.contrast)

emm.contrast.df <- 
  summary(emmeans(fit, ~Treatment.o))
ggplot(emm.contrast.df,
       aes(x = Treatment.o, y = emmean)) +
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge2(width = 0.7))+
  ylab("Relative hazard") +
  xlab("Treatment")

## Pairwise comparisons-- letter display ##
CLD.df <- CLD(emm.contrast)

## Multcomp alternative ##
### specify all pair-wise comparisons among levels of variable "tension"
tuk <- glht(fit, linfct = mcp(Treatment.o = "Tukey"))
### extract information
tuk.cld <- multcomp::cld(tuk) ##warnings
plot(tuk.cld)

CLD.df <- multcomp::cld(emm.contrast) %>%
  mutate(group = mapvalues(.group, 
                           from = c(" 1 ", "  2"),
                           to = c("a", "b")
                           )
         )##OK

CLD.df

## Merge this to original data
## Format plot to match order of consumption plot
Df.n <- 
  Df.impute %>% ##CONTINUE TO MERGE TO ORIGINAL DATA
  ## THEN MERGE TO CLD for graphic annotation
  group_by(Parasite.o, Treatment.o, Compound.o, Conc.o) %>%
  dplyr::count() ## next merge to CLD and emmeans
  
## Hitch CLD & emmeans to other variables
### Merge to emmeans ##
cld.plus <-
  left_join(CLD.df, Df.n, by = "Treatment.o")

Compounds <- c("None", "Carvacrol", "Eugenol", "Cinnamaldehyde", "Thymol")
Xlab <- Compounds
Ylab <- "Relative hazard"

ggplot(cld.plus,
       aes(x = Compound.o, y = response, 
       color = Parasite.o, 
       shape = Conc.o)) +
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge2(width = 0.75))+
  ylab("Relative hazard") +
  xlab("Treatment") +
  scale_color_discrete(name = "Inoculation",
                       labels = c("Sham", "Parasite"),
                       direction = -1)+
  scale_shape_discrete(name = "Concentration", 
                       labels = c("0", "Low", "High"))+
  #ylim(0, 1.05*max(Cons.summ$Mean + Cons.summ$SE))+
  scale_x_discrete(labels = Compounds) +
  ylab(Ylab) +
  xlab("Compound") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1,
                                   vjust = 1, color = "black")) +
  geom_text(aes(label = group, y = asymp.UCL + 0.1),
            position = position_dodge2(width = 0.75),
            vjust = 0)

# ggsave("./Figures/Hazards_v1.2_bycompound.pdf", height = 4, width = 5)

################################
## Survival curves ##
###############################
Fit <- surv_fit(Surv(End_time, Dead) ~ Treatment.o,
                data= Df.impute)
str(Fit)
names(Fit)
length(Fit)
ggsurvplot(Fit,
           short.panel.labs =TRUE,
           short.legend.labs = TRUE,
           legend.labs = 1:10,
           legend.title = "Treatment",
           surv.median.line = "v",
           facet.by = "Treatment") 

##Facet to make labels visible
Fit <- surv_fit(Surv(End_time, Dead) ~ Treatment.o,
                data= Df.impute)
Fit

Fit2plot<- surv_summary(Fit, data= Df.impute) %>%
  mutate(Treatment.o = factor(Treatment.o,
                              levels = c("0-S", "0-P",
                                         "C-L", "C-H",
                                         "E-L", "E-H",
                                         "N-L", "N-H",
                                         "T-L", "T-H"))
  )

##New column for concentration
##Facet by concentration:
#Top row: none
#Second row-- low
#Third row-- high

P <- ggsurvplot(Fit2plot,
                data= Df.impute, 
                conf.int = TRUE,
           short.panel.labs =TRUE,
           short.legend.labs = TRUE,
           #legend.labs = 1:10,
           legend.title = "Treatment",
           surv.median.line = "v",
           palette = rep("black", 10)
      )
P
##Adjust labels
##labels for Treatment
#levels(Fit2plot$Treatment.o)
Mylabs <- c(
  "No compound, sham",
  "No compound, parasite",
  "Carvacrol-Low",
  "Carvacrol-High",
  "Eugenol-Low",
  "Eugenol-High",
  "Cinnamaldehyde-Low",
  "Cinnamaldehyde-High",
  "Thymol-Low",
  "Thymol-High"
  )

# New facet label names for dose variable
dose.labs <- Mylabs
names(dose.labs) <- levels(Fit2plot$Treatment.o)
names(dose.labs)

# Create the plot
P +
  facet_wrap(Treatment.o ~.,
             labeller = labeller(Treatment.o = dose.labs),
             ncol = 2) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black")) +
  xlab("Time (d)")
ggsave("./Figures/Survival-curves-v0-facet.pdf",
       height = 7, width = 7)


##################################

## Prep raw data for export ##
str(Df.impute)

#Optional: rename treatments #
(Originals = unique(Df.impute$Compound))
Full_name = c( "Carvacrol", "Eugenol", "Cinnamaldehyde", "Thymol", "None")
Treatment.short <- c("P", "S")
Treatment.long <- c("Parasite", "Sham")
Concentration.short <- c("L", "H", "0")
Concentration.long <- c("Low", "High", "0")

Df.exp <- 
  Df.impute %>%
  mutate(Compound.name = mapvalues(Compound,
                                   from = Originals,
                                   to = Full_name),
         Parasite.treatment = mapvalues(Parasite,
                                        from = Treatment.short,
                                        to = Treatment.long),
         Concentration = mapvalues(Conc,
                                   from = Concentration.short,
                                   to = Concentration.long))
Df.exp %>% ##Change dead-before to dead (new imputed column)
  dplyr::select(Parasite:Dead, End_time:Dead) %>%
  write_excel_csv("./Data/Processed/Survival_data.csv")

Df.exp %>%
  dplyr::select(Parasite.treatment,
                Compound.name,
                Concentration,
                Rep,
                Number,
                End_time:Dead) %>%
  write_excel_csv("./Data/Processed/Survival_data_named.csv")

##############################################
#################################################




#######################
### OLD SCRIPT #####
######################
 
# let's set up the variables
time<-TimeToDeath 
event<-Death
X<-PollenDiet
time<-as.numeric(paste(time))
X<-as.factor(X)
BeeID<-as.factor(Data$BeeID)
CageID<-as.factor(Data$CageID)

library(coxme)


str(Data)
Data$BeeID<-as.factor(Data$BeeID)
Data$CageID<-as.factor(Data$CageID)
str(Data) #Looks good

####################### Cox-mixed effects model ###############################
###############################################################################


cox1 <- coxme(Surv(TimeToDeath,Death) ~ PollenDiet +  (1|CageID) , data=Data)
summary(cox1) #no error, no treat effects
# Cox mixed-effects model fit by maximum likelihood
# Data: Data
# events, n = 202, 1700
# Iterations= 6 39 
# NULL Integrated    Fitted
# Log-likelihood -1490.11  -1463.933 -1452.019
# 
# Chisq    df          p   AIC   BIC
# Integrated loglik 52.35  3.00 2.5173e-11 46.35 36.43
# Penalized loglik 76.18 12.01 2.2076e-11 52.17 12.45
# 
# Model:  Surv(TimeToDeath, Death) ~ PollenDiet + (1 | CageID) 
# Fixed coefficients
# coef exp(coef)  se(coef)    z       p
# PollenDietNoPollen 1.194782  3.302838 0.2585345 4.62 3.8e-06
# PollenDietSun      1.340166  3.819679 0.2589473 5.18 2.3e-07
# 
# Random effects
# Group  Variable  Std Dev    Variance  
# CageID Intercept 0.29480565 0.08691037


library(car)
Anova(cox1)
# Response: Surv(TimeToDeath, Death)
# Df  Chisq Pr(>Chisq)    
# PollenDiet  2 28.497  6.487e-07 ***


######CUSTOM CONTRASTS###################
contrasts<-c("contr.sum", "contr.poly")
##OPTION 1:
#MOST obviously, change reference level to NoPollen (control group)
Data$PollenDiet2<-factor(Data$PollenDiet, levels= c("NoPollen", "Buck", "Sun"))
#and re-run the model
Coxme.releveled<-coxme(Surv(TimeToDeath,Death) ~ PollenDiet2 +  (1|CageID) , data=Data)
summary(Coxme.releveled)

#Change reference to buckwheat
Data$PollenDiet3<-factor(Data$PollenDiet, levels= c("Buck","NoPollen","Sun"))
Coxme.buckref<-coxme(Surv(TimeToDeath,Death) ~ PollenDiet3 +  (1|CageID) , data=Data)
summary(Coxme.buckref)
# Fixed coefficients
# coef exp(coef)  se(coef)    z       p
# PollenDiet3NoPollen 1.194782  3.302838 0.2585345 4.62 3.8e-06
# PollenDiet3Sun      1.340166  3.819679 0.2589473 5.18 2.3e-07



#OPTION 2: 
#USE PACKAGE 'lsmeans' to define contrasts
library(lsmeans)
library("multcomp")
summary(glht(cox1, lsm(pairwise ~ PollenDiet)))
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Buck - NoPollen == 0  -1.1948     0.2585  -4.621 1.09e-05 ***
#   Buck - Sun == 0       -1.3402     0.2589  -5.175  < 1e-05 ***
#   NoPollen - Sun == 0   -0.1454     0.1948  -0.746    0.733    

#And with lsmeans:
lsm.cox1<-lsmeans(cox1, ~PollenDiet)
lsm.cox1
# PollenDiet     lsmean        SE df  asymp.LCL  asymp.UCL
# NoPollen    0.3395110 0.1221311 NA  0.1001384  0.5788836
# Buck       -0.8552710 0.1621775 NA -1.1731330 -0.5374090
# Sun         0.4848954 0.1268896 NA  0.2361964  0.7335943

#these are on log scale
lsm.cox.df<-data.frame(summary(lsm.cox1))
exp(lsm.cox.df$lsmean)
lsm.cox.df$Hazrate<-exp(lsm.cox.df$lsmean)
lsm.cox.df$upper<-exp(lsm.cox.df$lsmean + lsm.cox.df$SE)
lsm.cox.df$lower<-exp(lsm.cox.df$lsmean - lsm.cox.df$SE)
lsm.cox.df
# PollenDiet     lsmean        SE df  asymp.LCL  asymp.UCL   Hazrate     upper     lower
# 1   NoPollen  0.3395110 0.1221311 NA  0.1001384  0.5788836 1.4042608 1.5866774 1.2428161
# 2       Buck -0.8552710 0.1621775 NA -1.1731330 -0.5374090 0.4251679 0.5000268 0.3615162
# 3        Sun  0.4848954 0.1268896 NA  0.2361964  0.7335943 1.6240051 1.8437193 1.4304739

#significance letters
Lets<-cld(lsm.cox1, Letters = letters, sort = FALSE)
Lets
Lets$.group<-gsub(" ", "", Lets$.group)
Lets$.group
Lets.slim<-dplyr::select(Lets, PollenDiet, .group)
Lsm2<-join(lsm.cox.df, Lets, by = "PollenDiet")
Lsm2


Contrasts<-contrast(lsm.cox1, "pairwise")
Contrasts
# contrast          estimate        SE df z.ratio p.value
# NoPollen - Buck  1.1947820 0.2585345 NA   4.621  <.0001
# NoPollen - Sun  -0.1453843 0.1947780 NA  -0.746  0.7358
# Buck - Sun      -1.3401664 0.2589473 NA  -5.175  <.0001

#Get odds ratios:
Contrasts.df<-as.data.frame(summary(Contrasts))
Contrasts.df$OR<-exp(abs(Contrasts.df$estimate))
Contrasts.df$ORhi<-exp(abs(Contrasts.df$estimate) )
Contrasts.df
# contrast   estimate        SE df    z.ratio      p.value       OR
# 1 NoPollen - Buck  1.1947820 0.2585345 NA  4.6213642 1.137239e-05 3.302838
# 2  NoPollen - Sun -0.1453843 0.1947780 NA -0.7464106 7.357962e-01 1.156484
# 3      Buck - Sun -1.3401664 0.2589473 NA -5.1754409 6.807034e-07 3.819679
#3.3 and 3.8-fold elevated risk of death!


#Type = response
lsm.resp<-lsmeans(cox1, ~PollenDiet, type = "response")
lsm.resp
# PollenDiet  response         SE df asymp.LCL asymp.UCL
# NoPollen   1.4042608 0.17150395 NA 1.0681192 1.7404023
# Buck       0.4251679 0.06895266 NA 0.2900232 0.5603127
# Sun        1.6240051 0.20606928 NA 1.2201167 2.0278934


lsm.coxsum<-as.data.frame(summary(lsm.cox1))
lsm.coxsum$hazard<-exp(lsm.coxsum$lsmean)

lsm.coxsum


#OPTION 3: (ADVANCED), use only if you need to compare one group to combo of other groups, 
#ie you want to compare diet x to pooled subjects of diet y and z
#can do a custom set of contrasts
levels(Data$PollenDiet2)
# "NoPollen" "Buck"     "Sun"   
#Make the contrast matrix: each number refers to one level of the factor
contr <- rbind("nopollen vs. pollen" = c(-2,1,1),
               "sunflower vs other" = c(-1,-1,2),
               "buckwheat vs other"=c(-1,2,-1),
               "sunflower vs buckwheat"=c(0,-1,1),
               "nopollen vs sunflower"=c(-1,0,1))
contr
summary(glht(Coxme.releveled, linfct = mcp(PollenDiet2 = contr)))
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# nopollen vs. pollen == 0     -1.0494     0.3775  -2.780   0.0183 *  
#   sunflower vs other == 0       1.4856     0.3783   3.926   <0.001 ***
#   buckwheat vs other == 0      -2.5349     0.4794  -5.287   <0.001 ***
#   sunflower vs buckwheat == 0   1.3402     0.2589   5.175   <0.001 ***
#   nopollen vs sunflower == 0    0.1454     0.1948   0.746   0.7479    
# ---
#NOTE higher death hazard rate in sunflower vs others and sunflower vs buckwheat


#### Plot of exponentiated hazard rates #####################3
#nb: I think hazard is in proportion expected to die per unit time
library(cowplot)
Plabs<-c("No\npollen", "Buck", "Sun")
Ylabel<-"Death hazard rate"
p.Haz<- ggplot(Lsm2, aes(x=PollenDiet, y=Hazrate))+
  geom_bar(aes(fill = PollenDiet), color = "black", 
           stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymax = upper, ymin = lower), color = "black", 
                width = 0.4, size = 1)+
  scale_fill_manual(values = c("white", "blue", "orange"))+
  theme(legend.position = "none")+
  scale_x_discrete(labels = Plabs)+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.2 * max(Lsm2$upper)))+
  ylab (Ylabel) +
  xlab(NULL) +
  theme(axis.line = element_line(size = 2),
        text = element_text (face = "bold"))
Ph2<-p.Haz +geom_text(aes(y = upper+ 0.1* max(upper), label = .group), 
                      color = "black",
                      size = 8) 
Ph2
#ggsave("Nosema.survival.hazardrate.pdf", height = 3, width = 3)
getwd()
?scale_x_discrete


##plot_grid for 2 panel figure: 
##LSA request: align x axes

#run script Apis.nosema.infection.v3.barplot.R
#"C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/nosema.apis"
Ph2
p.let
Fig.nosema<- plot_grid(p.let, Ph2, labels = LETTERS,
                       label_size = 25,
                       rel_widths = c(2,1),
                       hjust = 0,
                       vjust = 1.1,
                       align = "h")
Fig.nosema
?plot_grid
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Nosema_survival_Apis")
# ggsave("Nosema.infection.and.survival.v3.reletter.pdf",
#        height = 4, width = 8)
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Manuscript Drafts/Figures")
# ggsave("Nosema.infection.and.survival.v3.reletter.pdf",
#        height = 4, width = 8)

############################### SURVIVAL PLOT (supplement)##############################################
############################################################################################

kmsurvival1 <- survfit(Surv(Data$TimeToDeath, Data$Death) ~ Data$PollenDiet)
summary(kmsurvival1)

Xlab<-expression(bold(Time~(d)))
Ylab<-expression(bold(Proportion~surviving))

plot(kmsurvival1, xlab=Xlab, ylab=Ylab)
levels(X)
#[1] "Buck"     "NoPollen" "Sun"  
Marks<- seq(0,14,2)
Marks
Colvec<- c("black", "blue", "orange")

setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Nosema_survival_Apis")
#pdf("Nosema.survival.base.v2.pdf", height =5, width = 5)
par(mfrow=c(1,1))
par(mar=c(5, 5.5, 1, 1) + 0.2) #bottom, left, top, right
par(mgp=c(3,1,0))
plot(kmsurvival1, frame=FALSE,
     ylim=c(0.83, 1.02),
     xlab=Xlab, ylab=Ylab, cex.lab=1.5, 
     #######careful here: line types and symbols must match up in the PLOT and the LEGEND
     lty=c(2,1,1), col = Colvec,
     mark.time=c(seq(0,16,1)), mark=c(15,16,17), cex=1, las=1
)
box(bty="l", lwd=2) 
legend("bottomleft", 
       legend=c("No Pollen", "Buckwheat",  "Sunflower"),
       #levels(X), #defined custom names instead
       cex=1.5, 
       lty=c(2,1,1), col = Colvec,
       pch=c(15,16,17), title = "", bty = "n")
#dev.off()











##### Second plot: ggsurvplot #################################


##Plot with ggsurvplot to match colors, etc
#try it with ggsurv
#https://cran.r-project.org/web/packages/survminer/survminer.pdf
Leg.tit.long<- c("No Pollen", "Buckwheat",  "Sunflower")
ggsurvplot(cox1) #doesn't like coxme

cox.fixed<-coxph(Surv(Data$TimeToDeath, Data$Death) ~ Data$PollenDiet)
ggsurvplot(cox.fixed) #can't handle coxph

?ggsurvplot
Fit<-survfit(Surv(TimeToDeath, Death)~PollenDiet, data = Data)
ggsurvplot(Fit, data = Data)

p1<-ggsurvplot(fit=Fit, data = Data,
               linetype = c("dashed", "dotted", "solid"),
               break.x.by=3, break.y.by=0.04,
               xlim = c(0,20), ylim=c(0.8, 1.01),
               legend.labs = Leg.tit.long,
               legend.title = "",
               font.main=c(20, "plain", "black"),
               xlab = "Time (d)", ylab = "Proportion alive",
               font.x = c(20, "plain", "black"), font.y= c(20, "plain", "black"),
               break.time.by=1,
               censor.shape=c(16,17,15), censor.size = 2 #ONLY WORKS IN DEVEL. VERSION OFF GGSURV
)
p1


?direct.label
p.lab<- direct.label(p1$plot, method = "last.bumpup")

#export
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Nosema_survival_Apis")
# ggsave("Nosema.survival.ggsurv.pdf", height = 5, width = 5)

#2-panel plot:



Nosema.inf
Nosema.inf.big<-Nosema.inf + theme(axis.title.x = element_text(size = 20),
                                   axis.title.y = element_text(size=16))
Nosema.inf.big
p.lab 

Fig.nosema<- plot_grid(Nosema.inf.big, p.lab, labels = LETTERS,
                       label_size = 25)
Fig.nosema

# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Nosema_survival_Apis")

# ggsave("Nosema.infection.and.survival.pdf",
#        height = 5, width = 10)




################################################################
################### Survival (proportional hazards)#############################
###################     Simple survival model  -- from Thermo Bee          #################

##Set contrasts: for interactions
contrasts <- c("contr.sum", "contr.poly")

str(D_real)
D_real <- 
  D_real %>%
  mutate(End_time = End_date - Inoc_date
        )


str(D_real$End_time)

library(lubridate)

D_real$End_time_days <- 
  D_real$End_time / 60 / 60 / 24
D_real$End_time_days

D_real$End_time_days <- as.numeric(D_real$End_time_days)
D_real$End_time_days

#make inoc date factor?

str(D_real)

#Change contrasts
contrasts <- c("contr.sum", "contr.poly")

Nested.Plus<-coxme(Surv(End_time_days, Dead_before) ~ 
                     Infect*Temp_factor + (1|Colony) + (1|Inoc_date), data=D_real)
library(car)
summary(Nested.Plus)
Anova(Nested.Plus)

# Df  Chisq Pr(>Chisq)
# Infect              1 0.8216     0.3647
# Temp_factor         4 5.9409     0.2036
# Infect:Temp_factor  4 1.5860     0.8113

EMM_coxme <-
  emmeans(Nested.Plus, ~Infect|Temp_factor)

summary(EMM_coxme, type = "response")

Rel_risks <- summary(EMM_coxme, type = "response")

Xlab <- expression(Temperature~(degree*C))
Ylab_risk <- expression("Relative risk")

#Plot relative risk of death:
ggplot(Rel_risks, aes(x = Temp_factor, y = response, color = Infect, shape = Infect)) +
  geom_pointrange(aes(ymin = response - SE, ymax = response + SE), position = position_dodge(width = 0.2)) +
  xlab(Xlab) + ylab(Ylab_risk) +
  scale_color_discrete(name = "Infection\n Treatment",
                       labels = c("Parasite", "Sham")) +
  scale_shape_discrete(name = "Infection\n Treatment",
                       labels = c("Parasite", "Sham")) 
getwd()
#ggsave("Thermo_mortality_v0.pdf", height = 5, width = 5)  



Plot_props_se <-
  ggplot(Emf_grid, aes(x = Temp_factor, y = prob, color = Infect, shape = Infect)) +
  #geom_line(aes(group = Infect)) +
  geom_pointrange(aes(ymin = prob - SE, ymax = prob + SE),
                  position = position_dodge(width = 0.2)) +
  xlab(Xlab) +
  ylab("Probability of death")


contrast(EMM_coxme, method = "pairwise", adjust = "Tukey")
contrast(EMM_coxme, method = "pairwise", adjust = "bonferroni")



summary(EMM_coxme)


#Plot odds ratios:
Cox_grid <-
  cld(EMM_coxme, type = "response", sort = FALSE, Letters = letters)
Cox_grid #why is there significant contrast with "cld" but not with "contrast"
contrast(EMM_coxme, "pairwise")
contrast(EMM_coxme, "pairwise", type = "response") #now contrast significant... no correction?
