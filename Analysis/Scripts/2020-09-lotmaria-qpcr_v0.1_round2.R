### Script to analyze infection intensity ###
### E Palmer-Young & L Markowitz, 2021-11 ###

## Steps ##
## Read in qPCR data 
## Merge to sample metadata
## Normalize infection intensity
## Compare groups
## Bonus: check for infection of mites

## Update 2022-03-12 Add second set of samples ##
 ## Kyle's qPCR for trypanosomatids, 8 more plates
####################

#update.packages(checkBuilt=TRUE, ask=FALSE)
#install.packages('TMB', type = 'source')
library(tidyverse)
library(qdap) ## for bracketX
library(readxl)
library(openxlsx)
library(cowplot)
library(directlabels)


library(glmmTMB)
library(lme4)
library(DHARMa)
library(MuMIn) #for r-squared
library(rsq)
#install.packages('see')
library(performance)
library(emmeans)

### Read in qPCR data ###
#library(readxl)

getwd()
setwd("C:/Users/Evan.Palmer-Young/OneDrive/2020-09-lotmaria-thymol-livebee/Analysis")
setwd("C:/Users/Evan/OneDrive/2020-09-lotmaria-thymol-livebee/Analysis")

files <- 
  list.files(path = "./qpcr_exports/Lotmaria", 
             pattern = "*.xlsx", full.names = T)
files.actin = 
  list.files(path = "./qpcr_exports/Actin", 
             pattern = "*.xlsx", full.names = T)
files

files[1]

prueba = read_excel(files[1], sheet = 1)

## looks like file names are too long ##
## Try a new folder 
files <- 
  list.files(
    path = "./qpcr_exports/Lotmaria/nicknames", 
                    pattern = "*.xlsx", full.names = T)

files.actin <- 
  list.files(
    path = "./qpcr_exports/Actin/nicknames", 
    pattern = "*.xlsx", full.names = T)
files.actin

read_excel(files[1], sheet = 1) 
read_excel(files.actin[1], sheet = 1)
read_excel(files[11], sheet = 1)
## need to open and save the file on this computer...

all.files = c(files, files.actin)
length(all.files) ## total number of plates

tbl <- 
  sapply(all.files, read_excel, simplify=FALSE) %>% 
  bind_rows(.id = "id") 


## MAKE A NOTE OF SAMPLES TO EXCLUDE ##
# !!! Make note of exclusions/exceptions!!
#   !Exclude TUBE 484 (T-L-6-8) 
# Note: tubes 340 & 581 are both labeled "P"


Samples = 
  tbl %>%
  rename(Run = id) %>%
  rename(N.copies = `Starting Quantity (SQ)`) %>%
  #filter(Content == "Unkn") %>% ## change this for round 2-- *contains* "Unk"
  filter(grepl('Unk', Content)) %>%
  filter(!is.na(Sample)) %>%
  dplyr::select(Target, Run, Well, Sample, Cq, N.copies) %>%
  ## impute zero copies or Cq of 40 for NA:
  mutate(Cq.impute = ifelse(is.na(Cq), 40, Cq),
         N.copies.impute = ifelse(is.na(N.copies), 0, N.copies))

Samples ## 1440 samples after round 1

length(unique(Samples$Sample))

ggplot(Samples, aes(x = log10(as.numeric(N.copies)))) +
  geom_histogram() +
  facet_grid(Target ~. )

unique(Samples$Run[Samples$Target == "Lotmaria"])

## Recode "Trypanosome" --> "Lotmaria
Samples = 
  Samples %>%
  mutate(Target = plyr::mapvalues(Target, 
                                  from = c("Trypanosome"),
                                  to = c("Lotmaria")))

ggplot(Samples, aes(x = log10(as.numeric(N.copies)))) +
  geom_histogram() +
  facet_grid(Target ~. )

hist(log10(Samples$N.copies))

Samples.summ =
  Samples %>%
  group_by(Target, Sample) %>%
  summarize(Cq.median = median(Cq.impute),
            Copies.median = median(N.copies.impute)*200,
            Cq.min = min(Cq.impute),
            Cq.max = max(Cq.impute),
            Cq.span = Cq.max - Cq.min)

#Two-thirds of samples span <1 cq (as of 2021-11-12)
ggplot(Samples.summ) +
  geom_histogram(aes(x = Cq.span))
#For those with a big span, usually one point is the outlier

row.names(Samples.summ)

ggplot(Samples.summ, aes(x = row.names(Samples.summ),
                         y = Cq.median)) +
  geom_point(color = "orange") +
  geom_point(aes(y = Cq.min), color = "red") +
  geom_point(aes(y = Cq.max), color = "blue") +
  facet_grid(Target ~ Cq.median > 30) +
  theme_bw()

ggsave("./Figures/amplification_scatter_target_cq30.pdf",
       height = 10, width = 30)

#most of the samples with a broad range 
#have either a low copy number (Cq.median >35)
#or one sample that didn't amplify

## Spread the sample names for merge ##
#?pivot_wider
Samples.summ.wide = 
  Samples.summ %>%
  pivot_wider(names_from = Target,
              values_from = Cq.median:Cq.span) 

names(Samples.summ.wide)

## Rename to keep old sample names
strip.lotmaria = 
  function(x){gsub(pattern = "_Lotmaria", 
                   replacement = "",
                   x)}

Samples.summ.rename = 
 Samples.summ.wide %>%
 rename_with(strip.lotmaria)

names(Samples.summ.rename)

## Where are the sample data ##
Samples.meta = 
  read_excel("../2020-09-lotmaria-thymol-livebee_data_v1_gdrive_dissections.xlsx",
             sheet = "Dissections")
glimpse(Samples.meta)

Samples.trim = 
  Samples.meta %>%
  dplyr::select(Tube, Label, N_MITES) %>%
  rename(Sample = Tube)%>%
  mutate(Sample = as.character(Sample))

str(Samples.trim)
str(Samples.summ)

Samples.join = 
  Samples.trim %>%
  right_join(Samples.summ.rename, by = "Sample")

## Get the proper metadata ##
Inocs.meta = 
  read_excel(
    "../2020-09-lotmaria-thymol-livebee_data_v1_gdrive_dissections.xlsx", 
    sheet = "Filtered_Data") %>% 
  ##  sheet "Filtered_Data"has data from bees frozen alive at 7d
  dplyr::select(Bee_ID, Label:Number) %>%
  #rename(Sample = Label) %>%
  filter(Label %in% unique(Samples.join$Label))

## check out the sample data ##
Samples.join.parse = 
  Samples.join %>%
  tidyr::separate(Label, 
                  into = c("Compound", "Conc", "Cage", "Beenumber"),
                  remove = FALSE) %>%
  mutate(Treatment = paste(Compound, Conc, sep = "_"))

## Join it properly ##
Samples.full =
  Samples.join %>%
  left_join(Inocs.meta, by = "Label") %>%
  # mutate(Treatment = paste(Parasite, Compound,Conc,
  #                          sep = "_"),
  #        .after = Label) %>%
  mutate(Treatment = ifelse(Compound == "X", 
                            paste(Conc, Parasite, sep = "-"),
                            paste(Compound, Conc, sep = "-")))

str(Samples.full$Treatment)

Highs <- unique(Samples.full$Treatment[str_detect(Samples.full$Treatment, "H")])
Lows <- unique(Samples.full$Treatment[str_detect(Samples.full$Treatment, "L")])

Treatment.order = c("0-S", "0-P", Lows, Highs)

Samples.full =
  Samples.full %>%
  mutate(Treatment.o = factor(Treatment, levels = Treatment.order),
         Conc.o = factor(Conc, levels = c("0", "L", "H")),
         Parasite.o = factor(Parasite, levels = c("S", "P")),
         Compound.o = factor(Compound, levels = c("X", "C", "E", "N", "T"))
  ) %>%
  mutate(Cage = paste(Compound, Conc, Rep, sep = "_"), .after = "Label")

Samples.summ =
  Samples.full %>%
  #group_by(Compound, Conc) %>%
  group_by(Treatment.o) %>%
  count()
Samples.summ

Distrib = 
  Samples.full %>%
  summarize(N.tot = length(Cq.median),
            Median = median(Cq.median, na.rm = TRUE),
            Zeroes = length(Cq.median [Cq.median == 40]),
            Prop.zeroes = Zeroes/N.tot)

Distrib  

ggplot(Samples.full %>% filter(!is.na(Treatment)), 
       aes(x = Treatment, y = Cq.median))+
  geom_boxplot() +
  geom_point(position = "jitter") +
  ylab("Median Cq") +
  xlab("Treatment (Parasite_Compound_Concentration)") +
  theme_bw()

ggsave("./Figures/exploratory_boxplot_lotmaria_v0.1_round2.pdf",
       height = 3, width = 5)

## Good news-- sham-inoculated have no parasites
## Good news-- high-dose carvacrol and eugenol ~8x reductions

### Can we gain anything by using actin to correct for extraction efficiency?
df.corr = 
  Samples.full %>%
  filter(!is.na(Treatment)) %>%
  filter(!is.na(Cq.median_Actin))

df.corr

ggplot(df.corr, aes(x = Cq.median_Actin, y = Cq.median)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Cq.median Lotmaria")
  
## doesn't look helpful-- negative rather than positive correlation

ggsave(
  "./Figures/Lotmaria.vs.actin.v1.round2.pdf", 
  height = 5, width = 5
)

## what proportion of variance explained by actin?
lm.correl = 
  lm(Cq.median ~ Cq.median_Actin, data = df.corr)
summary(lm.correl)
#Round 1:
# Multiple R-squared:  0.02317,	Adjusted R-squared:  0.01819 
# F-statistic: 4.649 on 1 and 196 DF,  p-value: 0.03228
#Both rounds: not significant
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      27.3865     4.4964   6.091 2.81e-09 ***
#   Cq.median_Actin   0.1625     0.1842   0.882    0.378    

model_performance(lm.correl) #pft, r-squared 0.002


##
## Let's check whether actin normalization cuts down on variance ##
glimpse(Samples.full)
Mean.actin = 
  Samples.full %>%
  filter(!is.na(Treatment)) %>%
  summarize(Mean.act = mean(Copies.median_Actin, na.rm = TRUE)) %>%
  pull(Mean.act)
Mean.actin

## Apply normalization ##
Samples.corrxn = 
  Samples.full %>%
  filter(!is.na(Treatment)) %>%
  mutate(grand.mean.actin = Mean.actin[1]) %>%
  mutate(correction.factor = Copies.median_Actin/grand.mean.actin,
         Copies.corrected = Copies.median/correction.factor)

## Melt ##
Corrected.long = 
  Samples.corrxn %>%
  filter(!is.na(Copies.corrected)) %>%
  dplyr::select(Sample, Label, Bee_ID:Compound.o,
                Copies.median, Copies.corrected) %>%
  pivot_longer(cols = c(Copies.median, Copies.corrected),
                        names_to = "Variable")

## Check for changes ##
ggplot(Corrected.long, aes(x = Treatment.o, 
                           y = value, 
                           group = forcats::fct_rev(Variable))) +
  geom_boxplot(aes(fill = forcats::fct_rev(Variable),
                   group = interaction(Treatment.o, forcats::fct_rev(Variable))), 
               position = position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(aes(shape = forcats::fct_rev(Variable)),
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_y_continuous(trans = "log1p", breaks = 10^seq(1, 7, by = 1)) +
  xlab("Treatment") +
  ylab(expression(italic(Lotmaria)~Copy~number~(bee^-1))) + 
  scale_fill_discrete(name = element_blank(), 
                      labels = rev(c("Normalized", "Raw"))) + 
  scale_shape_discrete(name = element_blank(), 
                      labels = rev(c("Normalized", "Raw"))) +
  theme_bw()

ggsave("./Figures/raw.vs.normalized.v0.pdf",
       height = 3.5, width = 6) ## 

## Check variance ##
SD.comparison = 
  Corrected.long %>%
  group_by(Treatment.o, Variable) %>%
  summarize(SD = sd(value, na.rm = TRUE)) 

ggplot(SD.comparison, aes(x = Treatment.o,
                          y = SD, 
                          fill = 
                            forcats::fct_rev(Variable))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "Copy number",
                      labels = rev(c("Corrected", 
                                     "Raw"))) +
  xlab("Treatment") +
  ylab("Standard deviation") +
  theme_bw()

ggsave("./Figures/Raw.vs.corrected.SD.v1.2round.pdf", height = 3, width = 5)
## Looks like correction is increasing variance, not decreasing it
## And only explains 2% of the variance in round 1, less in rounds 1-2!!
## Thinking to go with raw copy numbers rather than normalized

####### Look at the oddball samples ######
## Have a look at the mites ##
unique(Samples.full$Label)

Mites = 
  Samples.full %>%
  filter(grepl('MITE', Label))

ggplot(Mites, aes(x = Cq.median)) +
  geom_histogram()
Mites$Cq.median
## Shucks, no parasites deteted

## How about the controls ##
T0 = 
  Samples.full %>%
  filter(grepl("TIME", Label))

ggplot(T0, aes(x = Cq.median)) +
  geom_histogram() +
  facet_grid(.~Label)

t0.ct = 
  ggplot(T0, aes(x = Label, y = Cq.median)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = "jitter") +
  theme_bw()
 # facet_grid(.~Label)
t0.ct

t0.copies = 
  ggplot(T0, aes(x = Label, y = Copies.median)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = "jitter") +
  theme_bw()
# facet_grid(.~Label)

t0.copies
## a little frustrating to have 
## a 7 ct difference at baseline
2^7 #128-fold

cowplot::plot_grid(t0.ct, t0.copies, nrow = 1)

ggsave("./Figures/Time_0_lotmaria.pdf",
       height = 3, width = 6)

t0.summ = 
  T0 %>%
  group_by(Label) %>%
  summarize(Median = median(Copies.median),
            Q25 = quantile(Copies.median, 0.25),
            Q75 = quantile(Copies.median, 0.75),
            Mean = mean(Copies.median),
            SD = sd(Copies.median)
  )

t0.summ ## Initial count only 41K

T0.median = 
  T0 %>%
  filter(Label == "TIME 0 PAR") 
Median.t0 = median(T0.median$Copies.median)

Median.t0

## Check the inoculum ##
## Expected: 40K cells/uL * 150 uL in 200 uL extract
Expected = 40*10^3 * 150
Expected
Observed = 
  Samples.full %>%
  filter(Sample == "INOC150") %>%
  summarize(Median = median(Copies.median)) %>%
  pull(Median)
Observed

Ratio = Observed/Expected
Ratio

## ~87 copies per cell? 
## This is in range of 80-200 copies reported for Leishmania

## However, median time 0 is only ~40K copies, which equates to
Copies.per.cell = Ratio
Cells.t0 = Median.t0/Copies.per.cell
Cells.t0
#477.591, quite a bit below the goal of 100K cells inoculated

Expected.t0 = 10^5*Copies.per.cell
Fold.loss.vs.predicted = Expected.t0/Median.t0 
Fold.loss.vs.predicted
## 209.3842

## There must be a sizeable loss of cells during inoculation
## or failure to recover and/or amplify trypanosomatid DNA from the abdomen
## T0 samples probably give better indicator of parasite replication
## vs comparing to Cq of 100K cultured cells

####################################

##########  Modeling     #########

#####################################
names(Samples.full)
glimpse(Samples.full)

model.test = 
  lm(Cq.median ~ Treatment, data = Samples.full)

summary(model.test)  

car::Anova(model.test)

##############
#Analysis on quantities #
p = 
  ggplot(Samples.full %>% filter(!is.na(Treatment)), 
              aes(x = Treatment, y = Copies.median))+
         geom_boxplot() +
         geom_point(position = "jitter") +
         ylab("Copy number") +
         xlab("Treatment (Parasite_Compound_Concentration)") +
  scale_y_continuous(trans = "log1p", 
                     breaks = 10^seq(2,7,by=1)) +
  theme_bw()
  #scale_x_discrete(labels = c())
p 

ggsave("./Figures/exploratory_boxplot_lotmaria_v0.2_2round.pdf",
       height = 4, width = 6)

## Round 2 update: >90% reduction in each of high-concentration treats
## 90% reduction in low-dose eugenol and 
## 2 orders of magnitude reduction in eugenol-high!
## Clean it up ##

#Regroup - cluster by treatment
levels(Samples.full$Treatment.o)
levels(Samples.full$Compound.o)

## Temporary: for first 4 plates:
#Compounds <- c("Carvacrol", "Eugenol", "Cinnamaldehyde")

Compounds <- c("None", "Carvacrol", "Eugenol", 
               "Cinnamaldehyde", "Thymol")

Xlab <- Compounds

## Only graph relevant samples
Samples.graph = 
  Samples.full %>% 
  filter(!is.na(Treatment), Parasite == "P")

p = 
  ggplot(Samples.graph %>% 
           filter(!is.na(Treatment), Parasite == "P"), 
         aes(x = Compound.o, y = Copies.median,
             #color = Parasite.o, 
             shape = Conc.o, fill = Conc.o))+
  geom_boxplot(aes(group = interaction(Compound.o, Conc.o)),
               outlier.shape = NA) +
  #geom_boxplot(position = "dodge") +
 # geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.3), 
             aes(shape = Conc.o),
             size = 1, alpha = 2/3) +
  # scale_color_discrete(name = "Inoculation\ntreatment",
  #                      labels = c("Sham", "Parasite"),
  #                      direction = -1)+
  scale_fill_discrete(name = "Chemical\nconcentration", 
                       labels = c("0", "Low", "High")) +
                       #labels = c("Low", "High")) +
  scale_shape_discrete(name = "Chemical\nconcentration", 
                       labels = c("0", "Low", "High")) +
                       #labels = c("Low", "High")) +
  xlab("Compound") +
  scale_x_discrete(labels = Compounds)+
  ylab(expression(Infection~intensity~(copies~bee^-1))) +
  scale_y_continuous(trans = "log1p", 
                     breaks = 10^seq(2,7,by=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5,
                                   vjust = 0.5, color = "black"))
 
p

## this is pretty cool-- 2 log units (100x) diff's!

## Add reference lines:

Pooled.median = median(Samples.graph$Copies.median, na.rm = TRUE)
Pooled.median

quantile(Samples.graph$Copies.median, c(25, 50, 75)/100)

T0.median = 
  T0 %>%
  filter(Label == "TIME 0 PAR") 
Median.t0 = median(T0.median$Copies.median)

Median.t0

notes.df = 
  tibble(Var = c(#"Quantity\ninoculated", 
                 "Time 0\nmedian", 
                 "Pooled\nmedian"),
         Copies.median = c(#10^5, # Removing this for now
           Median.t0, Pooled.median)) %>%
  mutate(Var.f = factor(Var, levels = unique(Var)))

levels(notes.df$Var.f)

p + geom_hline(data = notes.df, aes(yintercept = Copies.median,
                                    linetype = Var.f)) +
                                    #linetype = fct_rev(Var))) +
                 scale_linetype_discrete(name = NULL) +
  guides(fill = guide_legend(order = 1),
         shape = guide_legend(order = 1),
         linetype = guide_legend(order = 2))

ggsave("./Figures/lotmaria.boxplot.v2.1.two.round.pdf",
       height = 3, width = 5)

## Pooled median seems low-- lower than 8 of 10 individual medians
## Fixed by focusing on inoculated samples (Treatment == "P")

m.quant =
  glmmTMB:: glmmTMB(round(Copies.median) ~ Treatment, 
                    data = Samples.graph,
          family = 'nbinom2')

model_performance(m.quant)

?performance

## Simulate residuals ##
m.quant_simres <- simulateResiduals(m.quant)

#plot the results:
plot(m.quant_simres)

summary(m.quant)

plot(residuals(m.quant), fitted(m.quant))

## General checks ##
performance::check_model(m.quant)

testDispersion(m.quant) 

## not looking good

## Try poisson ##
m.p =  
  glmmTMB:: glmmTMB(round(Copies.median) ~ Treatment, 
                    data = Samples.graph,
                    family = 'poisson')
plot(simulateResiduals(m.p)) ## ???

check_model(m.p)

glimpse(Samples.full)
str(Samples.full$Treatment.o)

unique(Samples.full$Treatment.o)

Samples.full.model = 
  Samples.full %>%
  filter(!(Treatment == "0-S"))

m.lnorm =  
  glmmTMB:: glmmTMB(log1p(Copies.median) ~ Treatment,
                    data = Samples.full.model,
                    family = 'gaussian')
plot(simulateResiduals(m.lnorm)) ## much better

performance::check_model(m.lnorm) #Good
## Log-normal model much better than negative binomial
## and much much better than Poisson!

model_performance(m.lnorm, metrics = "all")

## Try adding actin as covariate ##
m.lnorm.actin =  
  glmmTMB:: glmmTMB(log1p(Copies.median) ~ Treatment + 
                      Cq.median_Actin,
                    data = Samples.full.model,
                    family = 'gaussian')
summary(m.lnorm.actin)
car::Anova(m.lnorm.actin) ## not significant 

## Add random effect of cage
names(Samples.full.model)

m.lnorm.actin.re = 
  update(m.lnorm.actin, ~. + (1|Cage))

m.lnorm.re = update(m.lnorm.actin.re, ~. -Cq.median_Actin )

m.lnorm.noreml = update(m.lnorm, REML = FALSE)

anova(m.lnorm, m.lnorm.actin)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# m.lnorm       10 1915.9 1954.6 -947.97   1895.9                         
# m.lnorm.actin 11 1917.5 1960.1 -947.77   1895.5 0.4103      1     0.5218

anova(m.lnorm.re, m.lnorm.noreml)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# m.lnorm.noreml 10 1915.9 1954.6 -947.97   1895.9                        
# m.lnorm.re     11 1917.9 1960.5 -947.97   1895.9     0      1          1

## Compare with lme4
library(lme4)
m.lnorm.lme4 =  
  glm(log1p(Copies.median) ~ Treatment,
                    data = Samples.full.model)

m.lnorm.lme4.re =  
  lmer(log1p(Copies.median) ~ Treatment + (1|Cage),
      data = Samples.full.model,
      REML = FALSE)
## Inclusion of cage as random effect led to convergence warnings

anova(m.lnorm.lme4.re, m.lnorm.lme4)

glimpse(Samples.full.model$Cage)
unique(Samples.full.model$Cage)
compare_performance(m.lnorm, m.lnorm.actin,
                    m.lnorm.re, m.lnorm.actin.re )

## neither actin nor random effects helpful
car::Anova(m.lnorm.actin.re) ## Actin p = 0.52
anova(m.lnorm.actin.re, m.lnorm.re)

summary(m.lnorm.actin.re)

summary(m.lnorm)

10^seq(0, 3, by = 1)

contrast(emmeans(m.lnorm, ~Treatment), "trt.vs.ctrl")

contrast(emmeans(m.lnorm.re, ~Treatment), "trt.vs.ctrl")

## Decide on model here ##

Emm = emmeans(m.lnorm, ~Treatment, type = "response") 
Emm.df = 
  Emm %>%
  as_tibble()
ggplot(Emm.df, aes(x = Treatment, y = exp(response))) +
  geom_pointrange(aes(ymin = exp(lower.CL), ymax = exp(upper.CL))) +
  scale_y_continuous(trans = "log1p", breaks = 10^seq(0, 5, by = 1)) +
  ylab(expression(Infection~intensity~(copies~bee^-1))) +
  theme_bw()

## How to get this in proper order?
Treatments.meta =
  tibble(Treatment = unique(Emm.df$Treatment),
         Compound = c("X", 
                      rep("C", 2),
                      rep("E", 2),
                      rep("N", 2),
                      rep("T", 2)),
         Conc = c(0, rep(c("High", "Low"), 4))) %>%
  mutate(Compound.o = factor(Compound,
                             levels = unique(Compound)),
         Conc.o = factor(Conc, levels = c(0, "Low", "High")))

Emm.df.merged = 
  Emm.df %>%
  left_join(Treatments.meta, by = "Treatment")

p.ci = 
  ggplot(Emm.df.merged, 
       aes(x = Compound.o, y = exp(response),
           shape = Conc.o)) +
  geom_pointrange(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                  position = position_dodge(width = 0.7)) +
  scale_y_continuous(trans = "log1p", breaks = 10^seq(0, 5, by = 1)) +
  ylab(expression(Infection~intensity~(copies~bee^-1))) +
  theme_bw() +
  scale_x_discrete(labels = Compounds) +
  scale_shape_discrete(name = "Concentration") +
  xlab("Compound") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 30, 
                                   hjust = 1,
                                   vjust = 1))
                                  # vjust = 0.5, color = "black"))
p.ci
ggsave("./Figures/lotmaria.emmeans.v1.0_2round.pdf", height = 3.5, width = 5)

# Emm.df.ordered = 
#   Emm.df %>%
#   mutate(Conc = 
#            ifelse(
#              str_detect(Treatment, "H"), "High",
#              ifelse(str_detect(Treatment, "L"), "Low",
#                                NA))) 
  ## Could go on to extra Compond with str_split, 2?
?str_trunc

############# Add annotations #############

## pairwise comparisons between each treatment and control ##
#?`contrast-methods`

treat.vs.control = 
  as_tibble(contrast(Emm, "trt.vs.ctrl"))
treat.vs.control

pairwise.starry =
  treat.vs.control %>%
  mutate(
    Stars = 
      ifelse(`p.value`< 0.001, "***",
             ifelse(`p.value` < 0.01, "**",
                    ifelse(`p.value` < 0.05, "*", "")))
  ) %>%
  rename(estimate.difference = estimate,
         SE.difference = SE,
         df.difference = df)

pairwise.starry

#mutate(COL1 = str_replace_all(COL1, "\\*|\\(|\\)", ""))

## Merge this back to emmeans data?
## Need to get rid of the parentheses
pairwise.starry.treat = 
  pairwise.starry %>%
  tidyr::separate(contrast, 
                  into = c("Treat", "Baseline"), 
                  sep = " - ") %>%
  ##clutzy code to replace parentheses
  mutate(Treatment = str_replace_all(Treat, "\\(|\\)", ""), 
         .before = Treat) %>%
  right_join(Emm.df.merged, by = "Treatment")

pairwise.starry.treat

## add sample size ##
names(Samples.model.graph)
names(pairwise.starry.treat)

Df.n = 
  Samples.model.graph %>%
  group_by(Treatment) %>%
  summarize(N = n())

Pairwise.starry.N = 
  pairwise.starry.treat %>%
  left_join(Df.n)

p.ci.text = 
  p.ci +
  geom_text(data = pairwise.starry.treat,
            aes(y = exp(upper.CL), label = Stars),
            size = 8,
            position = position_dodge(width = 0.7)) +
  geom_text(data = Pairwise.starry.N,
            aes(y = exp(lower.CL), 
              #y = 1, 
              label = N),
            #size = 8,
            position = position_dodge(width = 0.7),
            vjust = 1) 
p.ci.text

ggsave("./Figures/lotmaria.emmeans.v1.1a.stars.vjust1.pdf", 
       height = 3.5, width = 5)

## Adding raw data ##
# aes(x = Compound.o, y = exp(response),
#     shape = Conc.o))

Samples.model.graph = 
  Samples.full.model %>%
  mutate(Conc.o = 
           plyr::mapvalues(Conc.o,
                           from = c("L", "H"),
                           to = c("Low", "High")))
p.ci.raw =
  p.ci.text + 
  geom_point(data = Samples.model.graph,
             aes(y = (1 + Copies.median)),
             position = 
               position_jitterdodge(
                 dodge.width = 0.7,
                 jitter.width = 0.2),
             size = 1,
             alpha = 0.3)

p.ci.raw

ggsave("./Figures/lotmaria.emmeans.v1.2_addraw.pdf", 
       height = 3.5, width = 5)

## Let's show a table for this
## Columns for compound, concentration, mean, SE (or CI)
## Difference (vs 0-P), 
 ## optional: fold change  (i.e., exp(difference))
  ##SE(difference), t, p
pairwise.starry.treat
 ##next steps: work on table!