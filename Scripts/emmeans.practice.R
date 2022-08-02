## Practice with emmeans ##
#'faking' transformations using 'regrid'
#https://mran.microsoft.com/snapshot/2017-12-11/web/packages/emmeans/vignettes/transformations.html#after
#'This allows us to get ratios from log-normal model

library(emmeans)
library(tidyverse)

data(pigs)
?pigs
glimpse(pigs)

pigroot.lm <- lm(sqrt(conc) ~ source + factor(percent), 
                 data = pigs)
summary(pigroot.lm)

piglog.emm.0 = 
  emmeans(pigroot.lm, "source")

piglog.emm.0

plot(piglog.emm.0) ## on square-root scale

piglog.emm.tf = 
  emmeans(pigroot.lm, "source", type = "response")
 ## clever emmeans figures out that response variable 
  ## was square-rooted 

plot(piglog.emm.tf)
pairs(piglog.emm.tf, type = "response")
#Use 'contrast(regrid(object), ...)' 
#to obtain contrasts of back-transformed estimates

piglog.emm.s <- regrid(emmeans(pigroot.lm, "source"), 
                       transform = "log")

confint(piglog.emm.s, type = "response")

plot(piglog.emm.s, type = "response")

pairs(piglog.emm.s, type = "response")
 ## all pairs
contrast(piglog.emm.s, "trt.vs.ctrl", type = "response")
  ## select pairs-- treatment vs reference

## this is clever, gives ratios ##
## we should use this for reporting changes in infection

## start with log-normal model
## use 'regrid' to fake a log transformation

p = ggplot(pigs, aes(x = source, y = conc)) +
  geom_point() + ggtitle("Pigs")
## hjust and vjust

## hjust = 0 to left-align (start at x = 0)
p + theme(axis.text.x = element_text(hjust = 0))
## left-justified

p +
  theme(plot.title = element_text(hjust = 0)) 
## aligned on left margin

## hjust = 0 to left-align (start at x = 0)
p + 
  theme(plot.title = element_text(hjust = 1)) 

p + 
  theme(plot.title = element_text(hjust = 1)) 

##vjust = 0 to bottom-align, margin at y = 0
p + 
  theme(axis.text.x = element_text(vjust = 0)) 

## vjust = 1 to top-align
p + 
  theme(axis.text.x = element_text(vjust = 1)) +
  ## hjust = 0 to left-align
  theme(axis.title.y = element_text(hjust = 0)) +
  theme(axis.text.y = element_text(vjust = -5))
  ## this puts top margin at lower value of y

###
## Try another model with log link ##
m.log = lm(
  #log1p(conc) ## not recognized
  log(conc + 1) ~ source + percent,
           data = pigs)

summary(m.log)

Emm.model = 
  emmeans(m.log, ~source, type = "response")

Emm.model ## auto-detects transformation

pairs(Emm.model) ## gives us the ratios

## What if we need ratios from model w/o log link?
m.lin = 
  lm(conc  ~ source + percent,
     data = pigs)

emm.m.lin = 
  emmeans(m.lin, ~source)
pairs(emm.m.lin)
# contrast    estimate   SE df t.ratio p.value
# fish - soy     -9.48 2.27 25 -4.182  0.0009 
# fish - skim   -15.62 2.34 25 -6.672  <.0001 
# soy - skim     -6.14 2.30 25 -2.668  0.0341

## type = "response" does not change pairwise tests
emm.m.lin.tf = 
  emmeans(m.lin, ~source, type = "response")
pairs(emm.m.lin.tf)
# contrast    estimate   SE df t.ratio p.value
# fish - soy     -9.48 2.27 25 -4.182  0.0009 
# fish - skim   -15.62 2.34 25 -6.672  <.0001 
# soy - skim     -6.14 2.30 25 -2.668  0.0341

## use regrid() to get ratios

## First use link scale emmeans as input
emm.m.lin

m.lin.emm.s <- 
  regrid(emmeans(emm.m.lin, "source"), transform = "log")

pairs(m.lin.emm.s)

# contrast    estimate     SE df t.ratio p.value
# fish - soy    -0.280 0.0687 25 -4.083  0.0011 
# fish - skim   -0.428 0.0668 25 -6.396  <.0001 
# soy - skim    -0.147 0.0552 25 -2.662  0.0345 

## close but not identical to tests without regrid()


## Does it change if I use backtransformed emmeans as input?

m.lin.emm.s.tf <- 
  regrid(emmeans(emm.m.lin.tf, "source"), transform = "log")
m.lin.emm.s.tf
pairs(m.lin.emm.s.tf)

# contrast    ratio     SE df t.ratio p.value
# fish / soy  0.755 0.0519 25 -4.083  0.0011 
# fish / skim 0.652 0.0436 25 -6.396  <.0001 
# soy / skim  0.863 0.0477 25 -2.662  0.0345
## slightly different from tests without regrid()

library(lme4)
