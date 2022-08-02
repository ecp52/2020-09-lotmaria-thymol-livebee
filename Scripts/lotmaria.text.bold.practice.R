library(tidyverse)
library(ggplot2)
data("diamonds")
summary(diamonds)
m1 = lm(price~color+cut, data = diamonds)
emm1 = emmeans::emmeans(m1, ~color + cut)
emm1.df = 
  as_tibble(emm1) %>%
  mutate(My.labs = 
           sample(letters, 
                  size = length(as_tibble(emm1)$emmean),
                            replace = TRUE))
pd = position_dodge(width = 0.8)

p.base = 
  ggplot(emm1.df, aes(x = color, y = emmean,
                    shape = cut)) +
         geom_pointrange(aes(ymin = lower.CL,
                             ymax = upper.CL),
                         position = pd) 
p.base + geom_text(aes(label = My.labs,
                       y = lower.CL),
                   position = pd,
                   fontface = "bold")
