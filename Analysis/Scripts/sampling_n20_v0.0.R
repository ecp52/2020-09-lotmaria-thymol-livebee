### 2021-10 ##
## SCript to randomly sample 20 bees per treatment group ##
## update 2021-11 additional 20 bees/group
setwd("C:/Users/Evan.Palmer-Young/OneDrive/2020-09-lotmaria-thymol-livebee/Analysis")

library(tidyverse)
library(readxl)
library(openxlsx)

## Data in ##
#2020-09-lotmaria-thymol-livebee_data_v1_gdrive_dissections.xlsx #
D = 
  read_excel(path = "../2020-09-lotmaria-thymol-livebee_data_v1_gdrive_dissections.xlsx",
             sheet = "Filtered_Data")

glimpse(D)

# Sampling ##
#https://stackoverflow.com/questions/64587692/sample-n-random-rows-per-group-in-a-dataframe-with-dplyr-when-some-observations

D.n =
  D %>%
  group_by(Parasite, Compound, Conc) %>%
  count()
D.n

# Parasite Compound Conc      n
# <chr>    <chr>    <chr> <int>
#   1 P        C        H        38
# 2 P        C        L        40
# 3 P        E        H        41
# 4 P        E        L        42
# 5 P        N        H        38
# 6 P        N        L        40
# 7 P        T        H        40
# 8 P        T        L        39
# 9 P        X        0        41
# 10 S        X        0        40

Nsize <- 20

Sample.20 =
  D %>% 
  group_by(Parasite, Compound, Conc) %>%
  slice_sample(n = Nsize)

# Samples.to.dissect =
#   write_excel_csv(Sample.20, file = "../Dissection_list.csv") 
  

## Get bee numbers ##
D.numbers = 
  read_excel(path = "../2020-09-lotmaria-thymol-livebee_data_v1_gdrive_dissections.xlsx",
             sheet = "Dissections")

glimpse(D.numbers)

D.numbers.trim = 
  D.numbers %>%
  dplyr::select(Tube, Label)

D.merged = 
  Samples.to.dissect %>%
  left_join(D.numbers.trim) %>%
  arrange(Tube) %>%
  ungroup %>%
  dplyr::select(Tube, Label)

glimpse(D.merged)

# for writing a data.frame or list of data.frames to an xlsx file
#write.xlsx(D.merged, file = "../Samples.to.extract.xlsx")

############
## Samples for round 2 ##

## Workflow:
## Get list of samples from round 1
## exclude these from full data frame of bees frozen alive at 8 d  ("D")

## Get samples from round 1 ##
df.r1 = 
  read_excel("../Samples.to.extract.xlsx")
glimpse(df.r1)

## Lindsey added some stuff here, so just take the relevant samples
df.r1.original = 
  df.r1[1:200,] %>%
  dplyr::select(Tube, Label)

## Get list of all samples ##
glimpse(D)

##
d.remaining = 
  D %>%
  filter(!(Label %in% unique(df.r1.original$Label)))

## Get corresponding tube numbers ##
length(unique(d.remaining$Label))

D.merged.r2 = 
  d.remaining %>%
  left_join(D.numbers.trim) %>%
  arrange(Tube) %>%
  ungroup %>%
  dplyr::select(Tube, Label) %>%
  filter(Label %in% unique(d.remaining$Label)) 

## Why the duplicate? ##
duplicates = 
  D.merged.r2 %>%
  group_by(Label) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n > 1)

# Label        n
# <chr>    <int>
#   1 P-X-1-4      2

duplicate.check = 
  D.merged.r2 %>% 
  filter(Label %in% unique(duplicates$Label)) %>%
  dplyr::select(Tube, Label)

duplicate.check
# Tube Label  
# <dbl> <chr>  
#   1   340 P-X-1-4
# 2   581 P-X-1-4

## Same label for tubes 340, 581 ##
## Must be 2 bees (2 tubes) with this label? ##

openxlsx::write.xlsx(D.merged.r2,
                     file = "../Samples.to.extract.round2.xlsx")

## We didn't find any parasites in the sham-inoculated samples, so
 ## write out another data frame without shams
d.r2.no.sham = 
  D.merged.r2 %>%
  filter(!is.na(Tube)) %>%
  filter(!grepl('S', Label))

openxlsx::write.xlsx(d.r2.no.sham,
                     file = "../Samples.to.extract.round2.no.sham.xlsx")

## 




  
