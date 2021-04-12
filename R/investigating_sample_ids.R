library(tidyverse)
library(readxl)
'%ni%' <- Negate('%in%')

gbs1 <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 1,skip = 1) %>% 
  select(1:3) %>% 
  mutate(country="USA")
names(gbs1) <- names(gbs1) %>% janitor::make_clean_names() %>% str_split("_") %>% map_chr(1)

gbs2 <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 1,skip = 1) %>% 
  select(5:7) %>% 
  mutate(country="China")
names(gbs2) <- names(gbs1) %>% janitor::make_clean_names() %>% str_split("_") %>% map_chr(1)

gbs3 <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 1,skip = 1) %>% 
  select(9:11) %>% 
  mutate(country="Japan")
names(gbs3) <- names(gbs1) %>% janitor::make_clean_names() %>% str_split("_") %>% map_chr(1)

gbs4 <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 1,skip = 1) %>% 
  select(13:15) %>% 
  mutate(country="Thailand")
names(gbs4) <- names(gbs1) %>% janitor::make_clean_names() %>% str_split("_") %>% map_chr(1)

gbs <- rbind(gbs1,gbs2,gbs3,gbs4)
rm(gbs1,gbs2,gbs3,gbs4)


loc <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 2) %>% 
  select(1:2)

writeLines(bams[bams %in% loc$`DNA ID`],"~/Desktop/DNA_IDs_IN_Spreadsheet.txt")
writeLines(bams[bams %ni% loc$`DNA ID`],"~/Desktop/DNA_IDs_NOTIN_Spreadsheet.txt")


s3 <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 3)
s4 <- read_xlsx("~/Desktop/GBS_Kudzu_Everything.xlsx",sheet = 4)



bams <- read_lines("~/Desktop/existant_bamfiles.txt")
