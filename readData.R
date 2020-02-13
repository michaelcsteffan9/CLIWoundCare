#########################################################
## Libraries
#########################################################

library(tidyverse)
library(readxl)
library(lubridate)
library(here)
here <- here::here
#########################################################
## Data
#########################################################

# You will need to change this to your path
# So that we don't have to comment out other people's directories, I'm putting some conditional logic
cur.directory <- here()
if(str_detect(string = cur.directory, pattern = "jsperger")){
  data.dir.path <- "/Users/jsperger/Dropbox/Consulting Project/Data/"
}
if(str_detect(string = cur.directory, pattern = "nikki")){
  data.dir.path <- "Fix me/Dropbox/Consulting Project/Data/"
}
if(str_detect(string = cur.directory, pattern = "haley")){
  data.dir.path <- "Fix me/Dropbox/Consulting Project/Data/"
}
if(str_detect(string = cur.directory, pattern = "michael")){
  data.dir.path <- "Fix me/Dropbox/Consulting Project/Data/"
}
if(exists("data.dir.path") == FALSE) stop("Data directory not specified")

disease <- read_xlsx(path = paste0(data.dir.path, "grid2NoMRN.xlsx"))
demo <- read_xlsx(path = paste0(data.dir.path, "masterPtFileNoMRN.xlsx"))
mortality <- read_xlsx(path = paste0(data.dir.path, "mortNoMRN.xlsx"))
procedures <- read_xlsx(path = paste0(data.dir.path, "procGroupedNoMRN.xlsx")) 
ep1.procedures <- read_xlsx(path = paste0(data.dir.path, "procLatNoMRN.xlsx")) 

# Exclude people who have wounds in both legs
pts.to.exclude <- study %>% pull(ptNum) %>% subset(duplicated(.))

study.data <- left_join(x = demo %>% filter(!ptNum %in% pts.to.exclude),
                        y = disease,
                        by = c("ptNum")) %>% 
  left_join(., 
            y = mortality %>% select("ptNum", "deathDate", "timeToDeathFromStudyStart","mort6mo", "mort1yr"), 
            by = c("ptNum"), suffix = c("_demo","_mort")) %>% select(ptNum, everything())

# Define 2 year mortality; using a 365 day year to define  
study.data$mort2yr <- if_else(study.data$timeToDeathFromStudyStart <= 700, 
                              TRUE, FALSE, missing = FALSE)
