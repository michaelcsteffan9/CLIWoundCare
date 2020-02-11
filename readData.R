#########################################################
## Libraries
#########################################################

library(tidyverse)
library(readxl)
library(here)

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

study <- read_xlsx(path = paste0(data.dir.path, "grid2NoMRN.xlsx"))
demo <- read_xlsx(path = paste0(data.dir.path, "masterPtFileNoMRN.xlsx"))
mortality <- read_xlsx(path = paste0(data.dir.path, "mortNoMRN.xlsx"))
procedures <- read_xlsx(path = paste0(data.dir.path, "procGroupedNoMRN.xlsx")) 
ep1.procedures <- read_xlsx(path = paste0(data.dir.path, "procLatNoMRN.xlsx")) 
