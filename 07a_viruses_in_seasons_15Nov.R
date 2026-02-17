## This script will check the frequency of the presence of the rhinovirus in the 2 seasons
library(dplyr)
library(stringr)
library(config)

# location 
conf <- config::get()

# load data
master.Table.ATLANTIS <- read.csv(file.path(conf$data_path, "Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"))

# save files from CZID
files <- list.files(path = file.path(conf$data_path, "Microbes/Output/nasal_brushes_score_check"), pattern="*.csv", full.names=TRUE, recursive=FALSE)

###### save samples which fits filtering conditions: 
##bpresence of enterovirus with z_score > 3 ##########
## presence of rhinovirus in enterovirus genus
samples <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE))%>%
    filter(max_z_score > 3)%>%
    filter(genus_tax_id == '12059') ## enterovirus 
  if (any(startsWith(filtered$name, 'Rhino'))){ #save only samples with Rhinoviruses present
    file_ID <- str_split(file, '/')[[1]][5]
    sample_ID <- str_split(file_ID, '_')[[1]][1]
    samples <- append (samples, sample_ID)
    print(filtered$name)
  }
}

## add data to master table ####
master.Table.ATLANTIS <- master.Table.ATLANTIS%>%
  mutate(Rhinovirus_presence = if_else(original_id %in% samples, 'yes', 'no'))

master.Table.ATLANTIS%>%
  group_by(new_seasons, Rhinovirus_presence) %>%
  summarise(n = n())

## chi-square 2 seasons - for Rhinovirus
table(master.Table.ATLANTIS$new_seasons, master.Table.ATLANTIS$Rhinovirus_presence)
chisq.test(master.Table.ATLANTIS$new_seasons, master.Table.ATLANTIS$Rhinovirus_presence, correct=FALSE)
## chi-square 4 seasons - for Rhinovirus
table(master.Table.ATLANTIS$season, master.Table.ATLANTIS$Rhinovirus_presence)
chisq.test(master.Table.ATLANTIS$season, master.Table.ATLANTIS$Rhinovirus_presence, correct=FALSE)


####### check all the viruses names 
## Check this out
## https://purrr.tidyverse.org/reference/map.html
## it's not just a for-loop under the hood ;p
## Also very easily paralizable with 
## https://furrr.futureverse.org/
names <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE))%>%
    filter(category == "viruses") 
  viral_names <- filtered$name 
  names <- append (names, viral_names)
}


## just use the console or store in an object if you want to use this. This clutters the script.
unique_names <- unique(names)
unique_names[str_detect(unique_names, "inf")]

unique_names[str_detect(unique_names, "adeno")]

unique_names[str_detect(unique_names, "boca")]

unique_names[str_detect(unique_names, "meta")]

#### check influenza ( winter from the papaer)####
samples_infl <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE)) %>%
    filter(max_z_score > 3) %>%
    filter(category == "viruses")
  #filter(str_detect(name, "influenzavirus"))
  if (any(str_detect(filtered$name, "influenzavirus"))){ #save
    file_ID <- str_split(file, '/')[[1]][5]
    sample_ID <- str_split(file_ID, '_')[[1]][1]
    samples_infl <- append (samples_infl, sample_ID)
    #print(filtered)
  }
}

## This is a high amount of repetition until line 140. This would be a perfect place for a function.
#### check adenovirus (all-year) ####
samples_adeno <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE)) %>%
    filter(max_z_score > 3) %>%  ## filter statements can be grouped with: &, |, !, xor() 
    filter(category == "viruses")
  #filter(str_detect(name, "influenzavirus")) ## remove dead code
  if (any(str_detect(filtered$name, "adenovirus"))){ #save   ## inline comments are discouraged because it is considered untidy
    file_ID <- str_split(file, '/')[[1]][5]
    sample_ID <- str_split(file_ID, '_')[[1]][1]
    samples_adeno <- append (samples_adeno, sample_ID)
  }
}

#### check bocavirus (all-year) ####
samples_boca <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE)) %>%
    filter(max_z_score > 3) %>%
    filter(category == "viruses")
  #filter(str_detect(name, "influenzavirus"))
  if (any(str_detect(filtered$name, "bocaparvovirus"))){ #save 
    file_ID <- str_split(file, '/')[[1]][5]
    sample_ID <- str_split(file_ID, '_')[[1]][1]
    samples_boca <- append (samples_boca, sample_ID)
  }
}

#### check metapneumovirus (all-year) ####
samples_metapneumovirus <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE)) %>%
    filter(max_z_score > 3) %>%
    filter(category == "viruses")
  #filter(str_detect(name, "influenzavirus"))
  if (any(str_detect(filtered$name, "metapneumovirus"))){ #save 
    file_ID <- str_split(file, '/')[[1]][5]
    sample_ID <- str_split(file_ID, '_')[[1]][1]
    samples_metapneumovirus <- append (samples_metapneumovirus, sample_ID)
  }
}

## add to master table

master.Table.ATLANTIS <- master.Table.ATLANTIS%>%
  mutate(influenza_presence = if_else(original_id %in% samples_infl, 'yes', 'no')) %>%
  mutate(adenovirus_presence = if_else(original_id %in% samples_adeno, 'yes', 'no')) %>%
  mutate(bocavirus_presence = if_else(original_id %in% samples_boca, 'yes', 'no')) %>%
  mutate(metapneumovirus_presence = if_else(original_id %in% samples_metapneumovirus, 'yes', 'no')) 

### to check all the viruses check the script: all_viruses_in_seasons_15Nov.py
