# This script will perform DE analysis comparing 2 seasons in PIAMA cohort
library(dplyr)
library(edgeR)
library(haven)
library(lubridate)
library(config)

# location 
conf <- config::get()

## upload PIAMA expression ##
load(file.path(conf$data_path, "Season/replication_PIAMA/PIAMA_expression_voom.Rdata"))
dim(data.piama)

PIAMA_dates <- read_sav(file.path(conf$data_path, "Season/replication_PIAMA/datesage16.sav")) %>%
  mutate(ID = as.character(ID))

pheno_PIAMA <- read.csv(file.path(conf$data_path, "Season/replication_PIAMA/piama_rnaseq_subjects.csv"), sep = ';') %>%
  mutate(ID_numb = substring(subject_ID, 1,5))

# nice, use of code reduces the need to repeat yourself. 
# other pieces of code might have benefited of it more.

### update IDs 
piama_ts_newid<- function(sampleid)
  ## transform old_sampleid to new_sampleid
{
  key<- read.csv(file.path(conf$data_path, "Season/replication_PIAMA/Keyid.csv"))
  nindex<-seq(1,length(sampleid))
  index<- match(sampleid,key[,1])
  index.c<- cbind(index,nindex); 
  index.c2<- na.omit(index.c)
  newid<- sampleid
  newid[index.c2[,2]]<- key[index.c2[,1],2] 
  return (newid)
}

#Functions
select.columns.in.order <- function(dataframe, columns) {
  dataframe[, columns]
}


IDs <- data.frame(RNA_IDs = colnames(data.piama))
IDs <- IDs%>%
  mutate(RNA_IDs_numb = substring(RNA_IDs,2,6))
IDs$RNA_IDs_new <- piama_ts_newid(IDs$RNA_IDs_numb)

length(intersect(PIAMA_dates$ID, IDs$RNA_IDs_new))

overlaped_samples <- (intersect(PIAMA_dates$ID, IDs$RNA_IDs_new))

PIAMA_dates<- PIAMA_dates%>%
  filter(ID %in% overlaped_samples)%>%
  left_join(IDs%>%
              dplyr::select(c(RNA_IDs, RNA_IDs_new)), by = c('ID'='RNA_IDs_new')) %>%
  mutate(date_yday = yday(pdate)) %>%
  mutate(new_seasons = as.factor(if_else(((date_yday >= 135) & (date_yday < 319)), "summer_autumn", "winter_spring")),
         week_day = weekdays(pdate))  

PIAMA_dates <- PIAMA_dates%>%
  left_join(pheno_PIAMA, by = c('ID'='ID_numb')) %>%
  mutate( asthma.status = if_else(asthma.status == "no asthma", "no_asthma", asthma.status))

#### differential expression ####
data.piama <- data.piama%>%
  as.data.frame()%>%
  select.columns.in.order(PIAMA_dates$RNA_IDs) %>%
  as.matrix()

mm <- model.matrix(~0 + new_seasons + gender + res.center + asthma.status, data = PIAMA_dates)
fit <- lmFit(data.piama, mm)
head(coef(fit))

contr <- makeContrasts(new_seasonswinter_spring - new_seasonssummer_autumn, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

write.csv(top.table, file.path(conf$data_path, "Season/Season_new_date/PIAMA_DE_seasons_sex_center_asthma.csv"))

## clinical characteristics 
## I understand, you too
library(tableone)

PIAMA_dates <- PIAMA_dates%>%
  mutate(gender = factor(gender, levels = c("girl", "boy")),
         asthma.status = factor(asthma.status, levels = c("no_asthma", "asthma")))

var_for_table <- c("asthma.status", "smoking.status", "gender")

tabtotal <- CreateTableOne(vars = var_for_table, strata = "new_seasons" , data = PIAMA_dates)

Final_statistics <- print(tabtotal, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)


