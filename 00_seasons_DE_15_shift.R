library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)
library(pheatmap)
library(ggrepel)
library(lubridate)
library(config)

# location 
conf <- config::get()

# load files
master.Table <- read.csv(file.path(conf$data_path,"/Season/Season_new_date/ATLANTIS_master_table_seasons_new.csv"))

# expression data

## The new standard pipe for R is: |> 
# This is built into R and is supported by every function and does not need dplyr to be imported
# Before any pipe place a spacebar. Not a hard rule, but the tidyverse team suggests this as a way 
# to keep the code tidy. At leas pick one or the other. Either use the space or dont.
expression.data <- read.csv(file.path(conf$data_path,"/Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header = TRUE) %>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

## define first day of each month 
st <- as.Date("2015-01-01")
en <- as.Date("2015-12-01")
first_month <- seq(st + 1, en + 1, by = "1 month") - 1
yday(first_month)

## add 15th of each month
breaks_df <- c()
for (day_n in first_month) {
  print(yday(as.Date(day_n)))
  breaks_df <- c(breaks_df, yday(as.Date(day_n)), yday(as.Date(day_n)) + 14)
}

# shift groups by 1 step - either 1st of the month or 15th
n_deg <- c()
results_df <- data.frame(
  iteration = integer(),
  n_down = integer(),
  n_up = integer(),
  n_total = integer(),
  stringsAsFactors = FALSE
)

# huge block of code in a for loop can benefit from the use of functions to keep your code neat and maintainable.
for (i in 1:(length(breaks_df)/2)) {
  group_1 = c(breaks_df[i],breaks_df[(i+12)]) # Space after every operating character x + y, space after comma 
  print(group_1)
  master.Table <- master.Table %>%
    mutate(groups = if_else(((DayYear.Maaike >= group_1[1]) & (DayYear.Maaike < group_1[2])), "group1", "group2"))
  print(table(master.Table$groups))
  design <- model.matrix(~0 + groups + age + gender + smoking.status + asthma.status, data = master.Table)
  DGEL <- edgeR::DGEList(expression.data)
  keep <- edgeR::filterByExpr(DGEL, design) 
  DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
  DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
  DGEL <- edgeR::estimateDisp(DGEL, design)
  fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE)
  #Conduct genewise statistical tests for a given coefficient or contrast.
  contrasts <- limma::makeContrasts(season_binary = groupsgroup1 - groupsgroup2,
                                    levels = design)
  qlf <- edgeR::glmQLFTest(fit, contrast = contrasts[,"season_binary"])
  summary <- summary(decideTests(qlf))
  n_down <- summary[1,1] # number of downregulated genes 
  n_up <- summary[3,1] ## number of upregulated genes 
  n_total <- n_down + n_up # total 
  results_df[i, ] <- list(i, n_down, n_up, n_total)
  n_deg[i] <- n_total
  print(total_n_diff_genes)
  }

# Nice use of portable paths
save(n_deg, results_df, file = file.path(conf$data_path,"Season/Season_new_date/15days_shift.Rdata"))

# plot N of genes vs date 

load(file.path(conf$data_path,"Season/Season_new_date/15days_shift"))
dates_genes_df <- data.frame(Date = as.Date(breaks_df[1:length(n_deg)]-1, origin = "2015-01-01"),
                             n_gene = n_deg)


ggplot(dates_genes_df, aes(Date, n_gene)) +
  geom_point() +
  geom_line() 
  

# make a circular plot 
dates_genes_df_polar <- data.frame(
  Date = as.Date(breaks_df -1, origin = "2015-01-01"),
  n_gene = n_deg)

dates_genes_df_polar <- rbind(dates_genes_df_polar, 
                              data.frame(Date = as.Date("2016-01-01"), n_gene = 45)) %>%
  mutate(date_month = format(Date, "%d-%b"))



ggplot(dates_genes_df_polar, aes(Date, n_gene)) +
  geom_point() +
  geom_line() +
  geom_rect(aes(xmin = as.Date("2015-05-15"), xmax = as.Date("2015-11-15"),
                ymin = 0, ymax = 1800), fill = alpha("#FFA559", 0.01), color = "black") +
  geom_rect(aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2015-05-15"),
                ymin = 0, ymax = 1800), fill = alpha("#BFCCB5", 0.05), color = "black") +
  geom_rect(aes(xmin = as.Date("2015-11-15"), xmax = as.Date("2015-12-31"),
                ymin = 0, ymax = 1800), fill = alpha("#BFCCB5", 0.05), color = "black") 
  

# Make sure to remove dead code before committing or keep it in with a conditional.  

round_date_genes_plt <- ggplot(dates_genes_df_polar, aes(Date, n_gene)) +
  # geom_rect(aes(xmin = as.Date("2015-05-15"), xmax = as.Date("2015-11-15"),
  #               ymin = 0, ymax = 1800), fill = alpha("#FFA559", 0.03), color = alpha("#FFA559", 0.0)) +
  # geom_rect(aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2015-05-15"),
  #               ymin = 0, ymax = 1800), fill = alpha("#BFCCB5", 0.03), color = alpha("#BFCCB5", 0.0)) +
  # geom_rect(aes(xmin = as.Date("2015-11-15"), xmax = as.Date("2016-01-01"),
  #               ymin = 0, ymax = 1800), fill = alpha("#BFCCB5", 0.03), color = alpha("#BFCCB5", 0.0)) +
  geom_rect(aes(xmin = as.Date("2015-05-15"), xmax = as.Date("2015-11-15"),
                ymin = 0, ymax = 1800), fill = alpha("#FFA559", 0.03), color = alpha("#FFA559", 0.0)) +
  geom_rect(aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2015-05-15"),
                ymin = 0, ymax = 1800), fill = alpha("#BFCCB5", 0.03), color = alpha("#BFCCB5", 0.0)) +
  geom_rect(aes(xmin = as.Date("2015-11-15"), xmax = as.Date("2016-01-01"),
                ymin = 0, ymax = 1800), fill = alpha("#BFCCB5", 0.03), color = alpha("#BFCCB5", 0.0)) +
  
  #scale_y_continuous(limits = c(2000, NA)) +
  geom_point() +
  geom_line() +
  ylim(0,2000) +
  scale_x_continuous(breaks = dates_genes_df_polar$Date[1:length(dates_genes_df_polar$Date)-1], 
                     labels = dates_genes_df_polar$date_month[1:length(dates_genes_df_polar$date_month)-1])+
  coord_polar(start = pi/4 , direction = 1) +
  theme_light() +
  ylab("N of significantly differentially expresed genes") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=8),
        legend.key = element_rect(fill = "white")) 

png(file.path(conf$data_path,"Season/Season_new_date/plots/circled_date_genes.png"), res = 150,  width = 1200, height = 800)
print(round_date_genes_plt)
dev.off()

fake_legend_df <- data.frame(x = c(rnorm(300, -3, 1.5),
                       rnorm(300, 0, 1)),
                 group = c(rep("winter_spring", 300),
                           rep("summer-autumn", 300)))
legend <- ggplot(fake_legend_df, aes(x = x, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("winter_spring" = "#BFCCB5", "summer-autumn" = "#FFA559")) +
  guides(fill = guide_legend(title = "season"))

png(file.path(conf$data_path, "./Season/Season_new_date/plots/legend_circled_date_genes.png"), res = 150,  width = 1200, height = 800)
print(legend)
dev.off()

