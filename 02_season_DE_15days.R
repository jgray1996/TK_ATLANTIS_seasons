library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)
library(pheatmap)
library(ggrepel)
library(biomaRt)
library(config)


# location 
conf <- config::get()

# load data
master.Table <- read.csv(file.path(conf$data_path, "Season/Season_new_date/ATLANTIS_master_table_seasons_new.csv"), header = TRUE) %>%
  dplyr::select(-c("VISDAT_c", "month", "season",  "season_2", "Date", "X", "Customer.ID",
            "PAT_NO.y", "Requestion.y", "RQN", "DayYear.Maarten", "difference", "season.Maaike", "season_2_Maaike"))

# 1. create 2 groups 
master.Table <- master.Table %>%
  mutate(new_seasons = as.factor(if_else(((DayYear.Maaike >= 135) & (DayYear.Maaike < 319)), "summer_autumn", "winter_spring")))
table(master.Table$new_seasons)
table(master.Table$new_seasons, master.Table$asthma.status)

# change the order of the factors
master.Table$asthma.status <- factor(master.Table$asthma.status, levels = c("H","A"))

#write.csv(master.Table, "./Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv", row.names = F)

# 2. expression data
expression.data <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header = TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

# 3. DE analysis season + asthma - main model 
design <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status, data = master.Table)
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) ## reproduce previous version of edgeR
# Conduct genewise statistical tests for a given coefficient or contrast.
# contrasts <- limma::makeContrasts(season_binary = new_seasonswinter_spring - new_seasonssummer_autumn,
#                                  levels = design)
qlf_seasons <- edgeR::glmQLFTest(fit, coef = 2)
summary(decideTests(qlf_seasons))

# new_seasonswinter_spring
# Down                        918
# NotSig                    16229
# Up                          813

qlf_asthma <- edgeR::glmQLFTest(fit, coef = 7)
summary(decideTests(qlf_asthma))

# asthma.statusA
# Down              128
# NotSig          17689
# Up                143

# add gene names
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_ids <- rownames(DGEL$counts)
all_new_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

de_result_seasons <- topTags(qlf_seasons, n = nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

# write.csv(de_result_seasons, file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"), row.names = FALSE)

de_result_asthma <-  topTags(qlf_asthma, n = nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

# write.csv(de_result_asthma, file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma.csv"), row.names = FALSE)


## Volcano plot 
de_result_seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))

plt <- ggplot(data = de_result_seasons, aes(x=logFC, y = -log10(FDR))) +
  geom_point(aes(color= ifelse((FDR < 0.05)&(logFC>0), "#B4251A", ifelse((FDR < 0.05)&(logFC<0),"#163A7D", 'gray')))) + 
  geom_text_repel(aes(label=ifelse(((FDR < 0.05) & (logFC > 0.8)) | ((FDR < 0.05) & (logFC < -0.8)), 
                                   (ifelse((!is.na(external_gene_name)), external_gene_name, Gene)), ''),
                      hjust= 0.4, vjust= 0.5),
                  max.overlaps =50, 
                  direction = 'both',
                  alpha = 0.6, 
                  box.padding=0.3,
                  point.padding=0.5,
                  size = 4) +
  xlab (expression (log[2]~fold~change)) +
  ylab (expression(-log[10]~FDR))+
  scale_color_identity(name = '', breaks= c("#B4251A", "#163A7D",'gray') ,
                       labels = c('Higher expressed in winter-spring','Lower expressed in winter-spring', 'not significant'), guide = "legend")+
  geom_hline(yintercept = -log10(0.05), col="black", linewidth = 0.3, linetype = 2) + #look at table
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=12),
        legend.key = element_rect(fill = "white"))+
  xlim(-1.8,2.5)

# nice use of config object!
png(file.path(conf$data_path, "Season/Season_new_date/plots/Volcano_seasons_15Nov.png"), width = 1200, height = 800, res = 150)
print(plt)
dev.off()

saveRDS(plt, file.path(conf$data_path, "Season/Season_new_date/plots/Volcano_seasons_15Nov.rds"))

## 4. DE analysis - Check model corrected for season AND ciliated cells 

# add cell types:
master.Table.cell.types <- master.Table %>%
  left_join(cell_types <- read.csv(file.path(conf$data_path, "Season/cell_types/CIBERSORT.proportions.csv"), row.names = "X") %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column("Sample") , by = c("GenomeScan_ID" = "Sample"))


design <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status + `Multiciliated.lineage`, data = master.Table.cell.types)
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) ## reproduce previous version of edgeR


# The way you name your objects and variables is inconsistent. Use either _ or . but not different naming conventions.
qlf_seasons_ciliated_corrected <- edgeR::glmQLFTest(fit, coef = 2)
summary(decideTests(qlf_seasons_ciliated_corrected))
# new_seasonswinter_spring
# Down                        102
# NotSig                    17694
# Up                          164

qlf_ciliated <- edgeR::glmQLFTest(fit, coef = 8)
summary(decideTests(qlf_ciliated))

# Multiciliated.lineage
# Down                    7473
# NotSig                  4459
# Up                      6028

qlf_asthma_season_ciliated_corrected <- edgeR::glmQLFTest(fit, coef = 7)
summary(decideTests(qlf_asthma_season_ciliated_corrected))

# asthma.statusA
# Down              227
# NotSig          17388
# Up                345

de_result_asthma_season_ciliated_corrected <- topTags(qlf_asthma_season_ciliated_corrected, n = nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id")) %>%
  write.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_season_ciliated_cor.csv"), row.names = FALSE)
  

## 5. DE analysis - Check model corrected ONLY for ciliated cells 

# add cell types:
master.Table.cell.types <- master.Table %>%
  left_join(cell_types <- read.csv(file.path(conf$data_path, "Season/cell_types/CIBERSORT.proportions.csv"), row.names = "X") %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column("Sample") , by = c("GenomeScan_ID" = "Sample"))

# add cell types:
master.Table.cell.types <- master.Table %>%
  left_join(cell_types <- read.csv(file.path(conf$data_path, "Season/cell_types/CIBERSORT.proportions.csv"), row.names = "X") %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column("Sample") , by = c("GenomeScan_ID" = "Sample"))

# This code is becoming very repetative and can easily be called inside a function.
# Try no to repeat yourself too much.
design <- model.matrix(~`Multiciliated.lineage` + age + gender + smoking.status + asthma.status, data = master.Table.cell.types)
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) ## reproduce previous version of edgeR

qlf_asthma_ciliated_corrected <- edgeR::glmQLFTest(fit, coef = 7)
summary(decideTests(qlf_asthma_ciliated_corrected))


# Again, keep data and code seperate. You have published, but it's bad practice and can lead to data leaks.

# asthma.statusA
# Down              148
# NotSig          17438
# Up                218

de_result_asthma_ciliated_corrected <- topTags(qlf_asthma_ciliated_corrected, n = nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id")) %>%
  write.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_ciliated_cor.csv"), row.names = FALSE)


