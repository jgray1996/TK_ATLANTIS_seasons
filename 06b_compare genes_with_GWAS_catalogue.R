## This script will overlap asthma genes from GWAS with genes from DE analysis

library(dplyr)
library(config)

# location 
conf <- config::get()

# 1. load GWAS catalogue - only asthma related
GWAS_catalogue <- read.csv(file.path(conf$data_path, "Multiomics/mofa_results/GWAS_catalogue/gwas_catalog_asthma.csv")) %>%
  filter(P.VALUE < 1e-08)
reported_gene_id <- unique(GWAS_catalogue$REPORTED.GENE.S.)
cleaned_reported_gene_id <- unlist(strsplit(reported_gene_id, split = ","))
cleaned_reported_gene_id <- unique(trimws(cleaned_reported_gene_id))

# 2. load DE genes
de.result.asthma <- read.csv(file.path(conf$data_path, "Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv"))
de.result.asthma.seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma.csv"))
de.result.asthma.ciliated <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_ciliated_cor.csv"))
de.result.asthma.seasons.ciliated <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_season_ciliated_cor.csv"))

# Gijs' comment: If you don't want to use the output anywhere else in your script, use a printing statement or save it in a variable if you do.
overlap.a <- length(intersect(cleaned_reported_gene_id, de.result.asthma[de.result.asthma$FDR < 0.05, "hgnc_symbol"]))
overlap.a/length(de.result.asthma[de.result.asthma$FDR < 0.05, "hgnc_symbol"])

overlap.a.s <- length(intersect(cleaned_reported_gene_id, de.result.asthma.seasons[de.result.asthma.seasons$FDR < 0.05, "external_gene_name"]))
overlap.a.s/length(de.result.asthma.seasons[de.result.asthma.seasons$FDR < 0.05, "external_gene_name"])

overlap.a.c <- length(intersect(cleaned_reported_gene_id, de.result.asthma.ciliated[de.result.asthma.ciliated$FDR < 0.05, "external_gene_name"]))
overlap.a.c/length(de.result.asthma.ciliated[de.result.asthma.ciliated$FDR < 0.05, "external_gene_name"])

overlap.a.s.c <- length(intersect(cleaned_reported_gene_id, de.result.asthma.seasons.ciliated[de.result.asthma.seasons.ciliated$FDR < 0.05, "external_gene_name"]))
overlap.a.s.c/length(de.result.asthma.seasons.ciliated[de.result.asthma.seasons.ciliated$FDR < 0.05, "external_gene_name"])

