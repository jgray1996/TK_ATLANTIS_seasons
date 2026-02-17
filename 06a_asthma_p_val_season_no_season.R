## This script will compare p-values which were obtained 
## 1. For asthma-associated genes in the model corrected and not corrected for season


library(duckplyr)
library(grid)
library(ggplot2)
library(config)

# location 
conf <- config::get()

######## Compare p-values ###########

# Create one table with DE statistics out of 2 
# Arguments: 2 tables resulting from DE analysis


# Nicely used through the script.
create_compare_df <- function(df_one, df_two, name_one, name_two) {
  overlap_genes <- intersect(df_one$Gene, df_two$Gene)
  comparison_df <- df_one %>%
    filter(Gene %in% overlap_genes) %>%
    dplyr::select(c(Gene, logFC, PValue, FDR)) %>%
    rename_with(~ paste0(name_one, "_", .x), all_of(c("logFC", "PValue", "FDR"))) %>%
    left_join(df_two %>%
                dplyr::select(c(Gene, logFC, PValue, FDR)), by = c("Gene" = "Gene")) %>%
    rename_with(~ paste0(name_two, "_", .x), all_of(c("logFC", "PValue", "FDR"))) %>%
    mutate(
      significant = case_when(
        .data[[paste0(name_one, "_FDR")]] < 0.05 & .data[[paste0(name_two, "_FDR")]] < 0.05 ~ "both models",
        .data[[paste0(name_one, "_FDR")]] < 0.05 & .data[[paste0(name_two, "_FDR")]] >= 0.05 ~ paste(name_one, "in model"),
        .data[[paste0(name_one, "_FDR")]] >= 0.05 & .data[[paste0(name_two, "_FDR")]] < 0.05 ~ paste(name_two, "in model"),
        TRUE ~ "none"
      )
    )
      
  return(comparison_df)
}
        
## 1. Compare models with and without seasons 

# inconsistent use of naming conventions 
de.result.asthma.seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma.csv"))
de.result.asthma <- read.csv(file.path(conf$data_path, "Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv"))

comparison_df <- create_compare_df(de.result.asthma, de.result.asthma.seasons, "asthma", "season_asthma")

# Check whether all genes from asthma design are in asthma+seasons design
comparison_only_asthma <- comparison_df%>%
  filter(asthma_FDR < 0.05)

# check number of different types of genes: 
comparison_df_sign <- comparison_df %>%
  filter((asthma_FDR < 0.05) | (season_asthma_FDR < 0.05)) %>%
  mutate(FDR_direction = (if_else(season_asthma_FDR < asthma_FDR, "smaller in SA", "bigger in SA")))

plt <- ggplot(data=comparison_df%>%
                filter(season_asthma_FDR < 0.05 | asthma_FDR<0.05), aes(x=-log10(season_asthma_FDR), y=-log10(asthma_FDR))) +
  #geom_point(aes(color = ifelse((Gene %in% de.results.asthma[de.results.asthma$FDR < 0.05,]$Gene), 'red', 'black')))+
  geom_point(aes(color = significant), size = 1.2) +
  scale_color_manual(values=c("#163A7D", "#B4251A", "#AFABAB"))+
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlab(expression(-log[10]~"FDR DE genes in asthma with season as a confounder")) +
  ylab(expression(-log[10]~"FDR DE genes in asthma")) +
  theme_classic() +
  #scale_x_break(breaks = c(2.8, 5.1)) +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text = element_text(size = rel(1.5)), 
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)))
png(file.path(conf$data_path, "Season/Season_new_date/plots/FDR_asthma_vs_asthma+season_dot_plot.png"), width = 1400, height = 800, res = 150)
print(plt)
dev.off()

## plot for the paper
comparison_df_sign <- comparison_df %>%
  filter((asthma_FDR < 0.05) | (season_asthma_FDR < 0.05)) %>%
  mutate(FDR_direction = (if_else(season_asthma_FDR < asthma_FDR, "smaller in SA", "bigger in SA")))

comparison_df_sign_longer <- comparison_df_sign %>%
  dplyr::select(c('Gene', 'asthma_FDR', 'season_asthma_FDR')) %>%
  tidyr::gather(model, FDR, -Gene) %>%
  mutate(model = factor(model, levels = c("asthma_FDR", "season_asthma_FDR"))) %>%
  #mutate(model = as.factor(model)) %>%
  mutate(model = if_else(model == "asthma_FDR", "model not accounting for seasons\nn = 126",
                         paste0("model accounting for seasons\nn = 271"))) %>%
  mutate(model = factor(model, levels = c("model not accounting for seasons\nn = 126", 
                                          "model accounting for seasons\nn = 271")))

plt <- ggplot(comparison_df_sign_longer, aes(x = model, y = FDR)) +
  geom_line(aes(group = Gene), col = "gray", linewidth = 0.1) +
  geom_point(aes(color = model), size = 1.5) +
  geom_hline(yintercept = 0.05, linetype='dashed') +
  #facet_wrap(~ ARM) +
  scale_color_manual(values=c("#163A7D","#B4251A"), name ='Model for DE analysis',
                     labels = c('~0 + asthma status + age + gender + smoking status',
                                '~0 + asthma status + age + gender + smoking status + season'))+
  xlab("asthma-associated DEGs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size =0.6),
        #axis.title.x=element_blank(),
        axis.text = element_text(face = "plain", size = 10, colour = 'black'),
        axis.title = element_text(face = "plain", size = 10, colour = 'black'),
        #plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        legend.position = "none")+
  guides(color = guide_legend(nrow=2,byrow=TRUE))

png(file.path(conf$data_path, "Season/Season_new_date/plots/FDR_asthma_vs_asthma+season_paired_plot.png"), width = 1000, height = 800, res = 150)
print(plt)
dev.off()

save(plt, file = file.path(conf$data_path, "Season/Season_new_date/plots/FDR_asthma_vs_asthma+season_paired_plot.Rdata"))

## 2. Compare models with seasons AND ciliated vs without seasons and ciliated

de.result.asthma.seasons.ciliated <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_season_ciliated_cor.csv"))
de.result.asthma <- read.csv(file.path(conf$data_path, "./Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv"))
comparison_df <- create_compare_df(de.result.asthma, de.result.asthma.seasons.ciliated, "asthma", "season_ciliated_asthma")

#plot 
plt <- ggplot(comparison_df %>%
                filter(asthma_FDR < 0.05 | season_ciliated_asthma_FDR<0.05), 
              aes(x = -log10(season_ciliated_asthma_FDR), y=-log10(asthma_FDR))) +
  #geom_point(aes(color = ifelse((Gene %in% de.results.asthma[de.results.asthma$FDR < 0.05,]$Gene), 'red', 'black')))+
  geom_point(aes(color = significant), size = 1.2) +
  scale_color_manual(values=c("#163A7D", "#B4251A", "#AFABAB"))+
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlab(expression(-log[10]~"FDR DE genes in asthma corrected for season and ciliated")) +
  ylab(expression(-log[10]~"FDR DE genes in asthma")) +
  theme_classic() +
  #scale_x_break(breaks = c(2.8, 5.1)) +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text = element_text(size = rel(1.5)), 
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)))
png(file.path(conf$data_path, "Season/Season_new_date/plots/FDR_asthma_vs_asthma+season+ciliated_dot_plot.png"), width = 1400, height = 800, res = 150)
print(plt)
dev.off()

## 3. Compare models with seasons vs with seasons AND ciliated 
comparison_df <- create_compare_df(de.result.asthma.seasons, de.result.asthma.seasons.ciliated, "asthma_seasons", "season_ciliated_asthma")

#plot <- just call the variable plot. No need to abbreviate if a self explaining variable name is enough. if you want to save a character, then don't comment on it.
plt <- ggplot(comparison_df %>%
                filter(asthma_seasons_FDR < 0.05 | season_ciliated_asthma_FDR < 0.05), 
              aes(x = -log10(season_ciliated_asthma_FDR), y=-log10(asthma_seasons_FDR))) +
  #geom_point(aes(color = ifelse((Gene %in% de.results.asthma[de.results.asthma$FDR < 0.05,]$Gene), 'red', 'black')))+
  geom_point(aes(color = significant), size = 1.2) +
  scale_color_manual(values=c("#163A7D", "#B4251A", "#AFABAB"))+
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlab(expression(-log[10]~"FDR DE genes in asthma corrected for season and ciliated")) +
  ylab(expression(-log[10]~"FDR DE genes in asthma corrected for season")) +
  theme_classic() +
  #scale_x_break(breaks = c(2.8, 5.1)) +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text = element_text(size = rel(1.5)), 
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)))
png(file.path(conf$data_path, "Season/Season_new_date/plots/FDR_asthma+season_vs_asthma+season+ciliated_dot_plot.png"), width = 1400, height = 800, res = 150)
print(plt)
dev.off()

## 4. Compare models with ciliated vs with seasons AND ciliated 
de.result.asthma.ciliated <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_ciliated_cor.csv"))
comparison_df <- create_compare_df( de.result.asthma.ciliated, de.result.asthma.seasons.ciliated, "asthma_cilited", "season_ciliated_asthma") %>%
  left_join(de.result.asthma.ciliated %>%
              dplyr::select(c(Gene, external_gene_name)), by = c("Gene" = "Gene"))

# plot 
plt <- ggplot(comparison_df %>%
                filter(asthma_cilited_FDR < 0.05 | season_ciliated_asthma_FDR < 0.05), 
              aes(x = -log10(season_ciliated_asthma_FDR), y=-log10(asthma_cilited_FDR))) +
  #geom_point(aes(color = ifelse((Gene %in% de.results.asthma[de.results.asthma$FDR < 0.05,]$Gene), 'red', 'black')))+
  geom_point(aes(color = significant), size = 1.2) +
  scale_color_manual(values=c("#163A7D", "#B4251A", "#AFABAB"))+
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = -log10(0.05), linetype="dashed") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlab(expression(-log[10]~"FDR DE genes in asthma corrected for season and ciliated")) +
  ylab(expression(-log[10]~"FDR DE genes in asthma corrected for ciliated")) +
  theme_classic() +
  #scale_x_break(breaks = c(2.8, 5.1)) +
  theme(legend.position="bottom", 
        legend.title=element_blank(), 
        legend.text = element_text(size = rel(1.5)), 
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)))

# good use of comments through the script. Made it very clear. 