library(here)
library(readr)
library(pheatmap)
library(dplyr)
library(tibble)

# reading in data
results_ashr <- readRDS(here("output", "results_ashr"))
results_wald <- readRDS(here("output", "results_wald"))
reslfc <- readRDS(here("output", "reslfc"))
deseq <- readRDS(here("output", "deseq"))
dds <- readRDS(here("output", "dds"))

counts <- counts(deseq, normalized = TRUE)
vst <- vst(deseq, blind = FALSE)


# reading in gene sets
# adhesion molecules
adhesion <- read_tsv(here("support", "GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES.v2024.1.Hs.tsv")) 
adhesion <- data.frame(strsplit(adhesion$GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES[17], ","))
adhesion <- data.frame(adhesion[!(adhesion$c..CDH2....CDH3....CDH4....CDH5....CDH6....CDH7....CDH8....CDH9... ==""),])
colnames(adhesion)[1] <- "genes"

# NOTCH pathway genes
notch <- read_tsv(here("support", "KEGG_NOTCH_SIGNALING_PATHWAY.v2024.1.Hs.tsv")) 
notch <- data.frame(strsplit(notch$KEGG_NOTCH_SIGNALING_PATHWAY[17], ","))
notch <- data.frame(notch[!(notch$c..DLL3....RBPJL....DTX2....CREBBP....CTBP1....CTBP2....DTX3L... ==""),])
colnames(notch)[1] <- "genes"



# subsetting adhesion molecules
adh_list <- list()
adh_assay <- list()
gene_set <- notch
for (i in names(results_ashr)){
  #i <- "V_VC_day19_VC"
  # gene_set molecules from deseq results
  adh_list[[i]] <- results_ashr[[i]] %>% as.data.frame() %>%
    filter(padj < 0.05) %>%
    filter(external_gene_name %in% gene_set$genes) %>%
    group_by(external_gene_name) %>%
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    arrange(desc(log2FoldChange))
  
  # vst assay for secreted proteins
  # b <- adh_assay[[i]] <- counts %>%
  #   merge(y = rowData(dds), by=0) %>%
  #   semi_join(y = adh_list[[i]], by=c("row")) %>%
  #   column_to_rownames(var = "external_gene_name") %>%
  #   dplyr::select(-c(Row.names, row))
  
  adh_assay[[i]] <- vst %>% assay %>%
    merge(y = rowData(dds), by=0) %>%
    semi_join(y = adh_list[[i]], by=c("row")) %>%
    column_to_rownames(var = "external_gene_name") %>%
    dplyr::select(-c(Row.names, row))
  
  adh_assay[[i]] <- adh_assay[[i]][match(adh_list[[i]]$external_gene_name, rownames(adh_assay[[i]])),]
  # annotation df
  df <- as.data.frame(colData(deseq)[,c("sort","rep.n")])
  
  plot_assay <- na.omit(adh_assay[[i]][1:70,])
  # plotting heatmaps
  if(nrow(plot_assay) > 0){
    print(pheatmap(plot_assay, cluster_rows=FALSE, show_rownames=TRUE,
    cluster_cols=TRUE, 
    #annotation_col=df, 
    main = i, fontsize_row = 10, fontsize_col = 10, angle_col = 90))
  }else{
    print("sorry there are no genes")
  }
  # print(EnhancedVolcano(sf_list[[i]],
  #                       lab = sf_list[[i]]$external_gene_name,
  #                       x = 'log2FoldChange',
  #                       y = 'padj',
  #                       pCutoff = 0.05,
  #                       FCcutoff = 1, 
  #                       title = paste0("Differentially expressed secreted factors for ", i),
  #                       drawConnectors = TRUE,
  #                       caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
  #                       max.overlaps = 20,
  #                       labSize = 4,
  #                       subtitle = ""))
}

EnhancedVolcano(sf_list[[i]],
                lab = sf_list[[i]]$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                main = paste0("Differentially expressed secreted factors for ", i)) #+
#ggrepel::geom_label_repel(mapping = ggplot2::aes(label = external_gene_name),max.overlaps = Inf)
