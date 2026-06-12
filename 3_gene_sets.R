library(here)
library(readr)
library(pheatmap)
library(dplyr)
library(tibble)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

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




######## heatmps for grant 



wanted_samples <- c(
  "V_VC_CD45_1", "V_VC_CD45_2",
  #"d19_1_VC_CD45", "d19_2_VC_CD45",
  "d12_1_VC_CD45", "d12_2_VC_CD45"
  #"V_VC_1", "V_VC_2",
  #"d19_1_VC", "d19_2_VC",
  #"d12_1_VC", "d12_2_VC"
)

sample_idx <- match(wanted_samples, colData(dds)$title)

mat <- assay(vsd)[genes, sample_idx]

# replace ENSG IDs with gene symbols
rownames(mat) <- c("KLF2", "EBF1", "EGR3")

# z-score by gene
mat_scaled <- t(scale(t(mat)))

# sample annotation
ann_col <- data.frame(
  sort = colData(dds)$sort[sample_idx],
  rep.n = colData(dds)$rep.n[sample_idx]
)

rownames(ann_col) <- colnames(mat_scaled)

# heatmap
library(pheatmap)

pheatmap(
  mat_scaled,
  annotation_col = ann_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 9,
  border_color = NA,
  main = "KLF2, EBF1 and EGR3", angle_col = 45
)

## compleex heatmap

vst_count <- list()
gene_list <- list() 
genes_to_label <- c("EBF1", "KLF2", "EGR3", "SOX17", "SNAI1", "YAP1", "WWTR1")
a <- gene_list[[i]]
col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c(
    "#313695",
    "#74ADD1",
    "#FFFFBF",
    "#FDAE61",
    "#A50026"
  )
)


for (i in names(results_ashr)){
  i <- "V_VC_CD45_d12_VC_CD45"
  gene_list[[i]] <- results_ashr[[i]] %>% as.data.frame() %>%
    filter(padj < 0.999) %>%
    group_by(external_gene_name) %>%
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    arrange(desc(log2FoldChange))
  
  
  vst_count[[i]] <- vst %>% assay %>%
    merge(y = rowData(dds), by=0) %>%
    semi_join(y = gene_list[[i]], by=c("row")) %>%
    column_to_rownames(var = "external_gene_name") %>%
    dplyr::select(-c(Row.names, row))
  
  labels <- ifelse(
    rownames(vst_count[[i]]) %in% genes_to_label,
    rownames(vst_count[[i]]),
    ""
  )

  wanted_samples <- c(
    "V_VC_CD45_1", "V_VC_CD45_2",
    #"d19_1_VC_CD45", "d19_2_VC_CD45",
    "d12_1_VC_CD45", "d12_2_VC_CD45"
    #"V_VC_1", "V_VC_2",
    #"d19_1_VC", "d19_2_VC",
    #"d12_1_VC", "d12_2_VC"
  )
  
  sample_idx <- match(wanted_samples, colData(dds)$title)
  
  vst_count[[i]] <- vst_count[[i]][, sample_idx]
  
  vst_count[[i]] <- vst_count[[i]][match(gene_list[[i]]$external_gene_name, rownames(vst_count[[i]])),]
  colnames(vst_count[[i]]) <- c("Embryo 1 (DP)", "Embryo 2 (DP)", "Day12 ESc 1 (DP)", "Day12 ESc 1 (DP)")
  
  # annotation df
  df <- as.data.frame(colData(deseq)[,c("sort","rep.n")])
  
  plot_assay <- na.omit(vst_count[[i]])
  # z-score by gene
  plot_assay <- t(scale(t(plot_assay)))
  
  idx <- which(rownames(plot_assay) %in% genes_to_label)
  
  # plotting heatmaps
  if(nrow(plot_assay) > 0){
    # print(pheatmap(plot_assay, cluster_rows=TRUE, show_rownames=TRUE,
    #                cluster_cols=TRUE, labels_row = labels,
    #                #annotation_col=df, 
    #                main = i, fontsize_row = 10, fontsize_col = 10, angle_col = 45))
    
    ComplexHeatmap::Heatmap(
      plot_assay,
      #col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      column_names_rot = 45,
      heatmap_legend_param = list(
        title = "Z-score"
      ),
      right_annotation = rowAnnotation(
        mark = anno_mark(
          at = idx,
          labels = rownames(plot_assay)[idx]
        )
      )
    )
  }else{
    print("sorry there are no genes")
  }
  
}
