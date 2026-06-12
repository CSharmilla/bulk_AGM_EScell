library(here)
library(DESeq2)
library(edgeR)
library(ashr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)


# readin in data
merged <- readRDS(here("output", "merged"))
meta <- readRDS(here("output", "meta"))
genes_list <- readRDS(here("support", "genes_list"))

# adding another variable to meta
meta <- mutate(meta, domain_day = gsub("_.*", "", sort))
meta <- mutate(meta, sort_pop = sub("^[^_]*_", "", sort)) # extract everything after first underscore
meta <- mutate(meta, rep.n = sub("_.*", "", rep))

## pre-processing
# coldata rownames and txi column names matching
meta <- meta[match(colnames(merged$counts), rownames(meta)), ] 

# check if coldata rownames and txi column names match and are in order
print(all(rownames(meta) == colnames(merged$counts)))
print(all(rownames(meta) %in% colnames(merged$counts)))

# getting gene list - make sure it is in the same order as txi rows
rowdata <- genes_list[which(genes_list$ensembl_gene_id %in% rownames(merged$counts)),]
rowdata <- rowdata[match(rownames(merged$counts), rowdata$ensembl_gene_id),] 
rownames(rowdata) <- 1:nrow(rowdata) 

# check if rowdata and txi rownames  match and are in order
print(all(rowdata$ensembl_gene_id == rownames(merged$counts)))
print(all(rowdata$ensembl_gene_id %in% rownames(merged$counts)))

# factor levels
meta$sort <- factor(meta$sort)
meta$rep.n <- factor(meta$rep.n)
meta$domain_day <- factor(meta$domain_day)
meta$sort_pop <- factor(meta$sort_pop)


# creating DESeqDataSet objects

if (!(file.exists(here("output", "dds")))) {
  
  # converting to deseq2 object
  dds <- dds_unfil <- DESeqDataSetFromTximport(merged, colData = meta, 
                                               design = ~rep.n + sort, #~domain_day + domain_day:rep.n + domain_day:sort_pop, 
                                               rowData = rowdata)
  
  # changing this to make merge easier
  colnames(rowData(dds))[1] <- colnames(rowData(dds_unfil))[1] <- "row"
  
  # basic pre-filtering
  #smallestGroupSize <- 2
  #keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  keep <- filterByExpr(dds, group = dds$sort)
  dds <- dds[keep,]
  
  # estimate size factors - another assay added with normalisation factors
  dds <- estimateSizeFactors(dds)
  dds_unfil <- estimateSizeFactors(dds_unfil)
  
  # saving these
  saveRDS(dds, file = here("output", "dds"))
  saveRDS(dds_unfil, file = here("output", "dds_unfil"))
  
  } else{
  dds <- readRDS(here("output", "dds"))
  dds_unfil <- readRDS(here("output", "dds_unfil"))
}


# transforming and plotting to see the data

# blind=TRUE should be used for comparing samples in a manner unbiased by prior 
# information on samples, for example to perform sample QA (quality assurance)
vst_blind <- vst(dds, blind=TRUE)

# plotting pca
plotPCA(vst_blind, intgroup=c("sort")) #+ geom_text(aes(label = vst_blind$title))


# running deseq
if (!(file.exists(here("output", "deseq")))) {
  deseq <- DESeq(dds)
  saveRDS(deseq, file = here("output", "deseq"))
} else{
  deseq <- readRDS(here("output", "deseq"))
}

resultsNames(deseq)

# the name of a factor in the design formula, 
# the name of the numerator level for the fold change, and the name of the denominator level for the fold change
#  log2 fold change of 1.5 for a specific gene in the “WT vs KO comparison” means that the 
# expression of that gene is increased in WT relative to KO by a multiplicative factor of 2^1.5 ≈ 2.82

# getting results and transforming
# creating contrasts
contrasts <- list(V_VC_day19_VC = c("sort", "V_VC", "d19_VC"),
                  V_CD45_d19_CD45 = c("sort", "V_CD45", "d19_CD45"),
                  V_VC_CD45_d19_VC_CD45 = c("sort", "V_VC_CD45", "d19_VC_CD45"),
                  
                  V_VC_day12_VC = c("sort", "V_VC", "d12_VC"),
                  V_CD45_d12_CD45 = c("sort", "V_CD45", "d12_CD45"),
                  V_VC_CD45_d12_VC_CD45 = c("sort", "V_VC_CD45", "d12_VC_CD45"), 
                  
                  V_VC_CD45_V_VC = c("sort", "V_VC_CD45", "V_VC"))

results_ashr <- list() # dataframes with shrunken lfc
results_wald <- list() # wald test dataframes
reslfc <- list() # list of results for plots

if (!file.exists(here("output", "results_ashr"))){
  for (i in names(contrasts)){
    
    name1 <- paste0(i, "_ashr_results")
    name2 <- paste0(i, "_results")
    
    res <- results(deseq, contrast = contrasts[[i]], alpha = 0.05) # cant use tidy = TRUE here since using lfc shrink
    
    # LFC shrinkage
    reslfc[[i]] <- lfcShrink(deseq, res = res, type = "ashr")
    ashr <- data.frame(reslfc[[i]])
    ashr$row <- rownames(ashr)
    ashr <- as.data.frame(merge(ashr, rowData(dds), all.x=TRUE))
    
    assign(name1, ashr)
    results_ashr[[i]] <- ashr
    
    results_wald[[i]] <- results(deseq, contrast = contrasts[[i]], alpha = 0.05, tidy = TRUE)
    results_wald[[i]] <- merge(results_wald[[i]], rowData(dds), all.x=TRUE)
    
    assign(name2, data.frame(results_wald[[i]]))
    
    
  }
  saveRDS(results_ashr, file = here("output", "results_ashr"))
  saveRDS(results_wald, file = here("output", "results_wald"))
  saveRDS(reslfc, file = here("output", "reslfc"))
  
}else{
  results_ashr <- readRDS(here("output", "results_ashr"))
  results_wald <- readRDS(here("output", "results_wald"))
  reslfc <- readRDS(here("output", "reslfc"))
}

# separate loop for this
for (i in names(contrasts)){
  
  name1 <- paste0(i, "_ashr_results")
  name2 <- paste0(i, "_results")
  
  assign(name1, results_ashr[[i]])
  assign(name2, data.frame(results_wald[[i]]))
  
}




# sanity check
plotCounts(dds, gene="ENSG00000179388", intgroup="sort", normalized = TRUE)

counts <- counts(deseq, normalized = TRUE)

# KLF2, EBF1, EGR3


plotCounts(
  dds,
  gene = "ENSG00000127528",
  intgroup = c("sort", "rep.n"),
  returnData = TRUE
)

results(
  deseq,
  contrast = c("sort", "V_VC_CD45", "d19_VC_CD45")
)["ENSG00000127528", ]

results_ashr[["V_VC_CD45_d19_VC_CD45"]] %>%
  filter(row == "ENSG00000127528") %>%
  select(log2FoldChange, lfcSE)


# KLF2
# Expression is higher in V_VC_CD45 than d19_VC_CD45 (MLE log2FC = 2.83, shrunken log2FC = 1.34), with nominal evidence of differential expression (p = 0.012), but it does not pass the chosen multiple-testing threshold (FDR = 0.105).

sum(results(
  deseq,
  contrast = c("sort", "V_VC_CD45", "d19_VC_CD45")
)$padj < 0.05, na.rm = TRUE)



# KLF2, EBF1, EGR3
genes <- c(
  "ENSG00000127528",  # KLF2
  "ENSG00000164330",  # EBF1
  "ENSG00000179388"   # EGR3
)

# variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# # choose samples to display
# keep_samples <- colData(dds)$title %in%
#   c(
#     "V_VC_CD45_1", "V_VC_CD45_2",
#     "d19_1_VC_CD45", "d19_2_VC_CD45",
#     "d12_1_VC_CD45", "d12_2_VC_CD45",
#     "V_VC_1", "V_VC_2",
#     "d19_1_VC", "d19_2_VC",
#     "d12_1_VC", "d12_2_VC"
#   )
# 
# # extract expression matrix
# mat <- assay(vsd)[genes, keep_samples]


wanted_samples <- c(
  "V_VC_CD45_1", "V_VC_CD45_2",
  "d19_1_VC_CD45", "d19_2_VC_CD45",
  "d12_1_VC_CD45", "d12_2_VC_CD45",
  "V_VC_1", "V_VC_2",
  "d19_1_VC", "d19_2_VC",
  "d12_1_VC", "d12_2_VC"
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

# genes to highlight
genes_of_interest <- 
#   c(
#   "ENSG00000127528",  # KLF2
#   "ENSG00000164330",  # EBF1
#   "ENSG00000179388"   # EGR3
# )
  
  c(
  "KLF2",
  "EBF1",
  "EGR3"
)

# contrasts to plot
contrasts_to_plot <- c(
  "V_VC_CD45_d19_VC_CD45",
  "V_VC_CD45_d12_VC_CD45",
  "V_VC_CD45_V_VC"
)


volcano_plots <- list()

# looping through contrasts
for (i in contrasts_to_plot) {
  
  res <- as.data.frame(results_ashr[[i]])
  
  # clean gene labels 
  res$gene_label <- res$external_gene_name
  
  # fallback safety
  res$gene_label[is.na(res$gene_label) | res$gene_label == ""] <- res$row
  
  # fix NA padj
  res$padj[is.na(res$padj)] <- 1
  
  # confirm genes exist
  print(intersect(genes_of_interest, res$gene_label))
  
  p <- EnhancedVolcano(
    res,
    subtitle = NULL,
    lab = res$gene_label,
    
    x = "log2FoldChange",
    y = "padj",
    
    selectLab = genes_of_interest,
    
    pCutoff = 0.05,
    FCcutoff = 1,
    
    drawConnectors = TRUE,
    widthConnectors = 0.8,
    
    arrowheads = TRUE,
    boxedLabels = TRUE,
    labSize = 7,
    maxoverlapsConnectors = Inf,
    title = i
  ) 
  
  volcano_plots[[i]] <- p
}


ggsave(
  "volcano_plot_1.pdf",
  volcano_plots[[1]],
  width = 12,
  height = 8
)


# excel sheets
library(openxlsx)

for (nm in names(results_wald)) {
  write.xlsx(
    results_wald[[nm]],
    file = here("output", paste0(nm, ".xlsx")),
    overwrite = TRUE
  )
}

a <- data.frame(results_wald[[1]])



#######

genes_of_interest <- c("KLF2", "EBF1", "EGR3")

contrasts_to_plot <- c(
  "V_VC_CD45_d19_VC_CD45",
  "V_VC_CD45_d12_VC_CD45",
  "V_VC_CD45_V_VC"
)

volcano_plots <- list()

for (i in contrasts_to_plot) {
  
  res <- as.data.frame(results_ashr[[i]])
  
  #ensuring ENSG is row
  res$row <- rownames(res)
  
  #force character
  res$external_gene_name <- as.character(res$external_gene_name)
  
  #fallback mapping 
  res$gene_label <- ifelse(
    is.na(res$external_gene_name) | res$external_gene_name == "",
    res$row,
    res$external_gene_name
  )
  
  #clean whitespace 
  res$gene_label <- trimws(res$gene_label)
  
  res$padj[is.na(res$padj)] <- 1
  
  
  #DEBUG 
 
  cat("\nContrast:", i, "\n")
  print(genes_of_interest %in% res$gene_label)
  print(intersect(genes_of_interest, res$gene_label))
  
  #label genes
  label_genes <- unique(c(
    genes_of_interest,
    res$gene_label[order(res$padj)][1:10]
  ))
  
  p <- EnhancedVolcano(
    res,
    
    lab = res$gene_label,
    
    x = "log2FoldChange",
    y = "padj",
    
    selectLab = label_genes,
    
    pCutoff = 0.05,
    FCcutoff = 1,
    
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.6,
    maxoverlapsConnectors = 5,
    
    labSize = 3,
    title = i
  )
  
  volcano_plots[[i]] <- p
}

volcano_plots[[1]]
