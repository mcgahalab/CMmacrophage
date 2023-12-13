renv::load("/cluster/home/quever/downloads/renvs/")

library(tidyverse)
library(ComplexHeatmap)
library(ggrastr)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(umap)
library(cowplot)
library(ggrepel)
library(GSVA)
library(org.Hs.eg.db)
library(icellnet)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(gridExtra)
library(jetset)
library(RColorBrewer)
library(clusterProfiler)


colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
species <- 'Homo sapiens'
pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/jalal_receptor_ligand/results/'

gprofiler_dir <- '/cluster/projects/mcgahalab/ref/gprofiler'
gprofiler_f <- file.path(gprofiler_dir, 'gprofiler_full_hsapiens.ENSG.gmt')

barcodes_f="~/git/mini_projects/ref/barcodes.tsv"

#### Functions ####
# code cloned from https://github.com/quevedor2/mini_projects
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/makeLoupe.R")
gm <- geneMap()

###################
#### Functions ####
ssGseaFun <- function(msig_ds, lfc_v, ss_method='ssgsea'){
  require(GSVA)
  ssgsea <- tryCatch({
    sig_ens_gs <- split(setNames(msig_ds$entrez_gene, msig_ds$entrez_gene), 
                        f=msig_ds$gs_name)
    gsva(lfc_v, sig_ens_gs, verbose=FALSE, method=ss_method)
  }, error=function(e){NULL})
  return(ssgsea)
}

getDEGedgeRwt <- function(cts, meta, group){
  idx <- match(rownames(meta), colnames(cts))
  if(any(is.na(idx))){
    na_idx <- which(is.na(idx))
    # print(paste0("NA idx: ", paste(na_idx, collapse=",")))
    idx <- idx[-na_idx]
    meta <- meta[-na_idx,]
  }
  
  ## Get EdgeR Results
  se <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(cts[,idx])), 
                             colData = meta)
  
  
  edgey <- tryCatch({
    edgey <- DGEList(counts=cts[,idx],
                     samples=rownames(meta),
                     group=meta[,group])
    # design <- with(meta, model.matrix(as.formula(formula)))
    design <- model.matrix(~ meta[,group])
    keep <- filterByExpr(edgey, design)
    edgey <- edgey[keep, , keep.lib.sizes=FALSE]
    edgey <- calcNormFactors(edgey, method = "TMM")
    edgey <- estimateDisp(edgey)
    edgey
  }, error=function(e){NULL})
  
  # Calculate CPM/TMM values
  edge_tmm <- cpm(edgey)
  edge_tmm_spl <- edge_tmm %>% t %>% as.data.frame %>%
    split(., meta[,group]) 
  
  # Differential testing
  er_dat <- exactTest(edgey)
  et_res <- er_dat$table  %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("ensemble") %>% 
    mutate(padj=p.adjust(PValue, method='BH')) %>%
    mutate(sig=ifelse(padj < 0.05, "er_sig", "er_ns"))
  
  glvl <- levels(factor(meta[,group]))
  wilcox_res <- sapply(seq_along(edge_tmm_spl[[1]]), function(idx){  
    wt <- wilcox.test(edge_tmm_spl[[glvl[1]]][,idx], edge_tmm_spl[[glvl[2]]][,idx])
    fc <- mean(edge_tmm_spl[[glvl[1]]][,idx],na.rm=T) /  mean(edge_tmm_spl[[glvl[2]]][,idx], na.rm=T)
    c('W'=wt$statistic, 'pval'=wt$p.value, "FC"=fc, 'Log2FC'=log2(fc+1),
      'ensemble'=colnames(edge_tmm_spl[[1]][idx]))
  }) %>% 
    t %>% as.data.frame %>%
    mutate(padj=p.adjust(pval, method='BH'),
           FC=as.numeric(FC),
           Log2FC=as.numeric(Log2FC)) %>%
    mutate(sig=ifelse(padj < 0.05, "wt_sig", "wt_ns"))
  
  # pdf("~/xfer/dds2.pdf") 
  # # ggplot(dds_et, aes(x=log2FoldChange, y=logFC)) +
  # ggplot(dds_et, aes(x=log2FoldChange, y=logFC, color=sig, group=sig)) +
  #   geom_point()
  # ggplot(wt_et, aes(x=Log2FC, y=logFC, color=sig, group=sig)) +
  #   geom_point()
  # dev.off()
  return(list("edger"=et_res,  "wilcox"=wilcox_res))
}

rmVstBatchEffects <- function(vsd, dds, samples_meta, 
                              condition='condition', batchcol='batch'){
  mat <- assay(vsd)
  mm <- model.matrix(as.formula(paste0("~", condition)), colData(dds))
  mat <- limma::removeBatchEffect(mat, design=mm, 
                                  batch=samples_meta[colnames(dds),batchcol])
  assay(vsd) <- mat
  return(vsd)
}


##############
#### Main ####
dir.create(file.path(pdir, "manual", "objs"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "pca"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "loupe"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "deg"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "gsea"), showWarnings = F, recursive = T)

deganaldir <- file.path(pdir, "manual", "differential_expression")
file <- 'all.tsv'
min_expr <- 3
expr_frac_samples <- 0.2

# Load bulk RNAseq data
data=read.table(file.path(pdir, "counts", file), header = T, check.names=FALSE, 
                stringsAsFactors = FALSE, na.strings = "")
genes <- gm$ENSEMBL$SYMBOL[data$gene]
genes[is.na(genes)] <- 'NA'
idx <- which(!duplicated(genes))
data <- as.data.frame(data[idx,]) %>%
  dplyr::select(-gene) %>% 
  magrittr::set_rownames(genes[idx]) %>%
  rename_with(., ~gsub("genesresults$", "", .)) 

# Remove genes where [expr_frac_sample] of the dataset has fewer than [min_expr] 
# reads linked to that gene
# e.g. Remove: LOC101929116 is expressed at higher than 3 counts in only 2/42 samples
low_expr_idx <- which(rowSums(data < min_expr) > (ncol(data) * expr_frac_samples))
data <- data[-low_expr_idx,]

#### 1. Create metadata ####
# Create metadata data frame
metadata <- colnames(data)
metad <- data.frame("sample"=colnames(data)) %>%
  mutate(HDsample=gsub(".*(HD[0-9]*|Donor[0-9]*).*", "\\1", sample),
         Cellline=gsub(".*(H[0-9]+|SHP[0-9]+|Donor).*$", "\\1", sample),
         Celltype=gsub(".*_(M[0-9IL]*|TAMlike)(acro)?.*", "\\1", sample),
         CM=grepl("_[HDSP0-9]*CM_", sample),
         Rep=gsub(".*(REP|Donor)([0-9]*).*", "\\2", sample),
         Dataset=ifelse(grepl("Donor", sample), "External", "Internal"))
saveRDS(metad, file.path(pdir, "manual", "objs", "metad.rds"))

#### 2. PCA reduction ####
max_pc <- 10
dds_all <- DESeqDataSetFromMatrix(countData=data,
                                   colData=metad,
                                   design=as.formula("~Cellline"))
saveRDS(dds_all, file=file.path(pdir, "manual", "objs", "deseq_all.rds"))

# Run the main PCA analysis on vst-counts
# obtain normalized counts
counts <- vst(dds_all, blind=T)
# counts <- rmVstBatchEffects(counts, dds_all, metad, condition='Cellline', batch='Dataset')

tcnts <- as.data.frame(t(assay(counts)))
pca <- prcomp(tcnts, scale=F)
pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc))],
                             "condition"=as.character(counts$condition)))
if(max_pc > ncol(pca_x)) max_pc <- ncol(pca_x)
for(id in paste0("PC", c(1:max_pc))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}
saveRDS(pca_x, file=file.path(pdir, "manual", "pca", "pca_dat.all.rds"))

#### 3. Loupe visualization ####
dds_all <- readRDS(file=file.path(pdir, "manual", "objs", "deseq_all.rds"))
pca_x <- readRDS(file=file.path(pdir, "manual", "pca", "pca_dat.all.rds"))
metad <- readRDS(file=file.path(pdir, "manual", "objs", "metad.rds"))
metad2 <- metad %>% 
  tibble::column_to_rownames("sample")
makeLoupe(mat=assay(dds_all), meta=metad2, projections=pca_x,
          output_dir=file.path(pdir, "manual", "loupe"),
          output_name="pca_all",
          barcodes_f=barcodes_f)
file.copy(file.path(pdir, "manual", "loupe", "pca_all.cloupe"), to="~/xfer", overwrite = T)


#### 4. Differential analysis ####
# Comparing H1048-CM samples to H841-CM samples in HD143 and HD151
idx <- intersect(grep("HD143|HD151", metad$HDsample),
                 grep("H1048|H841", metad$Cellline))
dds <- DESeqDataSetFromMatrix(countData=data[,idx],
                              colData=metad[idx,],
                              design=as.formula("~Cellline"))
# Setting H841-CM treated as the baseline (i.e. a +ve LFC will indicate that H1048 is larger than H841 samples)
dds$Cellline <- relevel(dds$Cellline, ref = "H841")
levels(dds$Cellline) # ensure proper levelling

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- lfcShrink(dds, coef="Cellline_H1048_vs_H841", type="apeglm") %>% # Cellline_H1048_vs_Donor = resultsNames(dds)[2]
  as.data.frame %>%
  arrange(padj)
saveRDS(res, file=file.path(pdir, "manual", "deg", "Cellline_H1048_vs_H841.rds"))

#### 5. GSEA ####
res <- readRDS(file=file.path(pdir, "manual", "deg", "Cellline_H1048_vs_H841.rds"))
res <- res %>% 
  tibble::rownames_to_column("symbol") %>%
  mutate(entrez = gm$SYMBOL$ENSEMBL[symbol])

# Read in geneset database
gmt <- GSA::GSA.read.gmt(gprofiler_f)
gprof_ds <-setNames(gmt$genesets, 
                    paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]


msig_ds <- lapply(names(gprof_ds), function(sublvl){
  data.frame("gs_name"=sublvl,
             "entrez_gene"=gprof_ds[[sublvl]])
}) %>% 
  do.call(rbind,.) %>% 
  mutate(classification =  gsub(":.*", "", gs_name))

# run gsea on the geneset
lfc_v <- setNames(res$log2FoldChange,
                  res$entrez)
gsea <- clusterProfiler::GSEA(sort(na.omit(lfc_v), decreasing = T), 
             TERM2GENE = msig_ds, pvalueCutoff = 1)
gsea_short <- as.data.frame(gsea) %>% 
  dplyr::select(-c(ID, Description, core_enrichment)) %>%
  arrange(p.adjust)
saveRDS(gsea, file=file.path(pdir, "manual", "gsea", "Cellline_H1048_vs_H841.rds"))

