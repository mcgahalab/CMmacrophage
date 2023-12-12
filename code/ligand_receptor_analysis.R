library(icellnet)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(gridExtra)
library(jetset)
library(RColorBrewer)

colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))


icellnetdb <- '/cluster/projects/mcgahalab/ref/icellnet/ICELLNETdb.tsv'
db <- read.table(icellnetdb, sep="\t", check.names = F,header=T,
                 stringsAsFactors = F, na.strings="")
db.name.couple=name.lr.couple(db, type="Subfamily")

#### Functions ####
source("~/git/mini_projects/mini_functions/geneMap.R")
gm <- geneMap()

#### Load in Partner Data ####
#download PC.data.all and PC.target.all objects from the github and open them on your Rstudio session - adapt path if needed
partnerdir <- '/cluster/home/quever/git/ICELLNET/data'

PC.data.all=as.data.frame(read.csv(file.path(partnerdir, "PC.data.all.csv"), 
                                   sep=",", header = T, check.names=FALSE, 
                                   stringsAsFactors = FALSE, na.strings = ""))
rownames(PC.data.all)=PC.data.all$ID
PC.target.all=as.data.frame(read.csv(file.path(partnerdir, "PC.target.all.csv"), sep=",",
                                     header = T, check.names=FALSE, 
                                     stringsAsFactors = FALSE, na.strings = ""))

my.selection=c("Epith", "Fblast_B", "Endoth","Mono", "Macroph", "pDC", "DC2", "DC1", "NK", "Neutrop","CD4 T cell","CD8 T cell", "Treg","B cell")
PC.target = PC.target.all[
  which(PC.target.all$Class %in% my.selection | PC.target.all$Class%in%my.selection),
  c("ID","Class","Cell_type")
]
PC.data = PC.data.all[,PC.target$ID]

### Convert the gene symbol to affy ID 
PC.affy.probes = as.data.frame(PC.data[,c(1,2)])
PC.affy.probes$ID = rownames(PC.affy.probes) # for format purpose
transform = db.hgu133plus2(db,PC.affy.probes) # creation of a new db2 database with AffyID instead of gene symbol

##Gene scaling of the partner cell dataset
PC.data=gene.scaling(data = PC.data, n=18, db = transform) 


#### Load in sample data ####
dir2='/cluster/home/quever/git/ICELLNET/data_CAF'
# Central cell data file (processed gene expression matrix)
demo_data=as.data.frame(read.table(file.path(dir2, "data_CAF.txt"), sep="\t", header = T))
rownames(demo_data)=demo_data$SYMBOL
demo_data=dplyr::select(demo_data, -SYMBOL)
CC.data= gene.scaling(data = demo_data, n=4, db = db) 

#Target central cell file (description of the different samples)
CC.target = as.data.frame(read.table(file.path(dir2, "target_CAF.txt"),sep = "\t",header=T))
head(CC.target)

#Rescale the data
CC.selection.S1 = CC.target[which(CC.target$Type=="T"&CC.target$subset=="S1"&CC.target$Cancer.subtype=="TN"),"Sample.Name"] # CAF-S1 in TNBC samples
CC.selection.S2 = CC.target[which(CC.target$Type=="T"&CC.target$subset=="S4"&CC.target$Cancer.subtype=="TN"),"Sample.Name"] # CAF-S4 in TNBC samples

CC.data.selection.S1 = CC.data[,which(colnames(CC.data)%in%CC.selection.S1)]
CC.data.selection.S2 = CC.data[,which(colnames(CC.data)%in%CC.selection.S2)]

#--- Computer communication scores ----
score.computation.1= icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S1,  
                                    PC.target = PC.target, PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db)
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]]


score.computation.2= icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S2,  
                                    PC.target = PC.target,PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db)
score2=as.data.frame(score.computation.2[[1]])
lr2=score.computation.2[[2]]

Scores=cbind(score1,score2)
colnames(Scores)=c("CAF-S1","CAF-S4")
Scores

#### Load Jalal data ####
dir <- '/cluster/projects/mcgahalab/data/mcgahalab/jalal_receptor_ligand/results/counts/'
lranaldir <- file.path(dir, "..", "manual", "ligand_receptor")
file <- 'all_tpm.tsv'
min_expr <- 1
expr_frac_samples <- 0.2

# Load bulk RNAseq data
data=read.table(file.path(dir, file), header = T, check.names=FALSE, 
                stringsAsFactors = FALSE, na.strings = "")
genes <- gm$ENSEMBL$SYMBOL[data$gene]
genes[is.na(genes)] <- 'NA'
idx <- which(!duplicated(genes))
data <- as.data.frame(data[idx,]) %>%
  dplyr::select(-gene) %>% 
  magrittr::set_rownames(genes[idx]) %>%
  rename_with(., ~gsub(".genes.results$", "", .)) 
low_expr_idx <- which(rowSums(data < min_expr) > (ncol(data) * expr_frac_samples))
data <- data[-low_expr_idx,]

# Data scaling
# n <- apply(data, 2, quantile, probs=seq(0.95, 1, by=0.01))
data.scaled=gene.scaling(data = data, n=1, db = db)

# data selection
ids <- gsub("_(REP[0-9]).*$", "", colnames(data.scaled)) %>% table
rmids_regex <- "HD150|HD96"
ids <- ids[grep(rmids_regex, names(ids), ignore.case = T, invert=T)]
sample_ids <- names(ids[ids > 1])
sample_ids <- lapply(sample_ids, function(i) grep(i, colnames(data.scaled), value=T)) %>%
  setNames(., sample_ids)
ctrl_ids <- grep("Donor", colnames(data.scaled), value=T) %>%
  gsub("^.*?_", "", .) %>% unique 
ctrl_ids <- lapply(ctrl_ids, function(i) grep(paste0("Donor.*", i), colnames(data.scaled), value=T)) %>%
  setNames(., ctrl_ids)

#compute communication score
score_lr <- lapply(c(sample_ids, ctrl_ids), function(id_i){
  print(paste(id_i, collapse=","))
  score.computation <- icellnet.score(direction="out", 
                                      PC.data=PC.data, 
                                      CC.data= data.scaled[,id_i],  
                                      PC.target = PC.target, 
                                      PC=my.selection, 
                                      CC.type = "RNAseq", 
                                      PC.type = "Microarray",  
                                      db = db)
  score=as.data.frame(score.computation[[1]])
  lr=score.computation[[2]]
  return(list("score"=score, "lr"=lr))
})
saveRDS(score_lr, file=file.path(lranaldir, "scores_lr.filter.hdrm.rds"))


#### Visualization ####
score_lr <- readRDS(file=file.path(lranaldir, "scores_lr.filter.hdrm.rds"))

Scores <- lapply(score_lr, function(i)i$score) %>%
  do.call(cbind, .) %>%
  magrittr::set_colnames(names(score_lr))
targetcells <- colnames(score_lr[[1]]$lr)
Lrs <- lapply(targetcells, function(tcellid){
  lapply(score_lr, function(i) i$lr[,tcellid]) %>%
    do.call(cbind, .) %>%
    magrittr::set_colnames(names(score_lr))
}) %>% 
  setNames(., targetcells)


# intercellular communication plot
PC.col = c("Epith"="#C37B90", "Muscle_cell"="#c100b9","Fblast_B"="#88b04b", "Fblast"="#88b04b","Endoth"="#88b04b",
           "Mono"="#ff962c","Macroph"="#ff962c","moDC"="#ff962c","DC1"="#ff962c","DC2"="#ff962c","pDC"="#ff962c","NK"="#ff962c","Neutrop"="#ff962c",
           "CD4 T cell"="#5EA9C3","CD8 T cell"="#5EA9C3","Treg"="#5EA9C3","B cell"="#5EA9C3")
network_plots <- lapply(seq_along(colnames(Scores)), function(scores_i) {
  network.create(icn.score = Scores[,scores_i, drop=F] %>%
                   rename_with(., ~stringr::str_wrap(gsub("_", " ", .),10)), 
                 scale = c(round(min(Scores)-1),round(max(Scores))+1), 
                 direction = "out", 
                 PC.col)
}) %>%
  setNames(., colnames(Scores))

pdf(file.path(lranaldir, "intercellular_comm.pdf"), width = 20, height = 20)
cowplot::plot_grid(plotlist = network_plots, ncol = ceiling(sqrt(ncol(Scores))))
dev.off()
file.copy(from = file.path(lranaldir, "intercellular_comm.pdf"), to="~/xfer", overwrite = T)


# individual communication scores
family.col = c( "CC"= "#AECBE3", "COL"= "#66ABDF", "CX3C"= "#1D1D18"  ,
            "EGFR family"="#156399", "Semaphorin" ="#676766", "TNF" = "#12A039",  "NA"="#908F90")
subfamilies <- unique(db.name.couple[,2])
family.col = setNames(sample(col_vec, length(subfamilies)), subfamilies)

balloon_plots <- lapply(names(Lrs), function(lr_id) {
  lr_ind <- Lrs[[lr_id]] %>%
    as.data.frame %>%
    rename_with(., ~stringr::str_wrap(gsub("_", " ", .),10))
  icellnet::LR.balloon.plot(lr = lr_ind, sort.by="var", 
                            thresh = 25, topn=15, 
                            db.name.couple=as.data.frame(db.name.couple), 
                            title=lr_id, 
                            family.col=family.col)
})
pdf(file.path(lranaldir, "balloon_plots.pdf"), width = 15)
balloon_plots
dev.off()
file.copy(from = file.path(lranaldir, "balloon_plots.pdf"), to="~/xfer", overwrite = T)









