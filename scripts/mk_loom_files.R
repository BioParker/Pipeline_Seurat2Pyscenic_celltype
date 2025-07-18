library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(SeuratDisk)
library(optparse)

option_list = list(
  make_option(c("-f", "--seurat"), type="character", default=NULL, 
              help="Path to Seurat Object", metavar="character"),
    make_option(c("-c", "--celltype"), type="character", default="NULL", 
              help="Cluster Identity", metavar="character"),
    make_option(c("-m", "--cellcol"), type="character", default="NULL",    
              help="Name of cell type column in seurat metadata", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sop <- opt$seurat

ct <- opt$celltype

cname <- opt$cellcol

so <- readRDS(sop)

colnames(so@meta.data)[which(colnames(so@meta.data) == cname)] <- "anno"

ido <- subset(so, subset = anno == ct)

SaveLoom(ido, paste0(ct,".loom"))

