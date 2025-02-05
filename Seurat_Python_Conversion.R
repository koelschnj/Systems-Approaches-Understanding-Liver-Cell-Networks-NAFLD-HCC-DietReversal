library(Seurat)
#!remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
#!remotes::install_github("PMBio/MuDataSeurat")
library(MuDataSeurat)
#!install.packages("anndata")
library(anndata)
#!remotes::install_github("ilia-kats/MuData")
library(MuData)
library(MultiAssayExperiment)
library(rhdf5)

sc.all.cells <- readRDS(file = "~/RD_sex_differences/Saved_Rfiles/RD.sex.diff.6group.5clusters.Sv4.all.cells.prelim.names.RDS")
table(sc.all.cells@meta.data$Prelim.annotation)

Idents(sc.all.cells) <- sc.all.cells@meta.data$new.group
table(sc.all.cells@meta.data$new.group)

##try to create AnnData object with SeuratDisk and save on cluster
SaveH5Seurat(sc.all.cells, filename = "RDsexdiff6groups.h5Seurat")
Convert("RDsexdiff6groups.h5Seurat", dest = "h5ad")
##these save directly to koelschnj directory so need to know that to load into python
##should reorganize these files once we are using python all into one directory

##try to create MuData object with MuDataSeurat and save on cluster as well
counts <- sc.all.cells[["RNA"]]@counts
rownames(counts) <- paste(rownames(counts), "RNA", sep = "-")
data <- sc.all.cells[["RNA"]]@data
rownames(data) <- rownames(counts)

all.cells <- CreateAssayObject(counts = counts)
all.cells@data <- data

all.cells.new <- CreateSeuratObject(sc.all.cells[["RNA"]])
all.cells.new[["RNA"]] <- all.cells
DefaultAssay(all.cells.new) <- "RNA"
all.cells.new

WriteH5MU(all.cells.new, "RDsexdiff6groups.h5mu")

##go back to check if we can load in h5 files and if they have the necessary meta.data
#!!H5.processed.data <- LoadH5Seurat("~/RDsexdiff6groups.h5Seurat")
#!h5ad <- ReadH5AD("RDsexdiff6groups.h5ad")
h5ad <- read_h5ad("RDsexdiff6groups.h5ad")
str(h5ad)
head(h5ad$raw, 5)
colnames(h5ad$raw)
head(h5ad$uns, 5)
head(h5ad$shape, 5)
head(h5ad$varp, 5)
head(h5ad$varm, 5)
head(h5ad$var_names, 5)
head(h5ad$var, 5)
head(h5ad$n_vars, 5)
head(h5ad$obsp, 5)
head(h5ad$obsm, 5)
head(h5ad$obs_names, 5)
head(h5ad$obs, 5)
head(h5ad$n_obs, 5)
head(h5ad$isbacked, 5)
head(h5ad$is_view, 5)
head(h5ad$T, 5)
head(h5ad$layers, 5)
head(h5ad$filename, 5)
head(h5ad$X, 5)

#!h5 <- H5.processed.data$new("RDsexdiff6groups.h5mu", mode = "r")
#!h5mu <- ReadH5MU("RDsexdiff6groups.h5mu")
h5mu <- readH5MU("RDsexdiff6groups.h5mu")
head(h5mu@colData, 5)
head(h5mu@sampleMap, 5)
head(h5mu@ExperimentList, 5)

###############################################################################################
##load in data from 2024 HCCp1 and extract MCD group only for test on LIANA after making files
HCCp1.2024.cells <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.new.abbreviated.group.names.rds")
table(HCCp1.2024.cells@meta.data$Abbreviated)

Idents(HCCp1.2024.cells) <- HCCp1.2024.cells@meta.data$new.group
table(HCCp1.2024.cells@meta.data$new.group)

CD <- subset(HCCp1.2024.cells, idents = "CD", invert = FALSE)
table(CD@meta.data$new.group)
table(CD@meta.data$Abbreviated)

WDt <- subset(HCCp1.2024.cells, idents = "WD.t", invert = FALSE)
table(WDt@meta.data$new.group)
table(WDt@meta.data$Abbreviated)

RDt <- subset(HCCp1.2024.cells, idents = "RD.t", invert = FALSE)
table(RDt@meta.data$new.group)
table(RDt@meta.data$Abbreviated)

RDn <- subset(HCCp1.2024.cells, idents = "RD.n", invert = FALSE)
table(RDn@meta.data$new.group)
table(RDn@meta.data$Abbreviated)

##try to create AnnData object with SeuratDisk and save on cluster
SaveH5Seurat(CD, filename = "MCD2024.h5Seurat")
Convert("MCD2024.h5Seurat", dest = "h5ad")

SaveH5Seurat(WDt, filename = "WDt2024.h5Seurat")
Convert("WDt2024.h5Seurat", dest = "h5ad")

SaveH5Seurat(RDt, filename = "RDt2024.h5Seurat")
Convert("RDt2024.h5Seurat", dest = "h5ad")

SaveH5Seurat(RDn, filename = "RDn2024.h5Seurat")
Convert("RDn2024.h5Seurat", dest = "h5ad")

##these save directly to koelschnj directory so need to know that to load into python
##should reorganize these files once we are using python all into one directory

##try to create MuData object with MuDataSeurat and save on cluster as well
counts <- CD[["RNA"]]@counts
rownames(counts) <- paste(rownames(counts), "RNA", sep = "-")
data <- CD[["RNA"]]@data
rownames(data) <- rownames(counts)

CD.cells <- CreateAssayObject(counts = counts)
CD.cells@data <- data

CD.cells.new <- CreateSeuratObject(CD[["RNA"]])
CD.cells.new[["RNA"]] <- CD.cells
DefaultAssay(CD.cells.new) <- "RNA"
CD.cells.new

WriteH5MU(CD.cells.new, "MCD2024.h5mu")

##can make MU data files for other groups later...should only need h5ad for basic liana analysis
