library("ggplot2")
library("Seurat")
library("scSorter")

HCCp1.sorted.individually <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.individually.w.scSorter5.RDS")
table(HCCp1.sorted.individually@meta.data$experimental.group)
table(HCCp1.sorted.individually@meta.data$group)
table(HCCp1.sorted.individually@meta.data$replicate)
table(HCCp1.sorted.individually@meta.data$Predicted_Type)

Idents(HCCp1.sorted.individually) <- HCCp1.sorted.individually@meta.data$group

##subset cells from each group to sort each individually on Macrophages
CD.cells <- subset(x = HCCp1.sorted.individually, idents = "CD", invert = FALSE)
table(CD.cells@meta.data$group)

PreT.cells <- subset(x = HCCp1.sorted.individually, idents = "PreT", invert = FALSE)
table(PreT.cells@meta.data$group)

MWD.T.cells <- subset(x = HCCp1.sorted.individually, idents = "MWD.T", invert = FALSE)
table(MWD.T.cells@meta.data$group)

MRD.T.cells <- subset(x = HCCp1.sorted.individually, idents = "MRD.T", invert = FALSE)
table(MRD.T.cells@meta.data$group)

MRD.NT.cells <- subset(x = HCCp1.sorted.individually, idents = "MRD.NT", invert = FALSE)
table(MRD.NT.cells@meta.data$group)

##first we need to sort T cell subsets in each group with GL14 
##sort CD T cell with GL14
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
table(CD.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList14.csv", sep = ",")

Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.Tcells <- subset (CD.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.Tcells <- NormalizeData(CD.Tcells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.Tcells <- FindVariableFeatures(CD.Tcells, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.Tcells <- ScaleData(CD.Tcells)
CD.Tcells <- RunPCA(CD.Tcells, features = VariableFeatures(CD.Tcells))

topgenes <- head(VariableFeatures(CD.Tcells), 3000)
expr = GetAssayData(CD.Tcells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
CD.Tcells$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.Tcells@meta.data$Predicted_subType
head(names(CD.Tcells@active.ident))
names(CellType) <- names(CD.Tcells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.Tcells) <- CD.Tcells@meta.data$group

warnings()

saveRDS(CD.Tcells, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.Tcells.sorted.w.scSorter5.GL14.RDS")

##sort PreT T cell with GL14
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
table(PreT.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList14.csv", sep = ",")

Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.Tcells <- subset (PreT.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.Tcells <- NormalizeData(PreT.Tcells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.Tcells <- FindVariableFeatures(PreT.Tcells, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.Tcells <- ScaleData(PreT.Tcells)
PreT.Tcells <- RunPCA(PreT.Tcells, features = VariableFeatures(PreT.Tcells))

topgenes <- head(VariableFeatures(PreT.Tcells), 3000)
expr = GetAssayData(PreT.Tcells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
PreT.Tcells$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.Tcells@meta.data$Predicted_subType
head(names(PreT.Tcells@active.ident))
names(CellType) <- names(PreT.Tcells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.Tcells) <- PreT.Tcells@meta.data$group

warnings()

saveRDS(PreT.Tcells, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.Tcells.sorted.w.scSorter5.GL14.RDS")

##sort MWD.T T cell with GL14
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
table(MWD.T.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList14.csv", sep = ",")

Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.Tcells <- subset (MWD.T.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.Tcells <- NormalizeData(MWD.T.Tcells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.Tcells <- FindVariableFeatures(MWD.T.Tcells, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.Tcells <- ScaleData(MWD.T.Tcells)
MWD.T.Tcells <- RunPCA(MWD.T.Tcells, features = VariableFeatures(MWD.T.Tcells))

topgenes <- head(VariableFeatures(MWD.T.Tcells), 3000)
expr = GetAssayData(MWD.T.Tcells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MWD.T.Tcells$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.Tcells@meta.data$Predicted_subType
head(names(MWD.T.Tcells@active.ident))
names(CellType) <- names(MWD.T.Tcells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.Tcells) <- MWD.T.Tcells@meta.data$group

warnings()

saveRDS(MWD.T.Tcells, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.Tcells.sorted.w.scSorter5.GL14.RDS")

##sort MRD.T T cell with GL14
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
table(MRD.T.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList14.csv", sep = ",")

Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.Tcells <- subset (MRD.T.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.Tcells <- NormalizeData(MRD.T.Tcells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.Tcells <- FindVariableFeatures(MRD.T.Tcells, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.Tcells <- ScaleData(MRD.T.Tcells)
MRD.T.Tcells <- RunPCA(MRD.T.Tcells, features = VariableFeatures(MRD.T.Tcells))

topgenes <- head(VariableFeatures(MRD.T.Tcells), 3000)
expr = GetAssayData(MRD.T.Tcells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MRD.T.Tcells$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.Tcells@meta.data$Predicted_subType
head(names(MRD.T.Tcells@active.ident))
names(CellType) <- names(MRD.T.Tcells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.Tcells) <- MRD.T.Tcells@meta.data$group

warnings()

saveRDS(MRD.T.Tcells, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.Tcells.sorted.w.scSorter5.GL14.RDS")

##sort MRD.NT T cell with GL14
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
table(MRD.NT.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList14.csv", sep = ",")

Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.Tcells <- subset (MRD.NT.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.Tcells <- NormalizeData(MRD.NT.Tcells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.Tcells <- FindVariableFeatures(MRD.NT.Tcells, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.Tcells <- ScaleData(MRD.NT.Tcells)
MRD.NT.Tcells <- RunPCA(MRD.NT.Tcells, features = VariableFeatures(MRD.NT.Tcells))

topgenes <- head(VariableFeatures(MRD.NT.Tcells), 3000)
expr = GetAssayData(MRD.NT.Tcells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MRD.NT.Tcells$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.Tcells@meta.data$Predicted_subType
head(names(MRD.NT.Tcells@active.ident))
names(CellType) <- names(MRD.NT.Tcells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.Tcells) <- MRD.NT.Tcells@meta.data$group

warnings()

saveRDS(MRD.NT.Tcells, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.Tcells.sorted.w.scSorter5.GL14.RDS")

##merge and save all results for GL14 in each group for reference
Tcells.annotated.individually <- merge(CD.Tcells, y = c(PreT.Tcells, MWD.T.Tcells, MRD.T.Tcells, MRD.NT.Tcells))
table(Tcells.annotated.individually@meta.data$group)
##now save RDS of the annotated data merged for GL14
saveRDS(Tcells.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Tcells.annotated.individually.w.scSorter5.GL14.RDS")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
##sort again with GL18 for each group 
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
table(CD.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList18.csv", sep = ",")

Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.Tcell <- subset (CD.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.Tcell <- NormalizeData(CD.Tcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.Tcell <- FindVariableFeatures(CD.Tcell, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.Tcell <- ScaleData(CD.Tcell)
CD.Tcell <- RunPCA(CD.Tcell, features = VariableFeatures(CD.Tcell))

topgenes <- head(VariableFeatures(CD.Tcell), 3000)
expr = GetAssayData(CD.Tcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
CD.Tcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.Tcell@meta.data$Predicted_subType
head(names(CD.Tcell@active.ident))
names(CellType) <- names(CD.Tcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.Tcell) <- CD.Tcell@meta.data$group

warnings()

saveRDS(CD.Tcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.Tcell.sorted.w.scSorter5.GL18.RDS")

##sort PreT T cell with GL18
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
table(PreT.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList18.csv", sep = ",")

Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.Tcell <- subset (PreT.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.Tcell <- NormalizeData(PreT.Tcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.Tcell <- FindVariableFeatures(PreT.Tcell, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.Tcell <- ScaleData(PreT.Tcell)
PreT.Tcell <- RunPCA(PreT.Tcell, features = VariableFeatures(PreT.Tcell))

topgenes <- head(VariableFeatures(PreT.Tcell), 3000)
expr = GetAssayData(PreT.Tcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
PreT.Tcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.Tcell@meta.data$Predicted_subType
head(names(PreT.Tcell@active.ident))
names(CellType) <- names(PreT.Tcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.Tcell) <- PreT.Tcell@meta.data$group

warnings()

saveRDS(PreT.Tcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.Tcell.sorted.w.scSorter5.GL18.RDS")

##sort MWD.T T cell with GL18
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
table(MWD.T.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList18.csv", sep = ",")

Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.Tcell <- subset (MWD.T.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.Tcell <- NormalizeData(MWD.T.Tcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.Tcell <- FindVariableFeatures(MWD.T.Tcell, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.Tcell <- ScaleData(MWD.T.Tcell)
MWD.T.Tcell <- RunPCA(MWD.T.Tcell, features = VariableFeatures(MWD.T.Tcell))

topgenes <- head(VariableFeatures(MWD.T.Tcell), 3000)
expr = GetAssayData(MWD.T.Tcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MWD.T.Tcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.Tcell@meta.data$Predicted_subType
head(names(MWD.T.Tcell@active.ident))
names(CellType) <- names(MWD.T.Tcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.Tcell) <- MWD.T.Tcell@meta.data$group

warnings()

saveRDS(MWD.T.Tcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.Tcell.sorted.w.scSorter5.GL18.RDS")

##sort MRD.T T cell with GL18
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
table(MRD.T.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList18.csv", sep = ",")

Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.Tcell <- subset (MRD.T.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.Tcell <- NormalizeData(MRD.T.Tcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.Tcell <- FindVariableFeatures(MRD.T.Tcell, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.Tcell <- ScaleData(MRD.T.Tcell)
MRD.T.Tcell <- RunPCA(MRD.T.Tcell, features = VariableFeatures(MRD.T.Tcell))

topgenes <- head(VariableFeatures(MRD.T.Tcell), 3000)
expr = GetAssayData(MRD.T.Tcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MRD.T.Tcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.Tcell@meta.data$Predicted_subType
head(names(MRD.T.Tcell@active.ident))
names(CellType) <- names(MRD.T.Tcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.Tcell) <- MRD.T.Tcell@meta.data$group

warnings()

saveRDS(MRD.T.Tcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.Tcell.sorted.w.scSorter5.GL18.RDS")

##sort MRD.NT T cell with GL18
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
table(MRD.NT.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList18.csv", sep = ",")

Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.Tcell <- subset (MRD.NT.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.Tcell <- NormalizeData(MRD.NT.Tcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.Tcell <- FindVariableFeatures(MRD.NT.Tcell, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.Tcell <- ScaleData(MRD.NT.Tcell)
MRD.NT.Tcell <- RunPCA(MRD.NT.Tcell, features = VariableFeatures(MRD.NT.Tcell))

topgenes <- head(VariableFeatures(MRD.NT.Tcell), 3000)
expr = GetAssayData(MRD.NT.Tcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MRD.NT.Tcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.Tcell@meta.data$Predicted_subType
head(names(MRD.NT.Tcell@active.ident))
names(CellType) <- names(MRD.NT.Tcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.Tcell) <- MRD.NT.Tcell@meta.data$group

warnings()

saveRDS(MRD.NT.Tcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.Tcell.sorted.w.scSorter5.GL18.RDS")

##merge and save all results for GL18 in each group for reference
Tcell.annotated.individually <- merge(CD.Tcell, y = c(PreT.Tcell, MWD.T.Tcell, MRD.T.Tcell, MRD.NT.Tcell))
table(Tcell.annotated.individually@meta.data$group)
##now save RDS of the annotated data merged for GL18
saveRDS(Tcell.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Tcell.annotated.individually.w.scSorter5.GL18.RDS")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
##sort again with GL19 for each group 
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
table(CD.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList19.csv", sep = ",")

Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.Tcel <- subset (CD.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.Tcel <- NormalizeData(CD.Tcel, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.Tcel <- FindVariableFeatures(CD.Tcel, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.Tcel <- ScaleData(CD.Tcel)
CD.Tcel <- RunPCA(CD.Tcel, features = VariableFeatures(CD.Tcel))

topgenes <- head(VariableFeatures(CD.Tcel), 3000)
expr = GetAssayData(CD.Tcel)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
CD.Tcel$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.Tcel@meta.data$Predicted_subType
head(names(CD.Tcel@active.ident))
names(CellType) <- names(CD.Tcel@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.Tcel) <- CD.Tcel@meta.data$group

warnings()

saveRDS(CD.Tcel, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.Tcel.sorted.w.scSorter5.GL19.RDS")

##sort PreT T cell with GL19
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
table(PreT.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList19.csv", sep = ",")

Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.Tcel <- subset (PreT.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.Tcel <- NormalizeData(PreT.Tcel, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.Tcel <- FindVariableFeatures(PreT.Tcel, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.Tcel <- ScaleData(PreT.Tcel)
PreT.Tcel <- RunPCA(PreT.Tcel, features = VariableFeatures(PreT.Tcel))

topgenes <- head(VariableFeatures(PreT.Tcel), 3000)
expr = GetAssayData(PreT.Tcel)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
PreT.Tcel$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.Tcel@meta.data$Predicted_subType
head(names(PreT.Tcel@active.ident))
names(CellType) <- names(PreT.Tcel@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.Tcel) <- PreT.Tcel@meta.data$group

warnings()

saveRDS(PreT.Tcel, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.Tcel.sorted.w.scSorter5.GL19.RDS")

##sort MWD.T T cell with GL19
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
table(MWD.T.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList19.csv", sep = ",")

Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.Tcel <- subset (MWD.T.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.Tcel <- NormalizeData(MWD.T.Tcel, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.Tcel <- FindVariableFeatures(MWD.T.Tcel, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.Tcel <- ScaleData(MWD.T.Tcel)
MWD.T.Tcel <- RunPCA(MWD.T.Tcel, features = VariableFeatures(MWD.T.Tcel))

topgenes <- head(VariableFeatures(MWD.T.Tcel), 3000)
expr = GetAssayData(MWD.T.Tcel)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MWD.T.Tcel$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.Tcel@meta.data$Predicted_subType
head(names(MWD.T.Tcel@active.ident))
names(CellType) <- names(MWD.T.Tcel@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.Tcel) <- MWD.T.Tcel@meta.data$group

warnings()

saveRDS(MWD.T.Tcel, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.Tcel.sorted.w.scSorter5.GL19.RDS")

##sort MRD.T T cell with GL19
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
table(MRD.T.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList19.csv", sep = ",")

Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.Tcel <- subset (MRD.T.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.Tcel <- NormalizeData(MRD.T.Tcel, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.Tcel <- FindVariableFeatures(MRD.T.Tcel, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.Tcel <- ScaleData(MRD.T.Tcel)
MRD.T.Tcel <- RunPCA(MRD.T.Tcel, features = VariableFeatures(MRD.T.Tcel))

topgenes <- head(VariableFeatures(MRD.T.Tcel), 3000)
expr = GetAssayData(MRD.T.Tcel)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MRD.T.Tcel$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.Tcel@meta.data$Predicted_subType
head(names(MRD.T.Tcel@active.ident))
names(CellType) <- names(MRD.T.Tcel@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.Tcel) <- MRD.T.Tcel@meta.data$group

warnings()

saveRDS(MRD.T.Tcel, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.Tcel.sorted.w.scSorter5.GL19.RDS")

##sort MRD.NT T cell with GL19
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
table(MRD.NT.cells@meta.data$Predicted_Type)

markers.T <- read.delim("~/Reverse_Diet/GeneLists/GeneList19.csv", sep = ",")

Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.Tcel <- subset (MRD.NT.cells, idents = "T cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.Tcel <- NormalizeData(MRD.NT.Tcel, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.Tcel <- FindVariableFeatures(MRD.NT.Tcel, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.Tcel <- ScaleData(MRD.NT.Tcel)
MRD.NT.Tcel <- RunPCA(MRD.NT.Tcel, features = VariableFeatures(MRD.NT.Tcel))

topgenes <- head(VariableFeatures(MRD.NT.Tcel), 3000)
expr = GetAssayData(MRD.NT.Tcel)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.T$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.T)
MRD.NT.Tcel$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.Tcel@meta.data$Predicted_subType
head(names(MRD.NT.Tcel@active.ident))
names(CellType) <- names(MRD.NT.Tcel@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.Tcel) <- MRD.NT.Tcel@meta.data$group

warnings()

saveRDS(MRD.NT.Tcel, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.Tcel.sorted.w.scSorter5.GL19.RDS")

##merge and save all results for GL19 in each group for reference
Tcel.annotated.individually <- merge(CD.Tcel, y = c(PreT.Tcel, MWD.T.Tcel, MRD.T.Tcel, MRD.NT.Tcel))
table(Tcel.annotated.individually@meta.data$group)
##now save RDS of the annotated data merged for GL14
saveRDS(Tcel.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Tcel.annotated.individually.w.scSorter5.GL19.RDS")

