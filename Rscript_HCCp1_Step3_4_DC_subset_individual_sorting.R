library("ggplot2")
library("Seurat")
library("scSorter")

HCCp1.sorted.individually <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.individually.w.scSorter5.RDS")
table(HCCp1.sorted.individually@meta.data$experimental.group)
table(HCCp1.sorted.individually@meta.data$group)
table(HCCp1.sorted.individually@meta.data$replicate)
table(HCCp1.sorted.individually@meta.data$Predicted_Type)

Idents(HCCp1.sorted.individually) <- HCCp1.sorted.individually@meta.data$group

##subset cells from each group to sort each individually on DCs
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

##first subset CD group DCs and sort these for subsets
markers.DC <- read.delim("~/Reverse_Diet/GeneLists/GeneList10.csv", sep = ",")

Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.DCs <- subset (CD.cells, idents = "Dendritic cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.DCs <- NormalizeData(CD.DCs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.DCs <- FindVariableFeatures(CD.DCs, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.DCs <- ScaleData(CD.DCs)
CD.DCs <- RunPCA(CD.DCs, features = VariableFeatures(CD.DCs))

topgenes <- head(VariableFeatures(CD.DCs), 3000)
expr = GetAssayData(CD.DCs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.DC$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.DC)
CD.DCs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.DCs@meta.data$Predicted_subType
head(names(CD.DCs@active.ident))
names(CellType) <- names(CD.DCs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.DCs) <- CD.DCs@meta.data$group

warnings()

saveRDS(CD.DCs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.DCs.sorted.w.scSorter5.GL10.RDS")

##second subset PreT group DCs and sort these for subsets
markers.DC <- read.delim("~/Reverse_Diet/GeneLists/GeneList10.csv", sep = ",")

Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.DCs <- subset (PreT.cells, idents = "Dendritic cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.DCs <- NormalizeData(PreT.DCs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.DCs <- FindVariableFeatures(PreT.DCs, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.DCs <- ScaleData(PreT.DCs)
PreT.DCs <- RunPCA(PreT.DCs, features = VariableFeatures(PreT.DCs))

topgenes <- head(VariableFeatures(PreT.DCs), 3000)
expr = GetAssayData(PreT.DCs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.DC$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.DC)
PreT.DCs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.DCs@meta.data$Predicted_subType
head(names(PreT.DCs@active.ident))
names(CellType) <- names(PreT.DCs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.DCs) <- PreT.DCs@meta.data$group

warnings()

saveRDS(PreT.DCs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.DCs.sorted.w.scSorter5.GL10.RDS")

##third subset MWD.T group DCs and sort these for subsets
markers.DC <- read.delim("~/Reverse_Diet/GeneLists/GeneList10.csv", sep = ",")

Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.DCs <- subset (MWD.T.cells, idents = "Dendritic cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.DCs <- NormalizeData(MWD.T.DCs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.DCs <- FindVariableFeatures(MWD.T.DCs, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.DCs <- ScaleData(MWD.T.DCs)
MWD.T.DCs <- RunPCA(MWD.T.DCs, features = VariableFeatures(MWD.T.DCs))

topgenes <- head(VariableFeatures(MWD.T.DCs), 3000)
expr = GetAssayData(MWD.T.DCs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.DC$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.DC)
MWD.T.DCs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.DCs@meta.data$Predicted_subType
head(names(MWD.T.DCs@active.ident))
names(CellType) <- names(MWD.T.DCs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.DCs) <- MWD.T.DCs@meta.data$group

warnings()

saveRDS(MWD.T.DCs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.DCs.sorted.w.scSorter5.GL10.RDS")

##fourth subset MRD.T group DCs and sort these for subsets
markers.DC <- read.delim("~/Reverse_Diet/GeneLists/GeneList10.csv", sep = ",")

Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.DCs <- subset (MRD.T.cells, idents = "Dendritic cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.DCs <- NormalizeData(MRD.T.DCs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.DCs <- FindVariableFeatures(MRD.T.DCs, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.DCs <- ScaleData(MRD.T.DCs)
MRD.T.DCs <- RunPCA(MRD.T.DCs, features = VariableFeatures(MRD.T.DCs))

topgenes <- head(VariableFeatures(MRD.T.DCs), 3000)
expr = GetAssayData(MRD.T.DCs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.DC$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.DC)
MRD.T.DCs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.DCs@meta.data$Predicted_subType
head(names(MRD.T.DCs@active.ident))
names(CellType) <- names(MRD.T.DCs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.DCs) <- MRD.T.DCs@meta.data$group

warnings()

saveRDS(MRD.T.DCs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.DCs.sorted.w.scSorter5.GL10.RDS")

##finally subset MRD.NT group DCs and sort these for subsets
markers.DC <- read.delim("~/Reverse_Diet/GeneLists/GeneList10.csv", sep = ",")

Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.DCs <- subset (MRD.NT.cells, idents = "Dendritic cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.DCs <- NormalizeData(MRD.NT.DCs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.DCs <- FindVariableFeatures(MRD.NT.DCs, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.DCs <- ScaleData(MRD.NT.DCs)
MRD.NT.DCs <- RunPCA(MRD.NT.DCs, features = VariableFeatures(MRD.NT.DCs))

topgenes <- head(VariableFeatures(MRD.NT.DCs), 3000)
expr = GetAssayData(MRD.NT.DCs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.DC$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.DC)
MRD.NT.DCs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.DCs@meta.data$Predicted_subType
head(names(MRD.NT.DCs@active.ident))
names(CellType) <- names(MRD.NT.DCs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.DCs) <- MRD.NT.DCs@meta.data$group

warnings()

saveRDS(MRD.NT.DCs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.DCs.sorted.w.scSorter5.GL10.RDS")

##combine all sorted subsets and save.RDS for re-combining all sorted subtypes later
DCs.annotated.individually <- merge(CD.DCs, y = c(PreT.DCs, MWD.T.DCs, MRD.T.DCs, MRD.NT.DCs))
table(DCs.annotated.individually@meta.data$group)
##now save RDS of the annotated data to try and visualize with UMAP downstream
saveRDS(DCs.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/DCs.annotated.individually.w.scSorter5.GL10.RDS")
