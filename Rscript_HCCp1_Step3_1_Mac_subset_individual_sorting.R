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

##first subset CD group Macrophages (MPs) and sort these for subsets
markers.MP <- read.delim("~/Reverse_Diet/GeneLists/GeneList8.csv", sep = ",")

##run scSorter to based on marker genes provided for macrophage subsets
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.MPs <- subset (CD.cells, idents = "Macrophage", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.MPs <- NormalizeData(CD.MPs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.MPs <- FindVariableFeatures(CD.MPs, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.MPs <- ScaleData(CD.MPs)
CD.MPs <- RunPCA(CD.MPs, features = VariableFeatures(CD.MPs))

topgenes <- head(VariableFeatures(CD.MPs), 3000)
expr = GetAssayData(CD.MPs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.MP$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.MP)
CD.MPs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.MPs@meta.data$Predicted_subType
head(names(CD.MPs@active.ident))
names(CellType) <- names(CD.MPs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.MPs) <- CD.MPs@meta.data$group

warnings()

saveRDS(CD.MPs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.MPs.sorted.w.scSorter5.RDS")

##second for PreT MPs and sort subsets
markers.MP <- read.delim("~/Reverse_Diet/GeneLists/GeneList8.csv", sep = ",")

##run scSorter to based on marker genes provided for macrophage subsets
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.MPs <- subset (PreT.cells, idents = "Macrophage", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.MPs <- NormalizeData(PreT.MPs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.MPs <- FindVariableFeatures(PreT.MPs, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.MPs <- ScaleData(PreT.MPs)
PreT.MPs <- RunPCA(PreT.MPs, features = VariableFeatures(PreT.MPs))

topgenes <- head(VariableFeatures(PreT.MPs), 3000)
expr = GetAssayData(PreT.MPs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.MP$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.MP)
PreT.MPs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.MPs@meta.data$Predicted_subType
head(names(PreT.MPs@active.ident))
names(CellType) <- names(PreT.MPs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.MPs) <- PreT.MPs@meta.data$group

warnings()

saveRDS(PreT.MPs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.MPs.sorted.w.scSorter5.RDS")

##third for MWD.T MPs and sort subsets
markers.MP <- read.delim("~/Reverse_Diet/GeneLists/GeneList8.csv", sep = ",")

##run scSorter to based on marker genes provided for macrophage subsets
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.MPs <- subset (MWD.T.cells, idents = "Macrophage", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.MPs <- NormalizeData(MWD.T.MPs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.MPs <- FindVariableFeatures(MWD.T.MPs, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.MPs <- ScaleData(MWD.T.MPs)
MWD.T.MPs <- RunPCA(MWD.T.MPs, features = VariableFeatures(MWD.T.MPs))

topgenes <- head(VariableFeatures(MWD.T.MPs), 3000)
expr = GetAssayData(MWD.T.MPs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.MP$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.MP)
MWD.T.MPs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.MPs@meta.data$Predicted_subType
head(names(MWD.T.MPs@active.ident))
names(CellType) <- names(MWD.T.MPs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.MPs) <- MWD.T.MPs@meta.data$group

warnings()

saveRDS(MWD.T.MPs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.MPs.sorted.w.scSorter5.RDS")

##fourth for MRD.T MPs and sort subsets
markers.MP <- read.delim("~/Reverse_Diet/GeneLists/GeneList8.csv", sep = ",")

##run scSorter to based on marker genes provided for macrophage subsets
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.MPs <- subset (MRD.T.cells, idents = "Macrophage", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.MPs <- NormalizeData(MRD.T.MPs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.MPs <- FindVariableFeatures(MRD.T.MPs, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.MPs <- ScaleData(MRD.T.MPs)
MRD.T.MPs <- RunPCA(MRD.T.MPs, features = VariableFeatures(MRD.T.MPs))

topgenes <- head(VariableFeatures(MRD.T.MPs), 3000)
expr = GetAssayData(MRD.T.MPs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.MP$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.MP)
MRD.T.MPs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.MPs@meta.data$Predicted_subType
head(names(MRD.T.MPs@active.ident))
names(CellType) <- names(MRD.T.MPs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.MPs) <- MRD.T.MPs@meta.data$group

warnings()

saveRDS(MRD.T.MPs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.MPs.sorted.w.scSorter5.RDS")

##finally perform for MRD.NT and sort subsets
markers.MP <- read.delim("~/Reverse_Diet/GeneLists/GeneList8.csv", sep = ",")

##run scSorter to based on marker genes provided for macrophage subsets
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.MPs <- subset (MRD.NT.cells, idents = "Macrophage", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.MPs <- NormalizeData(MRD.NT.MPs, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.MPs <- FindVariableFeatures(MRD.NT.MPs, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.MPs <- ScaleData(MRD.NT.MPs)
MRD.NT.MPs <- RunPCA(MRD.NT.MPs, features = VariableFeatures(MRD.NT.MPs))

topgenes <- head(VariableFeatures(MRD.NT.MPs), 3000)
expr = GetAssayData(MRD.NT.MPs)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.MP$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.MP)
MRD.NT.MPs$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.MPs@meta.data$Predicted_subType
head(names(MRD.NT.MPs@active.ident))
names(CellType) <- names(MRD.NT.MPs@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.MPs) <- MRD.NT.MPs@meta.data$group

warnings()

saveRDS(MRD.NT.MPs, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.MPs.sorted.w.scSorter5.RDS")

##combine all sorted subsets and save.RDS for re-combining all sorted subtypes later
MPs.annotated.individually <- merge(CD.MPs, y = c(PreT.MPs, MWD.T.MPs, MRD.T.MPs, MRD.NT.MPs))
table(MPs.annotated.individually@meta.data$group)
##now save RDS of the annotated data to try and visualize with UMAP downstream
saveRDS(MPs.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MPs.annotated.individually.w.scSorter5.RDS")
