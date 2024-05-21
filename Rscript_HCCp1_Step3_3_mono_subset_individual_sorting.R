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

##load in marker gene list for monocytes/mMDSCs and sort for CD first
markers.mono <- read.delim("~/Reverse_Diet/GeneLists/GeneList17.csv", sep = ",")

##run scSorter to based on GL17 marker genes provided for monocytic subtypes
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.Mono <- subset (CD.cells, idents = "Monocyte", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.Mono <- NormalizeData(CD.Mono, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.Mono <- FindVariableFeatures(CD.Mono, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.Mono <- ScaleData(CD.Mono)
CD.Mono <- RunPCA(CD.Mono, features = VariableFeatures(CD.Mono))

topgenes <- head(VariableFeatures(CD.Mono), 3000)
expr = GetAssayData(CD.Mono)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.mono$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.mono)
CD.Mono$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.Mono@meta.data$Predicted_subType
head(names(CD.Mono@active.ident))
names(CellType) <- names(CD.Mono@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.Mono) <- CD.Mono@meta.data$group

warnings()

saveRDS(CD.Mono, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.Mono.sorted.w.scSorter5.RDS")

##second subset PreT monocytes and sort 
markers.mono <- read.delim("~/Reverse_Diet/GeneLists/GeneList17.csv", sep = ",")

##run scSorter to based on GL17 marker genes provided for monocytic subtypes
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.Mono <- subset (PreT.cells, idents = "Monocyte", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.Mono <- NormalizeData(PreT.Mono, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.Mono <- FindVariableFeatures(PreT.Mono, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.Mono <- ScaleData(PreT.Mono)
PreT.Mono <- RunPCA(PreT.Mono, features = VariableFeatures(PreT.Mono))

topgenes <- head(VariableFeatures(PreT.Mono), 3000)
expr = GetAssayData(PreT.Mono)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.mono$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.mono)
PreT.Mono$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.Mono@meta.data$Predicted_subType
head(names(PreT.Mono@active.ident))
names(CellType) <- names(PreT.Mono@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.Mono) <- PreT.Mono@meta.data$group

warnings()

saveRDS(PreT.Mono, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.Mono.sorted.w.scSorter5.RDS")

##third subset MWD.T monocytes and sort
markers.mono <- read.delim("~/Reverse_Diet/GeneLists/GeneList17.csv", sep = ",")

##run scSorter to based on GL17 marker genes provided for monocytic subtypes
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.Mono <- subset (MWD.T.cells, idents = "Monocyte", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.Mono <- NormalizeData(MWD.T.Mono, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.Mono <- FindVariableFeatures(MWD.T.Mono, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.Mono <- ScaleData(MWD.T.Mono)
MWD.T.Mono <- RunPCA(MWD.T.Mono, features = VariableFeatures(MWD.T.Mono))

topgenes <- head(VariableFeatures(MWD.T.Mono), 3000)
expr = GetAssayData(MWD.T.Mono)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.mono$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.mono)
MWD.T.Mono$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.Mono@meta.data$Predicted_subType
head(names(MWD.T.Mono@active.ident))
names(CellType) <- names(MWD.T.Mono@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.Mono) <- MWD.T.Mono@meta.data$group

warnings()

saveRDS(MWD.T.Mono, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.Mono.sorted.w.scSorter5.RDS")

##fourth sort MRD.T monocytes
markers.mono <- read.delim("~/Reverse_Diet/GeneLists/GeneList17.csv", sep = ",")

##run scSorter to based on GL17 marker genes provided for monocytic subtypes
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.Mono <- subset (MRD.T.cells, idents = "Monocyte", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.Mono <- NormalizeData(MRD.T.Mono, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.Mono <- FindVariableFeatures(MRD.T.Mono, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.Mono <- ScaleData(MRD.T.Mono)
MRD.T.Mono <- RunPCA(MRD.T.Mono, features = VariableFeatures(MRD.T.Mono))

topgenes <- head(VariableFeatures(MRD.T.Mono), 3000)
expr = GetAssayData(MRD.T.Mono)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.mono$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.mono)
MRD.T.Mono$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.Mono@meta.data$Predicted_subType
head(names(MRD.T.Mono@active.ident))
names(CellType) <- names(MRD.T.Mono@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.Mono) <- MRD.T.Mono@meta.data$group

warnings()

saveRDS(MRD.T.Mono, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.Mono.sorted.w.scSorter5.RDS")

##finally sort MRD.NT monocytes and merge all to save
markers.mono <- read.delim("~/Reverse_Diet/GeneLists/GeneList17.csv", sep = ",")

##run scSorter to based on GL17 marker genes provided for monocytic subtypes
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.Mono <- subset (MRD.NT.cells, idents = "Monocyte", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.Mono <- NormalizeData(MRD.NT.Mono, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.Mono <- FindVariableFeatures(MRD.NT.Mono, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.Mono <- ScaleData(MRD.NT.Mono)
MRD.NT.Mono <- RunPCA(MRD.NT.Mono, features = VariableFeatures(MRD.NT.Mono))

topgenes <- head(VariableFeatures(MRD.NT.Mono), 3000)
expr = GetAssayData(MRD.NT.Mono)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.mono$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.mono)
MRD.NT.Mono$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.Mono@meta.data$Predicted_subType
head(names(MRD.NT.Mono@active.ident))
names(CellType) <- names(MRD.NT.Mono@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.Mono) <- MRD.NT.Mono@meta.data$group

warnings()

saveRDS(MRD.NT.Mono, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.Mono.sorted.w.scSorter5.RDS")

##merge all and combine to save
Mono.annotated.individually <- merge(CD.Mono, y = c(PreT.Mono, MWD.T.Mono, MRD.T.Mono, MRD.NT.Mono))
table(Mono.annotated.individually@meta.data$group)
##now save RDS of the annotated data to try and visualize with UMAP downstream
saveRDS(Mono.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Mono.annotated.individually.w.scSorter5.RDS")
