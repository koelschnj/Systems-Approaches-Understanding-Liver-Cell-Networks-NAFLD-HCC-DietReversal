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

##first sort CD Bcells 
markers.B <- read.delim("~/Reverse_Diet/GeneLists/GeneList16.csv", sep = ",")

##run scSorter to based on GL16 marker genes provided for B cell subsets 
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
CD.Bcell <- subset (CD.cells, idents = "B cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
CD.Bcell <- NormalizeData(CD.Bcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.Bcell <- FindVariableFeatures(CD.Bcell, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.Bcell <- ScaleData(CD.Bcell)
CD.Bcell <- RunPCA(CD.Bcell, features = VariableFeatures(CD.Bcell))

topgenes <- head(VariableFeatures(CD.Bcell), 3000)
expr = GetAssayData(CD.Bcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.B$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.B)
CD.Bcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.Bcell@meta.data$Predicted_subType
head(names(CD.Bcell@active.ident))
names(CellType) <- names(CD.Bcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.Bcell) <- CD.Bcell@meta.data$group

warnings()

saveRDS(CD.Bcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/CD.Bcell.sorted.w.scSorter5.RDS")

##second sort PreT B cells
markers.B <- read.delim("~/Reverse_Diet/GeneLists/GeneList16.csv", sep = ",")

##run scSorter to based on GL16 marker genes provided for B cell subsets 
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
PreT.Bcell <- subset (PreT.cells, idents = "B cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
PreT.Bcell <- NormalizeData(PreT.Bcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.Bcell <- FindVariableFeatures(PreT.Bcell, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.Bcell <- ScaleData(PreT.Bcell)
PreT.Bcell <- RunPCA(PreT.Bcell, features = VariableFeatures(PreT.Bcell))

topgenes <- head(VariableFeatures(PreT.Bcell), 3000)
expr = GetAssayData(PreT.Bcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.B$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.B)
PreT.Bcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.Bcell@meta.data$Predicted_subType
head(names(PreT.Bcell@active.ident))
names(CellType) <- names(PreT.Bcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.Bcell) <- PreT.Bcell@meta.data$group

warnings()

saveRDS(PreT.Bcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/PreT.Bcell.sorted.w.scSorter5.RDS")

##third sort MWD.T B cells
markers.B <- read.delim("~/Reverse_Diet/GeneLists/GeneList16.csv", sep = ",")

##run scSorter to based on GL16 marker genes provided for B cell subsets 
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
MWD.T.Bcell <- subset (MWD.T.cells, idents = "B cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MWD.T.Bcell <- NormalizeData(MWD.T.Bcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.Bcell <- FindVariableFeatures(MWD.T.Bcell, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.Bcell <- ScaleData(MWD.T.Bcell)
MWD.T.Bcell <- RunPCA(MWD.T.Bcell, features = VariableFeatures(MWD.T.Bcell), npcs = 30)

topgenes <- head(VariableFeatures(MWD.T.Bcell), 3000)
expr = GetAssayData(MWD.T.Bcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.B$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.B)
MWD.T.Bcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.Bcell@meta.data$Predicted_subType
head(names(MWD.T.Bcell@active.ident))
names(CellType) <- names(MWD.T.Bcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.Bcell) <- MWD.T.Bcell@meta.data$group

warnings()

saveRDS(MWD.T.Bcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MWD.T.Bcell.sorted.w.scSorter5.RDS")

##fourth sort MRD.T B cells
markers.B <- read.delim("~/Reverse_Diet/GeneLists/GeneList16.csv", sep = ",")

##run scSorter to based on GL16 marker genes provided for B cell subsets 
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
MRD.T.Bcell <- subset (MRD.T.cells, idents = "B cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.T.Bcell <- NormalizeData(MRD.T.Bcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.Bcell <- FindVariableFeatures(MRD.T.Bcell, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.Bcell <- ScaleData(MRD.T.Bcell)
MRD.T.Bcell <- RunPCA(MRD.T.Bcell, features = VariableFeatures(MRD.T.Bcell))

topgenes <- head(VariableFeatures(MRD.T.Bcell), 3000)
expr = GetAssayData(MRD.T.Bcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.B$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.B)
MRD.T.Bcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.Bcell@meta.data$Predicted_subType
head(names(MRD.T.Bcell@active.ident))
names(CellType) <- names(MRD.T.Bcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.Bcell) <- MRD.T.Bcell@meta.data$group

warnings()

saveRDS(MRD.T.Bcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.T.Bcell.sorted.w.scSorter5.RDS")

##finally sort MRD.NT B cells
markers.B <- read.delim("~/Reverse_Diet/GeneLists/GeneList16.csv", sep = ",")

##run scSorter to based on GL16 marker genes provided for B cell subsets 
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
MRD.NT.Bcell <- subset (MRD.NT.cells, idents = "B cell", invert = FALSE)

##find variable features again for Mac population then run scSorter
MRD.NT.Bcell <- NormalizeData(MRD.NT.Bcell, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.Bcell <- FindVariableFeatures(MRD.NT.Bcell, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.Bcell <- ScaleData(MRD.NT.Bcell)
MRD.NT.Bcell <- RunPCA(MRD.NT.Bcell, features = VariableFeatures(MRD.NT.Bcell))

topgenes <- head(VariableFeatures(MRD.NT.Bcell), 3000)
expr = GetAssayData(MRD.NT.Bcell)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers.B$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers.B)
MRD.NT.Bcell$Predicted_subType <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.Bcell@meta.data$Predicted_subType
head(names(MRD.NT.Bcell@active.ident))
names(CellType) <- names(MRD.NT.Bcell@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.Bcell) <- MRD.NT.Bcell@meta.data$group

warnings()

saveRDS(MRD.NT.Bcell, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MRD.NT.Bcell.sorted.w.scSorter5.RDS")

##combine all sorted subsets and save.RDS for re-combining all sorted subtypes later
Bcells.annotated.individually <- merge(CD.Bcell, y = c(PreT.Bcell, MWD.T.Bcell, MRD.T.Bcell, MRD.NT.Bcell))
table(Bcells.annotated.individually@meta.data$group)
##now save RDS of the annotated data to try and visualize with UMAP downstream
saveRDS(Bcells.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Bcells.annotated.individually.w.scSorter5.RDS")
