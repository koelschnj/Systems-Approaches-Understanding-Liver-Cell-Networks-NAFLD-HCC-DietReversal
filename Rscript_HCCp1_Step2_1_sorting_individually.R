library("ggplot2")
library("Seurat")
library("scSorter")

##load in QC and filtered cells from Rscript_Step1_QC_Filtering
HCCp1.cells <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.group.combo.QCandFiltered.RDS")

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$group
table(HCCp1.cells@meta.data$group)

##going to sort each group individually to see how that changes the results
##separate groups first and then run scSorter on each individually
CD.cells <- subset(x = HCCp1.cells, idents = "CD", invert = FALSE)
table(CD.cells@meta.data$group)

PreT.cells <- subset(x = HCCp1.cells, idents = "PreT", invert = FALSE)
table(PreT.cells@meta.data$group)

MWD.T.cells <- subset(x = HCCp1.cells, idents = "MWD.T", invert = FALSE)
table(MWD.T.cells@meta.data$group)

MRD.T.cells <- subset(x = HCCp1.cells, idents = "MRD.T", invert = FALSE)
table(MRD.T.cells@meta.data$group)

MRD.NT.cells <- subset(x = HCCp1.cells, idents = "MRD.NT", invert = FALSE)
table(MRD.NT.cells@meta.data$group)

##load in custom list of marker genes from CellMarker2.0 database from liver specific experiments
##annotated CD first, then PreT, MWD.T, MRD.T, and finally MRD.NT
markers <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter6.csv", sep = ",")

Idents(CD.cells) <- CD.cells@meta.data$group
table(CD.cells@meta.data$group)

CD.cells <- NormalizeData(CD.cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
CD.cells <- FindVariableFeatures(CD.cells, selection.method = "vst", nfeatures = 3000, verbose = F)
CD.cells <- ScaleData(CD.cells)
CD.cells <- RunPCA(CD.cells, features = VariableFeatures(CD.cells))

##run scSorter to based on marker genes from CellMarker2.0
topgenes <- head(VariableFeatures(CD.cells), 3000)
expr = GetAssayData(CD.cells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers)
CD.cells$Predicted_Type <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- CD.cells@meta.data$Predicted_Type
head(names(CD.cells@active.ident))
names(CellType) <- names(CD.cells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(CD.cells) <- CD.cells@meta.data$group

##these will all be seen in log files and can then be recorded in excel
Idents(CD.cells) <- CD.cells@meta.data$Predicted_Type
table(CD.cells@meta.data$Predicted_Type)

Idents(CD.cells) <- CD.cells@meta.data$group

warnings()

saveRDS(CD.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.CD.cells.sorted.w.scSorter6.RDS")
##save annotated CD.cells and move on to sort PreT
markers1 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter6.csv", sep = ",")

Idents(PreT.cells) <- PreT.cells@meta.data$group
table(PreT.cells@meta.data$group)

PreT.cells <- NormalizeData(PreT.cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
PreT.cells <- FindVariableFeatures(PreT.cells, selection.method = "vst", nfeatures = 3000, verbose = F)
PreT.cells <- ScaleData(PreT.cells)
PreT.cells <- RunPCA(PreT.cells, features = VariableFeatures(PreT.cells))

##run scSorter to based on marker genes from CellMarker2.0
topgenes <- head(VariableFeatures(PreT.cells), 3000)
expr = GetAssayData(PreT.cells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers1$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers1)
PreT.cells$Predicted_Type <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- PreT.cells@meta.data$Predicted_Type
head(names(PreT.cells@active.ident))
names(CellType) <- names(PreT.cells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(PreT.cells) <- PreT.cells@meta.data$group

##these will all be seen in log files and can then be recorded in excel
Idents(PreT.cells) <- PreT.cells@meta.data$Predicted_Type
table(PreT.cells@meta.data$Predicted_Type)

Idents(PreT.cells) <- PreT.cells@meta.data$group

warnings()

saveRDS(PreT.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.PreT.cells.sorted.w.scSorter6.RDS")
##save annotated PreT.cells and move on to sort MWD.T
markers2 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter6.csv", sep = ",")

Idents(MWD.T.cells) <- MWD.T.cells@meta.data$group
table(MWD.T.cells@meta.data$group)

MWD.T.cells <- NormalizeData(MWD.T.cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MWD.T.cells <- FindVariableFeatures(MWD.T.cells, selection.method = "vst", nfeatures = 3000, verbose = F)
MWD.T.cells <- ScaleData(MWD.T.cells)
MWD.T.cells <- RunPCA(MWD.T.cells, features = VariableFeatures(MWD.T.cells))

##run scSorter to based on marker genes from CellMarker2.0
topgenes <- head(VariableFeatures(MWD.T.cells), 3000)
expr = GetAssayData(MWD.T.cells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers2$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers2)
MWD.T.cells$Predicted_Type <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MWD.T.cells@meta.data$Predicted_Type
head(names(MWD.T.cells@active.ident))
names(CellType) <- names(MWD.T.cells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$group

##these will all be seen in log files and can then be recorded in excel
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Predicted_Type
table(MWD.T.cells@meta.data$Predicted_Type)

Idents(MWD.T.cells) <- MWD.T.cells@meta.data$group

warnings()

saveRDS(MWD.T.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.MWD.T.cells.sorted.w.scSorter6.RDS")
##save annotated MWD.T.cells and move on to sort MRD.T
markers3 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter6.csv", sep = ",")

Idents(MRD.T.cells) <- MRD.T.cells@meta.data$group
table(MRD.T.cells@meta.data$group)

MRD.T.cells <- NormalizeData(MRD.T.cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.T.cells <- FindVariableFeatures(MRD.T.cells, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.T.cells <- ScaleData(MRD.T.cells)
MRD.T.cells <- RunPCA(MRD.T.cells, features = VariableFeatures(MRD.T.cells))

##run scSorter to based on marker genes from CellMarker2.0
topgenes <- head(VariableFeatures(MRD.T.cells), 3000)
expr = GetAssayData(MRD.T.cells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers3$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers3)
MRD.T.cells$Predicted_Type <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.T.cells@meta.data$Predicted_Type
head(names(MRD.T.cells@active.ident))
names(CellType) <- names(MRD.T.cells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$group

##these will all be seen in log files and can then be recorded in excel
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Predicted_Type
table(MRD.T.cells@meta.data$Predicted_Type)

Idents(MRD.T.cells) <- MRD.T.cells@meta.data$group

warnings()

saveRDS(MRD.T.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.MRD.T.cells.sorted.w.scSorter6.RDS")
##save annotated MRD.T.cells and move on to sort MRD.NT
markers4 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter6.csv", sep = ",")

Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$group
table(MRD.NT.cells@meta.data$group)

MRD.NT.cells <- NormalizeData(MRD.NT.cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
MRD.NT.cells <- FindVariableFeatures(MRD.NT.cells, selection.method = "vst", nfeatures = 3000, verbose = F)
MRD.NT.cells <- ScaleData(MRD.NT.cells)
MRD.NT.cells <- RunPCA(MRD.NT.cells, features = VariableFeatures(MRD.NT.cells))

##run scSorter to based on marker genes from CellMarker2.0
topgenes <- head(VariableFeatures(MRD.NT.cells), 3000)
expr = GetAssayData(MRD.NT.cells)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(markers4$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, markers4)
MRD.NT.cells$Predicted_Type <- rts$Pred_Type

print(table(rts$Pred_Type))

CellType <- MRD.NT.cells@meta.data$Predicted_Type
head(names(MRD.NT.cells@active.ident))
names(CellType) <- names(MRD.NT.cells@active.ident)
CellType_Subset <- CellType[CellType == "Sorted Cells"]
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$group

##these will all be seen in log files and can then be recorded in excel
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Predicted_Type
table(MRD.NT.cells@meta.data$Predicted_Type)

Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$group

warnings()

saveRDS(MRD.NT.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.MRD.NT.cells.sorted.w.scSorter6.RDS")
##save annotated MRD.NT.cells and then merge all into one object for use later

HCCp1.cells.annotated.individually <- merge(CD.cells, y = c(PreT.cells, MWD.T.cells, MRD.T.cells, MRD.NT.cells))
table(HCCp1.cells.annotated.individually@meta.data$group)
##now save RDS of the annotated data to try and visualize with UMAP downstream
saveRDS(HCCp1.cells.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.individually.w.scSorter6.RDS")
