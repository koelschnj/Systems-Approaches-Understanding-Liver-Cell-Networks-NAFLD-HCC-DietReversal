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

##use this now to just sort MWD.T and MRD.T (we got errors after PreT sorting and are troubleshooting)
markers3 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter.csv", sep = ",")

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

saveRDS(MRD.T.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.MRD.T.cells.sorted.w.scSorter.Step2_2.RDS")
##save annotated MRD.T.cells and move on to sort MWD.T

markers2 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter.csv", sep = ",")

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

saveRDS(MWD.T.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.MWD.T.cells.sorted.w.scSorter.Step2_2.RDS")
