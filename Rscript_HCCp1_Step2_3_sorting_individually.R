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

##use this now to just sort MRD.NT (we got errors after PreT sorting and are troubleshooting)
markers4 <- read.delim("~/Reverse_Diet/GeneLists/CellMarker2.0_compiled_gene_list_scSorter.csv", sep = ",")

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

saveRDS(MRD.NT.cells, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.MRD.NT.cells.sorted.w.scSorter.Step2_3.RDS")
##save annotated MRD.NT.cells and then merge all into one object for use later
##need to load in data from other sorting results...
##do NOT run code below until we get all the sorting results together, so we can be sure they all got annotated individually first
#!CD.sorted.cells <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.CD.cells.sorted.w.scSorter.RDS")
#!PreT.sorted.cells <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.CD.cells.sorted.w.scSorter.RDS")

#!HCCp1.cells.annotated.individually <- merge(CD.cells, y = c(PreT.cells, MWD.T.cells, MRD.T.cells, MRD.NT.cells))
#!table(HCCp1.cells.annotated.individually@meta.data$group)
##now save RDS of the annotated data to try and visualize with UMAP downstream
#!saveRDS(HCCp1.cells.annotated.individually, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.individually.w.scSorter.RDS")
