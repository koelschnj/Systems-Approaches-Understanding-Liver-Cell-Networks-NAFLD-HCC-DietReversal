library("Seurat")
library("png")
library("optparse")
library("ggplot2")

CD1.data <- Read10X(data.dir = "~/HCC_paper1/Data/CD/CDm1/filtered_feature_bc_matrix.tar/filtered_feature_bc_matrix")
CD2.data <- Read10X(data.dir = "~/HCC_paper1/Data/CD/CDm2/filtered_feature_bc_matrix (1).tar/filtered_feature_bc_matrix (1)")

PreT1.data <- Read10X(data.dir = "~/HCC_paper1/Data/PreT/WDm1/filtered_feature_bc_matrix (2).tar/filtered_feature_bc_matrix (2)")
PreT2.data <- Read10X(data.dir = "~/HCC_paper1/Data/PreT/WDm2/filtered_feature_bc_matrix (3).tar/filtered_feature_bc_matrix (3)")

MWD4.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleWD4/filtered_feature_bc_matrix")
MWD5.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleWD5/filtered_feature_bc_matrix")
MWD7.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleWD7/filtered_feature_bc_matrix")

MRD1.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleR1/filtered_feature_bc_matrix")
MRD2.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleR2/filtered_feature_bc_matrix")
MRD3.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleR3/filtered_feature_bc_matrix")
MRD4.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleR4/filtered_feature_bc_matrix")
MRD5.data <- Read10X(data.dir = "~/Reverse_Diet/RD_data/48wksmaleR5/filtered_feature_bc_matrix")

##make data into seurat objects for QC and Filtering
CD1 <- CreateSeuratObject(counts = CD1.data, min.cells = 3, min.features = 200, project = "CD")
CD2 <- CreateSeuratObject(counts = CD2.data, min.cells = 3, min.features = 200, project = "CD")

PreT1 <- CreateSeuratObject(counts = PreT1.data, min.cells = 3, min.features = 200, project = "PreT")
PreT2 <- CreateSeuratObject(counts = PreT2.data, min.cells = 3, min.features = 200, project = "PreT")

MWD4 <- CreateSeuratObject(counts = MWD4.data, min.cells = 3, min.features = 200, project = "MWD")
MWD5 <- CreateSeuratObject(counts = MWD5.data, min.cells = 3, min.features = 200, project = "MWD")
MWD7 <- CreateSeuratObject(counts = MWD7.data, min.cells = 3, min.features = 200, project = "MWD")

MRD1 <- CreateSeuratObject(counts = MRD1.data, min.cells = 3, min.features = 200, project = "MRD")
MRD2 <- CreateSeuratObject(counts = MRD2.data, min.cells = 3, min.features = 200, project = "MRD")
MRD3 <- CreateSeuratObject(counts = MRD3.data, min.cells = 3, min.features = 200, project = "MRD")
MRD4 <- CreateSeuratObject(counts = MRD4.data, min.cells = 3, min.features = 200, project = "MRD")
MRD5 <- CreateSeuratObject(counts = MRD5.data, min.cells = 3, min.features = 200, project = "MRD")

##first add any meta.data needed to differentiate groups in visualization/analyses/replicates/etc.
##add new column to differentiate individual replicates
CD1.label <- "CD1"
CD2.label <- "CD2"
CD1 <- AddMetaData(CD1, metadata = CD1.label, col.name = "replicate")
CD2 <- AddMetaData(CD2, metadata = CD2.label, col.name = "replicate")

PreT1.label <- "PreT1"
PreT2.label <- "PreT2"
PreT1 <- AddMetaData(PreT1, metadata = PreT1.label, col.name = "replicate")
PreT2 <- AddMetaData(PreT2, metadata = PreT2.label, col.name = "replicate")

MWD.T1.label <- "MWD.T1"
MWD.T2.label <- "MWD.T2"
MWD.T3.label <- "MWD.T3"

MWD4 <- AddMetaData(MWD4, metadata = MWD.T1.label, col.name = "replicate")
MWD5 <- AddMetaData(MWD5, metadata = MWD.T2.label, col.name = "replicate")
MWD7 <- AddMetaData(MWD7, metadata = MWD.T3.label, col.name = "replicate")

MRD1.label <- "MRD.T1"
MRD2.label <- "MRD.NT1"
MRD3.label <- "MRD.NT2"
MRD4.label <- "MRD.NT3"
MRD5.label <- "MRD.T2"

MRD1 <- AddMetaData(MRD1, metadata = MRD1.label, col.name = "replicate")
MRD2 <- AddMetaData(MRD2, metadata = MRD2.label, col.name = "replicate")
MRD3 <- AddMetaData(MRD3, metadata = MRD3.label, col.name = "replicate")
MRD4 <- AddMetaData(MRD4, metadata = MRD4.label, col.name = "replicate")
MRD5 <- AddMetaData(MRD5, metadata = MRD5.label, col.name = "replicate")

##add new column to differentiate groupings/tumor bearing animals
CD.label <- "CD"
CD1 <- AddMetaData(CD1, metadata = CD.label, col.name = "group")
CD2 <- AddMetaData(CD2, metadata = CD.label, col.name = "group")

PreT.label <- "PreT"
PreT1 <- AddMetaData(PreT1, metadata = PreT.label, col.name = "group")
PreT2 <- AddMetaData(PreT2, metadata = PreT.label, col.name = "group")


MWD.tumor.label <- "MWD.T"
MWD4 <- AddMetaData(MWD4, metadata = MWD.tumor.label, col.name = "group")
MWD5 <- AddMetaData(MWD5, metadata = MWD.tumor.label, col.name = "group")
MWD7 <- AddMetaData(MWD7, metadata = MWD.tumor.label, col.name = "group")

MRD.no.tumor.label <- "MRD.NT"
MRD2 <- AddMetaData(MRD2, metadata = MRD.no.tumor.label, col.name = "group")
MRD3 <- AddMetaData(MRD3, metadata = MRD.no.tumor.label, col.name = "group")
MRD4 <- AddMetaData(MRD4, metadata = MRD.no.tumor.label, col.name = "group")

MRD.tumor.label <- "MRD.T"
MRD1 <- AddMetaData(MRD1, metadata = MRD.tumor.label, col.name = "group")
MRD5 <- AddMetaData(MRD5, metadata = MRD.tumor.label, col.name = "group")

##then merge all samples to give last identifying labels
CD.combo <- merge(CD1, y = CD2, add.cell.ids = c("CD1", "CD2"))
PreT.combo <- merge(PreT1, y = PreT2, add.cell.ids = c("PreT1", "PreT2"))
MWD.combo <- merge(MWD4, y = c(MWD5, MWD7), add.cell.ids = c("MWD4", "MWD5", "MWD7"))
MRD.combo <- merge(MRD1, y = c(MRD2, MRD3, MRD4, MRD5), add.cell.ids = c("MRD1", "MRD2", "MRD3", "MRD4", "MRD5"))

CD.label1 <- "CD"
CD.combo <- AddMetaData(CD.combo, metadata = CD.label1, col.name = 'experimental.group')

PreT.label1 <- "PreT"
PreT.combo <- AddMetaData(PreT.combo, metadata = PreT.label1, col.name = 'experimental.group')

MWD.label <- "MWD"
MWD.combo <- AddMetaData(MWD.combo, metadata = MWD.label, col.name = 'experimental.group')

MRD.label <- "MRD"
MRD.combo <- AddMetaData(MRD.combo, metadata = MRD.label, col.name = 'experimental.group')

M.group.combo <- merge(CD.combo, y = c(PreT.combo, MWD.combo, MRD.combo))
table(M.group.combo@meta.data$experimental.group)
table(M.group.combo@meta.data$group)
table(M.group.combo@meta.data$replicate)

##perform qc and filtering steps and try to save image of QC metrics 
M.group.combo[["percent.mt"]] <- PercentageFeatureSet(M.group.combo, pattern = "^mt-")

png(filename = "~/HCC_paper1/Output_Images/QC.HCCp1.group.combo.VlnPlot.png", width = 1500, height = 1000, res = 250) 
print(VlnPlot(M.group.combo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "group"))
dev.off()

##save the Rfile following filtering and initial visualization to re-load later for steps like annotation 
##normalization and QC of data done first and then run PCA and UMAP
M.group.combo <- subset(M.group.combo, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

M.group.combo <- NormalizeData(M.group.combo, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
M.group.combo <- FindVariableFeatures(M.group.combo, selection.method = "vst", nfeatures = 3000, verbose = F)
M.group.combo <- ScaleData(M.group.combo)
M.group.combo <- RunPCA(M.group.combo, features = VariableFeatures(M.group.combo))

Idents(M.group.combo) <- M.group.combo@meta.data$group
png(filename = "~/HCC_paper1/Output_Images/HCCp1.group.combo.PCAplot.png", width = 1500, height = 1000, res = 250) 
print(DimPlot(M.group.combo, reduction = "pca", split.by = "group"))
dev.off()

M.group.combo <- FindNeighbors(M.group.combo, dims = 1:10, verbose = FALSE)
M.group.combo <- FindClusters(M.group.combo, resolution = 0.4, dims = 1:10, verbose = FALSE)
M.group.combo <- RunUMAP(M.group.combo, dims = 1:10, n.neighbors = 15L)
options(ggrepel.max.overlaps = Inf)

new.colors <- DiscretePalette(n = 25, "polychrome")

png(filename = "~/HCC_paper1/Output_Images/HCCp1.group.combo.UMAP.w.legend.png", width = 1750, height = 1000, res = 250) 
print(DimPlot(M.group.combo, cols = new.colors, reduction = "umap", pt.size = 0.5, split.by = "group", repel = TRUE) + labs(title = NULL))
dev.off()

png(filename = "~/HCC_paper1/Output_Images/HCCp1.group.combo.UMAP.wo.leg.png", width = 1750, height = 1000, res = 250) 
print(DimPlot(M.group.combo, cols = new.colors, reduction = "umap", pt.size = 0.4, split.by = "group", label = TRUE, label.size = 4, repel = TRUE, ncol = 5) + NoLegend())
dev.off()

Idents(M.group.combo) <- M.group.combo@meta.data$group
table(M.group.combo@meta.data$group)

saveRDS(M.group.combo, file = "~/HCC_paper1/Saved_Rfiles/HCCp1.group.combo.QCandFiltered.RDS")
###################################################################
##end of step1