library(Seurat)
library(CellChat)
library(future)
library(patchwork)
options(stringsAsFactors = FALSE)

HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.labels.RDS")
#!HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.combo.stromal.labels.RDS")
##use this for more updated and uniform LSEC1 and LSEC2 labels
#!HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.LSEC1and2.all.new.labels.RDS")
##use this for abbreviated labels performed with analysis
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.labels.abbreviated.RDS")
##use this for select subsets (CD4/CD8/KC/Mac) w/ abbreviated labels
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.CD4.CD8.KC.Mac.select.subsets.RDS")

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$group
table(HCCp1.cells@meta.data$group)

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.cells@meta.data$group <- factor(x = HCCp1.cells@meta.data$group, levels = new.group.levels)
table(HCCp1.cells@meta.data$group)
##need new group labels for these five to be saved to rerun CellChat
CD <- subset(HCCp1.cells, idents = "CD", invert = FALSE)
table(CD@meta.data$group)
PreT <- subset(HCCp1.cells, idents = "PreT", invert = FALSE)
table(PreT@meta.data$group)
MWD.T <- subset(HCCp1.cells, idents = "MWD.T", invert = FALSE)
table(MWD.T@meta.data$group)
MRD.T <- subset(HCCp1.cells, idents = "MRD.T", invert = FALSE)
table(MRD.T@meta.data$group)
MRD.NT <- subset(HCCp1.cells, idents = "MRD.NT", invert = FALSE)
table(MRD.NT@meta.data$group)

CD.name <- "CD"
PreT.name <- "WD.nf"
MWD.T.name <- "WD.t"
MRD.T.name <- "RD.t"
MRD.NT.name <- "RD.n"

CD <- AddMetaData(CD, metadata = CD.name, col.name = "new.group")
PreT <- AddMetaData(PreT, metadata = PreT.name, col.name = "new.group")
MWD.T <- AddMetaData(MWD.T, metadata = MWD.T.name, col.name = "new.group")
MRD.T <- AddMetaData(MRD.T, metadata = MRD.T.name, col.name = "new.group")
MRD.NT <- AddMetaData(MRD.NT, metadata = MRD.NT.name, col.name = "new.group")

all.cells <- merge(CD, y = c(PreT, MWD.T, MRD.T, MRD.NT))
table(all.cells@meta.data$new.group)
saveRDS(all.cells, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.rds")
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
##now we can load from this new labeled data for UMAPs too
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.rds")
##use this to load in abbreviated cell tyep names for updated group names as well analysis
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.abbreviated.group.names.rds")
##use this to load in select subsets (CD4/CD8/KC/Mac) and abbreviated names for updated group names as well analysis
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.select.subsets.CD4.8.KC.Mac.new.group.labels.rds")

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$new.group
table(HCCp1.cells@meta.data$new.group)

new.group.levels <- c("CD", "WD.nf", "WD.t", "RD.t", "RD.n")
HCCp1.cells@meta.data$new.group <- factor(x = HCCp1.cells@meta.data$new.group, levels = new.group.levels)
table(HCCp1.cells@meta.data$new.group)

##make UMAPs of the cells and after this combine stromal population (endo/LSEC/stromal) to make another UMAP
new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")
#!new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")
##use this for LSEC1 and LSEC2 labels
#!new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC1", "LSEC2", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")
##use this for abbreviated cell type levels
new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutro", "Mono", "Mac", "Endo", "LSEC", "Stromal", "HSC", "Fibro", "Myofibro", "Cholangio", "Hep", "Cancer")
##use this for select subsets (CD4/CD8/KC/Mac)
new.cell.type.levels <- c("Bcell", "CD4", "CD8", "DC", "NKT", "NK", "Neutro", "Mono", "KC", "Mac", "Endo", "LSEC", "Stromal", "HSC", "Fibro", "Myofibro", "Cholangio", "Hep", "Cancer")

HCCp1.cells@meta.data$Cell.Type <- factor(x = HCCp1.cells@meta.data$Cell.Type, levels = new.cell.type.levels)
table(HCCp1.cells@meta.data$Cell.Type)
##again for LSEC1 and LSEC2 labels
#!HCCp1.cells@meta.data$LSEC1and2 <- factor(x = HCCp1.cells@meta.data$LSEC1and2, levels = new.cell.type.levels)
#!table(HCCp1.cells@meta.data$LSEC1and2)
##again use this for abbreviated cell type levels
HCCp1.cells@meta.data$Abbreviated <- factor(x = HCCp1.cells@meta.data$Abbreviated, levels = new.cell.type.levels)
table(HCCp1.cells@meta.data$Abbreviated)
##again for select subsets (CD4/CD8/KC/Mac)
HCCp1.cells@meta.data$Select.subsets <- factor(x = HCCp1.cells@meta.data$Select.subsets, levels = new.cell.type.levels)
table(HCCp1.cells@meta.data$Select.subsets)

new.colors <- DiscretePalette(n = 24, "polychrome")

HCCp1.cells <- FindVariableFeatures(HCCp1.cells, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(HCCp1.cells)
HCCp1.cells <- ScaleData(HCCp1.cells, features = all.genes)

HCCp1.cells <- RunPCA(HCCp1.cells, features = VariableFeatures(object = HCCp1.cells))

DimPlot(HCCp1.cells, reduction = "pca", split.by = "replicate")

ElbowPlot(HCCp1.cells)

HCCp1.cells <- FindNeighbors(HCCp1.cells, dims = 1:15)
HCCp1.cells <- FindClusters(HCCp1.cells)

HCCp1.cells <- RunUMAP(HCCp1.cells, dims = 1:15)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Cell.Type
DimPlot(HCCp1.cells, reduction = "umap", cols = new.colors)

DimPlot(HCCp1.cells, reduction = "umap", cols = new.colors, label = TRUE) + NoLegend()

new.group.levels <- c("CD", "WD.nf", "WD.t", "RD.t", "RD.n")
HCCp1.cells@meta.data$new.group <- factor(x = HCCp1.cells@meta.data$new.group, levels = new.group.levels)
table(HCCp1.cells@meta.data$new.group)

DimPlot(HCCp1.cells, reduction = "umap", cols = new.colors, label = TRUE, split.by = "new.group") + NoLegend()

DimPlot(HCCp1.cells, reduction = "umap", cols = new.colors, split.by = "new.group") + NoLegend()

HCCp1.ICs <- subset(HCCp1.cells, idents = c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Macrophage", "Monocyte"), invert = FALSE)
DimPlot(HCCp1.ICs, reduction = "umap", cols = new.colors)
DimPlot(HCCp1.ICs, reduction = "umap", cols = new.colors, label = TRUE) + NoLegend()
DimPlot(HCCp1.ICs, reduction = "umap", cols = new.colors, label = TRUE, split.by = "new.group") + NoLegend()
DimPlot(HCCp1.ICs, reduction = "umap", cols = new.colors, split.by = "new.group") + NoLegend()

HCCp1.nonICs <- subset(HCCp1.cells, idents = c("Cancer", "Cholangiocyte", "Endothelial", "Fibroblast", "HSC", "Hepatocyte", "LSEC", "Myofibroblast", "Stromal"), invert = FALSE)
#!HCCp1.nonICs <- subset(HCCp1.cells, idents = c("Cancer", "Cholangiocyte", "Fibroblast", "HSC", "Hepatocyte", "Myofibroblast", "Stromal"), invert = FALSE)

DimPlot(HCCp1.nonICs, reduction = "umap", cols = new.colors)
DimPlot(HCCp1.nonICs, reduction = "umap", cols = new.colors, label = TRUE) + NoLegend()
DimPlot(HCCp1.nonICs, reduction = "umap", cols = new.colors, label = TRUE, split.by = "group") + NoLegend()
DimPlot(HCCp1.nonICs, reduction = "umap", cols = new.colors, split.by = "new.group") + NoLegend()

DoHeatmap(HCCp1.cells, features = c("Ighm", "Cd79b", "Cd3e", "Cd247", "Trbc2", "Ccl5", "Klra4", "Gzma", "Ccr9", "H2-Aa", "Il1b", "S100a8", "S100a4", "Ly6g", "Itgam", "Adgre1", "Bmp2", "Dpp4", "Lyve1", "Stab2", "Gpc6", "Prkg1", "Col3a1", "Reln", "Pdgfra", "Col15a1", "Acta2", "Fbln1", "Spp1", "Krt19", "Alb", "Apoa1", "Hnf4a", "Afp"), group.by = "Cell.Type", size = 2)
DoHeatmap(HCCp1.cells, features = c("Cd19", "Cd22", "Cd79a", "Cd79b", "Ebf1", "Ighd", "Ighm", "Ms4a1", "Fcmr",
                                    "Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Thy1", "Il7r", "Prkcq", "Lef1", "Cd247",
                                    "Cst3", "Cd80", "Cd83", "H2-Ab1", "Ly6d", "Siglech", "Zbtb46", "Rnase6", "Itgax", "Irf8",
                                    "Il15ra", "Tox", "Cd36", "Ccl5",
                                    "Itga1", "Eomes", "Gzma", "Gzmb", "Klf2", "Klra4", "Klra9", "Klrb1c", "Klre1", "Lgals1",
                                    "S100a8", "S100a9", "Retnlg", "Mmp9", "Il1b", "Il1r2", "Csf3r", "Cxcr2",
                                    "Cd209a", "Ccr2", "Cd74", "Chil3", "Fn1", "Ly6c2", "S100a4", "Itgam",
                                    "Adgre1", "Apoe", "C1qa", "C1qb", "Ccnb2", "Cd68", "Clec4f", "Csf1r", "Ctss", "Lyz2", "Sdc3", "Stmn1", "Top2a", "Ube2c",
                                    "Dpp4", "Bmp2", "C1qtnf1", "F8", "Il1a", "Mmrn2", "Oit3", "Pcdh12", "Ushbp1",
                                    "Stab2", "Lyve1",
                                    "Prkg1", "Gpc6",
                                    "Acta2", "Bgn", "Col14a1", "Col3a1", "Dcn", "Ecm1", "Lrat", "Pdgfra", "Pdgfrb", "Reln", "Tcf21", "Tmem56", "Tpm2",
                                    "Cd34", "Thy1", "Clec3b", "Col15a1", "Entpd2", "Fbln2", "Gli1",
                                    "Nt5e", "Des", "Dpt", "Eln", "Fbln1",
                                    "Spp1", "Cftr", "Epcam", "Krt19", "Krt7", "Muc1", "Pkhd1", "Sox9", "St14",
                                    "Alb", "Acly", "Apoa1", "Asl", "Ass1", "Cyp2e1", "Cyp2f2", "G6pc", "Glul", "Mup3", "Pck1",
                                    "Hnf4a", "Afp", "Gpc3"), group.by = "Cell.Type", size = 3)
##IC heatmap alone 
DoHeatmap(HCCp1.ICs, features = c("Cd19", "Cd22", "Cd79a", "Cd79b", "Ebf1", "Ighd", "Ighm", "Ms4a1", "Fcmr",
                                    "Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Thy1", "Il7r", "Prkcq", "Lef1", "Cd247",
                                    "Cst3", "Cd80", "Cd83", "H2-Ab1", "Ly6d", "Siglech", "Zbtb46", "Rnase6", "Itgax", "Irf8",
                                    "Il15ra", "Tox", "Cd36", "Ccl5",
                                    "Itga1", "Eomes", "Gzma", "Gzmb", "Klf2", "Klra4", "Klra9", "Klrb1c", "Klre1", "Lgals1",
                                    "S100a8", "S100a9", "Retnlg", "Mmp9", "Il1b", "Il1r2", "Csf3r", "Cxcr2",
                                    "Cd209a", "Ccr2", "Cd74", "Chil3", "Fn1", "Ly6c2", "S100a4", "Itgam",
                                    "Adgre1", "Apoe", "C1qa", "C1qb", "Ccnb2", "Cd68", "Clec4f", "Csf1r", "Ctss", "Lyz2", "Sdc3", "Stmn1", "Top2a", "Ube2c"), size = 3)

########################################################
########################################################
##make heatmap of macrophage phagocytosis genes that may relate to entropy resolution information
Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Cell.Type
HCCp1.Mac <- subset(HCCp1.cells, idents = "Macrophage", invert = FALSE)
table(HCCp1.Mac@meta.data$group)

new.group.levels <- c("CD", "WD.nf", "WD.t", "RD.t", "RD.n")
HCCp1.Mac@meta.data$new.group <- factor(x = HCCp1.Mac@meta.data$new.group, levels = new.group.levels)
table(HCCp1.Mac@meta.data$new.group)

all.genes <- rownames(HCCp1.Mac)
HCCp1.Mac <- ScaleData(HCCp1.Mac, features = all.genes)

DoHeatmap(HCCp1.Mac, features = c("Msr1", "Cd36", "Cd163", "Olr1", "Stab1", "Stab2", "Cxcl16"), group.by = "new.group", size = 3)
##quantify some numbers for Mac phagocytosis related genes
Mac.Msr1 <- subset(HCCp1.Mac, subset = Msr1 > 0)
table(Mac.Msr1@meta.data$replicate)
Mac.Cd36 <- subset(HCCp1.Mac, subset = Cd36 > 0)
table(Mac.Cd36@meta.data$replicate)
Mac.Cd163 <- subset(HCCp1.Mac, subset = Cd163 > 0)
table(Mac.Cd163@meta.data$replicate)
Mac.Olr1 <- subset(HCCp1.Mac, subset = Olr1 > 0)
table(Mac.Olr1@meta.data$replicate)
Mac.Cxcl16 <- subset(HCCp1.Mac, subset = Cxcl16 > 0)
table(Mac.Cxcl16@meta.data$replicate)
Mac.Stab1 <- subset(HCCp1.Mac, subset = Stab1 > 0)
table(Mac.Stab1@meta.data$replicate)
Mac.Stab2 <- subset(HCCp1.Mac, subset = Stab2 > 0)
table(Mac.Stab2@meta.data$replicate)

AverageExpression(Hep, feature = "Hnf4a", group.by = "replicate")
AverageExpression(Hep, feature = "Afp", group.by = "replicate")
AverageExpression(Hep, feature = "Gpc3", group.by = "replicate")
###########################################################################
##nonIC heatmap alone
DoHeatmap(HCCp1.nonICs, features = c("Dpp4", "Bmp2", "C1qtnf1", "F8", "Il1a", "Mmrn2", "Oit3", "Pcdh12", "Ushbp1",
                                    "Stab2", "Lyve1",
                                    "Prkg1", "Gpc6",
                                    "Acta2", "Bgn", "Col14a1", "Col3a1", "Dcn", "Ecm1", "Lrat", "Pdgfra", "Pdgfrb", "Reln", "Tcf21", "Tmem56", "Tpm2",
                                    "Cd34", "Thy1", "Clec3b", "Col15a1", "Entpd2", "Fbln2", "Gli1",
                                    "Nt5e", "Des", "Dpt", "Eln", "Fbln1",
                                    "Spp1", "Cftr", "Epcam", "Krt19", "Krt7", "Muc1", "Pkhd1", "Sox9", "St14",
                                    "Alb", "Acly", "Apoa1", "Asl", "Ass1", "Cyp2e1", "Cyp2f2", "G6pc", "Glul", "Mup3", "Pck1",
                                    "Hnf4a", "Afp", "Gpc3"), group.by = "Cell.Type", size = 3)
##heatmap for endo/LSEC/stromal only
HCCp1.stromal <- subset(HCCp1.cells, idents = c("Endothelial", "LSEC", "Stromal"), invert = FALSE)

DoHeatmap(HCCp1.stromal, features = c("Dpp4", "Bmp2", "C1qtnf1", "F8", "Il1a", "Mmrn2", "Oit3", "Pcdh12", "Ushbp1",
                                     "Stab2", "Lyve1",
                                     "Prkg1", "Gpc6"), size = 3)
##same heatmap but now also including HSCs
Idents(HCCp1.cells) <- HCCp1.cells@meta.data$LSEC1and2

HCCp1.stromal <- subset(HCCp1.cells, idents = c("Endothelial", "HSC", "LSEC1", "LSEC2", "Stromal"), invert = FALSE)

DoHeatmap(HCCp1.stromal, features = c("Dpp4", "Bmp2", "C1qtnf1", "F8", "Il1a", "Mmrn2", "Oit3", "Pcdh12", "Ushbp1",
                                      "Acta2", "Bgn", "Col14a1", "Col3a1", "Dcn", "Ecm1", "Lrat", "Pdgfra", "Pdgfrb", "Reln", "Tcf21", "Tmem56", "Tpm2",
                                      "Stab2", "Lyve1",
                                      "Prkg1", "Gpc6"), size = 3)
##new heatmap for other nonIC species when using HSCs above 
HCCp1.nonICs <- subset(HCCp1.cells, idents = c("Cancer", "Cholangiocyte", "Fibroblast", "Hepatocyte", "Myofibroblast"), invert = FALSE)

DoHeatmap(HCCp1.nonICs, features = c("Cd34", "Clec3b", "Entpd2", "Fbln2", "Lrat", "Pdgfra",
                                     "Nt5e", "Des", "Dpt", "Eln", "Fbln1",
                                     "Spp1", "Cftr", "Epcam", "Krt19", "Krt7", "Muc1", "Pkhd1", "Sox9", "St14",
                                     "Alb", "Acly", "Apoa1", "Asl", "Ass1", "Cyp2e1", "Cyp2f2", "G6pc", "Glul", "Mup3", "Pck1",
                                     "Hnf4a", "Afp", "Gpc3"), size = 3)
##make heatmap to differentiate LSEC1 and LSEC2 functional identities based on identified markers
HCCp1.LSEC <- subset(HCCp1.cells, idents = c("LSEC1", "LSEC2"), invert = FALSE)

DoHeatmap(HCCp1.LSEC, features = c("Car3", "Fabp1", "Cyp3a25", "Scd1", "Sntb1", "Lama4", "Tek", "Slco1b2",
                                   "Lcn2", "Orm2", "Cxcl13", "Mt1", "Saa1", "Apcs", "Fgl1", "Saa2"), size = 3)

##use this section to merge stromal,LSEC, and endo populations as a collective stromal like cell signal for analysis
#!Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Cell.Type
StackedVlnPlot(HCCp1.cells, features = c("Cd19", "Cd22", "Cd79a", "Cd79b", "Ebf1", "Ighd", "Ighm", "Ms4a1", "Fcmr",
                                         "Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Thy1", "Il7r", "Prkcq", "Lef1", "Cd247",
                                         "Cst3", "Cd80", "Cd83", "H2-Ab1", "Ly6d", "Siglech", "Zbtb46", "Rnase6", "Itgax", "Irf8",
                                         "Fas", "Nkg7", "Trac", "Trbc2", "Il15ra", "Tox", "Cd36", "Ccl5",
                                         "Borcs7", "Itga1", "Eomes", "Gzma", "Gzmb", "Klf2", "Klra4", "Klra7", "Klra8", "Klra9", "Klrb1c", "Klre1", "Klrg1", "Lgals1",
                                         "S100a8", "S100a9", "Retnlg", "Mmp9", "Il1b", "Il1r2", "Csf3r", "Cxcr2",
                                         "Cd209a", "Ccr2", "Cd74", "Chil3", "Fn1", "Ly6c1", "Ly6c2", "Ly6g", "Mycl", "S100a4", "Itgam",
                                         "Adgre1", "Apoe", "C1qa", "C1qb", "Ccnb2", "Cd68", "Clec4f", "Csf1r", "Ctss", "Lyz2", "Sdc3", "Stmn1", "Top2a", "Ube2c",
                                         "Dpp4", "Bmp2", "C1qtnf1", "F8", "Il1a", "Mmrn2", "Oit3", "Pcdh12", "Ushbp1",
                                         "Stab2", "Lyve1",
                                         "Prkg1", "Gpc6",
                                         "Acta2", "Bgn", "Col14a1", "Col3a1", "Dcn", "Ecm1", "Lrat", "Pdgfra", "Pdgfrb", "Reln", "Tcf21", "Tmem56", "Tpm2",
                                         "Cd34", "Thy1", "Clec3b", "Col15a1", "Entpd2", "Fbln2", "Gli1",
                                         "Nt5e", "Des", "Dpt", "Eln", "Fbln1",
                                         "Spp1", "Cftr", "Epcam", "Krt19", "Krt7", "Muc1", "Pkhd1", "Sox9", "St14",
                                         "Alb", "Acly", "Apoa1", "Asl", "Ass1", "Cyp2e1", "Cyp2f2", "G6pc", "Glul", "Mup3", "Pck1",
                                         "Hnf4a", "Afp", "Gpc3"), group.by = "Cell.Type")

##########################################################################################################################################################################
##try to look at specific pathway expression of components in each group like MHCI/II
##because MHC-I and MHC-II are not detected even at 5%, need to manually plot expression across groups for all genes in non-cellchat processed data
all.genes <- rownames(CD)
CD <- ScaleData(CD, features = all.genes)
all.genes <- rownames(WD.nf)
WD.nf <- ScaleData(WD.nf, features = all.genes)
all.genes <- rownames(MWD.T)
MWD.T <- ScaleData(MWD.T, features = all.genes)
all.genes <- rownames(MRD.T)
MRD.T <- ScaleData(MRD.T, features = all.genes)
all.genes <- rownames(MRD.NT)
MRD.NT <- ScaleData(MRD.NT, features = all.genes)

Idents(CD) <- CD@meta.data$Select.subsets
Idents(WD.nf) <- WD.nf@meta.data$Select.subsets
Idents(MWD.T) <- MWD.T@meta.data$Select.subsets
Idents(MRD.T) <- MRD.T@meta.data$Select.subsets
Idents(MRD.NT) <- MRD.NT@meta.data$Select.subsets

##MHC-I ligands
StackedVlnPlot(CD, features = c("H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q10", "H2-T3", "H2-T22", "H2-T23", "H2-T24", "H2-M2", "H2-M3", "H2-M5"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(WD.nf, features = c("H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q10", "H2-T3", "H2-T22", "H2-T23", "H2-T24", "H2-M2", "H2-M3", "H2-M5"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MWD.T, features = c("H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q10", "H2-T3", "H2-T22", "H2-T23", "H2-T24", "H2-M2", "H2-M3", "H2-M5"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.T, features = c("H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q10", "H2-T3", "H2-T22", "H2-T23", "H2-T24", "H2-M2", "H2-M3", "H2-M5"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.NT, features = c("H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q10", "H2-T3", "H2-T22", "H2-T23", "H2-T24", "H2-M2", "H2-M3", "H2-M5"), group.by = "Select.subsets", color.use = new.colors)

##MHC-I receptors
StackedVlnPlot(CD, features = c("Klrd1", "Klrc1", "Klrc2", "Cd1d1", "Cd8a", "Cd8b1", "Kir3dl1", "Klrk1"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(WD.nf, features = c("Klrd1", "Klrc1", "Klrc2", "Cd1d1", "Cd8a", "Cd8b1", "Kir3dl1", "Klrk1"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MWD.T, features = c("Klrd1", "Klrc1", "Klrc2", "Cd1d1", "Cd8a", "Cd8b1", "Kir3dl1", "Klrk1"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.T, features = c("Klrd1", "Klrc1", "Klrc2", "Cd1d1", "Cd8a", "Cd8b1", "Kir3dl1", "Klrk1"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.NT, features = c("Klrd1", "Klrc1", "Klrc2", "Cd1d1", "Cd8a", "Cd8b1", "Kir3dl1", "Klrk1"), group.by = "Select.subsets", color.use = new.colors)

##MHC-II (much fewer so put receptor Cd4 in with ligands)
StackedVlnPlot(CD, features = c("Cd4", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Eb2", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Oa", "H2-Ob"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(WD.nf, features = c("Cd4", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Eb2", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Oa", "H2-Ob"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MWD.T, features = c("Cd4", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Eb2", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Oa", "H2-Ob"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.T, features = c("Cd4", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Eb2", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Oa", "H2-Ob"), group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.NT, features = c("Cd4", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Eb2", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Oa", "H2-Ob"), group.by = "Select.subsets", color.use = new.colors)

StackedVlnPlot(CD, features = "Cd1d1", group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(WD.nf, features = "Cd1d1", group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MWD.T, features = "Cd1d1", group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.T, features = "Cd1d1", group.by = "Select.subsets", color.use = new.colors)
StackedVlnPlot(MRD.NT, features = "Cd1d1", group.by = "Select.subsets", color.use = new.colors)

########################################################
########################################################
##try to select LSEC population alone on a UMAP to separate it into two groups
LSEC <- subset(HCCp1.cells, idents = "LSEC", invert = FALSE)

DimPlot(LSEC, reduction = "umap", pt.size = 0.4, cols = new.colors, repel = TRUE) 
DimPlot(LSEC, reduction = "umap", pt.size = 0.4, cols = c("red", "blue"), repel = TRUE, split.by = "group") 

plot <- DimPlot(LSEC, reduction = "umap", pt.size = 0.4, cols = new.colors, group.by = "Cell.Type", repel = TRUE) 
select.cells <- CellSelector(plot = plot)
head(select.cells)
tail(select.cells)

Idents(LSEC, cells = select.cells) <- "Selected.cells"
LSEC.selected <- subset(LSEC, cells = select.cells)

saveRDS(LSEC.selected, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/LSEC2.saved.rds")
# # # # # # # # # # # # ## # # # # # # # # # # # # # # #
LSEC1 <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/LSEC1.saved.rds")
LSEC2 <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/LSEC2.saved.rds")
table(LSEC1@meta.data$group)
table(LSEC2@meta.data$group)
table(LSEC1@meta.data$replicate)
table(LSEC2@meta.data$replicate)

LSEC1.label <- "LSEC1"
LSEC2.label <- "LSEC2"

new.LSEC1 <- AddMetaData(LSEC1, metadata = LSEC1.label, col.name = 'Cell.Type')
new.LSEC2 <- AddMetaData(LSEC2, metadata = LSEC2.label, col.name = 'Cell.Type')

new.LSEC.pop <- merge(new.LSEC1, y = new.LSEC2)
table(new.LSEC.pop@meta.data$Cell.Type)

saveRDS(new.LSEC.pop, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/new.LSEC.pop.1and2.saved.rds")

LSEC <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/new.LSEC.pop.1and2.saved.rds")
DimPlot(LSEC, reduction = "umap", pt.size = 0.4, cols = c("red", "blue"), repel = TRUE) 
DimPlot(LSEC, reduction = "umap", pt.size = 0.4, cols = c("red", "blue"), repel = TRUE, split.by = "group") 

LSEC1.markers <- FindMarkers(HCCp1.cells, ident.1 = "LSEC1", ident.2 = "LSEC2", only.pos = TRUE)
LSEC2.markers <- FindMarkers(HCCp1.cells, ident.1 = "LSEC2", ident.2 = "LSEC1", only.pos = TRUE)

write.csv(LSEC1.markers, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2 LSEC1 and LSEC2/UMAPs/LSEC1.markersvLSEC2.csv")
write.csv(LSEC2.markers, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2 LSEC1 and LSEC2/UMAPs/LSEC2.markersvLSEC1.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###try to do something similar with macrophages
Mac <- subset(HCCp1.cells, idents = "Macrophage", invert = FALSE)

DimPlot(Mac, reduction = "umap", pt.size = 0.4, cols = c("red", "blue"), repel = TRUE, split.by = "group") 
DimPlot(Mac, reduction = "umap", pt.size = 0.4, cols = c("red", "blue"), repel = TRUE, split.by = "group") + NoLegend()

DimPlot(Mac, reduction = "umap", pt.size = 0.4, cols = c("red", "blue"), repel = TRUE) + NoLegend()

########################################################
########################################################
##use this portion for another CellChat analysis with LSEC1 and LSEC2 populations after loading in cells at top
LSEC1.LSEC2 <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/new.LSEC.pop.1and2.saved.rds")

table(LSEC1.LSEC2@meta.data$Cell.Type)
table(HCCp1.cells@meta.data$Cell.Type)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Cell.Type
HCCp1.cells.noLSEC <- subset(HCCp1.cells, idents = c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer"), invert = FALSE)
table(HCCp1.cells.noLSEC@meta.data$Cell.Type)

HCCp1.cells.newLSEC <- merge(HCCp1.cells.noLSEC, y = LSEC1.LSEC2)
table(HCCp1.cells.newLSEC@meta.data$Cell.Type)
##now that we have the LSEC1/2 in this object, we can just adjust code below to pull cells from correct object

#!all.other.cells <- subset(HCCp1.cells, idents = c("Bcell", "Tcell", "NK", "NKT", "DC", "Macrophage", "Monocyte", "Neutrophil", "Cancer", "Cholangiocyte", "Fibroblast", "Myofibroblast", "HSC", "Hepatocyte"), invert = FALSE)
#!table(all.other.cells@meta.data$Cell.Type)
#!all.stromal <- subset(HCCp1.cells, idents = c("LSEC", "Endothelial", "Stromal"), invert = FALSE)
#!table(all.stromal@meta.data$Cell.Type)

#!all.stromal.label <- "Stromal"

#!new.all.stromal <- AddMetaData(all.stromal, metadata = all.stromal.label, col.name = "Cell.Type")

#!merged.cells <- merge(all.other.cells, y = new.all.stromal)
#!table(merged.cells@meta.data$Cell.Type)
#!saveRDS(merged.cells, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.combo.stromal.labels.RDS")
##after saving this updated object with labels reload it for analysis
##re-loading the object will be done at the top so we can also make new order levels below

#!new.cell.type.levels <- c("Bcell", "Tcell", "NK", "NKT", "DC", "Macrophage", "Monocyte", "Neutrophil", "Cancer", "Cholangiocyte", "Endothelial", "LSEC", "Fibroblast", "Myofibroblast", "HSC", "Hepatocyte", "Stromal")
#!new.cell.type.levels <- c("Bcell", "Tcell", "NK", "NKT", "DC", "Macrophage", "Monocyte", "Neutrophil", "Cancer", "Cholangiocyte", "Fibroblast", "Myofibroblast", "HSC", "Hepatocyte", "Stromal")
new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")
#!new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")
#!new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC1", "LSEC2", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$new.group
CD.cells <- subset(HCCp1.cells, idents = "CD", invert = FALSE)
table(CD.cells@meta.data$new.group)
CD.cells@meta.data$Select.subsets <- factor(x = CD.cells@meta.data$Select.subsets, levels = new.cell.type.levels)
table(CD.cells@meta.data$Select.subsets)
Idents(CD.cells) <- CD.cells@meta.data$Select.subsets

PreT.cells <- subset(HCCp1.cells, idents = "WD.nf", invert = FALSE)
table(PreT.cells@meta.data$new.group)
PreT.cells@meta.data$Select.subsets <- factor(x = PreT.cells@meta.data$Select.subsets, levels = new.cell.type.levels)
table(PreT.cells@meta.data$Select.subsets)
Idents(PreT.cells) <- PreT.cells@meta.data$Select.subsets

MWD.T.cells <- subset(HCCp1.cells, idents = "WD.t", invert = FALSE)
table(MWD.T.cells@meta.data$new.group)
MWD.T.cells@meta.data$Select.subsets <- factor(x = MWD.T.cells@meta.data$Select.subsets, levels = new.cell.type.levels)
table(MWD.T.cells@meta.data$Select.subsets)
Idents(MWD.T.cells) <- MWD.T.cells@meta.data$Select.subsets

MRD.T.cells <- subset(HCCp1.cells, idents = "RD.t", invert = FALSE)
table(MRD.T.cells@meta.data$new.group)
MRD.T.cells@meta.data$Select.subsets <- factor(x = MRD.T.cells@meta.data$Select.subsets, levels = new.cell.type.levels)
table(MRD.T.cells@meta.data$Select.subsets)
Idents(MRD.T.cells) <- MRD.T.cells@meta.data$Select.subsets

MRD.NT.cells <- subset(HCCp1.cells, idents = "RD.n", invert = FALSE)
table(MRD.NT.cells@meta.data$new.group)
MRD.NT.cells@meta.data$Select.subsets <- factor(x = MRD.NT.cells@meta.data$Select.subsets, levels = new.cell.type.levels)
table(MRD.NT.cells@meta.data$Select.subsets)
Idents(MRD.NT.cells) <- MRD.NT.cells@meta.data$Select.subsets

cellchat.CD <- createCellChat(object = CD.cells, group.by = "Select.subsets", assay = "RNA")
cellchat.PreT <- createCellChat(object = PreT.cells, group.by = "Select.subsets", assay = "RNA")
cellchat.MWD.T <- createCellChat(object = MWD.T.cells, group.by = "Select.subsets", assay = "RNA")
cellchat.MRD.T <- createCellChat(object = MRD.T.cells, group.by = "Select.subsets", assay = "RNA")
cellchat.MRD.NT <- createCellChat(object = MRD.NT.cells, group.by = "Select.subsets", assay = "RNA")

##after adjusting large section of code above use this for LSEC1 and LSEC2
#!cellchat.CD <- createCellChat(object = CD.cells, group.by = "LSEC1and2", assay = "RNA")
#!cellchat.PreT <- createCellChat(object = PreT.cells, group.by = "LSEC1and2", assay = "RNA")
#!cellchat.MWD.T <- createCellChat(object = MWD.T.cells, group.by = "LSEC1and2", assay = "RNA")
#!cellchat.MRD.T <- createCellChat(object = MRD.T.cells, group.by = "LSEC1and2", assay = "RNA")
#!cellchat.MRD.NT <- createCellChat(object = MRD.NT.cells, group.by = "LSEC1and2", assay = "RNA")

##load in reference databases and run cellchat on each condition
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Adhesion", "ECM-Receptor", "Non-protein Signaling"))

cellchat.CD@DB <- CellChatDB.use
cellchat.PreT@DB <- CellChatDB.use
cellchat.MWD.T@DB <- CellChatDB.use
cellchat.MRD.T@DB <- CellChatDB.use
cellchat.MRD.NT@DB <- CellChatDB.use

##start with CD 
cellchat.CD <- subsetData(cellchat.CD)
future::plan("multisession", workers = 2)

cellchat.CD <- identifyOverExpressedGenes(cellchat.CD)
cellchat.CD <- identifyOverExpressedInteractions(cellchat.CD)

unique(cellchat.CD@idents)

cellchat.CD@idents = droplevels(cellchat.CD@idents, exclude = setdiff(levels(cellchat.CD@idents),unique(cellchat.CD@idents)))

cellchat.CD <- computeCommunProb(cellchat.CD)
#!cellchat.CD <- computeCommunProb(cellchat.CD, type = "truncatedMean", trim = 0.05)
#!cellchat.CD <- computeCommunProb(cellchat.CD, type = "truncatedMean", trim = 0.75)

cellchat.CD <- filterCommunication(cellchat.CD, min.cells = 2)

cellchat.CD <- computeCommunProbPathway(cellchat.CD)

cellchat.CD <- aggregateNet(cellchat.CD)

groupSize <- as.numeric(table(cellchat.CD@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat.CD@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.CD@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat.CD@netP$pathways
pathways.show <- cellchat.CD@netP$pathways

cellchat.CD <- netAnalysis_computeCentrality(cellchat.CD, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2
##reordering
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                    "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                    "SEMA3", "AGRN", "HSPG", "IL1", "CXCL", "KIT", "RELN", "TENASCIN", "VWF"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                    "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                    "SEMA3", "AGRN", "HSPG", "IL1", "CXCL", "KIT", "RELN", "TENASCIN", "VWF"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##default heatmaps removing signaling detected in 75% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("BMP", "VTN", "VEGF", "DHEA", "Cholesterol", "NRG", "GALECTIN", "SPP1", "PROS", "HGF", "FGF", "TENASCIN", "PDGF", "IL1", "Androsterone", "RELN", "27HC", "DHEAS", "HSPG", "Testosterone", "AGRN", "CXCL", "ANGPTL", "CSF", "KIT", "VWF", "AGT"), pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("BMP", "VTN", "VEGF", "DHEA", "Cholesterol", "NRG", "GALECTIN", "SPP1", "PROS", "HGF", "FGF", "TENASCIN", "PDGF", "IL1", "Androsterone", "RELN", "27HC", "DHEAS", "HSPG", "Testosterone", "AGRN", "CXCL", "ANGPTL", "CSF", "KIT", "VWF", "AGT"), pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2

##stromal combined analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                    "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                    "ANGPT", "SEMA3", "AGRN", "HSPG", "IL1", "CXCL", "RELN", "TENASCIN", "VWF"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("AGT","ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                    "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                    "ANGPT", "SEMA3", "AGRN", "HSPG", "IL1", "CXCL", "RELN", "TENASCIN", "VWF"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##LSEC1 and LSEC2 analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                    "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                    "SEMA3", "AGRN", "HSPG", "IL1", "CXCL", "KIT", "RELN", "TENASCIN", "VWF"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                    "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                    "SEMA3", "AGRN", "HSPG", "IL1", "CXCL", "KIT", "RELN", "TENASCIN", "VWF"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% analysis 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "outgoing", width = 5, height = 18, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "incoming", width = 5, height = 18, font.size = 8, font.size.title = 9)
ht1 + ht2

##immune related in 5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("APRIL", "CD40", "EDN", "FASLG", "IL16", "IL2", "TRAIL"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("APRIL", "CD40", "EDN", "FASLG", "IL16", "IL2", "TRAIL"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (all even when not detected)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (only those detected in each group)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("IL1", "TGFb", "IGF", "CypA", "FLT3", "VISFATIN", "EPO", "27HC", "TRAIL", "FASLG", "IL2", "CD40", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "MIF", "AGT"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("IL1", "TGFb", "IGF", "CypA", "FLT3", "VISFATIN", "EPO", "27HC", "TRAIL", "FASLG", "IL2", "CD40", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "MIF", "AGT"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("GALECTIN", "NRG", "COMPLEMENT", "SPP1", "Glutamate", "TGFb", "IL1", "FLT3", "FASLG", "IL2", "CD40", "TRAIL"), pattern = "outgoing", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("GALECTIN", "NRG", "COMPLEMENT", "SPP1", "Glutamate", "TGFb", "IL1", "FLT3", "FASLG", "IL2", "CD40", "TRAIL"), pattern = "incoming", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% filtered pathways based on those detected in 25 and 75% analyses
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("COMPLEMENT", "Glutamate", "TGFb", "EGF", "SLIT", "CCL", "WNT", "GDF", "VISFATIN", "CypA", "CHEMERIN", "GRN", "THBS", "SerotoninDopamin", "EDN", "PERIOSTIN", "PROCR", "IL2", "FLT3", "MIF", "IGFBP", "ncWNT", "NT", "IL16", "PTPR", "Prostaglandin", "DHT", "2-AG", "BMP10", "NGF", "FASLG", "ACTIVIN", "NTS", "Histamine", "Desmosterol", "Adenosine", "CD40", "TRAIL", "EPO", "APRIL", "ANNEXIN"), pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, signaling = c("COMPLEMENT", "Glutamate", "TGFb", "EGF", "SLIT", "CCL", "WNT", "GDF", "VISFATIN", "CypA", "CHEMERIN", "GRN", "THBS", "SerotoninDopamin", "EDN", "PERIOSTIN", "PROCR", "IL2", "FLT3", "MIF", "IGFBP", "ncWNT", "NT", "IL16", "PTPR", "Prostaglandin", "DHT", "2-AG", "BMP10", "NGF", "FASLG", "ACTIVIN", "NTS", "Histamine", "Desmosterol", "Adenosine", "CD40", "TRAIL", "EPO", "APRIL", "ANNEXIN"), pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

##75% analysis 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "outgoing", width = 5, height = 4, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.CD, pattern = "incoming", width = 5, height = 4, font.size = 8, font.size.title = 9)
ht1 + ht2

CD.df.net <- subsetCommunication(cellchat.CD, slot.name = "net")
write.csv(CD.df.net, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/CD.scS5.v2.25pt.CD4CD8.df.net.default.csv")

saveRDS(cellchat.CD, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.CD.scS5.25pt.CD4CD8.default.rds")
##after saving file for later run on PreT
cellchat.PreT <- subsetData(cellchat.PreT)
future::plan("multisession", workers = 2)

cellchat.PreT <- identifyOverExpressedGenes(cellchat.PreT)
cellchat.PreT <- identifyOverExpressedInteractions(cellchat.PreT)

unique(cellchat.PreT@idents)

cellchat.PreT@idents = droplevels(cellchat.PreT@idents, exclude = setdiff(levels(cellchat.PreT@idents),unique(cellchat.PreT@idents)))

cellchat.PreT <- computeCommunProb(cellchat.PreT)
#!cellchat.PreT <- computeCommunProb(cellchat.PreT, type = "truncatedMean", trim = 0.05)
#!cellchat.PreT <- computeCommunProb(cellchat.PreT, type = "truncatedMean", trim = 0.75)

cellchat.PreT <- filterCommunication(cellchat.PreT, min.cells = 2)

cellchat.PreT <- computeCommunProbPathway(cellchat.PreT)

cellchat.PreT <- aggregateNet(cellchat.PreT)

groupSize <- as.numeric(table(cellchat.PreT@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat.PreT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.PreT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat.PreT@netP$pathways
pathways.show <- cellchat.PreT@netP$pathways

cellchat.PreT <- netAnalysis_computeCentrality(cellchat.PreT, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "outgoing", width = 5, height = 11, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "incoming", width = 5, height = 11, font.size = 8, font.size.title = 9)
ht1 + ht2
##reordering
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                      "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                      "CCL", "SEMA3", "AGRN", "HSPG", "IL1", "COMPLEMENT", "TGFB", "EGF", "GDF", "GRN", "KIT", "NT", "RELN", "SLIT", "VISFATIN", "VWF", "ncWNT", "PERIOSTIN", "PTN"), pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                      "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                      "CCL", "SEMA3", "AGRN", "HSPG", "IL1", "COMPLEMENT", "TGFB", "EGF", "GDF", "GRN", "KIT", "NT", "RELN", "SLIT", "VISFATIN", "VWF", "ncWNT", "PERIOSTIN", "PTN"), pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

##default heatmaps filtered based on 75% results
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("COMPLEMENT", "DHEA", "SEMA3", "NRG", "IL1", "SLIT", "GAS", "TGFb", "ANGPT", "CSF", "PTN", "HSPG", "GRN", "RELN", "AGT", "SerotoninDopamin", "AGRN", "GDF", "PERIOSTIN", "VISFATIN", "Glutamate", "KIT", "VWF", "NT", "ncWNT"), pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("COMPLEMENT", "DHEA", "SEMA3", "NRG", "IL1", "SLIT", "GAS", "TGFb", "ANGPT", "CSF", "PTN", "HSPG", "GRN", "RELN", "AGT", "SerotoninDopamin", "AGRN", "GDF", "PERIOSTIN", "VISFATIN", "Glutamate", "KIT", "VWF", "NT", "ncWNT"), pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2

##for stromal combined analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                      "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                      "ANGPT", "SEMA3", "CCL", "AGRN", "COMPLEMENT", "HSPG", "IL1",  "TGFB", "EGF", "GDF", "GRN", "NT", "RELN", "SLIT", "VISFATIN", "ncWNT", "PERIOSTIN", "PTN"), pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                      "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                      "ANGPT", "SEMA3", "CCL", "AGRN", "COMPLEMENT", "HSPG", "IL1",  "TGFB", "EGF", "GDF", "GRN", "NT", "RELN", "SLIT", "VISFATIN", "ncWNT", "PERIOSTIN", "PTN"), pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

##LSEC1 and LSEC2 analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                      "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                      "CCL", "SEMA3", "AGRN", "COMPLEMENT", "HSPG", "IL1",  "TGFB", "VISFATIN", "EGF", "GDF", "GRN", "KIT", "NT", "RELN", "SLIT", "VWF", "ncWNT", "PERIOSTIN", "PTN"), pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                      "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                      "CCL", "SEMA3", "AGRN", "COMPLEMENT", "HSPG", "IL1",  "TGFB", "VISFATIN", "EGF", "GDF", "GRN", "KIT", "NT", "RELN", "SLIT", "VWF", "ncWNT", "PERIOSTIN", "PTN"), pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "outgoing", width = 5, height = 20, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "incoming", width = 5, height = 20, font.size = 8, font.size.title = 9)
ht1 + ht2

##immune related in 5% analysis 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("BAFF", "CX3C", "EDN", "FASLG", "IL12", "IL16", "IL2", "LIFR", "TNF", "TRAIL", "TWEAK"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("BAFF", "CX3C", "EDN", "FASLG", "IL12", "IL16", "IL2", "LIFR", "TNF", "TRAIL", "TWEAK"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (all even when not detected)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (only those detected in each group)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("IL1", "TGFb", "IGF", "FLT3", "VISFATIN", "27HC", "TRAIL", "TNF", "FASLG", "IL2", "IL12", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "AGT", "SAA", "TWEAK", "CX3C"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("IL1", "TGFb", "IGF", "FLT3", "VISFATIN", "27HC", "TRAIL", "TNF", "FASLG", "IL2", "IL12", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "AGT", "SAA", "TWEAK", "CX3C"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("GALECTIN", "COMPLEMENT", "SPP1", "TGFb", "NRG", "Glutamate", "IL1", "IL2", "FASLG", "FLT3","TRAIL"), pattern = "outgoing", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("GALECTIN", "COMPLEMENT", "SPP1", "TGFb", "NRG", "Glutamate", "IL1", "IL2", "FASLG", "FLT3","TRAIL"), pattern = "incoming", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% filtered pathways based on those detected already in 25 and 75% analyses
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("TENASCIN", "Prostaglandin", "DHT", "CHEMERIN", "WNT", "IGFBP", "SAA", "CXCL", "ACTIVIN", "ANNEXIN", "THBS", "IL2", "FASLG", "Adenosine", "EDN", "Desmosterol", "FLT3", "Histamine", "Androsterone", "NGF", "2-AG", "CysLTs", "EDA", "BMP10", "BRADYKININ", "Androstenedione", "LIFR", "TNF", "TRAIL", "IL12", "IL16", "HH", "TWEAK", "BAFF", "CX3C"), pattern = "outgoing", width = 5, height = 15, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, signaling = c("TENASCIN", "Prostaglandin", "DHT", "CHEMERIN", "WNT", "IGFBP", "SAA", "CXCL", "ACTIVIN", "ANNEXIN", "THBS", "IL2", "FASLG", "Adenosine", "EDN", "Desmosterol", "FLT3", "Histamine", "Androsterone", "NGF", "2-AG", "CysLTs", "EDA", "BMP10", "BRADYKININ", "Androstenedione", "LIFR", "TNF", "TRAIL", "IL12", "IL16", "HH", "TWEAK", "BAFF", "CX3C"), pattern = "incoming", width = 5, height = 15, font.size = 8, font.size.title = 9)
ht1 + ht2

##75% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PreT, pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

PreT.df.net <- subsetCommunication(cellchat.PreT, slot.name = "net")
write.csv(PreT.df.net, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/PreT.scS5.v2.25pt.CD4CD8.df.net.default.csv")

saveRDS(cellchat.PreT, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.PreT.scS5.25pt.CD4CD8.default.rds")
##after saving PreT move on to MWD.T
cellchat.MWD.T <- subsetData(cellchat.MWD.T)
future::plan("multisession", workers = 2)

cellchat.MWD.T <- identifyOverExpressedGenes(cellchat.MWD.T)
cellchat.MWD.T <- identifyOverExpressedInteractions(cellchat.MWD.T)

unique(cellchat.MWD.T@idents)

cellchat.MWD.T@idents = droplevels(cellchat.MWD.T@idents, exclude = setdiff(levels(cellchat.MWD.T@idents),unique(cellchat.MWD.T@idents)))

cellchat.MWD.T <- computeCommunProb(cellchat.MWD.T)
#!cellchat.MWD.T <- computeCommunProb(cellchat.MWD.T, type = "truncatedMean", trim = 0.05)
#!cellchat.MWD.T <- computeCommunProb(cellchat.MWD.T, type = "truncatedMean", trim = 0.75)

cellchat.MWD.T <- filterCommunication(cellchat.MWD.T, min.cells = 2)

cellchat.MWD.T <- computeCommunProbPathway(cellchat.MWD.T)

cellchat.MWD.T <- aggregateNet(cellchat.MWD.T)

groupSize <- as.numeric(table(cellchat.MWD.T@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat.MWD.T@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.MWD.T@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat.MWD.T@netP$pathways
pathways.show <- cellchat.MWD.T@netP$pathways

cellchat.MWD.T <- netAnalysis_computeCentrality(cellchat.MWD.T, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2
#reordering
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3", "AGRN", "HSPG", "COMPLEMENT", "TGFB", "CXCL", "EGF", "GDF", "MIF", "NT", "SLIT", "TENASCIN"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3", "AGRN", "HSPG", "COMPLEMENT", "TGFB", "CXCL", "EGF", "GDF", "MIF", "NT", "SLIT", "TENASCIN"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##default analysis for filtered pathways detected in 75% already
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("COLLAGEN", "VEGF", "DHEA", "COMPLEMENT", "EGF", "Testosterone", "TENASCIN", "SLIT", "SEMA3", "GAS", "CXCL", "Androstenedione", "CSF", "GDF", "SerotoninDopamin", "AGRN", "HSPG", "TGFb", "MIF", "NT", "Prostaglandin"), pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("COLLAGEN", "VEGF", "DHEA", "COMPLEMENT", "EGF", "Testosterone", "TENASCIN", "SLIT", "SEMA3", "GAS", "CXCL", "Androstenedione", "CSF", "GDF", "SerotoninDopamin", "AGRN", "HSPG", "TGFb", "MIF", "NT", "Prostaglandin"), pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2

##for combined stromal analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "ANGPT", "CCL", "SEMA3", "AGRN", "COMPLEMENT", "HSPG", "TGFB", "CXCL", "EGF", "GDF", "MIF", "NT", "SLIT", "TENASCIN"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "ANGPT", "CCL", "SEMA3", "AGRN", "COMPLEMENT", "HSPG", "TGFB", "CXCL", "EGF", "GDF", "MIF", "NT", "SLIT", "TENASCIN"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##LSEC1 and LSEC2 analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3", "AGRN", "COMPLEMENT", "HSPG", "TGFB", "VISFATIN", "CXCL", "EGF", "GDF", "MIF", "NT", "SLIT", "TENASCIN"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3", "AGRN", "COMPLEMENT", "HSPG", "TGFB", "VISFATIN", "CXCL", "EGF", "GDF", "MIF", "NT", "SLIT", "TENASCIN"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2
##5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "outgoing", width = 5, height = 20, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "incoming", width = 5, height = 20, font.size = 8, font.size.title = 9)
ht1 + ht2

##immune related in 5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("BAFF", "EDN", "FASLG", "IFN-I", "IL10", "IL16", "IL2", "LIFR", "TNF", "TRAIL", "TWEAK"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("BAFF", "EDN", "FASLG", "IFN-I", "IL10", "IL16", "IL2", "LIFR", "TNF", "TRAIL", "TWEAK"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (all even when not detected)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (only those detected in each group)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("IL1", "TGFb", "IGF", "CypA", "VISFATIN", "EPO", "27HC", "TRAIL", "TNF", "FASLG", "IL2", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "MIF", "AGT", "SAA", "TWEAK", "IFN-I", "IL10"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("IL1", "TGFb", "IGF", "CypA", "VISFATIN", "EPO", "27HC", "TRAIL", "TNF", "FASLG", "IL2", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "MIF", "AGT", "SAA", "TWEAK", "IFN-I", "IL10"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("SPP1", "COMPLEMENT", "GALECTIN", "NRG", "TGFb", "IL1", "Glutamate", "IL2", "FASLG", "TRAIL"), pattern = "outgoing", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("SPP1", "COMPLEMENT", "GALECTIN", "NRG", "TGFb", "IL1", "Glutamate", "IL2", "FASLG", "TRAIL"), pattern = "incoming", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% filtered pathways based on those detected in 25 and 75% analyses
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("IL1", "THBS", "SAA", "GRN", "RELN", "VISFATIN", "ACTIVIN", "Glutamate", "CHEMERIN", "ncWNT", "IL2", "IGFBP", "Adenosine", "PERIOSTIN", "Androsterone", "WNT", "VWF", "2-AG", "KIT", "FASLG", "DHT", "IFN-I", "LIFR", "TNF", "Desmosterol", "ANNEXIN", "CysLTs", "PTPR", "EDN", "EDA", "BMP10", "TWEAK", "HH", "IL16", "TRAIL", "IL10", "BAFF", "Histamine"), pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, signaling = c("IL1", "THBS", "SAA", "GRN", "RELN", "VISFATIN", "ACTIVIN", "Glutamate", "CHEMERIN", "ncWNT", "IL2", "IGFBP", "Adenosine", "PERIOSTIN", "Androsterone", "WNT", "VWF", "2-AG", "KIT", "FASLG", "DHT", "IFN-I", "LIFR", "TNF", "Desmosterol", "ANNEXIN", "CysLTs", "PTPR", "EDN", "EDA", "BMP10", "TWEAK", "HH", "IL16", "TRAIL", "IL10", "BAFF", "Histamine"), pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

##75% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MWD.T, pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

MWD.T.df.net <- subsetCommunication(cellchat.MWD.T, slot.name = "net")
write.csv(MWD.T.df.net, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/MWD.T.scS5.v2.25pt.CD4CD8.df.net.default.csv")

saveRDS(cellchat.MWD.T, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.MWD.T.scS5.25pt.CD4CD8.default.rds")
##after saving MWD.T perform for MRD.T
cellchat.MRD.T <- subsetData(cellchat.MRD.T)
future::plan("multisession", workers = 2)

cellchat.MRD.T <- identifyOverExpressedGenes(cellchat.MRD.T)
cellchat.MRD.T <- identifyOverExpressedInteractions(cellchat.MRD.T)

unique(cellchat.MRD.T@idents)

cellchat.MRD.T@idents = droplevels(cellchat.MRD.T@idents, exclude = setdiff(levels(cellchat.MRD.T@idents),unique(cellchat.MRD.T@idents)))

cellchat.MRD.T <- computeCommunProb(cellchat.MRD.T)
#!cellchat.MRD.T <- computeCommunProb(cellchat.MRD.T, type = "truncatedMean", trim = 0.05)
#!cellchat.MRD.T <- computeCommunProb(cellchat.MRD.T, type = "truncatedMean", trim = 0.75)

cellchat.MRD.T <- filterCommunication(cellchat.MRD.T, min.cells = 2)

cellchat.MRD.T <- computeCommunProbPathway(cellchat.MRD.T)

cellchat.MRD.T <- aggregateNet(cellchat.MRD.T)

groupSize <- as.numeric(table(cellchat.MRD.T@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat.MRD.T@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.MRD.T@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat.MRD.T@netP$pathways
pathways.show <- cellchat.MRD.T@netP$pathways

cellchat.MRD.T <- netAnalysis_computeCentrality(cellchat.MRD.T, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2
#reordering
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##default analysis filtering out pathways detected in 75% already
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("COLLAGEN", "BMP", "LAMININ", "Cholesterol", "SPP1", "FGF", "CypA", "AGT", "HGF", "GAS", "PDGF", "SEMA3", "Testosterone", "NRG", "CCL", "SerotoninDopamin", "ANGPT", "CSF"), pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("COLLAGEN", "BMP", "LAMININ", "Cholesterol", "SPP1", "FGF", "CypA", "AGT", "HGF", "GAS", "PDGF", "SEMA3", "Testosterone", "NRG", "CCL", "SerotoninDopamin", "ANGPT", "CSF"), pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2

##for stromal combined analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##LSEC1 and LSEC2 analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                       "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                       "CCL", "SEMA3"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2
##5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "outgoing", width = 5, height = 18, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "incoming", width = 5, height = 18, font.size = 8, font.size.title = 9)
ht1 + ht2

##immune related in 5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("EDN", "FASLG", "IFN-I", "IL16", "IL2", "LIFR", "TWEAK"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("EDN", "FASLG", "IFN-I", "IL16", "IL2", "LIFR", "TWEAK"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (all even when not detected)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (only those detected in each group)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("IL1", "TGFb", "IGF", "CypA", "VISFATIN", "27HC", "FASLG", "IL2", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "MIF", "AGT", "SAA", "TWEAK", "IFN-I"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("IL1", "TGFb", "IGF", "CypA", "VISFATIN", "27HC", "FASLG", "IL2", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "MIF", "AGT", "SAA", "TWEAK", "IFN-I"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("GALECTIN", "COMPLEMENT", "SPP1", "NRG", "TGFb", "IL1", "Glutamate", "FASLG", "IL2"), pattern = "outgoing", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("GALECTIN", "COMPLEMENT", "SPP1", "NRG", "TGFb", "IL1", "Glutamate", "FASLG", "IL2"), pattern = "incoming", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% filtered pathways based on those detected in 25 and 75% analyses
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("COMPLEMENT", "DHEA", "TENASCIN", "EGF", "SLIT", "MIF", "TGFb", "CXCL", "IL1", "VISFATIN", "Prostaglandin", "GRN", "Glutamate", "AGRN", "RELN", "SAA", "HSPG", "ANNEXIN", "PERIOSTIN", "THBS", "Androstenedione", "PTPR", "GDF", "CHEMERIN", "IGFBP", "EDN", "NT", "VWF", "Androsterone", "KIT", "FASLG", "ACTIVIN", "Adenosine", "2-AG", "IL2", "DHT", "ncWNT", "WNT", "BMP10", "LIFR", "IL16", "Desmosterol", "IFN-I", "TWEAK"), pattern = "outgoing", width = 5, height = 11, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, signaling = c("COMPLEMENT", "DHEA", "TENASCIN", "EGF", "SLIT", "MIF", "TGFb", "CXCL", "IL1", "VISFATIN", "Prostaglandin", "GRN", "Glutamate", "AGRN", "RELN", "SAA", "HSPG", "ANNEXIN", "PERIOSTIN", "THBS", "Androstenedione", "PTPR", "GDF", "CHEMERIN", "IGFBP", "EDN", "NT", "VWF", "Androsterone", "KIT", "FASLG", "ACTIVIN", "Adenosine", "2-AG", "IL2", "DHT", "ncWNT", "WNT", "BMP10", "LIFR", "IL16", "Desmosterol", "IFN-I", "TWEAK"), pattern = "incoming", width = 5, height = 11, font.size = 8, font.size.title = 9)
ht1 + ht2

##75% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "outgoing", width = 5, height = 4, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.T, pattern = "incoming", width = 5, height = 4, font.size = 8, font.size.title = 9)
ht1 + ht2

MRD.T.df.net <- subsetCommunication(cellchat.MRD.T, slot.name = "net")
write.csv(MRD.T.df.net, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/MRD.T.scS5.v2.25pt.CD4CD8.df.net.default.csv")

saveRDS(cellchat.MRD.T, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.MRD.T.scS5.25pt.CD4CD8.default.rds")
##after saving MRD.T results run on last group MRD.NT
cellchat.MRD.NT <- subsetData(cellchat.MRD.NT)
future::plan("multisession", workers = 2)

cellchat.MRD.NT <- identifyOverExpressedGenes(cellchat.MRD.NT)
cellchat.MRD.NT <- identifyOverExpressedInteractions(cellchat.MRD.NT)

unique(cellchat.MRD.NT@idents)

cellchat.MRD.NT@idents = droplevels(cellchat.MRD.NT@idents, exclude = setdiff(levels(cellchat.MRD.NT@idents),unique(cellchat.MRD.NT@idents)))

cellchat.MRD.NT <- computeCommunProb(cellchat.MRD.NT)
#!cellchat.MRD.NT <- computeCommunProb(cellchat.MRD.NT, type = "truncatedMean", trim = 0.05)
#!cellchat.MRD.NT <- computeCommunProb(cellchat.MRD.NT, type = "truncatedMean", trim = 0.75)

cellchat.MRD.NT <- filterCommunication(cellchat.MRD.NT, min.cells = 2)

cellchat.MRD.NT <- computeCommunProbPathway(cellchat.MRD.NT)

cellchat.MRD.NT <- aggregateNet(cellchat.MRD.NT)

groupSize <- as.numeric(table(cellchat.MRD.NT@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat.MRD.NT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.MRD.NT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat.MRD.NT@netP$pathways
pathways.show <- cellchat.MRD.NT@netP$pathways

cellchat.MRD.NT <- netAnalysis_computeCentrality(cellchat.MRD.NT, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2
#reordering
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                        "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                        "CCL", "IL1", "COMPLEMENT", "TGFB", "GRN", "MIF", "VISFATIN", "SAA"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                        "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                        "CCL", "IL1", "COMPLEMENT", "TGFB", "GRN", "MIF", "VISFATIN", "SAA"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##default analysis filtering pathways detected in 75% ALREADY
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("LAMININ", "IGF", "COLLAGEN", "DHEA", "SAA", "AGT", "COMPLEMENT", "CCL", "GAS", "PDGF", "IL1", "ANGPT", "CSF", "Androstenedione", "SerotoninDopamin", "TGFb", "Prostaglandin", "VISFATIN", "GRN", "MIF"), pattern = "outgoing", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("LAMININ", "IGF", "COLLAGEN", "DHEA", "SAA", "AGT", "COMPLEMENT", "CCL", "GAS", "PDGF", "IL1", "ANGPT", "CSF", "Androstenedione", "SerotoninDopamin", "TGFb", "Prostaglandin", "VISFATIN", "GRN", "MIF"), pattern = "incoming", width = 5, height = 9, font.size = 8, font.size.title = 9)
ht1 + ht2

##for stromal combined analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                        "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                        "ANGPT", "CCL", "COMPLEMENT", "IL1", "TGFB", "GRN", "MIF", "VISFATIN", "SAA"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("AGT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                        "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                        "ANGPT", "CCL", "COMPLEMENT", "IL1", "TGFB", "GRN", "MIF", "VISFATIN", "SAA"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##LSEC1 and LSEC2 analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                        "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                        "CCL", "COMPLEMENT", "IL1", "TGFB", "VISFATIN", "GRN", "MIF", "SAA"), pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("AGT", "ANGPT", "ANGPTL", "BMP", "COLLAGEN", "CSF", "FGF", "FN1", "GALECTIN", "GAS",
                                                                        "HGF", "IGF", "LAMININ", "NRG", "PARs", "PDGF", "PROS", "SPP1", "VEGF", "VTN",
                                                                        "CCL", "COMPLEMENT", "IL1", "TGFB", "VISFATIN", "GRN", "MIF", "SAA"), pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "outgoing", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "incoming", width = 5, height = 8, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "outgoing", width = 5, height = 20, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "incoming", width = 5, height = 20, font.size = 8, font.size.title = 9)
ht1 + ht2

##immune related in 5% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("BAFF", "EDN", "FASLG", "IFN-I", "IL12", "IL16", "IL2", "LIFR", "TNF", "TWEAK"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("BAFF", "EDN", "FASLG", "IFN-I", "IL12", "IL16", "IL2", "LIFR", "TNF", "TWEAK"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (all even when not detected)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "outgoing", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("IL1", "IL2", "IL12", "FASLG", "TRAIL", "CD40", "TNF", "IFN-I", "IL16", "FLT3", "CCL", "CXCL", "CSF", "MIF", "SAA", "CX3C", "TWEAK", "IGF", "AGT", "GRN", "VISFATIN", "LIFR", "EPO", "TGFb", "27HC", "CypA", "IL10"), pattern = "incoming", width = 5, height = 7, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% selected pathways (only those detected in each group)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("IL1", "TGFb", "IGF", "CypA", "VISFATIN", "27HC", "TNF", "FASLG", "IL2", "IL12", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "MIF", "AGT", "SAA", "TWEAK", "IFN-I"), pattern = "outgoing", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("IL1", "TGFb", "IGF", "CypA", "VISFATIN", "27HC", "TNF", "FASLG", "IL2", "IL12", "GRN", "IL16", "GALECTIN", "CCL", "CSF", "CXCL", "LIFR", "MIF", "AGT", "SAA", "TWEAK", "IFN-I"), pattern = "incoming", width = 5, height = 6, font.size = 8, font.size.title = 9)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("SPP1", "COMPLEMENT", "NRG", "IL1", "Glutamate", "TGFb", "IL2", "FASLG"), pattern = "outgoing", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("SPP1", "COMPLEMENT", "NRG", "IL1", "Glutamate", "TGFb", "IL2", "FASLG"), pattern = "incoming", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht1 + ht2

##5% filtered pathways based on those detected in 25 and 75% analyses
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("SEMA3", "TENASCIN", "EGF", "Glutamate", "CXCL", "CHEMERIN", "IGFBP", "THBS", "NT", "AGRN", "SLIT", "RELN", "PROCR", "ACTIVIN", "HSPG", "IL2", "GDF", "ncWNT", "ANNEXIN", "KIT", "Desmosterol", "LIFR", "WNT", "Adenosine", "DHT", "TNF", "VWF", "Androsterone", "BMP10", "2-AG", "FASLG", "CysLTs", "IFN-I", "PTPR", "CTSG", "IL12", "Histamine", "TWEAK", "IL16", "LXA4", "ENHO", "BAFF"), pattern = "outgoing", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, signaling = c("SEMA3", "TENASCIN", "EGF", "Glutamate", "CXCL", "CHEMERIN", "IGFBP", "THBS", "NT", "AGRN", "SLIT", "RELN", "PROCR", "ACTIVIN", "HSPG", "IL2", "GDF", "ncWNT", "ANNEXIN", "KIT", "Desmosterol", "LIFR", "WNT", "Adenosine", "DHT", "TNF", "VWF", "Androsterone", "BMP10", "2-AG", "FASLG", "CysLTs", "IFN-I", "PTPR", "CTSG", "IL12", "Histamine", "TWEAK", "IL16", "LXA4", "ENHO", "BAFF"), pattern = "incoming", width = 5, height = 10, font.size = 8, font.size.title = 9)
ht1 + ht2

##75% analysis
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "outgoing", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.MRD.NT, pattern = "incoming", width = 5, height = 5, font.size = 8, font.size.title = 9)
ht1 + ht2

MRD.NT.df.net <- subsetCommunication(cellchat.MRD.NT, slot.name = "net")
write.csv(MRD.NT.df.net, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/MRD.NT.scS5.v2.25pt.CD4CD8.df.net.default.csv")

saveRDS(cellchat.MRD.NT, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.MRD.NT.scS5.25pt.CD4CD8.default.rds")

##now with all the saved CellChat analyses we can try to run comparative cellchat analyses
##prototypical vignettes only do 2 conditions, but creator of CellChat says more is possible

##load these files for default CellChatv2 analyses
cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt/Saved_Rfiles/cellchatv2.CD.scS5.new.labels.default.rds")
cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt/Saved_Rfiles/cellchatv2.PreT.scS5.new.labels.default.rds")
cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt/Saved_Rfiles/cellchatv2.MWD.T.scS5.new.labels.default.rds")
cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt/Saved_Rfiles/cellchatv2.MRD.T.scS5.new.labels.default.rds")
cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt/Saved_Rfiles/cellchatv2.MRD.NT.scS5.new.labels.default.rds")

##load these files for 5% CellChatv2 analyses
cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 5pt/Saved_Rfiles/cellchatv2.CD.scS5.5pt.new.labels.default.rds")
cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 5pt/Saved_Rfiles/cellchatv2.PreT.scS5.5pt.new.labels.default.rds")
cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 5pt/Saved_Rfiles/cellchatv2.MWD.T.scS5.5pt.new.labels.default.rds")
cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 5pt/Saved_Rfiles/cellchatv2.MRD.T.scS5.5pt.new.labels.default.rds")
cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 5pt/Saved_Rfiles/cellchatv2.MRD.NT.scS5.5pt.new.labels.default.rds")

##load these files for 75% CellChatv2 analyses
cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.CD.scS5.75pt.new.labels.default.rds")
cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.PreT.scS5.75pt.new.labels.default.rds")
cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.MWD.T.scS5.75pt.new.labels.default.rds")
cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.MRD.T.scS5.75pt.new.labels.default.rds")
cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.MRD.NT.scS5.75pt.new.labels.default.rds")

##load these files for 75% abbreviated names CellChatv2 analyses
cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.CD.scS5.75pt.abbreviated.new.labels.default.rds")
cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.PreT.scS5.75pt.abbreviated.new.labels.default.rds")
cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.MWD.T.scS5.75pt.abbreviated.new.labels.default.rds")
cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.MRD.T.scS5.75pt.abbreviated.new.labels.default.rds")
cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 75pt/Saved_Rfiles/cellchatv2.MRD.NT.scS5.75pt.abbreviated.new.labels.default.rds")

##this the files for LSEC1 and LSEC2 analysis
#!cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 LSEC1 and LSEC2/Saved_Rfiles/cellchatv2.CD.scS5.LSEC1and2.default.rds")
#!cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 LSEC1 and LSEC2/Saved_Rfiles/cellchatv2.PreT.scS5.LSEC1and2.default.rds")
#!cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 LSEC1 and LSEC2/Saved_Rfiles/cellchatv2.MWD.T.scS5.LSEC1and2.default.rds")
#!cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 LSEC1 and LSEC2/Saved_Rfiles/cellchatv2.MRD.T.scS5.LSEC1and2.default.rds")
#!cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 LSEC1 and LSEC2/Saved_Rfiles/cellchatv2.MRD.NT.scS5.LSEC1and2.default.rds")

##use these files for 5% analysis including LSEC1 and LSEC2
#!cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 5pt/Saved_Rfiles/cellchatv2.CD.scS5.5pt.default.rds")
#!cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 5pt/Saved_Rfiles/cellchatv2.PreT.scS5.5pt.default.rds")
#!cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 5pt/Saved_Rfiles/cellchatv2.MWD.T.scS5.5pt.default.rds")
#!cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 5pt/Saved_Rfiles/cellchatv2.MRD.T.scS5.5pt.default.rds")
#!cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChat v2 5pt/Saved_Rfiles/cellchatv2.MRD.NT.scS5.5pt.default.rds")

##use these files for abbreviated cell type default analysis
cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt abbreviated.labels/Saved_Rfiles/cellchatv2.CD.scS5.abbreviated.new.labels.default.rds")
cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt abbreviated.labels/Saved_Rfiles/cellchatv2.PreT.scS5.abbreviated.new.labels.default.rds")
cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt abbreviated.labels/Saved_Rfiles/cellchatv2.MWD.T.scS5.abbreviated.new.labels.default.rds")
cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt abbreviated.labels/Saved_Rfiles/cellchatv2.MRD.T.scS5.abbreviated.new.labels.default.rds")
cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt abbreviated.labels/Saved_Rfiles/cellchatv2.MRD.NT.scS5.abbreviated.new.labels.default.rds")

##use these files for 25pt select.subsets analysis (CD4/CD8/KC/Mac) 
cellchat.CD <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.CD.scS5.25pt.CD4CD8.default.rds")
cellchat.PreT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.PreT.scS5.25pt.CD4CD8.default.rds")
cellchat.MWD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.MWD.T.scS5.25pt.CD4CD8.default.rds")
cellchat.MRD.T <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.MRD.T.scS5.25pt.CD4CD8.default.rds")
cellchat.MRD.NT <- readRDS("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/SciReports 2024/CellChat/CellChatv2 25pt CD4CD8/Saved_Rfiles/cellchatv2.MRD.NT.scS5.25pt.CD4CD8.default.rds")

object.list <- list(CD = cellchat.CD, WD.nf = cellchat.PreT, WD.t = cellchat.MWD.T, RD.t = cellchat.MRD.T, RD.n = cellchat.MRD.NT)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##this is for CDvPreT comparisons with rankNet and signalingrole Scatter ONLY
object.list.CDvPreT <- list(CD = cellchat.CD, PreT = cellchat.PreT)

cellchat.CDvPreT <- mergeCellChat(object.list.CDvPreT, add.names = names(object.list.CDvPreT))

##this is for PreTvMWD.T comparisons with rankNet and signalingrole Scatter ONLY
object.list.PreTvMWD.T <- list(PreT = cellchat.PreT, MWD.T = cellchat.MWD.T)

cellchat.PreTvMWD.T <- mergeCellChat(object.list.PreTvMWD.T, add.names = names(object.list.PreTvMWD.T))

##once all analyzed groups are merged together we can compare total number of interactions and strengths
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5), measure = "weight")
gg1 + gg2

##look for differential number of interactons with circos plots
##slightly too many interaxns with 5 groups so not very useful
#!par(mfrow = c(1,2), xpd=TRUE)
#!netVisual_diffInteraction(cellchat, weight.scale = T)
#!netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

##but employing this on the heatmap could help
##however with more than 2 datasets it is not that useful either
#!gg1 <- netVisual_heatmap(cellchat)
#!gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#!gg1 + gg2

##can show number of interaxns for any specific cell types across all groups
##this requires aggregation of cell populations into a new CellChat object
group.cellType <- c(rep("Stromal"), rep("Tcell"), rep("NKT"), rep("NK"), rep("Neutrophil"), rep("Monocyte"), rep("Macrophage"), rep("LSEC"), rep("HSC"), rep("Hepatocyte"), rep("Fibroblast"), rep("Myofibroblast"), rep("Endothelial"), rep("DC"), rep("Cholangiocyte"), rep("Cancer"), rep("Bcell"))
##use this for combined stromal analysis (stromal = stromal, LSEC, endo)
group.cellType <- c(rep("Stromal"), rep("Tcell"), rep("NKT"), rep("NK"), rep("Neutrophil"), rep("Monocyte"), rep("Macrophage"), rep("HSC"), rep("Hepatocyte"), rep("Fibroblast"), rep("Myofibroblast"), rep("DC"), rep("Cholangiocyte"), rep("Cancer"), rep("Bcell"))
##use this for LSEC1 and LSEC2 analysis
group.cellType <- c(rep("Stromal"), rep("Tcell"), rep("NKT"), rep("NK"), rep("Neutrophil"), rep("Monocyte"), rep("Macrophage"), rep("LSEC1"), rep("LSEC2"), rep("HSC"), rep("Hepatocyte"), rep("Fibroblast"), rep("Myofibroblast"), rep("Endothelial"), rep("DC"), rep("Cholangiocyte"), rep("Cancer"), rep("Bcell"))
##use for abbreviated groups
group.cellType <- c(rep("Stromal"), rep("Tcell"), rep("NKT"), rep("NK"), rep("Neutro"), rep("Mono"), rep("Mac"), rep("LSEC"), rep("HSC"), rep("Hep"), rep("Fibro"), rep("Myofibro"), rep("Endo"), rep("DC"), rep("Cholangio"), rep("Cancer"), rep("Bcell"))

#!group.cellType <- c(rep("Stromal", 4), rep("Tcell", 4), rep("NKT", 4), rep("NK", 4), rep("Neutrophil", 4), rep("Monocyte", 4), rep("Macrophage", 4), rep("LSEC", 4), rep("HSC", 4), rep("Hepatocyte", 4), rep("Fibro", 4), rep("Endothelial", 4), rep("DC", 4), rep("Cholangiocyte", 4), rep("Cancer", 4), rep("Bcell", 4))
group.cellType <- factor(group.cellType, levels = c("Stromal", "Tcell", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "LSEC", "HSC", "Hepatocyte", "Fibroblast", "Myofibroblast", "Endothelial", "DC", "Cholangiocyte", "Cancer", "Bcell"))
##same here for combined stromal analysis (stromal = stromal, LSEC, endo)
group.cellType <- factor(group.cellType, levels = c("Stromal", "Tcell", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "HSC", "Hepatocyte", "Fibroblast", "Myofibroblast", "DC", "Cholangiocyte", "Cancer", "Bcell"))
##same here for LSEC1 and LSEC2 analysis
group.cellType <- factor(group.cellType, levels = c("Stromal", "Tcell", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "LSEC1", "LSEC2", "HSC", "Hepatocyte", "Fibroblast", "Myofibroblast", "Endothelial", "DC", "Cholangiocyte", "Cancer", "Bcell"))
##same use for abbreviated groups
group.cellType <- factor(group.cellType, levels = c("Stromal", "Tcell", "NKT", "NK", "Neutro", "Mono", "Mac", "LSEC", "HSC", "Hep", "Fibro", "Myofibro", "Endo", "DC", "Cholangio", "Cancer", "Bcell"))

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))


weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, comparison = c(1,2), weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat.merged, comparison = c(2,3), weight.scale = T, measure = "count.merged", label.edge = T)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, comparison = c(1,2), sources.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvWD.nf")
netVisual_diffInteraction(cellchat.merged, comparison = c(2,3), sources.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "WD.nfvWD.t")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, comparison = c(3,4), sources.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "WD.tvRD.t")
netVisual_diffInteraction(cellchat.merged, comparison = c(4,5), sources.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "RD.tvRD.n")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, comparison = c(1,2), targets.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvWD.nf")
netVisual_diffInteraction(cellchat.merged, comparison = c(2,3), targets.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "WD.nfvWD.t")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, comparison = c(3,4), targets.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "WD.tvRD.t")
netVisual_diffInteraction(cellchat.merged, comparison = c(4,5), targets.use = "HSC", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "RD.tvRD.n")

##additional group comparisons to perform (CDvRDn, WDnfvRDn, CDvRDt, CDvWDt)
netVisual_diffInteraction(cellchat.merged, comparison = c(1,5), sources.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvRDn")
netVisual_diffInteraction(cellchat.merged, comparison = c(2,5), sources.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "WD.nfvRDn")
netVisual_diffInteraction(cellchat.merged, comparison = c(1,4), sources.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvRD.t")
netVisual_diffInteraction(cellchat.merged, comparison = c(1,3), sources.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvWD.t")

netVisual_diffInteraction(cellchat.merged, comparison = c(1,5), targets.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvRDn")
netVisual_diffInteraction(cellchat.merged, comparison = c(2,5), targets.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "WD.nfvRDn")
netVisual_diffInteraction(cellchat.merged, comparison = c(1,4), targets.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvRD.t")
netVisual_diffInteraction(cellchat.merged, comparison = c(1,3), targets.use = "Cancer", weight.scale = T, measure = "count", label.edge = T, vertex.size.max = 3, alpha.edge = 0.8, vertex.label.cex = 1.15, edge.label.cex = 1.15, margin = 0.1, title.name = "CDvWD.t")

##or we can show the number of differential interaxns and interaxn strength b/w cell types
##this is also mainly just for two dataset comparison though
#!par(mfrow = c(1,2), xpd=TRUE)
#!netVisual_diffInteraction(cellchat.S.T.HSC, weight.scale = T, measure = "count.merged", label.edge = T)
#!netVisual_diffInteraction(cellchat.S.T.HSC, weight.scale = T, measure = "weight.merged", label.edge = T)


##compare major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + scale_y_continuous(breaks = seq(0, 16, by = 2)) + scale_x_continuous(breaks = seq(0, 20, by = 2))
}
patchwork::wrap_plots(plots = gg)

##try to make same plots but only for stromal from each group
##try to identify specific changes of a cell type b/w conditions 
##this can only handle 2 comparisons at a time
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(2,3), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(3,4), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(4,5), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,3), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,4), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,5), idents.use = "Stromal") + scale_y_continuous(limits = c(-2.0,2.0)) + scale_x_continuous(limits = c(-2.0,2.0))

patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))
#!patchwork::wrap_plots(plots = list(gg1,gg2,gg3))
##make images smaller and custom titles
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), idents.use = "Stromal", xlabel = "Outgoing", ylabel = "Incoming") + ggtitle("Stromal CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), idents.use = "Stromal", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Stromal CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

patchwork::wrap_plots(plots = list(gg1))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(2,3), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(3,4), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(5,4), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))

patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,3), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,4), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,5), idents.use = "Cancer") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))

patchwork::wrap_plots(plots = list(gg1))

##LSEC1 and LSEC2 analysis
#!gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), idents.use = "LSEC1")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(2,3), idents.use = "LSEC1") + scale_y_continuous(limits = c(-3.0,3.0)) + scale_x_continuous(limits = c(-3.0,3.0))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(5,3), idents.use = "LSEC1") + scale_y_continuous(limits = c(-3.0,3.0)) + scale_x_continuous(limits = c(-3.0,3.0))
#!gg4 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,5), idents.use = "LSEC1")
patchwork::wrap_plots(plots = gg2)
patchwork::wrap_plots(plots = gg3)

#!patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))

#!gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), idents.use = "LSEC2")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(2,3), idents.use = "LSEC2") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
#!gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(3,4), idents.use = "LSEC2")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(5,3), idents.use = "LSEC2") + scale_y_continuous(limits = c(-4.0,4.0)) + scale_x_continuous(limits = c(-3.0,3.0))
patchwork::wrap_plots(plots = gg2)
patchwork::wrap_plots(plots = gg4)

#!patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##use this section to try to perform pairwise on all cell types for different analyses
##due to issues/inconsistencies in MRD.T cell type idents (LSEC split) liftCellChat workaround may not work so only do PreTvMWD.T this way
##when using analysis without LSEC split all should work correctly
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), signaling = c("FN1", "PARs", "LAMININ", "VTN", "BMP", "Cholesterol", "GALECTIN", "VEGF", "SPP1", "COLLAGEN", "COMPLEMENT", "IGF", "DHEAS", "SEMA3", "DHEA"), idents.use = "Bcell")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1,2), signaling = c("FN1", "PARs", "LAMININ", "VTN", "BMP", "Cholesterol", "GALECTIN", "VEGF", "SPP1", "COLLAGEN", "COMPLEMENT", "IGF", "DHEAS", "SEMA3", "DHEA"), idents.use = "Bcell", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Bcell CD vs. WD.nf")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(2,3), signaling = c("FN1", "PARs", "LAMININ", "VTN", "BMP", "Cholesterol", "GALECTIN", "VEGF", "SPP1", "COLLAGEN", "COMPLEMENT", "IGF", "DHEAS", "SEMA3", "DHEA", "FGF", "CypA", "27HC", "NRG", "ANGPTL", "AGT", "HGF", "EGF", "PROS", "Testosterone", "TENASCIN"), idents.use = "Bcell")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(2,3), signaling = c("FN1", "PARs", "LAMININ", "VTN", "BMP", "Cholesterol", "GALECTIN", "VEGF", "SPP1", "COLLAGEN", "COMPLEMENT", "IGF", "DHEAS", "SEMA3", "DHEA", "FGF", "CypA", "27HC", "NRG", "ANGPTL", "AGT", "HGF", "EGF", "PROS", "Testosterone", "TENASCIN"), idents.use = "Bcell", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Bcell WD.nf vs. WD.t")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(3,4), signaling = c("PARs", "FN1", "VTN", "DHEAS", "FGF", "Cholesterol", "LAMININ", "COLLAGEN", "VEGF", "CypA", "27HC", "BMP", "GALECTIN", "NRG", "SPP1", "IGF", "DHEA", "ANGPTL", "AGT", "HGF", "COMPLEMENT", "EGF", "PROS", "Testosterone", "TENASCIN"), idents.use = "Bcell")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(3,4), signaling = c("PARs", "FN1", "VTN", "DHEAS", "FGF", "Cholesterol", "LAMININ", "COLLAGEN", "VEGF", "CypA", "27HC", "BMP", "GALECTIN", "NRG", "SPP1", "IGF", "DHEA", "ANGPTL", "AGT", "HGF", "COMPLEMENT", "EGF", "PROS", "Testosterone", "TENASCIN"), idents.use = "Bcell", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Bcell WD.t vs. RD.t")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(4,5), signaling = c("PARs", "FN1", "VTN", "DHEAS", "GALECTIN", "COLLAGEN", "BMP", "PROS", "LAMININ", "Cholesterol", "27HC", "FGF", "CypA", "VEGF", "ANGPTL", "SPP1", "NRG", "IGF"), idents.use = "Bcell")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(4,5), signaling = c("PARs", "FN1", "VTN", "DHEAS", "GALECTIN", "COLLAGEN", "BMP", "PROS", "LAMININ", "Cholesterol", "27HC", "FGF", "CypA", "VEGF", "ANGPTL", "SPP1", "NRG", "IGF"), idents.use = "Bcell", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Bcell RD.t vs. RD.n")
patchwork::wrap_plots(plots = gg1)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

##this section just to compare CDvPreT analysis
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), idents.use = "Tcell")
patchwork::wrap_plots(plots = gg1)

##try and make images smaller 
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), idents.use = "LSEC1", show.legend = F)
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), idents.use = "LSEC1", show.legend = F, font.size = 8, font.size.title = 8, label.size = 2)
patchwork::wrap_plots(plots = gg1)

#!options(ggrepel.max.overlaps = Inf)
##just adjust this code for any cell types of interest to make smaller plots
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), idents.use = "LSEC1", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("LSEC1 CD vs. PreT")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), idents.use = "LSEC1", xlims = c(-0.5, 0.5), ylims = c(-1.5, 1.5), show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("LSEC1 CD vs. PreT")
patchwork::wrap_plots(plots = gg1)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##for CDvPreT analysis need to save all pathways that are enriched in flow info p-value < 0.01 analysis for every cell type
##save a larger and smaller size plot over the previous files so there is no confusion
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("IGF", "PARs", "RELN", "FN1", "VTN", "ANGPTL", "GRN", "COMPLEMENT", "LAMININ", "COLLAGEN", "GALECTIN", "SLIT", "EGF", "SEMA3", "SPP1"), idents.use = "Stromal")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("IGF", "PARs", "RELN", "FN1", "VTN", "ANGPTL", "GRN", "COMPLEMENT", "LAMININ", "COLLAGEN", "GALECTIN", "SLIT", "EGF", "SEMA3", "SPP1"), idents.use = "Stromal", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Stromal CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("SEMA3", "FN1", "VTN", "ANGPTL", "COMPLEMENT", "PARs", "HGF", "TENASCIN", "FGF", "EGF", "GALECTIN"), idents.use = "Cancer")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("SEMA3", "FN1", "VTN", "ANGPTL", "COMPLEMENT", "PARs", "HGF", "TENASCIN", "FGF", "EGF", "GALECTIN"), idents.use = "Cancer", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Cancer CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("CSF", "RELN", "LAMININ", "FN1", "COLLAGEN", "VTN", "GRN", "COMPLEMENT", "PARs", "BMP", "IL1", "NRG", "ANGPTL", "GALECTIN", "SLIT", "EGF"), idents.use = "HSC")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("CSF", "RELN", "LAMININ", "FN1", "COLLAGEN", "VTN", "GRN", "COMPLEMENT", "PARs", "BMP", "IL1", "NRG", "ANGPTL", "GALECTIN", "SLIT", "EGF"), idents.use = "HSC", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("HSC CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "VTN", "LAMININ", "RELN", "NT", "ANGPTL", "GRN", "PERIOSTIN", "COMPLEMENT", "PARs", "GALECTIN", "COLLAGEN", "ANGPT", "GDF", "VWF", "PTN", "EGF", "TGFb", "SLIT"), idents.use = "Fibroblast")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "VTN", "LAMININ", "RELN", "NT", "ANGPTL", "GRN", "PERIOSTIN", "COMPLEMENT", "PARs", "GALECTIN", "COLLAGEN", "ANGPT", "GDF", "VWF", "PTN", "EGF", "TGFb", "SLIT"), idents.use = "Fibroblast", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Fibroblast CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("COLLAGEN", "NRG", "FGF", "SEMA3", "VTN", "FN1", "KIT", "LAMININ", "PDGF", "VISFATIN", "COMPLEMENT", "EGF", "CCL", "PTN", "ANGPTL", "GALECTIN", "IGF", "AGRN", "HGF", "PARs", "ncWNT", "PERIOSTIN", "SLIT", "GRN", "SPP1"), idents.use = "Myofibroblast")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("COLLAGEN", "NRG", "FGF", "SEMA3", "VTN", "FN1", "KIT", "LAMININ", "PDGF", "VISFATIN", "COMPLEMENT", "EGF", "CCL", "PTN", "ANGPTL", "GALECTIN", "IGF", "AGRN", "HGF", "PARs", "ncWNT", "PERIOSTIN", "SLIT", "GRN", "SPP1"), idents.use = "Myofibroblast", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Myofibroblast CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("PROS", "PARs", "ANGPTL", "FN1", "GRN", "EGF", "COMPLEMENT", "SLIT", "COLLAGEN", "LAMININ", "GALECTIN", "PTN"), idents.use = "Hepatocyte")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("PROS", "PARs", "ANGPTL", "FN1", "GRN", "EGF", "COMPLEMENT", "SLIT", "COLLAGEN", "LAMININ", "GALECTIN", "PTN"), idents.use = "Hepatocyte", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Hepatocyte CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("GAS", "VTN", "FN1", "PARs", "LAMININ", "VISFATIN", "COMPLEMENT", "EGF", "CCL", "PTN", "ANGPTL", "GALECTIN", "IGF", "COLLAGEN", "ncWNT", "PERIOSTIN", "SLIT", "GRN"), idents.use = "Cholangiocyte")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("GAS", "VTN", "FN1", "PARs", "LAMININ", "VISFATIN", "COMPLEMENT", "EGF", "CCL", "PTN", "ANGPTL", "GALECTIN", "IGF", "COLLAGEN", "ncWNT", "PERIOSTIN", "SLIT", "GRN"), idents.use = "Cholangiocyte", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Cholangiocyte CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("LAMININ", "FN1", "VTN", "COLLAGEN", "GRN", "ANGPTL", "COMPLEMENT", "PARs", "GALECTIN", "NRG", "GDF", "PTN", "TGFb", "SLIT"), idents.use = "LSEC1")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("LAMININ", "FN1", "VTN", "COLLAGEN", "GRN", "ANGPTL", "COMPLEMENT", "PARs", "GALECTIN", "NRG", "GDF", "PTN", "TGFb", "SLIT"), idents.use = "LSEC1", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("LSEC1 CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("LAMININ", "FN1", "VTN", "VWF", "ANGPTL", "COMPLEMENT", "PARs", "GALECTIN", "COLLAGEN", "TENASCIN", "NRG", "GDF", "SLIT"), idents.use = "Endothelial")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("LAMININ", "FN1", "VTN", "VWF", "ANGPTL", "COMPLEMENT", "PARs", "GALECTIN", "COLLAGEN", "TENASCIN", "NRG", "GDF", "SLIT"), idents.use = "Endothelial", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Endothelial CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("VTN", "COMPLEMENT", "NRG", "GALECTIN", "GDF"), idents.use = "Bcell")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("VTN", "COMPLEMENT", "NRG", "GALECTIN", "GDF"), idents.use = "Bcell", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Bcell CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "VTN", "ANGPTL", "COMPLEMENT", "PARs", "NRG", "GALECTIN", "GDF"), idents.use = "Tcell")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "VTN", "ANGPTL", "COMPLEMENT", "PARs", "NRG", "GALECTIN", "GDF"), idents.use = "Tcell", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Tcell CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("IGF", "FN1", "VTN", "ANGPTL", "GRN", "TGFb", "COMPLEMENT", "PARs", "LAMININ", "RELN", "COLLAGEN", "GALECTIN", "EGF", "CCL"), idents.use = "DC")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("IGF", "FN1", "VTN", "ANGPTL", "GRN", "TGFb", "COMPLEMENT", "PARs", "LAMININ", "RELN", "COLLAGEN", "GALECTIN", "EGF", "CCL"), idents.use = "DC", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("DC CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "VTN", "PARs", "ANGPTL", "VISFATIN", "GRN", "TGFb", "COMPLEMENT", "TENASCIN", "LAMININ", "COLLAGEN", "GALECTIN", "PTN", "GDF", "PERIOSTIN", "CCL"), idents.use = "Macrophage")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "VTN", "PARs", "ANGPTL", "VISFATIN", "GRN", "TGFb", "COMPLEMENT", "TENASCIN", "LAMININ", "COLLAGEN", "GALECTIN", "PTN", "GDF", "PERIOSTIN", "CCL"), idents.use = "Macrophage", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Macrophage CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("VTN", "FN1", "TGFb", "COMPLEMENT", "PARs", "GALECTIN", "ANGPTL", "GDF", "PERIOSTIN", "COLLAGEN", "LAMININ"), idents.use = "Monocyte")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("VTN", "FN1", "TGFb", "COMPLEMENT", "PARs", "GALECTIN", "ANGPTL", "GDF", "PERIOSTIN", "COLLAGEN", "LAMININ"), idents.use = "Monocyte", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Monocyte CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "ncWNT", "GRN", "VISFATIN", "GDF", "CCL", "COMPLEMENT", "EGF", "LAMININ", "ANGPTL", "COLLAGEN", "VTN", "PARs", "GALECTIN", "PERIOSTIN", "PTN", "SLIT", "TGFb"), idents.use = "NKT")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "ncWNT", "GRN", "VISFATIN", "GDF", "CCL", "COMPLEMENT", "EGF", "LAMININ", "ANGPTL", "COLLAGEN", "VTN", "PARs", "GALECTIN", "PERIOSTIN", "PTN", "SLIT", "TGFb"), idents.use = "NKT", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("NKT CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "PARs", "VTN", "ANGPTL", "CCL", "COMPLEMENT", "COLLAGEN", "GALECTIN"), idents.use = "NK")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("FN1", "PARs", "VTN", "ANGPTL", "CCL", "COMPLEMENT", "COLLAGEN", "GALECTIN"), idents.use = "NK", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("NK CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("IGF", "VTN", "ANGPTL", "COMPLEMENT", "PARs", "FN1", "LAMININ", "COLLAGEN", "GALECTIN"), idents.use = "Neutrophil")
patchwork::wrap_plots(plots = gg1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.CDvPreT, comparison = c(1,2), signaling = c("IGF", "VTN", "ANGPTL", "COMPLEMENT", "PARs", "FN1", "LAMININ", "COLLAGEN", "GALECTIN"), idents.use = "Neutrophil", show.legend = FALSE, xlabel = "Outgoing", ylabel = "Incoming", font.size = 8, font.size.title = 8, label.size = 2) + ggtitle("Neutrophil CD vs. PreT")
patchwork::wrap_plots(plots = gg1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

##try to identify conserved and context-specific signaling pathways 
##functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")

Y <- methods::slot(cellchat, "net")$similarity[[type]]$dr[[comparison.name]]

pathways.ignore <- rownames( Y[rowSums(!is.finite(Y))>0, ] )
cellchat@options$pathways.ignore = pathways.ignore
Y <- Y[!rowSums(!is.finite(Y)),] # filter out rows with NaN, not working downstream
methods::slot(cellchat, "net")$similarity[[type]]$dr[[comparison.name]] <- Y
data.use <- Y

cellchat <- netClustering(cellchat, type = "functional", comparison = c(1,2,3,4,5))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
##structural similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")

##identify conserved and context-specific signaling pathways
##we can then compare overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), stacked = F, do.stat = TRUE)
gg1 + gg2

##altering this to include specific sources enables us to look at specific cell types
##ex. we could look at stromal cells as a source or as a receiver of signal
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), sources.use = "Stromal", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), sources.use = "Stromal", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), targets.use = "Stromal", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), targets.use = "Stromal", stacked = F, do.stat = TRUE)
gg1 + gg2

##for LSEC1 and LSEC2 cell types
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), sources.use = "LSEC1", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), sources.use = "LSEC1", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), targets.use = "LSEC1", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), targets.use = "LSEC1", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), sources.use = "LSEC1", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), sources.use = "LSEC1", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), targets.use = "LSEC1", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,5), targets.use = "LSEC1", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), sources.use = "LSEC2", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), sources.use = "LSEC2", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), targets.use = "LSEC2", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), targets.use = "LSEC2", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), sources.use = "LSEC2", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), sources.use = "LSEC2", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), targets.use = "LSEC2", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(3,4,5), targets.use = "LSEC2", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

##use these section for specific immune cells (just adjust name in code for cell types)
##cell types to assess are Tcell, Bcell, DC, NK, NKT, Neutro, Mono, and Macrophages
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), sources.use = "Macrophage", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), sources.use = "Macrophage", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), targets.use = "Macrophage", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), targets.use = "Macrophage", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), sources.use = "Macrophage", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), sources.use = "Macrophage", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), targets.use = "Macrophage", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), targets.use = "Macrophage", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

##perform only CDvPreT rankNET analyses of interest in default (25% trimean) analysis with LSEC1 and LSEC2
##this one piece of code below found from (https://github.com/sqjin/CellChat/issues/164) and seemed to fix issues with LSEC1 and LSEC2 idents that were occurring
cellchat.CDvPreT <- liftCellChat(cellchat.CDvPreT , group.new = levels(cellchat.CDvPreT@idents$joint))

gg1 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), sources.use = "Tcell", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), sources.use = "Tcell", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), targets.use = "Tcell", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), targets.use = "Tcell", stacked = F, do.stat = TRUE)
gg1 + gg2
##filter these on pvalues too
#!gg1 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), sources.use = "Tcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
#!gg2 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), sources.use = "Tcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
#!gg1 + gg2

#!gg1 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), targets.use = "Tcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
#!gg2 <- rankNet(cellchat.CDvPreT, mode = "comparison", comparison = c(1,2), targets.use = "Tcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
#!gg1 + gg2

##do PreTvMWD.T analyses for flow info filtered on pvalue
#!gg1 <- rankNet(cellchat.PreTvMWD.T, mode = "comparison", comparison = c(1,2), sources.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
#!gg2 <- rankNet(cellchat.PreTvMWD.T, mode = "comparison", comparison = c(1,2), sources.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
#!gg1 + gg2

#!gg1 <- rankNet(cellchat.PreTvMWD.T, mode = "comparison", comparison = c(1,2), targets.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
#!gg2 <- rankNet(cellchat.PreTvMWD.T, mode = "comparison", comparison = c(1,2), targets.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
#!gg1 + gg2

##perform only MRD.NTvMRD.T rankNET analyses of interest in default (25% trimean) analysis with LSEC1 and LSEC2
cellchat <- liftCellChat(cellchat , group.new = levels(cellchat@idents$joint))
##also utilize this for p-value filtered all cell types source/target in default no LSEC split due to levels issue
##filter these on pvalues too
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2), sources.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2), sources.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2), targets.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2), targets.use = "Bcell", cutoff.pvalue = 0.01, thresh = 0.01, stacked = F, do.stat = TRUE)
gg1 + gg2

##now we can compare outgoing or incoming signals associated with each cell population
library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, y = c(object.list[[i+1]]@netP$pathways, object.list[[i+2]]@netP$pathways, object.list[[i+3]]@netP$pathways, object.list[[i+4]]@netP$pathways))
#!ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
#!ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
#!draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

##incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 10, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 10, color.heatmap = "GnBu")
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+4], width = 5, height = 10, color.heatmap = "GnBu")

draw(ht1 + ht2 + ht3 + ht4 + ht5, ht_gap = unit(0.5, "cm"))

##outgoing
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 10, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 10, color.heatmap = "GnBu")
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+4], width = 5, height = 10, color.heatmap = "GnBu")

draw(ht1 + ht2 + ht3 + ht4 + ht5, ht_gap = unit(0.5, "cm"))

##overall
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 10, color.heatmap = "OrRd")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 10, color.heatmap = "OrRd")
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+4], width = 5, height = 10, color.heatmap = "OrRd")

draw(ht1 + ht2 + ht3 + ht4 + ht5, ht_gap = unit(0.5, "cm"))


##try to identify dysfunctional signaling 
##find upregulated and downregulated signals

netVisual_bubble(cellchat, signaling = "CSF",  comparison = c(1, 2, 3, 4, 5), angle.x = 45)

##try to visualize specific pathways of interest
pathways.show <- c("CSF") 

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

##or can visualize with heatmaps
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]] + ht[[4]] + ht[[5]], ht_gap = unit(0.5, "cm"))

##try to visualize specific signaling pathways for each group in merged cellchat object
group.cellType <- c(rep("Stromal"), rep("Tcell"), rep("NKT"), rep("NK"), rep("Neutrophil"), rep("Monocyte"), rep("Macrophage"), rep("LSEC"), rep("HSC"), rep("Hepatocyte"), rep("Fibroblast"), rep("Myofibroblast"), rep("Endothelial"), rep("DC"), rep("Cholangiocyte"), rep("Cancer"), rep("Bcell"))
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
##this is a bit difficult for those that are NOT shared by all groups
##likely easier to use each saved CellChat group object and pick pathways from each group one by one
##can also quantify all interactions and the cell type across groups with this
plotGeneExpression(cellchat, signaling = "FN1", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "HGF", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "IGF", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "PARs", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "VTN", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "FGF", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "NRG", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "ANGPTL", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "BMP", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "COLLAGEN", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "COMPLEMENT", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "LAMININ", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "GRN", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "EGF", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "SEMA3", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "VEGF", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "VISFATIN", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "CSF", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "CCL", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "GAS", split.by = "group", type = c("violin"), colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "TGFb", split.by = "group", type = c("violin"), colors.ggplot = T)

###############################################################################################################
##from each saved group we loaded in try to extract the required pathways for stromal cells sending and receiving signals

##first for CD cellchat analysis 
##stromal specific pathways first
netVisual_chord_cell(cellchat.CD, signaling = "FN1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "FN1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "LAMININ", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "LAMININ", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "COLLAGEN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "COLLAGEN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "HGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "HGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "PARs", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "PARs", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "RELN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "RELN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "VTN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "VTN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "Cholesterol", lab.cex = 1, thresh = 0.01, cell.order = c("Hepatocyte", "Cancer", "Cholangiocyte", "Tcell", "Myofibroblast", "DC", "Fibroblast", "NKT", "HSC", "NK", "Stromal", "Neutrophil", "LSEC", "Monocyte", "Macrophage", "Endothelial"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "Cholesterol", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "FGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "FGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "NRG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "KIT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "KIT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "HSPG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "HSPG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "PDGF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "PDGF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "BMP", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "BMP", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "PROS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "PROS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "VEGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "SEMA3", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "SEMA3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "DHEAS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "DHEAS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "DHEA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "DHEA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "PROS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "PROS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "ANGPTL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "ANGPTL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "Testosterone", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "Testosterone", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "TENASCIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "TENASCIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "GAS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "GAS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##CD immune related signals in 5% analysis
netVisual_chord_cell(cellchat.CD, signaling = "APRIL", cell.order = c("Fibroblast", "DC", "Bcell"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "APRIL", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "CD40", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "CD40", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "EDN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "EDN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "FASLG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "FASLG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "IL16", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "IL16", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "IL2", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "IL2", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "TRAIL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "TRAIL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "IL1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "IL1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "FLT3", thresh = 0.01, cell.order = c("NKT", "LSEC", "NK", "Stromal", "HSC", "Myofibroblast", "Tcell", "Hepatocyte", "DC"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "FLT3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "CCL", thresh = 0.01, cell.order = c("NK", "Cancer", "Neutrophil", "Stromal", "Macrophage", "Monocyte", "Hepatocyte", "Endothelial", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Bcell", "NKT", "Tcell", "LSEC", "DC", "HSC"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "CCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "CXCL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "CXCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "CSF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "CSF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "MIF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "MIF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "GRN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "GRN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "VISFATIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "VISFATIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "EPO", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "EPO", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "TGFb", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "TGFb", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "27HC", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "27HC", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "CypA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "CypA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "Glutamate", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "Glutamate", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "NRG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##second for PreT cellchat
##stromal specific pathways first
netVisual_chord_cell(cellchat.PreT, signaling = "ANGPTL", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "LSEC", "Stromal", "NK", "HSC", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "DC", "Endothelial", "Myofibroblast", "Monocyte", "Cholangiocyte", "Tcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "ANGPTL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "BMP", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "BMP", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "COLLAGEN", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "COLLAGEN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "COMPLEMENT", font.size = 16, thresh = 0.01, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "FN1", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "FN1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "GRN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "GRN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "HGF", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "HGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IGF", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "LAMININ", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "LAMININ", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PARs", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PARs", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "VTN", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "VTN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "EGF", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "EGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "Cholesterol", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "Cholesterol", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "FGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "FGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "VEGF", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "VISFATIN", lab.cex = 1, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "VISFATIN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "SEMA3", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "SEMA3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "LSEC", "Stromal", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "RELN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "RELN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "SPP1", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "SPP1", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "SLIT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "SLIT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "VWF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "VWF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "HSPG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "HSPG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PDGF", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PDGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "EGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "EGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "VEGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "NRG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PROS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PROS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "BMP", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "BMP", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "DHEAS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "DHEAS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "DHEA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "DHEA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "Testosterone", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "Testosterone", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##PreT immune related signals from 5% analysis
netVisual_chord_cell(cellchat.PreT, signaling = "BAFF", cell.order = c("HSC", "DC", "Bcell"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "BAFF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "CX3C", thresh = 0.01, cell.order = c("Stromal", "Endothelial", "LSEC", "Macrophage", "NK", "NKT", "Cancer", "Myofibroblast", "Cholangiocyte", "Monocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "CX3C", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "EDN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "EDN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "FASLG", thresh = 0.01, cell.order = c("NK", "Macrophage", "Neutrophil", "Monocyte", "Endothelial", "Stromal", "LSEC", "HSC", "Myofibroblast", "Cholangiocyte", "Cancer", "Hepatocyte", "Tcell", "NKT"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "FASLG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IL12", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IL12", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IL16", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IL16", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IL2", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IL2", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "LIFR", thresh = 0.01, cell.order = c("NKT", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "Bcell", "DC", "Tcell"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "LIFR", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "TNF", thresh = 0.01, cell.order = c("Monocyte", "Bcell", "DC", "Tcell", "Macrophage", "Hepatocyte", "Cancer", "Stromal", "LSEC", "Endothelial", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte", "NKT", "NK", "Neutrophil"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "TNF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "TRAIL", thresh = 0.01, cell.order = c("Fibroblast", "Cholangiocyte", "Endothelial", "LSEC", "Cancer", "NK", "Tcell", "Stromal"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "TRAIL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "TWEAK", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "TWEAK", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IL1", thresh = 0.01, cell.order = c("Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "DC", "NKT", "Tcell", "Neutrophil", "NK", "Monocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IL1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "FLT3", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "FLT3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "CCL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "CCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "CXCL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "CXCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "CSF", thresh = 0.01, cell.order = c("HSC", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "Monocyte", "Bcell", "Neutrophil", "DC", "Tcell", "Macrophage", "NK", "Stromal", "NKT", "LSEC", "Endothelial", "Fibroblast"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "CSF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "SAA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "SAA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "GRN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "GRN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "VISFATIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "VISFATIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "TGFb", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "TGFb", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "27HC", thresh = 0.01, cell.order = c("Cancer", "Stromal", "HSC", "LSEC", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "NKT"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "27HC", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "Glutamate", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "Glutamate", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "NRG", lab.cex = 1, cell.order = c("Bcell", "DC", "Tcell", "NKT", "NK", "Endothelial", "Neutrophil", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Macrophage", "Cholangiocyte", "Hepatocyte", "Monocyte", "Cancer"), thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##third for MWD.T cellchat
##stromal specific pathways first
netVisual_chord_cell(cellchat.MWD.T, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "ANGPT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "ANGPT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "ANGPTL", lab.cex = 1, thresh = 0.01, cell.order = c("Hepatocyte", "DC", "Cancer", "Bcell", "Tcell", "NKT", "NK", "Monocyte", "Neutrophil", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "ANGPTL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "BMP", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "BMP", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "COLLAGEN", thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Stromal", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "Endothelial", "Myofibroblast", "Monocyte", "Cholangiocyte", "Bcell"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "COLLAGEN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "CSF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CSF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "CXCL", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CXCL", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "EGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "EGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "FGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "FGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "FN1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "FN1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "GAS", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "GAS", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "GDF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "GDF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "HGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "HGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "LAMININ", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "LAMININ", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "MIF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "MIF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "NT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "NT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "PARs", lab.cex = 1, thresh = 0.01, cell.order = c("Hepatocyte", "Cancer", "Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "HSC", "LSEC", "Fibroblast", "Endothelial", "Myofibroblast", "Cholangiocyte", "Stromal"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "PARs", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "PDGF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "PDGF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "PROS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "PROS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "SEMA3", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "SEMA3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TENASCIN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TENASCIN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TGFb", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TGFb", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "VEGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "SLIT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "SLIT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "VTN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "VTN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "HSPG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "HSPG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "EGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "EGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "NRG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "Cholesterol", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "Cholesterol", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "DHEAS", lab.cex = 1, thresh = 0.01, cell.order = c("Cancer", "NKT", "Monocyte", "Cholangiocyte", "Macrophage", "Endothelial", "Fibroblast", "LSEC", "Myofibroblast", "Stromal", "Hepatocyte", "HSC"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "DHEAS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "DHEA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "DHEA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "Testosterone", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "Testosterone", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TENASCIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TENASCIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##MWD.T immune related in 5% analysis
netVisual_chord_cell(cellchat.MWD.T, signaling = "BAFF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "BAFF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "EDN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "EDN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "FASLG", thresh = 0.01, cell.order = c("Tcell", "HSC", "LSEC", "Endothelial", "Stromal", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "Macrophage", "Monocyte", "NKT"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "FASLG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IFN-I", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IFN-I", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IL10", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IL10", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IL16", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IL16", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IL2", thresh = 0.01, cell.order = c("Tcell", "NKT", "Hepatocyte", "Cancer", "Cholangiocyte", "HSC", "Fibroblast", "Stromal", "Myofibroblast", "LSEC", "Endothelial", "Monocyte", "Macrophage"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IL2", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "LIFR", thresh = 0.01, cell.order = c("HSC", "Endothelial", "LSEC", "Stromal", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Cancer", "NKT", "Monocyte", "Macrophage"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "LIFR", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TNF", thresh = 0.01, cell.order = c("Macrophage", "Cholangiocyte", "Hepatocyte", "Fibroblast", "Myofibroblast", "HSC", "Stromal", "LSEC", "Endothelial", "Cancer", "Bcell", "Tcell", "DC", "NK", "NKT", "Neutrophil", "Monocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TNF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TRAIL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TRAIL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TWEAK", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TWEAK", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IL1", thresh = 0.01, cell.order = c("Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "Tcell", "NKT", "NK", "Neutrophil", "DC", "Monocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IL1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "CCL", thresh = 0.01, cell.order = c("Hepatocyte", "Cholangiocyte", "Cancer", "Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Endothelial", "Macrophage", "LSEC", "Fibroblast", "Stromal", "Myofibroblast", "HSC"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "CXCL", thresh = 0.01, cell.order = c("Cancer", "Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Fibroblast", "Macrophage", "Endothelial", "LSEC", "Stromal", "Myofibroblast", "HSC", "Cholangiocyte", "Hepatocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CXCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "CSF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CSF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "MIF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "MIF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "SAA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "SAA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "GRN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "GRN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "VISFATIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "VISFATIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "EPO", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "EPO", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "TGFb", thresh = 0.01, cell.order = c("LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "Bcell", "Tcell", "NKT", "NK", "Monocyte", "DC", "Macrophage", "Endothelial"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "TGFb", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "27HC", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "27HC", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "CypA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CypA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "Glutamate", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "Glutamate", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "NRG", lab.cex = 1, cell.order = c("Cancer", "DC", "NKT", "Macrophage", "Monocyte", "Endothelial", "Fibroblast", "HSC", "Myofibroblast", "LSEC", "Hepatocyte", "Stromal", "Cholangiocyte"), thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##fourth for MRD.T cellchat
##stromal specific pathways first
netVisual_chord_cell(cellchat.MRD.T, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "ANGPTL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "ANGPTL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "BMP", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "BMP", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "CCL", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "CCL", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "COLLAGEN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "COLLAGEN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "FGF", thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Bcell", "LSEC", "Endothelial", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "NK", "Myofibroblast", "DC", "Monocyte", "Cholangiocyte", "Stromal"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "FGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "FN1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "FN1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "GALECTIN", thresh = 0.01, cell.order = c("NKT", "Hepatocyte", "Bcell", "LSEC", "Endothelial", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "DC", "NK", "Myofibroblast", "Monocyte", "Cholangiocyte", "Stromal"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "GAS", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "GAS", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "HGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "HGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "LAMININ", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "LAMININ", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "NRG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "PARs", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "PARs", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "PROS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "PROS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "SEMA3", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "SEMA3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "VEGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "VTN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "VTN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "Cholesterol", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "Cholesterol", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "DHEAS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "DHEAS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "DHEA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "DHEA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "Testosterone", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "Testosterone", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##MRD.T immune related pathways in 5% analysis
netVisual_chord_cell(cellchat.MRD.T, signaling = "EDN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "EDN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "FASLG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "FASLG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "IFN-I", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "IFN-I", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "IL16", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "IL16", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "IL2", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "IL2", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "LIFR", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "LIFR", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "TWEAK", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "TWEAK", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "IL1", thresh = 0.01, cell.order = c("Bcell", "DC", "Tcell", "Neutrophil", "NK", "Monocyte", "Macrophage", "Endothelial", "LSEC", "Stromal", "Cholangiocyte", "HSC", "Fibroblast", "Myofibroblast", "Hepatocyte", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "IL1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "CCL", thresh = 0.01, cell.order = c("Tcell", "Endothelial", "DC", "NK", "Neutrophil", "Monocyte", "Macrophage", "LSEC", "Stromal", "Fibroblast", "Cholangiocyte", "Cancer", "Hepatocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "CCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "CXCL", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "CXCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "CSF", thresh = 0.01, cell.order = c("HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "DC", "Bcell", "Monocyte", "Macrophage", "Endothelial", "Stromal"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "CSF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "MIF", lab.cex = 1, thresh = 0.01, cell.order = c("Hepatocyte", "Cancer", "Bcell", "Tcell", "DC", "NK", "Neutrophil", "Monocyte", "Macrophage", "HSC", "Endothelial", "LSEC", "Stromal", "Fibroblast", "Myofibroblast", "Cholangiocyte"), remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "MIF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "SAA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "SAA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "IGF", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "GRN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "GRN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "VISFATIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "VISFATIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "TGFb", thresh = 0.01, cell.order = c("Bcell", "Tcell", "DC", "Neutrophil", "Monocyte", "Macrophage", "Cancer", "Endothelial", "LSEC", "NK", "Stromal", "HSC", "Fibroblast", "Hepatocyte", "Myofibroblast", "Cholangiocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "TGFb", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "27HC", thresh = 0.01, cell.order = c("Bcell", "Tcell", "DC", "Neutrophil", "Monocyte", "NK", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "27HC", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "CypA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "CypA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "Glutamate", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "Glutamate", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "NRG", lab.cex = 1, cell.order = c("Cancer", "DC", "Macrophage", "Monocyte", "Endothelial", "LSEC", "Stromal", "HSC", "Myofibroblast", "Cholangiocyte", "Hepatocyte"), thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.T, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.T, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##lastly for MRD.NT cellchat
##stromal specific pathways first
netVisual_chord_cell(cellchat.MRD.NT, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "ANGPT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "ANGPT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "ANGPTL", thresh = 0.01, cell.order = c("NKT", "Bcell", "LSEC", "Endothelial", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "DC", "NK", "Myofibroblast", "Cholangiocyte", "Monocyte", "Stromal", "Hepatocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "ANGPTL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "BMP", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "BMP", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "COLLAGEN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "COLLAGEN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "FGF", thresh = 0.01, cell.order = c("NKT", "LSEC", "HSC", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "Macrophage", "DC", "NK", "Myofibroblast", "Cholangiocyte", "Endothelial", "Monocyte", "Stromal", "Bcell", "Hepatocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "FGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "FN1", thresh = 0.01, cell.order = c("LSEC", "NKT", "Tcell", "Fibroblast", "Cancer", "Neutrophil", "HSC", "Macrophage", "DC", "Cholangiocyte", "NK", "Myofibroblast", "Endothelial", "Monocyte", "Stromal", "Bcell", "Hepatocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "FN1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "GAS", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "GAS", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "HGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "HGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "LAMININ", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "LAMININ", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "PARs", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "PARs", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "PROS", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "PROS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "PDGF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "PDGF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "SAA", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "SAA", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "TGFb", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "TGFb", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "VEGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "VISFATIN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "VISFATIN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "VTN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "VTN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "VEGF", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "VEGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "NRG", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "SEMA3", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "SEMA3", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "Cholesterol", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "Cholesterol", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "DHEAS", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "DHEAS", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "DHEA", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "DHEA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "Testosterone", thresh = 0.01, lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "Testosterone", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##MRD.NT immune related pathways from 5% analysis
netVisual_chord_cell(cellchat.MRD.NT, signaling = "BAFF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "BAFF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "EDN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "EDN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "FASLG", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "FASLG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IFN-I", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IFN-I", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IL12", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IL12", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IL16", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IL16", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IL2", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IL2", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "LIFR", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "LIFR", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "TNF", thresh = 0.01, cell.order = c("Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer", "NKT", "Neutrophil", "Monocyte", "Macrophage", "Stromal", "LSEC", "Endothelial", "HSC", "Fibroblast"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "TNF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "TWEAK", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "TWEAK", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IL1", thresh = 0.01, cell.order = c("Tcell", "DC", "Neutrophil", "NKT", "Monocyte", "NK", "Macrophage", "Endothelial", "Stromal", "HSC", "Fibroblast", "LSEC", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IL1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "CCL", thresh = 0.01, cell.order = c("Tcell", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "LSEC", "Fibroblast", "Myofibroblast", "Stromal", "Hepatocyte", "HSC", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "CCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "CXCL", thresh = 0.01, cell.order = c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC", "HSC", "Stromal", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "CXCL", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "CSF", thresh = 0.01, cell.order = c("NKT", "Neutrophil", "NK", "Monocyte", "Macrophage", "Endothelial", "Stromal", "LSEC", "HSC", "Fibroblast", "Myofibroblast", "Cancer", "Cholangiocyte"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "CSF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "MIF", thresh = 0.01, cell.order = c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "MIF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "SAA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "SAA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IGF", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IGF", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "AGT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "AGT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "GRN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "GRN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "VISFATIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "VISFATIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "TGFb", thresh = 0.01, cell.order = c("Bcell", "Tcell", "NK", "Monocyte", "Neutrophil", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte", "Hepatocyte", "Cancer"), lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "TGFb", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "27HC", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "27HC", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "CypA", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "CypA", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "GALECTIN", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "GALECTIN", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "Glutamate", lab.cex = 1, cell.order = c("Cancer", "Bcell", "Monocyte", "DC", "Endothelial", "NKT", "Stromal", "Macrophage", "HSC", "Fibroblast", "Neutrophil", "Myofibroblast", "Cholangiocyte", "Hepatocyte"), thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "Glutamate", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "COMPLEMENT", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "COMPLEMENT", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "NRG", lab.cex = 1, cell.order = c("Endothelial", "Hepatocyte", "Neutrophil", "Cancer", "DC", "Monocyte", "Macrophage", "LSEC", "Fibroblast", "HSC", "Myofibroblast", "Stromal", "Cholangiocyte"), thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "NRG", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "SPP1", lab.cex = 1, thresh = 0.01, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "SPP1", thresh = 0.01, font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##from each saved group we loaded in extract the  pathways for cancer cells sending and receiving signals
##save all pathways that are specific for cancer that were not saved/shared with stromal cells

##first for CD cancer
netVisual_chord_cell(cellchat.CD, signaling = "COLLAGEN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "COLLAGEN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "LAMININ", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "LAMININ", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "SEMA3", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "SEMA3", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.CD, signaling = "TENASCIN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.CD, signaling = "TENASCIN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##second for PreT cancer
netVisual_chord_cell(cellchat.PreT, signaling = "AGT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "AGT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "CCL", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "CCL", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "GAS", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "GAS", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "GDF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "GDF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "ncWNT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "ncWNT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "NRG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "NRG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PDGF", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PDGF", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PROS", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PROS", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "AGRN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "AGRN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "ANGPT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "ANGPT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "HSPG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "HSPG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "IL1", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "IL1", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "NT", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "NT", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PERIOSTIN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PERIOSTIN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "PTN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "PTN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.PreT, signaling = "TGFb", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.PreT, signaling = "TGFb", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##third for MWD.T cancer
netVisual_chord_cell(cellchat.MWD.T, signaling = "CCL", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "CCL", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "NRG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "NRG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "AGRN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "AGRN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MWD.T, signaling = "HSPG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MWD.T, signaling = "HSPG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##fourth for MRD.T cancer
##ALL PATHS WERE SHARED/SAVED IN STROMAL SECTION

##lastly, for MRD.NT cancer
netVisual_chord_cell(cellchat.MRD.NT, signaling = "CCL", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "CCL", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "GALECTIN", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "GALECTIN", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "IL1", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "IL1", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

netVisual_chord_cell(cellchat.MRD.NT, signaling = "NRG", lab.cex = 1, remove.isolate = TRUE)
netAnalysis_contribution(cellchat.MRD.NT, signaling = "NRG", font.size = 16, font.size.title = 16) + theme(text = element_text(size=16))

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##set up analyses for looking for cells expressing genes of interest and afterwards perform DESEQ2 for IPA analysis
library(ggplot2)
library(Seurat)
library(scuttle)

HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.labels.RDS")
##use this file to get our unknown population from sorting results
#!HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.individually.w.scSorter5.RDS")

##use these files to directly look at LSEC1 and LSEC2 to get counts for IPA analyses
LSEC1 <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/LSEC1.saved.rds")
LSEC2 <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Output images/HCCp1/CellChat v2/Saved_Rfiles/LSEC2.saved.rds")
table(LSEC1@meta.data$group)
table(LSEC2@meta.data$group)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$group
table(HCCp1.cells@meta.data$group)

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.cells@meta.data$group <- factor(x = HCCp1.cells@meta.data$group, levels = new.group.levels)
table(HCCp1.cells@meta.data$group)
table(HCCp1.cells@meta.data$Cell.Type)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Cell.Type
#!Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Predicted_Type
#!table(HCCp1.cells@meta.data$Predicted_Type)
##VlnPlot(HCCp1.cells, features = "Il13", idents = "NKT")
##Il13 is not detected in our data

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##setup DESeq2 to extract individual files for all non-IC species individually and altogther
##this is to uplaod to IPA and trying to show non-ICs also exhibit different functions together or individually
##rational for this is some genes may have low level of expression in the individual population...
##but when evaluated as a whole, minor expression level may have a significant role overall in the non-IC population

##first extract cell type of interest only (just adjust code for each for subsequent analyses)
##non-ICs = subset(HCCp1.cells, idents = c(Cancer, Cholangiocyte, Endothelial, Fibroblasts, HSC, Hepatocyte, LSEC, Myofibroblast, and Stromal), invert = FALSE)
#!Non.ICs <- subset(HCCp1.cells, idents = c("Cancer", "Cholangiocyte", "Endothelial", "Fibroblast", "HSC", "Hepatocyte", "LSEC", "Myofibroblast", "Stromal"), invert = FALSE)
#!Cancer <- subset(HCCp1.cells, idents = "Cancer", invert = FALSE)
#!Cholangio <- subset(HCCp1.cells, idents = "Cholangiocyte", invert = FALSE)
#!Endo <- subset(HCCp1.cells, idents = "Endothelial", invert = FALSE)
#!Fibro <- subset(HCCp1.cells, idents = "Fibroblast", invert = FALSE)
#!HSC <- subset(HCCp1.cells, idents = "HSC", invert = FALSE)
#!Hep <- subset(HCCp1.cells, idents = "Hepatocyte", invert = FALSE)
#!LSEC <- subset(HCCp1.cells, idents = "LSEC", invert = FALSE)
#!Myofibro <- subset(HCCp1.cells, idents = "Myofibroblast", invert = FALSE)
#!Stromal <- subset(HCCp1.cells, idents = "Stromal", invert = FALSE)
all.Stromal <- subset(HCCp1.cells, idents = c("Stromal", "LSEC", "Endothelial"), invert = FALSE)

##now do the same for immune cells
#!ICs <- subset(HCCp1.cells, idents = c("Bcell", "Tcell", "DC", "Macrophage", "Monocyte", "NK", "NKT", "Neutrophil"), invert = FALSE)
#!Bcell <- subset(HCCp1.cells, idents = "Bcell", invert = FALSE)
#!Tcell <- subset(HCCp1.cells, idents = "Tcell", invert = FALSE)
#!NK.NKT <- subset(HCCp1.cells, idents = c("NK", "NKT"), invert = FALSE)
#!DC <- subset(HCCp1.cells, idents = "DC", invert = FALSE)
#!Mac <- subset(HCCp1.cells, idents = "Macrophage", invert = FALSE)
#!Mono <- subset(HCCp1.cells, idents = "Monocyte", invert = FALSE)
#!Neutro <- subset(HCCp1.cells, idents = "Neutrophil", invert = FALSE)
#!NK <- subset(HCCp1.cells, idents = "NK", invert = FALSE)
NKT <- subset(HCCp1.cells, idents = "NKT", invert = FALSE)

##last for any unknown cells
#!Unk <- subset(HCCp1.cells, idents = "Unknown", invert = FALSE)

table(NKT@meta.data$group)

##second separate each group after extracting the cell population of interest
Idents(NKT) <- NKT@meta.data$group
table(NKT@meta.data$group)

CD <- subset(NKT, idents = "CD", invert = FALSE)
PreT <- subset(NKT, idents = "PreT", invert = FALSE)
MWD.T <- subset(NKT, idents = "MWD.T", invert = FALSE)
MRD.T <- subset(NKT, idents = "MRD.T", invert = FALSE)
MRD.NT <- subset(NKT, idents = "MRD.NT", invert = FALSE)

table(CD@meta.data$Cell.Type)
table(PreT@meta.data$Cell.Type)
table(MWD.T@meta.data$Cell.Type)
table(MRD.T@meta.data$Cell.Type)
table(MRD.NT@meta.data$Cell.Type)

##first CD
CD$replicate <- unlist(lapply(names(CD$group), function(x) strsplit(x, split="_")[[1]][[1]]))
head(CD$replicate)
tail(CD$replicate)

CD.sce <- as.SingleCellExperiment(CD)
CD.sce <- aggregateAcrossCells(CD.sce, ids = colData(CD.sce)[,"replicate"])
head(assay(CD.sce))
colData(CD.sce)[,"replicate"]
write.csv(assay(CD.sce), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/CD/CD.NKT.counts.csv")
##second PreT
PreT$replicate <- unlist(lapply(names(PreT$group), function(x) strsplit(x, split="_")[[1]][[1]]))
head(PreT$replicate)
tail(PreT$replicate)

PreT.sce <- as.SingleCellExperiment(PreT)
PreT.sce <- aggregateAcrossCells(PreT.sce, ids = colData(PreT.sce)[,"replicate"])
head(assay(PreT.sce))
colData(PreT.sce)[,"replicate"]
write.csv(assay(PreT.sce), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/PreT/PreT.NKT.counts.csv")

##third MWD.T
MWD.T$replicate <- unlist(lapply(names(MWD.T$group), function(x) strsplit(x, split="_")[[1]][[1]]))
head(MWD.T$replicate)
tail(MWD.T$replicate)

MWD.T.sce <- as.SingleCellExperiment(MWD.T)
MWD.T.sce <- aggregateAcrossCells(MWD.T.sce, ids = colData(MWD.T.sce)[,"replicate"])
head(assay(MWD.T.sce))
colData(MWD.T.sce)[,"replicate"]
write.csv(assay(MWD.T.sce), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/MWD.T/MWD.T.NKT.counts.csv")

##fourth MRD.T
MRD.T$replicate <- unlist(lapply(names(MRD.T$group), function(x) strsplit(x, split="_")[[1]][[1]]))
head(MRD.T$replicate)
tail(MRD.T$replicate)

MRD.T.sce <- as.SingleCellExperiment(MRD.T)
MRD.T.sce <- aggregateAcrossCells(MRD.T.sce, ids = colData(MRD.T.sce)[,"replicate"])
head(assay(MRD.T.sce))
colData(MRD.T.sce)[,"replicate"]
write.csv(assay(MRD.T.sce), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/MRD.T/MRD.T.NKT.counts.csv")

##last MRD.NT
MRD.NT$replicate <- unlist(lapply(names(MRD.NT$group), function(x) strsplit(x, split="_")[[1]][[1]]))
head(MRD.NT$replicate)
tail(MRD.NT$replicate)

MRD.NT.sce <- as.SingleCellExperiment(MRD.NT)
MRD.NT.sce <- aggregateAcrossCells(MRD.NT.sce, ids = colData(MRD.NT.sce)[,"replicate"])
head(assay(MRD.NT.sce))
colData(MRD.NT.sce)[,"replicate"]
write.csv(assay(MRD.NT.sce), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/MRD.NT/MRD.NT.NKT.counts.csv")

##once we have desired files, we need to set up counts comparison and data.info files based on group comparisons 
##these is tedious and time consuming, but needs to be done for each analysis
##comparisons to make are to model different stages of disease progression
##first ones to do are CDvPreT, PreTvMWD.T, MWD.TvMRD.T, and MRD.TvMRD.NT
library(tidyverse)
library(DESeq2)

CD.PreT.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/CDvPreT.NKT.counts.csv", header = TRUE, sep = ",")
PreT.MWD.T.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/PreTvMWD.T.NKT.counts.csv", header = TRUE, sep = ",")
MWD.TvMRD.T.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/MWD.TvMRD.T.NKT.counts.csv", header = TRUE, sep = ",")
MRD.TvMRD.NT.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/MRD.TvMRD.NT.NKT.counts.csv", header = TRUE, sep = ",")

head(CD.PreT.counts)
head(PreT.MWD.T.counts)
head(MWD.TvMRD.T.counts)
head(MRD.TvMRD.NT.counts)

CD.PreT.counts$X[duplicated(CD.PreT.counts$X)]
PreT.MWD.T.counts$X[duplicated(PreT.MWD.T.counts$X)]
MWD.TvMRD.T.counts$X[duplicated(MWD.TvMRD.T.counts$X)]
MRD.TvMRD.NT.counts$X[duplicated(MRD.TvMRD.NT.counts$X)]

which(duplicated(CD.PreT.counts$X))
which(duplicated(PreT.MWD.T.counts$X))
which(duplicated(MWD.TvMRD.T.counts$X))
which(duplicated(MRD.TvMRD.NT.counts$X))

CD.PreT.counts$X[9519] <- "1-Mar_dup"
CD.PreT.counts$X[18441] <- "2-Mar_dup"
PreT.MWD.T.counts$X[9519] <- "1-Mar_dup"
PreT.MWD.T.counts$X[18441] <- "2-Mar_dup"
MWD.TvMRD.T.counts$X[9519] <- "1-Mar_dup"
MWD.TvMRD.T.counts$X[18441] <- "2-Mar_dup"
MRD.TvMRD.NT.counts$X[9519] <- "1-Mar_dup"
MRD.TvMRD.NT.counts$X[18441] <- "2-Mar_dup"

which(duplicated(CD.PreT.counts$X))
which(duplicated(PreT.MWD.T.counts$X))
which(duplicated(MWD.TvMRD.T.counts$X))
which(duplicated(MRD.TvMRD.NT.counts$X))

row.names(CD.PreT.counts) <- CD.PreT.counts$X
row.names(PreT.MWD.T.counts) <- PreT.MWD.T.counts$X
row.names(MWD.TvMRD.T.counts) <- MWD.TvMRD.T.counts$X
row.names(MRD.TvMRD.NT.counts) <- MRD.TvMRD.NT.counts$X

CD.PreT.counts <- CD.PreT.counts[,-1]
PreT.MWD.T.counts <- PreT.MWD.T.counts[,-1]
MWD.TvMRD.T.counts <- MWD.TvMRD.T.counts[,-1]
MRD.TvMRD.NT.counts <- MRD.TvMRD.NT.counts[,-1]

colData.CDvPreT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/CDvPreT.NKT.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.PreTvMWD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/PreTvMWD.T.NKT.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.MWD.TvMRD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/MWD.TvMRD.T.NKT.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.MRD.TvMRD.NT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/ICs/Comparisons/MRD.TvMRD.NT.NKT.data.info.csv", header = TRUE, sep = ",", row.names = 1)

View(colData.CDvPreT)
View(colData.PreTvMWD.T)
View(colData.MWD.TvMRD.T)
View(colData.MRD.TvMRD.NT)

all(names(CD.PreT.counts) %in% row.names(colData.CDvPreT))
all(names(PreT.MWD.T.counts) %in% row.names(colData.PreTvMWD.T))
all(names(MWD.TvMRD.T.counts) %in% row.names(colData.MWD.TvMRD.T))
all(names(MRD.TvMRD.NT.counts) %in% row.names(colData.MRD.TvMRD.NT))

all(names(CD.PreT.counts) == row.names(colData.CDvPreT))
all(names(PreT.MWD.T.counts) == row.names(colData.PreTvMWD.T))
all(names(MWD.TvMRD.T.counts) == row.names(colData.MWD.TvMRD.T))
all(names(MRD.TvMRD.NT.counts) == row.names(colData.MRD.TvMRD.NT))

##construct DESeqDataSet object
dds.CDvPreT <- DESeqDataSetFromMatrix(countData = CD.PreT.counts, colData = colData.CDvPreT, design = ~ condition)
dds.PreTvMWD.T <- DESeqDataSetFromMatrix(countData = PreT.MWD.T.counts, colData = colData.PreTvMWD.T, design = ~ condition)
dds.MWD.TvMRD.T <- DESeqDataSetFromMatrix(countData = MWD.TvMRD.T.counts, colData = colData.MWD.TvMRD.T, design = ~ condition)
dds.MRD.TvMRD.NT <- DESeqDataSetFromMatrix(countData = MRD.TvMRD.NT.counts, colData = colData.MRD.TvMRD.NT, design = ~ condition)

dds.CDvPreT
dds.PreTvMWD.T
dds.MWD.TvMRD.T
dds.MRD.TvMRD.NT

##pre-filtering to remove rows with low gene counts
keep.CDvPreT <- rowSums(counts(dds.CDvPreT)) >= 3
dds.CDvPreT <- dds.CDvPreT[keep.CDvPreT,]
dds.CDvPreT

keep.PreTvMWD.T <- rowSums(counts(dds.PreTvMWD.T)) >= 3
dds.PreTvMWD.T <- dds.PreTvMWD.T[keep.PreTvMWD.T,]
dds.PreTvMWD.T

keep.MWD.TvMRD.T <- rowSums(counts(dds.MWD.TvMRD.T)) >= 3
dds.MWD.TvMRD.T <- dds.MWD.TvMRD.T[keep.MWD.TvMRD.T,]
dds.MWD.TvMRD.T

keep.MRD.TvMRD.NT <- rowSums(counts(dds.MRD.TvMRD.NT)) >= 3
dds.MRD.TvMRD.NT <- dds.MRD.TvMRD.NT[keep.MRD.TvMRD.NT,]
dds.MRD.TvMRD.NT

##run DESeq
dds.CDvPreT <- DESeq(dds.CDvPreT)
dds.PreTvMWD.T <- DESeq(dds.PreTvMWD.T)
dds.MWD.TvMRD.T <- DESeq(dds.MWD.TvMRD.T)
dds.MRD.TvMRD.NT <- DESeq(dds.MRD.TvMRD.NT)

##Set up for the contrast function is the variable (condition), one to test (T1/high), and reference (C1/low)
res.CDvPreT <- results(dds.CDvPreT, contrast = c("condition", "PreT", "CD"))
res.PreTvMWD.T <- results(dds.PreTvMWD.T, contrast = c("condition", "MWD.T", "PreT"))
res.MWD.TvMRD.T <- results(dds.MWD.TvMRD.T, contrast = c("condition", "MRD.T", "MWD.T"))
res.MRD.TvMRD.NT <- results(dds.MRD.TvMRD.NT, contrast = c("condition", "MRD.NT", "MRD.T"))

res.CDvPreT
res.PreTvMWD.T
res.MWD.TvMRD.T
res.MRD.TvMRD.NT

write.csv(as.data.frame(res.CDvPreT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/ICs/CDvPreT.NKT.IPA.csv")
write.csv(as.data.frame(res.PreTvMWD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/ICs/PreTvMWD.T.NKT.IPA.csv")
write.csv(as.data.frame(res.MWD.TvMRD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/ICs/MWD.TvMRD.T.NKT.IPA.csv")
write.csv(as.data.frame(res.MRD.TvMRD.NT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/ICs/MRD.TvMRD.NT.NKT.IPA.csv")

#####################################################################################################################################################################################
#####################################################################################################################################################################################
##use this section to explore gene expression from CellChat analyses involving MRD.T group
#!HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.LSEC1and2.all.new.labels.RDS")

#!Idents(HCCp1.cells) <- HCCp1.cells@meta.data$group
#!table(HCCp1.cells@meta.data$group)

#!new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
#!HCCp1.cells@meta.data$group <- factor(x = HCCp1.cells@meta.data$group, levels = new.group.levels)
#!table(HCCp1.cells@meta.data$group)

##first need to look at gene expression per group...later can quantify # in replicates too for error bars/p-values
VlnPlot(HCCp1.cells, features = c("Lama1", "Lama2", "Lama4", "Lamb1", "Lamb3", "Lamc1", "Lamc3"), pt.size = 1, idents = c("PreT", "MWD.T"), split.by = "group", sort = "decreasing", group.by = "LSEC1and2", log = TRUE, stack = TRUE, same.y.lims = TRUE)
VlnPlot(HCCp1.cells, features = c("Cd44", "Dag1", "Itga1", "Itgb1", "Itga9", "Itgav", "Itgb8"), pt.size = 1, idents = c("PreT", "MWD.T"), split.by = "group", sort = "decreasing", group.by = "LSEC1and2", log = TRUE, stack = TRUE, same.y.lims = TRUE)

RidgePlot(HCCp1.cells, features = c("Lama1", "Lama2", "Lama4", "Lamb1", "Lamb3", "Lamc1", "Lamc3"), idents = "PreT", group.by = "LSEC1and2", sort = "decreasing", log = TRUE, stack = TRUE, same.y.lims = TRUE, fill.by = "ident")
RidgePlot(HCCp1.cells, features = c("Lama1", "Lama2", "Lama4", "Lamb1", "Lamb3", "Lamc1", "Lamc3"), idents = "MWD.T", group.by = "LSEC1and2", sort = "decreasing", log = TRUE, stack = TRUE, same.y.lims = TRUE, fill.by = "ident")

DotPlot(HCCp1.cells, features = c("Lama1", "Lama2", "Lama4", "Lamb1", "Lamb3", "Lamc1", "Lamc3"), cols = c("blue", "red"), idents = c("PreT", "MWD.T"), split.by = "group", group.by = "LSEC1and2")
DotPlot(HCCp1.cells, features = c("Igf1", "Igf1r", "Itgav", "Itgb3"), cols = c("blue", "red"), idents = c("PreT"), group.by = "LSEC1and2")
DotPlot(HCCp1.cells, features = c("Igf1", "Igf1r", "Itgav", "Itgb3"), cols = c("blue", "red"), idents = c("MWD.T"), group.by = "LSEC1and2")

VlnPlot(HCCp1.cells, features = c("Igf1", "Igf1r", "Itgav", "Itgb3"), pt.size = 1, idents = c("PreT", "MWD.T"), split.by = "group", sort = "decreasing", group.by = "LSEC1and2", log = TRUE, stack = TRUE, same.y.lims = TRUE)

##load in new.group.label data to investigate hepatocyte and cancer cell marker expression and extract replicate info
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.rds")
Idents(HCCp1.cells) <- HCCp1.cells@meta.data$new.group
table(HCCp1.cells@meta.data$new.group)

new.group.levels <- c("CD", "WD.nf", "WD.t", "RD.t", "RD.n")
HCCp1.cells@meta.data$new.group <- factor(x = HCCp1.cells@meta.data$new.group, levels = new.group.levels)
table(HCCp1.cells@meta.data$new.group)

new.cell.type.levels <- c("Bcell", "Tcell", "DC", "NKT", "NK", "Neutrophil", "Monocyte", "Macrophage", "Endothelial", "LSEC", "Stromal", "HSC", "Fibroblast", "Myofibroblast", "Cholangiocyte",  "Hepatocyte", "Cancer")

HCCp1.cells@meta.data$Cell.Type <- factor(x = HCCp1.cells@meta.data$Cell.Type, levels = new.cell.type.levels)
table(HCCp1.cells@meta.data$Cell.Type)

all.genes <- rownames(HCCp1.cells)
HCCp1.cells <- ScaleData(HCCp1.cells, features = all.genes)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Cell.Type
Hep <- subset(HCCp1.cells, idents = "Hepatocyte", invert = FALSE)
table(Hep@meta.data$replicate)

DoHeatmap(Hep, features = c("Hnf4a", "Afp", "Gpc3"), group.by = "new.group", size = 3)
DoHeatmap(Hep, features = c("Hnf4a", "Afp", "Gpc3"), group.by = "new.group", size = 3) + NoLegend()

Hep.Hnf4a <- subset(Hep, subset = Hnf4a > 0)
table(Hep.Hnf4a@meta.data$replicate)
Hep.Afp <- subset(Hep, subset = Afp > 0)
table(Hep.Afp@meta.data$replicate)
Hep.Gpc3 <- subset(Hep, subset = Gpc3 > 0)
table(Hep.Gpc3@meta.data$replicate)

Cancer <- subset(HCCp1.cells, idents = "Cancer", invert = FALSE)
table(Cancer@meta.data$replicate)

DoHeatmap(Cancer, features = c("Hnf4a", "Afp", "Gpc3"), group.by = "new.group", size = 3)
DoHeatmap(Cancer, features = c("Hnf4a", "Afp", "Gpc3"), group.by = "new.group", size = 3) + NoLegend()

Cancer.Hnf4a <- subset(Cancer, subset = Hnf4a > 0)
table(Cancer.Hnf4a@meta.data$replicate)
Cancer.Afp <- subset(Cancer, subset = Afp > 0)
table(Cancer.Afp@meta.data$replicate)
Cancer.Gpc3 <- subset(Cancer, subset = Gpc3 > 0)
table(Cancer.Gpc3@meta.data$replicate)

Hep.and.Cancer <- merge(Hep, y = Cancer)
table(Hep.and.Cancer@meta.data$new.group)

VlnPlot(Hep.and.Cancer, cols = c("red", "blue", "green", "black", "orange"), features = c("Hnf4a", "Afp", "Gpc3"), pt.size = 1, split.by = "new.group", group.by = "Cell.Type", log = TRUE, stack = TRUE, same.y.lims = TRUE)

AverageExpression(Hep, feature = "Hnf4a", group.by = "replicate")
AverageExpression(Hep, feature = "Afp", group.by = "replicate")
AverageExpression(Hep, feature = "Gpc3", group.by = "replicate")

AverageExpression(Cancer, feature = "Hnf4a", group.by = "replicate")
AverageExpression(Cancer, feature = "Afp", group.by = "replicate")
AverageExpression(Cancer, feature = "Gpc3", group.by = "replicate")

AverageExpression(Hep.Hnf4a, feature = "Hnf4a", group.by = "replicate")
AverageExpression(Hep.Afp, feature = "Afp", group.by = "replicate")
AverageExpression(Hep.Gpc3, feature = "Gpc3", group.by = "replicate")

AverageExpression(Cancer.Hnf4a, feature = "Hnf4a", group.by = "replicate")
AverageExpression(Cancer.Afp, feature = "Afp", group.by = "replicate")
AverageExpression(Cancer.Gpc3, feature = "Gpc3", group.by = "replicate")

##############################################################
##save and read in Hep and Cancer populations only 
saveRDS(Hep, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.Hep.only.rds")
saveRDS(Cancer, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.Cancer.only.rds")

Hep <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.Hep.only.rds")
Cancer <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.group.names.Cancer.only.rds")

Hep.and.Cancer.Hnf4apos <- merge(Hep.Hnf4a, y = Cancer.Hnf4a)
table(Hep.and.Cancer.Hnf4apos@meta.data$new.group)

VlnPlot(Hep.and.Cancer.Hnf4apos, cols = c("red", "blue", "green", "black", "orange"), features = c("Hnf4a", "Afp", "Gpc3"), pt.size = 1, split.by = "new.group", group.by = "Cell.Type", log = TRUE, stack = TRUE, same.y.lims = TRUE)

Hep.and.Cancer.Afppos <- merge(Hep.Afp, y = Cancer.Afp)
table(Hep.and.Cancer.Afppos@meta.data$new.group)

VlnPlot(Hep.and.Cancer.Afppos, cols = c("red", "blue", "green", "black", "orange"), features = c("Hnf4a", "Afp", "Gpc3"), pt.size = 1, split.by = "new.group", group.by = "Cell.Type", log = TRUE, stack = TRUE, same.y.lims = TRUE)

Hep.and.Cancer.Gpc3pos <- merge(Hep.Gpc3, y = Cancer.Gpc3)
table(Hep.and.Cancer.Gpc3pos@meta.data$new.group)

VlnPlot(Hep.and.Cancer.Gpc3pos, cols = c("red", "blue", "green", "black", "orange"), features = c("Hnf4a", "Afp", "Gpc3"), pt.size = 1, split.by = "new.group", group.by = "Cell.Type", log = TRUE, stack = TRUE, same.y.lims = TRUE)

#####################################################################################################################################################################################
#####################################################################################################################################################################################
##load in this data to separate CD4/CD8 Subsets and KC/Mac to do CellChat analyses for assessing MHC interactions
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.LSEC1and2.all.new.labels.RDS")
Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Predicted_subType

CD4 <- subset(HCCp1.cells, idents = "CD4", invert = FALSE)
CD8 <- subset(HCCp1.cells, idents = "CD8", invert = FALSE)

KC <- subset(HCCp1.cells, idents = "Kupffer", invert = FALSE)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Predicted_Type
Mac <- subset(HCCp1.cells, idents = "Macrophage")
table(Mac@meta.data$Predicted_subType)

Idents(Mac) <- Mac@meta.data$Predicted_subType
Mac.no.KC <- subset(Mac, idents = c("M1", "M2", "Unknown"), invert = FALSE)
table(Mac.no.KC@meta.data$Predicted_subType)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Predicted_Type

Bcell <- subset(HCCp1.cells, idents = "B cell", invert = FALSE)
DC <- subset(HCCp1.cells, idents = "Dendritic cell", invert = FALSE)
NKT <- subset(HCCp1.cells, idents = "NKT cell", invert = FALSE)
NK <- subset(HCCp1.cells, idents = "NK cell", invert = FALSE)
Neutro <- subset(HCCp1.cells, idents = "Neutrophil", invert = FALSE)
Mono <- subset(HCCp1.cells, idents = "Monocyte", invert = FALSE)
Endo <- subset(HCCp1.cells, idents = "Endothelial cell", invert = FALSE)
LSEC <- subset(HCCp1.cells, idents = "Liver sinusoid endothelial cell", invert = FALSE)
Stromal <- subset(HCCp1.cells, idents = "Stromal cell", invert = FALSE)
HSC <- subset(HCCp1.cells, idents = "Hepatic stellate cell", invert = FALSE)
Fibro <- subset(HCCp1.cells, idents = "Fibroblast", invert = FALSE)
Myofibro <- subset(HCCp1.cells, idents = "Myofibroblast", invert = FALSE)
Cholangio <- subset(HCCp1.cells, idents = "Cholangiocyte", invert = FALSE)
Hep <- subset(HCCp1.cells, idents = "Hepatocyte", invert = FALSE)
Cancer <- subset(HCCp1.cells, idents = "Cancer cell", invert = FALSE)

new.label.Bcell <- "Bcell"
new.label.CD4 <- "CD4"
new.label.CD8 <- "CD8"
new.label.Cancer <- "Cancer"
new.label.Cholangio <- "Cholangio"
new.label.DC <- "DC"
new.label.Endo <- "Endo"
new.label.Fibro <- "Fibro"
new.label.HSC <- "HSC"
new.label.Hep <- "Hep"
new.label.LSEC <- "LSEC"
new.label.KC <- "KC"
new.label.Mac <- "Mac"
new.label.Mono <- "Mono"
new.label.Myofibro <- "Myofibro"
new.label.Neutro <- "Neutro"
new.label.NK <- "NK"
new.label.NKT <- "NKT"
new.label.Stromal <- "Stromal"

new.Bcell <- AddMetaData(Bcell, metadata = new.label.Bcell, col.name = 'Select.subsets')
new.CD4 <- AddMetaData(CD4, metadata = new.label.CD4, col.name = 'Select.subsets')
new.CD8 <- AddMetaData(CD8, metadata = new.label.CD8, col.name = 'Select.subsets')
new.Cancer <- AddMetaData(Cancer, metadata = new.label.Cancer, col.name = 'Select.subsets')
new.Cholangio <- AddMetaData(Cholangio, metadata = new.label.Cholangio, col.name = 'Select.subsets')
new.DC <- AddMetaData(DC, metadata = new.label.DC, col.name = 'Select.subsets')
new.Endo <- AddMetaData(Endo, metadata = new.label.Endo, col.name = 'Select.subsets')
new.Fibro <- AddMetaData(Fibro, metadata = new.label.Fibro, col.name = 'Select.subsets')
new.HSC <- AddMetaData(HSC, metadata = new.label.HSC, col.name = 'Select.subsets')
new.Hep <- AddMetaData(Hep, metadata = new.label.Hep, col.name = 'Select.subsets')
new.LSEC <- AddMetaData(LSEC, metadata = new.label.LSEC, col.name = 'Select.subsets')
new.KC <- AddMetaData(KC, metadata = new.label.KC, col.name = 'Select.subsets')
new.Mac <- AddMetaData(Mac.no.KC, metadata = new.label.Mac, col.name = 'Select.subsets')
new.Mono <- AddMetaData(Mono, metadata = new.label.Mono, col.name = 'Select.subsets')
new.Myofibro <- AddMetaData(Myofibro, metadata = new.label.Myofibro, col.name = 'Select.subsets')
new.Neutro <- AddMetaData(Neutro, metadata = new.label.Neutro, col.name = 'Select.subsets')
new.NK <- AddMetaData(NK, metadata = new.label.NK, col.name = 'Select.subsets')
new.NKT <- AddMetaData(NKT, metadata = new.label.NKT, col.name = 'Select.subsets')
new.Stromal <- AddMetaData(Stromal, metadata = new.label.Stromal, col.name = 'Select.subsets')

HCCp1.CD4.CD8.KC.Mac <- merge(new.Bcell, y = c(new.CD4, new.CD8, new.Cancer, new.Cholangio, new.DC, new.Endo, new.Fibro, new.HSC, new.Hep, new.LSEC, new.KC, new.Mac, new.Mono, new.Myofibro, new.Neutro, new.NK, new.NKT, new.Stromal))
table(HCCp1.CD4.CD8.KC.Mac@meta.data$Select.subsets)

saveRDS(HCCp1.CD4.CD8.KC.Mac, file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.CD4.CD8.KC.Mac.select.subsets.rds")
