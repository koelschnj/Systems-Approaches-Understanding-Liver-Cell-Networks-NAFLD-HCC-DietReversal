library("ggplot2")
library("Seurat")

##need to quantify cellular patterns of everything we annotated for each replicate
##first load in collectively annotated files to get overall cell types for each replicate

##first look at cell types in each replicate bulk sorted approach
HCCp1.bulk.sorted <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.w.scSorter.RDS")

Idents(HCCp1.bulk.sorted) <- HCCp1.bulk.sorted@meta.data$replicate
table(HCCp1.bulk.sorted@meta.data$replicate)

CD1.cells <- subset(x = HCCp1.bulk.sorted, idents = "CD1", invert = FALSE)
table(CD1.cells@meta.data$Predicted_Type)

CD2.cells <- subset(x = HCCp1.bulk.sorted, idents = "CD2", invert = FALSE)
table(CD2.cells@meta.data$Predicted_Type)

PreT1.cells <- subset(x = HCCp1.bulk.sorted, idents = "PreT1", invert = FALSE)
table(PreT1.cells@meta.data$Predicted_Type)

PreT2.cells <- subset(x = HCCp1.bulk.sorted, idents = "PreT2", invert = FALSE)
table(PreT2.cells@meta.data$Predicted_Type)

MWD.T1.cells <- subset(x = HCCp1.bulk.sorted, idents = "MWD.T1", invert = FALSE)
table(MWD.T1.cells@meta.data$Predicted_Type)

MWD.T2.cells <- subset(x = HCCp1.bulk.sorted, idents = "MWD.T2", invert = FALSE)
table(MWD.T2.cells@meta.data$Predicted_Type)

MWD.T3.cells <- subset(x = HCCp1.bulk.sorted, idents = "MWD.T3", invert = FALSE)
table(MWD.T3.cells@meta.data$Predicted_Type)

MRD.T1.cells <- subset(x = HCCp1.bulk.sorted, idents = "MRD.T1", invert = FALSE)
table(MRD.T1.cells@meta.data$Predicted_Type)

MRD.T2.cells <- subset(x = HCCp1.bulk.sorted, idents = "MRD.T2", invert = FALSE)
table(MRD.T2.cells@meta.data$Predicted_Type)

MRD.NT1.cells <- subset(x = HCCp1.bulk.sorted, idents = "MRD.NT1", invert = FALSE)
table(MRD.NT1.cells@meta.data$Predicted_Type)

MRD.NT2.cells <- subset(x = HCCp1.bulk.sorted, idents = "MRD.NT2", invert = FALSE)
table(MRD.NT2.cells@meta.data$Predicted_Type)

MRD.NT3.cells <- subset(x = HCCp1.bulk.sorted, idents = "MRD.NT3", invert = FALSE)
table(MRD.NT3.cells@meta.data$Predicted_Type)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
##now look at cell types in each replicate for individual sorting approach
HCCp1.sorted.individually <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/HCCp1.cells.sorted.individually.w.scSorter.RDS")

Idents(HCCp1.sorted.individually) <- HCCp1.sorted.individually@meta.data$replicate
table(HCCp1.sorted.individually@meta.data$replicate)

CD1.cells.i <- subset(x = HCCp1.sorted.individually, idents = "CD1", invert = FALSE)
table(CD1.cells.i@meta.data$Predicted_Type)

CD2.cells.i <- subset(x = HCCp1.sorted.individually, idents = "CD2", invert = FALSE)
table(CD2.cells.i@meta.data$Predicted_Type)

PreT1.cells.i <- subset(x = HCCp1.sorted.individually, idents = "PreT1", invert = FALSE)
table(PreT1.cells.i@meta.data$Predicted_Type)

PreT2.cells.i <- subset(x = HCCp1.sorted.individually, idents = "PreT2", invert = FALSE)
table(PreT2.cells.i@meta.data$Predicted_Type)

MWD.T1.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MWD.T1", invert = FALSE)
table(MWD.T1.cells.i@meta.data$Predicted_Type)

MWD.T2.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MWD.T2", invert = FALSE)
table(MWD.T2.cells.i@meta.data$Predicted_Type)

MWD.T3.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MWD.T3", invert = FALSE)
table(MWD.T3.cells.i@meta.data$Predicted_Type)

MRD.T1.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MRD.T1", invert = FALSE)
table(MRD.T1.cells.i@meta.data$Predicted_Type)

MRD.T2.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MRD.T2", invert = FALSE)
table(MRD.T2.cells.i@meta.data$Predicted_Type)

MRD.NT1.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MRD.NT1", invert = FALSE)
table(MRD.NT1.cells.i@meta.data$Predicted_Type)

MRD.NT2.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MRD.NT2", invert = FALSE)
table(MRD.NT2.cells.i@meta.data$Predicted_Type)

MRD.NT3.cells.i <- subset(x = HCCp1.sorted.individually, idents = "MRD.NT3", invert = FALSE)
table(MRD.NT3.cells.i@meta.data$Predicted_Type)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
##second load in sorted subset annotated files to quantify for each replicate
##look at bulk sorted subsets first and merge all at the end and save for potential future use
HCCp1.bulk.sorted.MPs <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.MPs.bulk.sorted.w.scSorter.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.MPs@meta.data$group <- factor(x = HCCp1.bulk.sorted.MPs@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.MPs) <- HCCp1.bulk.sorted.MPs@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.MPs@meta.data$Predicted_subType)

bulk.KC <- subset(x = HCCp1.bulk.sorted.MPs, idents = "Kupffer", invert = FALSE)
table(bulk.KC@meta.data$replicate)

bulk.M1 <- subset(x = HCCp1.bulk.sorted.MPs, idents = "M1", invert = FALSE)
table(bulk.M1@meta.data$replicate)

bulk.M2 <- subset(x = HCCp1.bulk.sorted.MPs, idents = "M2", invert = FALSE)
table(bulk.M2@meta.data$replicate)

bulk.MP.Unk <- subset(x = HCCp1.bulk.sorted.MPs, idents = "Unknown", invert = FALSE)
table(bulk.MP.Unk@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.bulk.sorted.Tcells <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.Tcells.bulk.sorted.w.GL14.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.Tcells@meta.data$group <- factor(x = HCCp1.bulk.sorted.Tcells@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.Tcells) <- HCCp1.bulk.sorted.Tcells@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.Tcells@meta.data$Predicted_subType)

bulk.CD4 <- subset(x = HCCp1.bulk.sorted.Tcells, idents = "CD4", invert = FALSE)
table(bulk.CD4@meta.data$replicate)

bulk.CD8 <- subset(x = HCCp1.bulk.sorted.Tcells, idents = "CD8", invert = FALSE)
table(bulk.CD8@meta.data$replicate)

bulk.Tcell.Unk <- subset(x = HCCp1.bulk.sorted.Tcells, idents = "Unknown", invert = FALSE)
table(bulk.Tcell.Unk@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.bulk.sorted.Tcell <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.Tcell.bulk.sorted.w.GL18.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.Tcell@meta.data$group <- factor(x = HCCp1.bulk.sorted.Tcell@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.Tcell) <- HCCp1.bulk.sorted.Tcell@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.Tcell@meta.data$Predicted_subType)

bulk.CD4.GL18 <- subset(x = HCCp1.bulk.sorted.Tcell, idents = "CD4", invert = FALSE)
table(bulk.CD4.GL18@meta.data$replicate)

bulk.CD8.GL18 <- subset(x = HCCp1.bulk.sorted.Tcell, idents = "CD8", invert = FALSE)
table(bulk.CD8.GL18@meta.data$replicate)

bulk.Tcell.Unk.GL18 <- subset(x = HCCp1.bulk.sorted.Tcell, idents = "Unknown", invert = FALSE)
table(bulk.Tcell.Unk.GL18@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.bulk.sorted.Tcel <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.Tcells.bulk.sorted.w.GL19.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.Tcel@meta.data$group <- factor(x = HCCp1.bulk.sorted.Tcel@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.Tcel) <- HCCp1.bulk.sorted.Tcel@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.Tcel@meta.data$Predicted_subType)

bulk.CD4.GL19 <- subset(x = HCCp1.bulk.sorted.Tcel, idents = "CD4", invert = FALSE)
table(bulk.CD4.GL19@meta.data$replicate)

bulk.CD8.GL19 <- subset(x = HCCp1.bulk.sorted.Tcel, idents = "CD8", invert = FALSE)
table(bulk.CD8.GL19@meta.data$replicate)

bulk.Tcell.Unk.GL19 <- subset(x = HCCp1.bulk.sorted.Tcel, idents = "Unknown", invert = FALSE)
table(bulk.Tcell.Unk.GL19@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.bulk.sorted.DC <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.DC.bulk.sorted.w.GL10.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.DC@meta.data$group <- factor(x = HCCp1.bulk.sorted.DC@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.DC) <- HCCp1.bulk.sorted.DC@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.DC@meta.data$Predicted_subType)

bulk.DC1 <- subset(x = HCCp1.bulk.sorted.DC, idents = "cDC1", invert = FALSE)
table(bulk.DC1@meta.data$replicate)

bulk.DC2 <- subset(x = HCCp1.bulk.sorted.DC, idents = "cDC2", invert = FALSE)
table(bulk.DC2@meta.data$replicate)

bulk.pDC<- subset(x = HCCp1.bulk.sorted.DC, idents = "pDC", invert = FALSE)
table(bulk.pDC@meta.data$replicate)

bulk.Unk.DC<- subset(x = HCCp1.bulk.sorted.DC, idents = "Unknown", invert = FALSE)
table(bulk.Unk.DC@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.bulk.sorted.Bcell <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.Bcell.bulk.sorted.w.GL16.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.Bcell@meta.data$group <- factor(x = HCCp1.bulk.sorted.Bcell@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.Bcell) <- HCCp1.bulk.sorted.Bcell@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.Bcell@meta.data$Predicted_subType)

bulk.B1 <- subset(x = HCCp1.bulk.sorted.Bcell, idents = "B1", invert = FALSE)
table(bulk.B1@meta.data$replicate)

bulk.B2 <- subset(x = HCCp1.bulk.sorted.Bcell, idents = "B2", invert = FALSE)
table(bulk.B2@meta.data$replicate)

bulk.Unk.Bcell<- subset(x = HCCp1.bulk.sorted.Bcell, idents = "Unknown", invert = FALSE)
table(bulk.Unk.Bcell@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.bulk.sorted.Mono <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/HCCp1.Mono.bulk.sorted.w.GL17.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.bulk.sorted.Mono@meta.data$group <- factor(x = HCCp1.bulk.sorted.Mono@meta.data$group, levels = new.group.levels)

Idents(HCCp1.bulk.sorted.Mono) <- HCCp1.bulk.sorted.Mono@meta.data$Predicted_subType
table(HCCp1.bulk.sorted.Mono@meta.data$Predicted_subType)

bulk.mMDSC <- subset(x = HCCp1.bulk.sorted.Mono, idents = "MDSC", invert = FALSE)
table(bulk.mMDSC@meta.data$replicate)

bulk.mono <- subset(x = HCCp1.bulk.sorted.Mono, idents = "Monocyte", invert = FALSE)
table(bulk.mono@meta.data$replicate)

bulk.Unk.mono<- subset(x = HCCp1.bulk.sorted.Mono, idents = "Unknown", invert = FALSE)
table(bulk.Unk.mono@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
##finally combine all loaded bulk sorted subsets and recombine to save as a single object
bulk.annotated.subsets <- merge(HCCp1.bulk.sorted.MPs, y = c(HCCp1.bulk.sorted.Tcel, HCCp1.bulk.sorted.DC, HCCp1.bulk.sorted.Bcell, HCCp1.bulk.sorted.Mono))
table(bulk.annotated.subsets@meta.data$group)

saveRDS(bulk.annotated.subsets, file = "~/HCC_paper1/Saved_Rfiles/Bulk_subsets/Bulk.annotated.subsets.combined.RDS")

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
##now look at individually sorted subsets and merge all at to save for future use
HCCp1.indiv.sorted.MPs <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/MPs.annotated.individually.w.scSorter.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.MPs@meta.data$group <- factor(x = HCCp1.indiv.sorted.MPs@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.MPs) <- HCCp1.indiv.sorted.MPs@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.MPs@meta.data$Predicted_subType)

indiv.KC <- subset(x = HCCp1.indiv.sorted.MPs, idents = "Kupffer", invert = FALSE)
table(indiv.KC@meta.data$replicate)

indiv.M1 <- subset(x = HCCp1.indiv.sorted.MPs, idents = "M1", invert = FALSE)
table(indiv.M1@meta.data$replicate)

indiv.M2 <- subset(x = HCCp1.indiv.sorted.MPs, idents = "M2", invert = FALSE)
table(indiv.M2@meta.data$replicate)

indiv.MP.Unk <- subset(x = HCCp1.indiv.sorted.MPs, idents = "Unknown", invert = FALSE)
table(indiv.MP.Unk@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.indiv.sorted.Tcells <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Tcells.annotated.individually.w.scSorter.GL14.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.Tcells@meta.data$group <- factor(x = HCCp1.indiv.sorted.Tcells@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.Tcells) <- HCCp1.indiv.sorted.Tcells@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.Tcells@meta.data$Predicted_subType)

indiv.CD4 <- subset(x = HCCp1.indiv.sorted.Tcells, idents = "CD4", invert = FALSE)
table(indiv.CD4@meta.data$replicate)

indiv.CD8 <- subset(x = HCCp1.indiv.sorted.Tcells, idents = "CD8", invert = FALSE)
table(indiv.CD8@meta.data$replicate)

indiv.Tcell.Unk <- subset(x = HCCp1.indiv.sorted.Tcells, idents = "Unknown", invert = FALSE)
table(indiv.Tcell.Unk@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.indiv.sorted.Tcell <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Tcell.annotated.individually.w.scSorter.GL18.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.Tcell@meta.data$group <- factor(x = HCCp1.indiv.sorted.Tcell@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.Tcell) <- HCCp1.indiv.sorted.Tcell@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.Tcell@meta.data$Predicted_subType)

indiv.CD4.GL18 <- subset(x = HCCp1.indiv.sorted.Tcell, idents = "CD4", invert = FALSE)
table(indiv.CD4.GL18@meta.data$replicate)

indiv.CD8.GL18 <- subset(x = HCCp1.indiv.sorted.Tcell, idents = "CD8", invert = FALSE)
table(indiv.CD8.GL18@meta.data$replicate)

indiv.Tcell.Unk.GL18 <- subset(x = HCCp1.indiv.sorted.Tcell, idents = "Unknown", invert = FALSE)
table(indiv.Tcell.Unk.GL18@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.indiv.sorted.Tcel <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Tcel.annotated.individually.w.scSorter.GL19.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.Tcel@meta.data$group <- factor(x = HCCp1.indiv.sorted.Tcel@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.Tcel) <- HCCp1.indiv.sorted.Tcel@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.Tcel@meta.data$Predicted_subType)

indiv.CD4.GL19 <- subset(x = HCCp1.indiv.sorted.Tcel, idents = "CD4", invert = FALSE)
table(indiv.CD4.GL19@meta.data$replicate)

indiv.CD8.GL19 <- subset(x = HCCp1.indiv.sorted.Tcel, idents = "CD8", invert = FALSE)
table(indiv.CD8.GL19@meta.data$replicate)

indiv.Tcell.Unk.GL19 <- subset(x = HCCp1.indiv.sorted.Tcel, idents = "Unknown", invert = FALSE)
table(indiv.Tcell.Unk.GL19@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.indiv.sorted.DC <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/DCs.annotated.individually.w.scSorter.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.DC@meta.data$group <- factor(x = HCCp1.indiv.sorted.DC@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.DC) <- HCCp1.indiv.sorted.DC@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.DC@meta.data$Predicted_subType)

indiv.DC1 <- subset(x = HCCp1.indiv.sorted.DC, idents = "cDC1", invert = FALSE)
table(indiv.DC1@meta.data$replicate)

indiv.DC2 <- subset(x = HCCp1.indiv.sorted.DC, idents = "cDC2", invert = FALSE)
table(indiv.DC2@meta.data$replicate)

indiv.pDC <- subset(x = HCCp1.indiv.sorted.DC, idents = "pDC", invert = FALSE)
table(indiv.pDC@meta.data$replicate)

indiv.Unk.DC <- subset(x = HCCp1.indiv.sorted.DC, idents = "Unknown", invert = FALSE)
table(indiv.Unk.DC@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.indiv.sorted.Bcell <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Bcells.annotated.individually.w.scSorter.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.Bcell@meta.data$group <- factor(x = HCCp1.indiv.sorted.Bcell@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.Bcell) <- HCCp1.indiv.sorted.Bcell@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.Bcell@meta.data$Predicted_subType)

indiv.B1 <- subset(x = HCCp1.indiv.sorted.Bcell, idents = "B1", invert = FALSE)
table(indiv.B1@meta.data$replicate)

indiv.B2 <- subset(x = HCCp1.indiv.sorted.Bcell, idents = "B2", invert = FALSE)
table(indiv.B2@meta.data$replicate)

indiv.Unk.B <- subset(x = HCCp1.indiv.sorted.Bcell, idents = "Unknown", invert = FALSE)
table(indiv.Unk.B@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
HCCp1.indiv.sorted.Mono <- readRDS(file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/Mono.annotated.individually.w.scSorter.RDS")

new.group.levels <- c("CD", "PreT", "MWD.T", "MRD.T", "MRD.NT")
HCCp1.indiv.sorted.Mono@meta.data$group <- factor(x = HCCp1.indiv.sorted.Mono@meta.data$group, levels = new.group.levels)

Idents(HCCp1.indiv.sorted.Mono) <- HCCp1.indiv.sorted.Mono@meta.data$Predicted_subType
table(HCCp1.indiv.sorted.Mono@meta.data$Predicted_subType)

indiv.mMDSC <- subset(x = HCCp1.indiv.sorted.Mono, idents = "MDSC", invert = FALSE)
table(indiv.mMDSC@meta.data$replicate)

indiv.mono <- subset(x = HCCp1.indiv.sorted.Mono, idents = "Monocyte", invert = FALSE)
table(indiv.mono@meta.data$replicate)

indiv.Unk.mono <- subset(x = HCCp1.indiv.sorted.Mono, idents = "Unknown", invert = FALSE)
table(indiv.Unk.mono@meta.data$replicate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
##finally combine all loaded bulk sorted subsets and recombine to save as a single object
indiv.annotated.subsets <- merge(HCCp1.indiv.sorted.MPs, y = c(HCCp1.indiv.sorted.Tcel, HCCp1.indiv.sorted.DC, HCCp1.indiv.sorted.Bcell, HCCp1.indiv.sorted.Mono))
table(indiv.annotated.subsets@meta.data$group)

saveRDS(indiv.annotated.subsets, file = "~/HCC_paper1/Saved_Rfiles/individual_subsets/individually.annotated.subsets.combined.RDS")
