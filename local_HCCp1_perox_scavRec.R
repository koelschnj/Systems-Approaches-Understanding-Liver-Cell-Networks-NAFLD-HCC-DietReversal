library(Seurat)
library(CellChat)
library(future)
library(patchwork)
library(ggplot2)
options(stringsAsFactors = FALSE)

##use this to load in abbreviated cell tyep names for updated group names as well analysis
HCCp1.cells <- readRDS(file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/Reverse Diet Data/HCCp1.cells.sorted.new.abbreviated.group.names.rds")

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$new.group
table(HCCp1.cells@meta.data$new.group)

new.group.levels <- c("CD", "WD.nf", "WD.t", "RD.t", "RD.n")
HCCp1.cells@meta.data$new.group <- factor(x = HCCp1.cells@meta.data$new.group, levels = new.group.levels)
table(HCCp1.cells@meta.data$new.group)

Idents(HCCp1.cells) <- HCCp1.cells@meta.data$Abbreviated
table(HCCp1.cells@meta.data$Abbreviated)

##first generate dotplots of peroxidase and scavenger receptors for monocytes across all groups
Mono <- subset(HCCp1.cells, idents = "Mono", invert = FALSE)
table(Mono@meta.data$new.group)

DotPlot(Mono, features = c("Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx2", "Gpx3", "Gpx4"), group.by = "new.group", cols = c("blue", "red")) + coord_flip()
DotPlot(Mono, features = c("Msr1", "Marco", "Scarb1", "Cd36", "Cd68", "Olr1", "Scarf1", "Stab1", "Stab2", "Timd4", "Mertk"), group.by = "new.group", cols = c("blue", "red")) + coord_flip()

DotPlot(Mono, features = c("Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx2", "Gpx3", "Gpx4"), group.by = "new.group", cols = c("blue", "red"), scale.max = 60, scale.min = 0, col.min = -1.5, col.max = 2.0) + scale_color_gradientn(colors = c("blue","red"), limits=c(-1.5,2.0)) + coord_flip()
DotPlot(Mono, features = c("Msr1", "Marco", "Scarb1", "Cd36", "Cd68", "Olr1", "Scarf1", "Stab1", "Stab2", "Timd4", "Mertk"), group.by = "new.group", cols = c("blue", "red"), scale.max = 40, scale.min = 0, col.min = -1.5, col.max = 2.0) + scale_color_gradientn(colors = c("blue","red"), limits=c(-1.5,2.0)) + coord_flip()

##second generate dotplots of peroxidase and scavenger receptors for monocytes in groups of interest (WD.t, RD.t, RD.n)
Idents(Mono) <- Mono@meta.data$new.group
table(Mono@meta.data$new.group)

GoI.Mono <- subset(Mono, idents = c("WD.t", "RD.t", "RD.n"), invert = FALSE)
table(GoI.Mono@meta.data$new.group)

DotPlot(GoI.Mono, features = c("Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx2", "Gpx3", "Gpx4"), group.by = "new.group", cols = c("blue", "red")) + coord_flip()
DotPlot(GoI.Mono, features = c("Msr1", "Marco", "Scarb1", "Cd36", "Cd68", "Olr1", "Scarf1", "Stab1", "Stab2", "Timd4", "Mertk"), group.by = "new.group", cols = c("blue", "red")) + coord_flip()

DotPlot(GoI.Mono, features = c("Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx2", "Gpx3", "Gpx4"), group.by = "new.group", cols = c("blue", "red"), scale.max = 60, scale.min = 0, col.min = -1.5, col.max = 1.5) + scale_color_gradientn(colors = c("blue","red"), limits=c(-1.5,1.5)) + coord_flip()
DotPlot(GoI.Mono, features = c("Msr1", "Marco", "Scarb1", "Cd36", "Cd68", "Olr1", "Scarf1", "Stab1", "Stab2", "Timd4", "Mertk"), group.by = "new.group", cols = c("blue", "red"), scale.max = 40, scale.min = 0, col.min = -1.5, col.max = 1.5) + scale_color_gradientn(colors = c("blue","red"), limits=c(-1.5,1.5)) + coord_flip()
