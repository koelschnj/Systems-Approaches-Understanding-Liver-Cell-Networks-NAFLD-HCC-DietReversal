##script devised for comparing two cell types in one sample to detect malignant events in cancer cells compared to hepatocyte
library(tidyverse)
library(DESeq2)

CD.HepvCancer.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/CD.HepvCancer.counts.csv", header = TRUE, sep = ",")
PreT.HepvCancer.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/PreT.HepvCancer.counts.csv", header = TRUE, sep = ",")
MWD.T.HepvCancer.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/MWD.T.HepvCancer.counts.csv", header = TRUE, sep = ",")
MRD.T.HepvCancer.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/MRD.T.HepvCancer.counts.csv", header = TRUE, sep = ",")
MRD.NT.HepvCancer.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/MRD.NT.HepvCancer.counts.csv", header = TRUE, sep = ",")

head(CD.HepvCancer.counts)
head(PreT.HepvCancer.counts)
head(MWD.T.HepvCancer.counts)
head(MRD.T.HepvCancer.counts)
head(MRD.NT.HepvCancer.counts)

CD.HepvCancer.counts$X[duplicated(CD.HepvCancer.counts$X)]
PreT.HepvCancer.counts$X[duplicated(PreT.HepvCancer.counts$X)]
MWD.T.HepvCancer.counts$X[duplicated(MWD.T.HepvCancer.counts$X)]
MRD.T.HepvCancer.counts$X[duplicated(MRD.T.HepvCancer.counts$X)]
MRD.NT.HepvCancer.counts$X[duplicated(MRD.NT.HepvCancer.counts$X)]

which(duplicated(CD.HepvCancer.counts$X))
which(duplicated(PreT.HepvCancer.counts$X))
which(duplicated(MWD.T.HepvCancer.counts$X))
which(duplicated(MRD.T.HepvCancer.counts$X))
which(duplicated(MRD.NT.HepvCancer.counts$X))

CD.HepvCancer.counts$X[9519] <- "1-Mar_dup"
CD.HepvCancer.counts$X[18441] <- "2-Mar_dup"
PreT.HepvCancer.counts$X[9519] <- "1-Mar_dup"
PreT.HepvCancer.counts$X[18441] <- "2-Mar_dup"
MWD.T.HepvCancer.counts$X[9519] <- "1-Mar_dup"
MWD.T.HepvCancer.counts$X[18441] <- "2-Mar_dup"
MRD.T.HepvCancer.counts$X[9519] <- "1-Mar_dup"
MRD.T.HepvCancer.counts$X[18441] <- "2-Mar_dup"
MRD.NT.HepvCancer.counts$X[9519] <- "1-Mar_dup"
MRD.NT.HepvCancer.counts$X[18441] <- "2-Mar_dup"

which(duplicated(CD.HepvCancer.counts$X))
which(duplicated(PreT.HepvCancer.counts$X))
which(duplicated(MWD.T.HepvCancer.counts$X))
which(duplicated(MRD.T.HepvCancer.counts$X))
which(duplicated(MRD.NT.HepvCancer.counts$X))

row.names(CD.HepvCancer.counts) <- CD.HepvCancer.counts$X
row.names(PreT.HepvCancer.counts) <- PreT.HepvCancer.counts$X
row.names(MWD.T.HepvCancer.counts) <- MWD.T.HepvCancer.counts$X
row.names(MRD.T.HepvCancer.counts) <- MRD.T.HepvCancer.counts$X
row.names(MRD.NT.HepvCancer.counts) <- MRD.NT.HepvCancer.counts$X

CD.HepvCancer.counts <- CD.HepvCancer.counts[,-1]
PreT.HepvCancer.counts <- PreT.HepvCancer.counts[,-1]
MWD.T.HepvCancer.counts <- MWD.T.HepvCancer.counts[,-1]
MRD.T.HepvCancer.counts <- MRD.T.HepvCancer.counts[,-1]
MRD.NT.HepvCancer.counts <- MRD.NT.HepvCancer.counts[,-1]

colData.CD <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/CD.HepvCancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.PreT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/PreT.HepvCancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.MWD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/MWD.T.HepvCancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.MRD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/MRD.T.HepvCancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.MRD.NT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/HepvCancer/MRD.NT.HepvCancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)

View(colData.CD)
View(colData.PreT)
View(colData.MWD.T)
View(colData.MRD.T)
View(colData.MRD.NT)

all(names(CD.HepvCancer.counts) %in% row.names(colData.CD))
all(names(PreT.HepvCancer.counts) %in% row.names(colData.PreT))
all(names(MWD.T.HepvCancer.counts) %in% row.names(colData.MWD.T))
all(names(MRD.T.HepvCancer.counts) %in% row.names(colData.MRD.T))
all(names(MRD.NT.HepvCancer.counts) %in% row.names(colData.MRD.NT))

all(names(CD.HepvCancer.counts) == row.names(colData.CD))
all(names(PreT.HepvCancer.counts) == row.names(colData.PreT))
all(names(MWD.T.HepvCancer.counts) == row.names(colData.MWD.T))
all(names(MRD.T.HepvCancer.counts) == row.names(colData.MRD.T))
all(names(MRD.NT.HepvCancer.counts) == row.names(colData.MRD.NT))

##construct DESeqDataSet object
dds.CD <- DESeqDataSetFromMatrix(countData = CD.HepvCancer.counts, colData = colData.CD, design = ~ condition)
dds.PreT <- DESeqDataSetFromMatrix(countData = PreT.HepvCancer.counts, colData = colData.PreT, design = ~ condition)
dds.MWD.T <- DESeqDataSetFromMatrix(countData = MWD.T.HepvCancer.counts, colData = colData.MWD.T, design = ~ condition)
dds.MRD.T <- DESeqDataSetFromMatrix(countData = MRD.T.HepvCancer.counts, colData = colData.MRD.T, design = ~ condition)
dds.MRD.NT <- DESeqDataSetFromMatrix(countData = MRD.NT.HepvCancer.counts, colData = colData.MRD.NT, design = ~ condition)

dds.CD
dds.PreT
dds.MWD.T
dds.MRD.T
dds.MRD.NT

##pre-filtering to remove rows with low gene counts
keep.CD <- rowSums(counts(dds.CD)) >= 3
dds.CD <- dds.CD[keep.CD,]
dds.CD

keep.PreT <- rowSums(counts(dds.PreT)) >= 3
dds.PreT <- dds.PreT[keep.PreT,]
dds.PreT

keep.MWD.T <- rowSums(counts(dds.MWD.T)) >= 3
dds.MWD.T <- dds.MWD.T[keep.MWD.T,]
dds.MWD.T

keep.MRD.T <- rowSums(counts(dds.MRD.T)) >= 3
dds.MRD.T <- dds.MRD.T[keep.MRD.T,]
dds.MRD.T

keep.MRD.NT <- rowSums(counts(dds.MRD.NT)) >= 3
dds.MRD.NT <- dds.MRD.NT[keep.MRD.NT,]
dds.MRD.NT

##run DESeq
dds.CD <- DESeq(dds.CD)
dds.PreT <- DESeq(dds.PreT)
dds.MWD.T <- DESeq(dds.MWD.T)
dds.MRD.T <- DESeq(dds.MRD.T)
dds.MRD.NT <- DESeq(dds.MRD.NT)

##Set up for the contrast function is the variable (condition), one to test (T1/high), and reference (C1/low)
res.CD <- results(dds.CD, contrast = c("condition", "Cancer", "Hepatocyte"))
res.PreT <- results(dds.PreT, contrast = c("condition", "Cancer", "Hepatocyte"))
res.MWD.T <- results(dds.MWD.T, contrast = c("condition", "Cancer", "Hepatocyte"))
res.MRD.T <- results(dds.MRD.T, contrast = c("condition", "Cancer", "Hepatocyte"))
res.MRD.NT <- results(dds.MRD.NT, contrast = c("condition", "Cancer", "Hepatocyte"))

res.CD
res.PreT
res.MWD.T
res.MRD.T
res.MRD.NT

write.csv(as.data.frame(res.CD), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/HepvCancer/CD.HepvCancer.IPA.csv")
write.csv(as.data.frame(res.PreT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/HepvCancer/PreT.HepvCancer.IPA.csv")
write.csv(as.data.frame(res.MWD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/HepvCancer/MWD.T.HepvCancer.IPA.csv")
write.csv(as.data.frame(res.MRD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/HepvCancer/MRD.T.HepvCancer.IPA.csv")
write.csv(as.data.frame(res.MRD.NT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/HepvCancer/MRD.NT.HepvCancer.IPA.csv")

####################################################################################################################################################################################
##use this section for additional IPA analyses in different group comparisons than standard 
library(tidyverse)
library(DESeq2)

MRD.NTvMRD.T.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/MRD.NTvMRD.T.Cancer.counts.csv", header = TRUE, sep = ",")
PreTvMRD.NT.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/PreTvMRD.NT.Cancer.counts.csv", header = TRUE, sep = ",")
CDvMRD.NT.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMRD.NT.Cancer.counts.csv", header = TRUE, sep = ",")

head(MRD.NTvMRD.T.counts)
head(PreTvMRD.NT.counts)
head(CDvMRD.NT.counts)

MRD.NTvMRD.T.counts$X[duplicated(MRD.NTvMRD.T.counts$X)]
PreTvMRD.NT.counts$X[duplicated(PreTvMRD.NT.counts$X)]
CDvMRD.NT.counts$X[duplicated(CDvMRD.NT.counts$X)]

which(duplicated(MRD.NTvMRD.T.counts$X))
which(duplicated(PreTvMRD.NT.counts$X))
which(duplicated(CDvMRD.NT.counts$X))

MRD.NTvMRD.T.counts$X[9519] <- "1-Mar_dup"
MRD.NTvMRD.T.counts$X[18441] <- "2-Mar_dup"
PreTvMRD.NT.counts$X[9519] <- "1-Mar_dup"
PreTvMRD.NT.counts$X[18441] <- "2-Mar_dup"
CDvMRD.NT.counts$X[9519] <- "1-Mar_dup"
CDvMRD.NT.counts$X[18441] <- "2-Mar_dup"

which(duplicated(MRD.NTvMRD.T.counts$X))
which(duplicated(PreTvMRD.NT.counts$X))
which(duplicated(CDvMRD.NT.counts$X))

row.names(MRD.NTvMRD.T.counts) <- MRD.NTvMRD.T.counts$X
row.names(PreTvMRD.NT.counts) <- PreTvMRD.NT.counts$X
row.names(CDvMRD.NT.counts) <- CDvMRD.NT.counts$X

MRD.NTvMRD.T.counts <- MRD.NTvMRD.T.counts[,-1]
PreTvMRD.NT.counts <- PreTvMRD.NT.counts[,-1]
CDvMRD.NT.counts <- CDvMRD.NT.counts[,-1]

colData.MRD.NTvMRD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/MRD.NTvMRD.T.Cancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.PreTvMRD.NT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/PreTvMRD.NT.Cancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.CDvMRD.NT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMRD.NT.Cancer.data.info.csv", header = TRUE, sep = ",", row.names = 1)

View(colData.MRD.NTvMRD.T)
View(colData.PreTvMRD.NT)
View(colData.CDvMRD.NT)

all(names(MRD.NTvMRD.T.counts) %in% row.names(colData.MRD.NTvMRD.T))
all(names(PreTvMRD.NT.counts) %in% row.names(colData.PreTvMRD.NT))
all(names(CDvMRD.NT.counts) %in% row.names(colData.CDvMRD.NT))

all(names(MRD.NTvMRD.T.counts) == row.names(colData.MRD.NTvMRD.T))
all(names(PreTvMRD.NT.counts) == row.names(colData.PreTvMRD.NT))
all(names(CDvMRD.NT.counts) == row.names(colData.CDvMRD.NT))

##construct DESeqDataSet object
dds.MRD.NTvMRD.T <- DESeqDataSetFromMatrix(countData = MRD.NTvMRD.T.counts, colData = colData.MRD.NTvMRD.T, design = ~ condition)
dds.PreTvMRD.NT <- DESeqDataSetFromMatrix(countData = PreTvMRD.NT.counts, colData = colData.PreTvMRD.NT, design = ~ condition)
dds.CDvMRD.NT <- DESeqDataSetFromMatrix(countData = CDvMRD.NT.counts, colData = colData.CDvMRD.NT, design = ~ condition)

dds.MRD.NTvMRD.T
dds.PreTvMRD.NT
dds.CDvMRD.NT

##pre-filtering to remove rows with low gene counts
keep.MRD.NTvMRD.T <- rowSums(counts(dds.MRD.NTvMRD.T)) >= 3
dds.MRD.NTvMRD.T <- dds.MRD.NTvMRD.T[keep.MRD.NTvMRD.T,]
dds.MRD.NTvMRD.T

keep.PreTvMRD.NT <- rowSums(counts(dds.PreTvMRD.NT)) >= 3
dds.PreTvMRD.NT <- dds.PreTvMRD.NT[keep.PreTvMRD.NT,]
dds.PreTvMRD.NT

keep.CDvMRD.NT <- rowSums(counts(dds.CDvMRD.NT)) >= 3
dds.CDvMRD.NT <- dds.CDvMRD.NT[keep.CDvMRD.NT,]
dds.CDvMRD.NT

##run DESeq
dds.MRD.NTvMRD.T <- DESeq(dds.MRD.NTvMRD.T)
dds.PreTvMRD.NT <- DESeq(dds.PreTvMRD.NT)
dds.CDvMRD.NT <- DESeq(dds.CDvMRD.NT)

##Set up for the contrast function is the variable (condition), one to test (T1/high), and reference (C1/low)
res.MRD.NTvMRD.T <- results(dds.MRD.NTvMRD.T, contrast = c("condition", "MRD.T", "MRD.NT"))
res.PreTvMRD.NT <- results(dds.PreTvMRD.NT, contrast = c("condition", "MRD.NT", "PreT"))
res.CDvMRD.NT <- results(dds.CDvMRD.NT, contrast = c("condition", "MRD.NT", "CD"))

res.MRD.NTvMRD.T
res.PreTvMRD.NT
res.CDvMRD.NT

write.csv(as.data.frame(res.MRD.NTvMRD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/non.ICs/MRD.NTvMRD.T.Cancer.IPA.csv")
write.csv(as.data.frame(res.PreTvMRD.NT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/non.ICs/PreTvMRD.NT.Cancer.IPA.csv")
write.csv(as.data.frame(res.CDvMRD.NT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/non.ICs/CDvMRD.NT.Cancer.IPA.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##use this for additonal CD comparisons against other groups to detect malignancy
CDvMWD.T.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMWD.T.Hep.counts.csv", header = TRUE, sep = ",")
CDvMRD.T.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMRD.T.Hep.counts.csv", header = TRUE, sep = ",")
CDvMRD.NT.counts <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMRD.NT.Hep.counts.csv", header = TRUE, sep = ",")

head(CDvMWD.T.counts)
head(CDvMRD.T.counts)
head(CDvMRD.NT.counts)

CDvMWD.T.counts$X[duplicated(CDvMWD.T.counts$X)]
CDvMRD.T.counts$X[duplicated(CDvMRD.T.counts$X)]
CDvMRD.NT.counts$X[duplicated(CDvMRD.NT.counts$X)]

which(duplicated(CDvMWD.T.counts$X))
which(duplicated(CDvMRD.T.counts$X))
which(duplicated(CDvMRD.NT.counts$X))

CDvMWD.T.counts$X[9519] <- "1-Mar_dup"
CDvMWD.T.counts$X[18441] <- "2-Mar_dup"
CDvMRD.T.counts$X[9519] <- "1-Mar_dup"
CDvMRD.T.counts$X[18441] <- "2-Mar_dup"
CDvMRD.NT.counts$X[9519] <- "1-Mar_dup"
CDvMRD.NT.counts$X[18441] <- "2-Mar_dup"

which(duplicated(CDvMWD.T.counts$X))
which(duplicated(CDvMRD.T.counts$X))
which(duplicated(CDvMRD.NT.counts$X))

row.names(CDvMWD.T.counts) <- CDvMWD.T.counts$X
row.names(CDvMRD.T.counts) <- CDvMRD.T.counts$X
row.names(CDvMRD.NT.counts) <- CDvMRD.NT.counts$X

CDvMWD.T.counts <- CDvMWD.T.counts[,-1]
CDvMRD.T.counts <- CDvMRD.T.counts[,-1]
CDvMRD.NT.counts <- CDvMRD.NT.counts[,-1]

colData.CDvMWD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMWD.T.Hep.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.CDvMRD.T <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMRD.T.Hep.data.info.csv", header = TRUE, sep = ",", row.names = 1)
colData.CDvMRD.NT <- read.csv("T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/DESeq2/HCCp1/Non.ICs/Comparisons/CDvMRD.NT.Hep.data.info.csv", header = TRUE, sep = ",", row.names = 1)

View(colData.CDvMWD.T)
View(colData.CDvMRD.T)
View(colData.CDvMRD.NT)

all(names(CDvMWD.T.counts) %in% row.names(colData.CDvMWD.T))
all(names(CDvMRD.T.counts) %in% row.names(colData.CDvMRD.T))
all(names(CDvMRD.NT.counts) %in% row.names(colData.CDvMRD.NT))

all(names(CDvMWD.T.counts) == row.names(colData.CDvMWD.T))
all(names(CDvMRD.T.counts) == row.names(colData.CDvMRD.T))
all(names(CDvMRD.NT.counts) == row.names(colData.CDvMRD.NT))

##construct DESeqDataSet object
dds.CDvMWD.T <- DESeqDataSetFromMatrix(countData = CDvMWD.T.counts, colData = colData.CDvMWD.T, design = ~ condition)
dds.CDvMRD.T <- DESeqDataSetFromMatrix(countData = CDvMRD.T.counts, colData = colData.CDvMRD.T, design = ~ condition)
dds.CDvMRD.NT <- DESeqDataSetFromMatrix(countData = CDvMRD.NT.counts, colData = colData.CDvMRD.NT, design = ~ condition)

dds.CDvMWD.T
dds.CDvMRD.T
dds.CDvMRD.NT

##pre-filtering to remove rows with low gene counts
keep.CDvMWD.T <- rowSums(counts(dds.CDvMWD.T)) >= 3
dds.CDvMWD.T <- dds.CDvMWD.T[keep.CDvMWD.T,]
dds.CDvMWD.T

keep.CDvMRD.T <- rowSums(counts(dds.CDvMRD.T)) >= 3
dds.CDvMRD.T <- dds.CDvMRD.T[keep.CDvMRD.T,]
dds.CDvMRD.T

keep.CDvMRD.NT <- rowSums(counts(dds.CDvMRD.NT)) >= 3
dds.CDvMRD.NT <- dds.CDvMRD.NT[keep.CDvMRD.NT,]
dds.CDvMRD.NT

##run DESeq
dds.CDvMWD.T <- DESeq(dds.CDvMWD.T)
dds.CDvMRD.T <- DESeq(dds.CDvMRD.T)
dds.CDvMRD.NT <- DESeq(dds.CDvMRD.NT)

##Set up for the contrast function is the variable (condition), one to test (T1/high), and reference (C1/low)
res.CDvMWD.T <- results(dds.CDvMWD.T, contrast = c("condition", "MWD.T", "CD"))
res.CDvMRD.T <- results(dds.CDvMRD.T, contrast = c("condition", "MRD.T", "CD"))
res.CDvMRD.NT <- results(dds.CDvMRD.NT, contrast = c("condition", "MRD.NT", "CD"))

res.CDvMWD.T
res.CDvMRD.T
res.CDvMRD.NT

write.csv(as.data.frame(res.CDvMWD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/non.ICs/CDvMWD.T.Hep.IPA.csv")
write.csv(as.data.frame(res.CDvMRD.T), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/non.ICs/CDvMRD.T.Hep.IPA.csv")
write.csv(as.data.frame(res.CDvMRD.NT), file = "T:/Microbiology and Immunology/ManjiliLab/Nick Koelsch/HCC/PhD RNAseq/IPA Outfile/HCCp1/non.ICs/CDvMRD.NT.Hep.IPA.csv")

