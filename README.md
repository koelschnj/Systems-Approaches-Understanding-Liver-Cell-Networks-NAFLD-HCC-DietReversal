# Systems-Approaches-Understanding-Liver-Cell-Networks-NAFLD-HCC-DietReversal
snRNA-seq data from various male mice on standard chow diet (CD), western diet with NAFLD (WD.nf/PreT), HCC (WD.t), and those that underwent diet reversal that developed HCC (RD.t) or were rescued from it (RD.n) were studied through systems biology perspectives to understand the complex and dynamic networks of immune and non-immune cells interacting in the liver microenvironment.

Primarily utilize Seurat for data handling, but also scSorter in conjunction with markers extracted from the CellMarker2.0 database to annotate samples accordingly.  Specific populations of cells were subjected to DESeq2 analysis to upload results to Ingenuity Pathway Analysis (IPA) to interpret results.

Major systems biology perspectives were employed through the use of CellChat, that enables great visualization and data interpretations methods of all cell populations detected in each group, as well as comparative analyses across groups to understand underlying differences.

Sorting steps and initial cell annotation was performed via the VCU HPRC (high performance research computer) core (High Performance Computing resources provided by the High Performance Research Computing (HPRC) core facility at Virginia Commonwealth University (https://hprc.vcu.edu) were used for conducting the research reported in this work.) In which the files are titled "Rscript_HCCp1_Step_number_" corresponding to the shell scripts titled "HCCp1_Step_number_" for job submission to the HPRC core servers.
