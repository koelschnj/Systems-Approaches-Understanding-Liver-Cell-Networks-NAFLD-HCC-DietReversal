# Systems-Approaches-Understanding-Liver-Cell-Networks-NAFLD-HCC-DietReversal
snRNA-seq data from various male mice on standard chow diet (CD), western diet with NAFLD (WD.nf/PreT), HCC (WD.t), and those that underwent diet reversal that developed HCC (RD.t) or were rescued from it (RD.n) were studied through systems biology perspectives to understand the complex and dynamic networks of immune and non-immune cells interacting in the liver microenvironment.

Primarily utilize Seurat for data handling, but also scSorter in conjunction with markers extracted from the CellMarker2.0 database to annotate samples accordingly.  Specific populations of cells were subjected to DESeq2 analysis to upload results to Ingenuity Pathway Analysis (IPA) to interpret results.

Major systems biology perspectives were employed through the use of CellChat, that enables great visualization and data interpretations methods of all cell populations detected in each group, as well as comparative analyses across groups to understand underlying differences.

Sorting steps and initial cell annotation was performed via the VCU HPRC (high performance research computer) core (High Performance Computing resources provided by the High Performance Research Computing (HPRC) core facility at Virginia Commonwealth University (https://hprc.vcu.edu) were used for conducting the research reported in this work.) In which the files are titled "Rscript_HCCp1_Step_number_" corresponding to the shell scripts titled "HCCp1_Step_number_" for job submission to the HPRC core servers.The R script titled, "local_CellChat_scSorter5_results.R" is a conglomeration of different sections of code for major analyses in Seurat and CellChat; this requires a stop-and-go approach for highlighting the desired sections of code to perform the analyses desired.

The versions of the main programs utilized include: R 4.4.1, Seurat 4.4.0, CellChat 2.1.2, future 1.33.2, patchwork 1.2.0, ggplot2 3.5.1, ComplexHeatmap 2.20.0

Seurat_Python_Conversion.R was used to generate AnnData objects for use in LIANA throught the VCU HPRC core. Shell scripts for executing python on slurm include Py_li_WDt2024.sh, Py_li_RDt2024.sh, and Py_li_RDn2024.sh for the python scripts LIANA_SSLRI_WDt2024.py, LIANA_SSLRI_RDt2024.py, and LIANA_SSLRI_RDn2024.py
