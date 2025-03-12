# Figures 1, 2 and S1 to S5

Code and metadata to reproduce figures 1, 2 and S1 to S5. The data relates to the Human single cell RNA-seq experiment on 6 patient PDAC samples.

The raw fastq files are available upon request via the Swedish National Data service (SND) at <https://doi.org/10.5878/0ehq-1434>

The UMI count files for the 6 samples and the two Seurat objects with integrated data ( `FAP.Mesenchymal.Souped.Singlets.rPCA.integrated.RDS` and `CAFs.rPCA.integrated.RDS`) are available on SND at <https://doi.org/10.5878/0ehq-1434>

## Code

-   **Human_PDAC_FAP_cells_6s_fastq_to_read_count.Rmd**

    Generation of UMI count tables from fastq files.

-   **Figure_1\_&\_2_Pre_Processing.Rmd**

    Normalization, scaling and integration of Samples using Seurat. Also includes re-analysis of publically available datsets by Peng et al, 2019 and Olaniru et al, 2023. This analysis requires the UMI count files named like `S1_filtered_feature_bc_matrix.zip` and `S1_raw_feature_bc_matrix.zip` (12 files in total for samples S1 to S6).

-   **Figure_1_S1_S2_S3_data_analysis.Rmd**

    Characterization of single cell clusters at low resolution and comparison to published datasets. Code to reproduce figures 1 and S1 to S3. This analysis requires the Seurat objects `FAP.Mesenchymal.Souped.Singlets.rPCA.integrated.RDS` generated in the pre-processing step.

-   **Figure_2_S4_S5_data_analysis.Rmd**

    Characterization of single cell clusters at high resolution.

