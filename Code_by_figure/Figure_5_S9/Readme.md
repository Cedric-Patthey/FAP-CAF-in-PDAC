# Figure 5 and S9

Code and metadata to reproduce figures 5 and S9. The data relates to the mouse single cell RNA-seq experiment with 12 samples from in vivo MSA-2/vehicle-treated tumour-bearing mice.

The raw fastq files are available on ArrayExpress with acession number E-MTAB-14940.

The UMI count file `Orthotopic_in_vivo_MSA2_CR_count.h5` and the Seurat object with integrated data (`Ortho.Tumors.MSA2.integrated.RDS`) are available on ArrayExpress with acession number E-MTAB-14940.

## Code

-   **Orthotopic_in_vivo_MSA2_10x_scSeq_fastq_to_UMI_count.Rmd**

    Generation of UMI count tables from fastq files. Used the reference genome described in R markdown `Musmu_GRCm39_Ens109_26Feb2024.Rmd`.

-   **Figure_5_S9_data_analysis.Rmd**

    Normalization, scaling and integration of samples using Seurat. Characterization of cell types and analysys of MSA-2 effects on cell type composition and gene expression. This analysis requires the UMI count file `Orthotopic_in_vivo_MSA2_CR_count.h5` . The Seurat object generated in the first steps of this analysis can be loaded from the file `Ortho.Tumors.MSA2.integrated.RDS`.

## Accessory files (metadata)

-   **Orthotopic_in_vivo_MSA2_10x_scSeq_sample_metadata.txt**

    Metadata table with sample information compiled manually

