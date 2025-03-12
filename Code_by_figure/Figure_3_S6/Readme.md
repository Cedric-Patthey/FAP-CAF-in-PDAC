# Figures 3 and S6

Code and metadata to reproduce figures 3 and S6. The data relates to the mouse single cell RNA-seq experiment with 11 samples from co-cultures of PSC and Tumour cells (and mono-culture controls) at day 6.

The raw fastq files are available on ArrayExpress with acession number E-MTAB-14944

The UMI count file `PSC_cocul_day6_11s_Raw_UMI_Count.h5` and the Seurat object with integrated data ( `PSC.singlets.rPCA.integrated.RDS`) are available on ArrayExpress with acession number E-MTAB-14944

## Code

-   **PSC_cocul_day6_11s_scSeq_fastq_to_count.Rmd**

    Generation of UMI count tables from fastq files. This pipeline requires scripts `ddSeq_fastq_to_RUC.sh` and `make_ddSeq_RUC_table.py` available in the pipeline directory. The barcode whitelist from Biorad `Whitelist_96_bc.txt` is also required. Used the reference genome described in R markdown `Musmu_GRCm38_Ens100_10May2020.Rmd`.

-   **Figure_3_S6_data_analysis.Rmd**

    Normalization, scaling and integration of samples using Seurat. Characterization of cell types. This analysis requires the UMI count file `PSC_cocul_day6_11s_Raw_UMI_Count.h5`. The Seurat object generated in the first steps of this analysis can be loaded from the file `PSC.singlets.rPCA.integrated.RDS`.

## Accessory files (metadata)

-   **PSC_cocul_day6_11s_sample_metadata.txt**

    Metadata table with sample information compiled manually

-   **Whitelist_96_bc.txt**

    List of 96 white-listed barcodes from the Surecell library prep method, used for barcode deconvolution.

