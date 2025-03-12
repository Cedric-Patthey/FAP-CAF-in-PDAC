# Figure 4 and S7

Code and metadata to reproduce figures 4 and S7. The data relates to the mouse single cell RNA-seq experiment with 15 samples from co-cultures of PSC and Tumour cells (and mono-culture controls) over a time course.

The raw fastq files are available on ArrayExpress with acession number E-MTAB-14945

The UMI count file `PSC_timecourse_15s_Raw_UMI_Count.h5` and two Seurat objects with integrated data ( `Time.course.no.day.0.Clustered.RDS` and `Time.course.Clustered.RDS`) are available on ArrayExpress with acession number E-MTAB-14945

## Code

-   **PSC_timecourse_15s_ddSeq_fastq_to_UMI_count.Rmd**

    Generation of UMI count tables from fastq files. This pipeline requires scripts `ddSeq_fastq_to_RUC.sh` and `make_ddSeq_RUC_table.py` available in the pipeline directory. The barcode whitelist from Biorad `Whitelist_96_bc.txt` is also required. Used the reference genome described in R markdown `Musmu_GRCm38_Ens100_10May2020.Rmd`.

-   **Figure_4_S7_data_analysis.Rmd**

    Normalization, scaling and integration of samples using Seurat. Characterization of cell types and differentiation trajectories. This analysis requires the UMI count file `PSC_timecourse_15s_Raw_UMI_Count.h5` . The Seurat object generated in the first steps of this analysis can be loaded from the file `Time.course.no.day.0.Clustered.RDS`.

## Accessory files (metadata)

-   **PSC_timecourse_15s_sample_metadata.txt**

    Metadata table with sample information compiled manually

-   **Whitelist_96_bc.txt**

    List of 96 white-listed barcodes from the Surecell library prep method, used for barcode deconvolution.

