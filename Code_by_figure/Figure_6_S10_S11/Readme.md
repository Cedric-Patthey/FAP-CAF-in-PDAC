# Figure 6, S10 and S11

Code and metadata to reproduce figures 6 and S10 to S11. The data relates to the folloing datasets:

-   Mouse single cell RNA-seq experiment with 12 samples from in vivo MSA-2/vehicle-treated tumour-bearing mice, as decribed for figure 5.

    The raw fastq files are available on ArrayExpress with acession number E-MTAB-14940.

    The UMI count file `Orthotopic_in_vivo_MSA2_CR_count.h5` and the Seurat object with integrated data (`Ortho.Tumors.MSA2.integrated.RDS`) are available on ArrayExpress with acession number E-MTAB-14940.

-   Mouse bulk RNA-seq experiment with 167 samples from co-cultures of PSC and Tumour cells (and mono-culture controls) , treated with DMXAA or H151

    The raw fastq files are available on ArrayExpress with acession number E-MTAB-14946.

    The raw read count file `cocul_STING_modulation_167s_HT_counts.txt` is available on ArrayExpress with acession number E-MTAB-14946.

-   Mouse bulk RNA-seq experiment with 18 samples from neutrophil cultures, treated with ifCAF conditionned media.

    The raw fastq files are available on ArrayExpress with acession number E-MTAB-14934.

    The raw read count file `Nphi_CM_18s_HT_counts_17jun2024.txt` is available on ArrayExpress with acession number E-MTAB-14934.

## Code

-   **Figure_6_S10_MSA2_in_vivo_CellChat.Rmd**

    Analysis of cell-cell signalling using CellChat and NicheNet. This analysis requires the Seurat object `Ortho.Tumors.MSA2.integrated.RDS` generated in Figure 5.

-   **cocul_STING_modulation_167s_fastq_to_read_count.Rmd**

    Generation of raw read count from fastq files for the bulk RNA-seq experiment with 167 samples treated with DMSO, DMXAA or H151. Used the reference genome described in R markdown `Musmu_GRCm39_Ens109_26Feb2024.Rmd`. This pipeline is described in R markdown `SS2_SE_FASTQ2RRC.Rmd` under Accessory_code \> SmartSeq2_bulk_RNAseq_read_count_pipeline. It requires the scripts `fastq2RRC.sh`, `symplicity.py`, `get_cdna_pos.py` and `get_fastq_not_from_id_list.py` available in the pipeline directory.

-   **Figure_6_S11_cocul_STING_modulation_data_analysis.Rmd**

    Differential expression analysis. This analysis requires the raw read count file `cocul_STING_modulation_167s_HT_counts.txt`.

-   **Neutrophil_CM_18s_fastq_to_read_count.Rmd**

    Generation of raw read count from fastq files for the bulk RNA-seq experiment with 18 samples of Neutrophil cultures. Used the reference genome described in R markdown `Musmu_GRCm39_Ens109_26Feb2024.Rmd`. This pipeline is described in R markdown `SS2_SE_FASTQ2RRC.Rmd` under Accessory_code \> SmartSeq2_bulk_RNAseq_read_count_pipeline. It requires the scripts `fastq2RRC.sh`, `symplicity.py`, `get_cdna_pos.py` and `get_fastq_not_from_id_list.py` available in the pipeline directory.

-   **Figure_6_S10_ifCAF_Neutrophil_CM data_analysis.Rmd**

    -   Differential expression analysis. This analysis requires the raw read count file `Nphi_CM_18s_HT_counts_17jun2024.txt`.

## Accessory files (metadata)

-   **cocul_STING_modulation_167s_sample_metadata.txt**

    Metadata table with sample information compiled manually for the STING modulation experiment

-   **Neutrophil_CM_18s_sample_metadata.txt**

    Metadata table with sample information compiled manually for the Neutrophil experiment

-   **cocul_STING_modulation_167s_12Nov2024_parameters.csv**

    Parameters table used for the `fastq2RRC.sh` pipeline on rthe STING modulation data.

-   **Nphi_CM_18s_parameters_17jun2024.csv**

    Parameters table used for the `fastq2RRC.sh` pipeline on the neutrophil data.

-   **Musmu_GRCm39.109_protcod_fd_Gene_md.txt**

    Gene metadata table with Gene ID and information related to gene models

