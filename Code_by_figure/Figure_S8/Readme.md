# Figure S8

Code and metadata to reproduce figures S8. The data relates to the mouse bulk RNA-seq experiment with 68 samples from co-cultures of PSC and Tumour cells (and mono-culture controls) , treated with DMXAA or Interferon gamma.

The raw fastq files are available on ArrayExpress with acession number E-MTAB-14943.

The raw read count file `ifCAF_32s_apCAF_36s_HT_counts_long.txt` is available on ArrayExpress with acession number E-MTAB-14943.

## Code

-   **ifCAF_32s_apCAF_36s_fastq_to_read_count.Rmd**

    Generation of raw read count from fastq files for the bulk RNA-seq experiment with 68 samples treated with DMSO, DMXAA or Interferon gamma. Used the reference genome Musmu_GRCm39.104 generated as described in R markdown `Musmu_GRCm39_Ens109_26Feb2024.Rmd`. This pipeline is described in R markdown `SS2_SE_FASTQ2RRC.Rmd` under Accessory_code \> SmartSeq2_bulk_RNAseq_read_count_pipeline. It requires the scripts `fastq2RRC.sh`, `symplicity.py`, `get_cdna_pos.py` and `get_fastq_not_from_id_list.py` available in the pipeline directory.

-   **Figure_S8_data_analysis.Rmd**

    Differential expression analysis. This analysis requires the raw read count file `ifCAF_32s_apCAF_36s_HT_counts_long.txt`.

## Accessory files (metadata)

-   **ifCAF_32s_apCAF_36s_sample_metadata.txt**

    Metadata table with sample information compiled manually for the ifCAF/apCAF experiment

-   **ifCAF_32s_apCAF_36s_parameters.csv**

    Parameters table used for the `fastq2RRC.sh` pipeline on the ifCAF/apCAF data.

-   **Musmu_GRCm39.104_protcod_fd_Gene_md.txt**

    Gene metadata table with Gene ID and information related to gene models

-   **ISMARA_activity_table_Mono_PSC_DMXAA.txt**

    Transcription factor activity table output by the ISMARA software, for the comparison of DMXAA-treated vs DMSO control PSC mono-cultures

-   **ISMARA_activity_table_Mono_PSC_IFNg.txt**

    Transcription factor activity table output by the ISMARA software, for the comparison of IFNg-treated vs DMSO control PSC mono-cultures

-   **ISMARA_delta_table_Mono_PSC_DMXAA.txt**

    Transcription factor activity delta table output by the ISMARA software, for the comparison of DMXAA-treated vs DMSO control PSC mono-cultures

-   **ISMARA_delta_table_Mono_PSC_IFNg.txt**

    Transcription factor activity delta table output by the ISMARA software, for the comparison of IFNg-treated vs DMSO control PSC mono-cultures

