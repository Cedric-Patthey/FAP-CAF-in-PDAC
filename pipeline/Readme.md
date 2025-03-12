# **Scripts for read count generation**

## *Bulk RNA-seq pipeline for Smart-Seq2 data*

-   **fastq2RRC.sh**

    Main bash script processing fastq files to raw read count via quality filtering, mapping to reference genome and read count

-   **symplicity.py**

    Python script removing exact duplicate reads in a coverage-dependent manner

-   **get_cdna_pos.py**

    Python script calculating 5'-3' distribution based on an alignment bam file and annotation bed file

-   **get_fastq_not_from_id_list.py**

    Python script to remove fastq entries from a fastq file based on a list

## Single cell RNA-seq pipeline for ddSeq/Surecell data

-   **ddSeq_fastq_to_RUC.sh**

    Main bash script processing fastq files to raw UMI count via quality filtering, barcode deconvolution, mapping to reference genome and read count

-   **make_ddSeq_RUC_table.py**

    Python script operating barcode deconvolution and parsing BEDTools intersection output to compile a single cell UMI count table

