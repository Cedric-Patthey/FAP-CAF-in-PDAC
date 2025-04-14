# FAP-CAF-in-PDAC

This repository gathers the code and metadata to reproduce the results presented in the publication by Cumming et al, 2025 vailable at [10.1158/0008-5472.CAN-23-3252](https://doi.org/10.1158/0008-5472.CAN-23-3252)

**Dissecting FAP+ Cell Diversity in Pancreatic Cancer Uncovers an Interferon-Response Subtype of Cancer-Associated Fibroblasts with Tumor-Restraining Properties** 

Joshua Cumming, Parniyan Maneshi, Mitesh Dongre, Tala Alsaed, Mohammad Javad Dehghan-Nayeri, Agnes Ling, Kristian Pietras, Cedric Patthey & Daniel Öhlund

Abstract:

Within the stroma of pancreatic ductal adenocarcinoma (PDAC), mesenchymal cells differentiate into cancer-associated fibroblast (CAF) subtypes that differentially mediate disease progression. Defining the regulatory mechanism and diversity of CAF subtypes could identify potential therapeutic strategies to harness the tumor suppressive activities of CAFs. To address this, we utilized single-cell RNA sequencing to profile fibroblast activation protein-alpha (FAP) expressing mesenchymal cells in human PDAC. The mesenchymal subpopulations in PDAC reflected mesenchymal cell heterogeneity found in the normal developing pancreas. In addition to characterizing inflammatory CAF (iCAF) and myofibroblastic CAF (myCAF) subpopulations in detail, the analysis uncovered a previously undescribed interferon-response CAF (ifCAF) subtype. Tumor-derived signals induced specific CAF subtypes from pancreatic stellate cells (PSCs) in an organoid-based co-culture model, and time-course experiments revealed regulatory mechanisms that govern subtype formation. STING agonists promoted an ifCAF phenotype in vivo and in vitro. Importantly, induction of an ifCAF phenotype suppressed tumor cell invasiveness and induced an anti-tumor phenotype in tumor-associated neutrophils. Together, this study resolves FAP+ stromal cell heterogeneity in PDAC and identifies an ifCAF subtype that can be induced to suppress pro-tumorigenic features of PDAC.  

The code is presented as R markdowns for each figure. The metadata files required to reproduce the analysis are included in the same directories.

The code required to prepare the reference genome and generate gene annotation files is found in the directory accessory_code. The analysis pipeline to generate the raw read count for bulk RNA-seq data is found under accessory_code \>  SmartSeq2_bulk_RNAseq_read_count_pipeline.

The raw and processed data are available on ArrayExpress and Swedish National Data service (SND). The specific links to six mouse datasets and one human dataset can be found on Biostudies at  <https://doi.org/10.6019/S-BSST1896>.
