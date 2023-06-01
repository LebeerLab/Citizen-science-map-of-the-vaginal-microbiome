# The Isala cross-sectional study

This repository contains all scripts that were used for the data processing and analysis of the Isala cross-sectional study. This analysis and its results are described in the paper "Citizen-science map of the vaginal microbiome", a preprint of which you can find on [Research Square](https://doi.org/10.21203/rs.3.rs-1350465/v1). 

The scripts are divided into six workflows: three for the main data analysis of the project and three that are relevant to the main project but not directly part of it (like the creation and validation of Lactobacillus subgenera). 

If you want to run these scripts yourself, you'll have to change the paths to the necessary input files. In all scripts, these paths are listed at the top of the script. 

## Main workflows

These workflows can be found in `scripts-main`. 

* `01-process-amplicon-seqs`: **scripts for the processing of the raw amplicon sequencing data**
    * Output: read count tables for taxa (ASVs) in all samples. 
    * The scripts are given for one of the nine MiSeq runs of the project. The same parameter settings and reference database were used for the other runs. 
* `02-process-read-counts`: **scripts for the processing of the read count tables**
    * Output: preprocessed read count data ready to be analyzed. 
    * The main steps are: the integration of the read count tables from the nine sequencing runs; the adding of metadata on the wetlab processing and sequencing (e.g., volumes of samples that were pooled into the final library); the reclassification of Lactobacillaceae sequences to implement [the taxonomic update of 2020](https://doi.org/10.1099/ijsem.0.004107); the reclassification of Lactobacillus sequences to the subgenera that we defined for this project; the quality control and filtering of taxa and samples. 
* `03-analyze-data`: **scripts for the actual data analysis**
    * Output: the main figures and tables of the paper. 

## Side workflows

These workflows can be found in `scripts-side`. 

* `01-create-amplicon-dbs`: **scripts for the creation of three amplicon databases used by the project**
    * Outputs: 
        * Nine monophyletic subgenera of the genus Lactobacillus. 
        * Three custom 16S reference databases for various purposes. All three databases are contain 16S sequences from the GTDB and use the GTDB taxonomy. 
    * The databases: 
        1. All bacteria: used to classify the amplicon samples of the Isala pilot study to the GTDB taxonomy. 
        2. Lactobacillaceae: used to reclassify the Lactobacillaceae reads of the cross-sectional amplicon samples to the [new Lactobacillaceae taxonomy](https://doi.org/10.1099/ijsem.0.004107).
        3. Lactobacillus: used to classify the Lactobacillus reads of the cross-sectional and pilot study samples to the custom-defined Lactobacillus subgenera. 
* `02-create-metagenome-db`: **scripts for the creation of a reference database to classify the metagenome reads**
    * Output: a kraken2 database that can be used to classify metagenome samples of the Isala pilot study to the GTDB taxonomy. 
* `03-validate-subgenera`: **scripts for the validation of the Lactobacillus subgenera**
    * Output: figures and tables that show that V4 16S sequences can be classified to the Lactobacillus subgenera that were defined for this project. 

## Notes 

We also make the custom amplicon databases available; all of them are based on 16S sequences from the GTDB and on the GTDB taxonomy: 

* [All bacteria](https://github.com/SWittouck/datasets/raw/main/SSUrRNA_GTDB05-bac-rep_DADA2.fna.gz) 
* [Lactobacillaceae only](https://github.com/SWittouck/datasets/raw/main/SSUrRNA_GTDB05-lactobacillaceae-all_DADA2.fna.gz)
* [Lactobacillus only, with Lactobacillus replaced with subgenera](https://github.com/SWittouck/datasets/raw/main/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna.gz)

If you want to reclassify your own ASV data to the new Lactobacillus taxonomy and/or the Lactobacillus subgenera, you can find example code in `scripts-main/02-process-read-counts/02_reclassify_lactos.R`. 

These are the type species we chose to define the Lactobacillus subgenera: 

| subgenus                         | type species                |
|---------------------------------|-----------------------------|
| Lactobacillus crispatus group   | Lactobacillus crispatus     |
| Lactobacillus iners group       | Lactobacillus iners         |
| Lactobacillus jensenii group    | Lactobacillus jensenii      |
| Lactobacillus gasseri group     | Lactobacillus gasseri       |
| Lactobacillus apis group        | Lactobacillus apis          |
| Lactobacillus pasteurii group   | Lactobacillus pasteurii     |
| Lactobacillus sp002418055 group | Lactobacillus sp002418055   |
| Lactobacillus delbrueckii group | Lactobacillus delbrueckii   |
| Lactobacillus sp002417825 group | Lactobacillus sp002417825   |
