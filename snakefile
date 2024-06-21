rule download-files:
    input:
    output: "./data/GSE100866/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz"
    shell: "Rscript code/CITE-seq_data_fetch.R"
    