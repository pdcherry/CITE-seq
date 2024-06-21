import os

def get_accession_ids(txt_file):
    accession_file  = open(txt_file, 'r')
    accession_id_list = []

    accession_id_list = accession_file.readlines()
    
    count = 0
    for line in accession_id_list:
        count += 1
        print("Line{}: {}".format(count, line.strip()))

    accession_file.close()
    
    return(accession_id_list)

get_accession_ids('/Users/patrick/bfx git projects/CITE-seq/accession_id.txt')

rule NAME:
    input: "path/to/inputfile", "path/to/other/inputfile"
    output: "path/to/outputfile", "path/to/another/outputfile"
    shell: "somecommand {input} {output}"

rule download-files
    input: "accession_id.txt"
    output:  expand(["../data/inputfile{dataset}/*.{ext}"],
        dataset=DATASETS, ext=FORMATS)
    