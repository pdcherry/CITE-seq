rule NAME:
    input: "path/to/inputfile", "path/to/other/inputfile"
    output: "path/to/outputfile", "path/to/another/outputfile"
    shell: "somecommand {input} {output}"

rule download-files
    input:  expand(["../data/inputfile{dataset}/*.{ext}"],
        dataset=DATASETS, ext=FORMATS)