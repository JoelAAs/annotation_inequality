include: "src/BaitUsage.smk"
include: "src/GetAnnotationData.smk"

rule all:
    input:
        "work_folder/data/disgenet/anotation_per_uniprot.csv"