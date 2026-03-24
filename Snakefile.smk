include: "src/BaitUsage.smk"
include: "src/GetAnnotationData.smk"
include: "src/DfJoin.smk"
include: "src/PrintCorrelationValuesPlots.smk"
include: "src/ComputeCorrelationValues.smk"
include: "src/FeatureMatrix.smk"
include: "src/ElasticNet.smk"
include: "src/AnnotationDepth.smk"

rule all:
    input:
        "work_folder/data/disgenet/anotation_per_uniprot.csv",
        