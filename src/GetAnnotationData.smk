
rule query_disgenet:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df = "work_folder/data/disgenet/uniprot_to_disgenet.csv"
    script: "query_disgenet.R"


rule get_annotations_per_id:
    input:
        annotation_df = "work_folder/data/disgenet/uniprot_to_disgenet.csv"
    output:
        annotations_per_id = "work_folder/data/disgenet/annotation_per_uniprot.csv"
    run:
        current_id = "uniprotids"
        current_annotations =  []
    
        with open(output.annotations_per_id, "w") as w:
            with open(input.annotation_df, "r") as f:
                header=True
                for line in f:
                    values = line.strip().split("\t")

                    if current_id != [values1]:
                        if not header:
                            w.write(f"{current_id}\t{len(set(current_annotations))}\n")
                        else:
                            w.write("uniprotid\tannotation_count\n")
                            header=False
                        current_id = values[1]
                        current_annotations = [values[3],]
                    else:
                        current_annotations.append(values[3]) # TODO last protein is not written


rule get_HDO:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df = "work_folder/data/HDO/entrez_to_HDO.csv"
    script:
        "query_HDO.R"

rule get_HDO_per_entrez:
    input:
        b_count = "work_folder/data/intact/bait_count.csv",
        annotation_df = "work_folder/data/HDO/entrez_to_HDO.csv"
    output:
        annotations_per_id = "work_folder/data/HDO/annotation_per_entrez.csv"
    run:
        df_freq = pd.read_csv(input.b_count, sep="\t")
        df_annot = pd.read_csv(input.annotation_df, sep="\t")

        df = pd.merge(df_freq, df_annot, left_on="entrez_id_bait", right_on="entrez_id", how="right", suffixes=("_degree", "_annot"))
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id, sep="\t", index=False) 
        # TODO compute the correlation and pvalue, also save plots