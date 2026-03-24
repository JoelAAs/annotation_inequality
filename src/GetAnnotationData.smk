import pandas as pd

rule query_disgenet:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_uniprot_ids.pq"
    output:
        annotation_df = "work_folder/data/disgenet/uniprot_to_disgenet.csv"
    script: "query_disgenet.R"


rule get_annotations_per_entrez_DISGENET:
    input:
        b_count = "work_folder/data/intact/bait_count.csv",
        p_count = "work_folder/data/intact/prey_count.csv",
        annotation_df = "work_folder/data/disgenet/uniprot_to_disgenet.csv"
    output:
        annotations_per_id_baits = "work_folder/data/disgenet/annotation_per_entrez_baits.csv",
        annotations_per_id_preys = "work_folder/data/disgenet/annotation_per_entrez_preys.csv"
    run:
        '''
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
        '''

        df_annot = pd.read_csv(input.annotation_df, sep ='\t')
        df_studies_baits = pd.read_csv(input.b_count, sep="\t")
        df_studies_preys = pd.read_csv(input.p_count, sep="\t")

        df_annot = df_annot[['geneid', 'diseaseUMLSCUI']]
        df_annot = df_annot.drop_duplicates(subset = ['geneid', 'diseaseUMLSCUI'])
        annot_count = df_annot.groupby('geneid').size().reset_index(name = 'count')
        annot_count = annot_count.dropna()
        annot_count = annot_count.rename(columns = {'geneid': 'entrez_id'})

        '''
        Bait annotation counts
        '''
        df = pd.merge(df_studies_baits, annot_count, left_on="entrez_id_bait", right_on="entrez_id", how="right", suffixes=("_studies", "_annot"))
        df = df[df['entrez_id_bait'].notna()]
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id_baits, sep="\t", index=False)

        '''
        Prey annotation counts
        '''
        df = pd.merge(df_studies_preys, annot_count, left_on="entrez_id_prey", right_on="entrez_id", how="right", suffixes=("_studies", "_annot"))
        df = df[df['entrez_id_prey'].notna()]
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id_preys, sep="\t", index=False)     


rule get_HDO:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df = "work_folder/data/HDO/entrez_to_HDO.csv"
    script:
        "query_HDO.R"

rule get_annotations_per_entrez_HDO:
    input:
        b_count = "work_folder/data/intact/bait_count.csv",
        p_count = "work_folder/data/intact/prey_count.csv",
        annotation_df = "work_folder/data/HDO/entrez_to_HDO.csv"
    output:
        annotations_per_id_baits = "work_folder/data/HDO/annotation_per_entrez_baits.csv",
        annotations_per_id_preys = "work_folder/data/HDO/annotation_per_entrez_preys.csv"
    run:
        df_annot = pd.read_csv(input.annotation_df, sep="\t")
        df_studies_baits = pd.read_csv(input.b_count, sep="\t")
        df_studies_preys = pd.read_csv(input.p_count, sep="\t")

        '''
        Bait annotation counts
        '''
        df = pd.merge(df_studies_baits, df_annot, left_on="entrez_id_bait", right_on="entrez_id", how="right", suffixes=("_studies", "_annot"))
        df = df[df['entrez_id_bait'].notna()]
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id_baits, sep="\t", index=False)

        '''
        Prey annotation counts
        '''
        df = pd.merge(df_studies_preys, df_annot, left_on="entrez_id_prey", right_on="entrez_id", how="right", suffixes=("_studies", "_annot"))
        df = df[df['entrez_id_prey'].notna()]
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id_preys, sep="\t", index=False)

rule get_GO:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df_BP = "work_folder/data/GO/entrez_to_GO_BP.csv",
        annotation_df_MF = "work_folder/data/GO/entrez_to_GO_MF.csv",
        annotation_df_CC = "work_folder/data/GO/entrez_to_GO_CC.csv"
    script: 
        "query_GO.py"

rule get_GO_with_given_depth:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df_BP = "work_folder/data/GO/entrez_to_GO_BP_depth_{depth}.csv",
        annotation_df_MF = "work_folder/data/GO/entrez_to_GO_MF_depth_{depth}.csv",
        annotation_df_CC = "work_folder/data/GO/entrez_to_GO_CC_depth_{depth}.csv"
    params:
        depth = "{depth}"
    script: 
        "query_GO_with_depth.py"
        
rule get_annotations_per_entrez_GO:
    input:
        b_count = "work_folder/data/intact/bait_count.csv",
        p_count = "work_folder/data/intact/prey_count.csv",
        annotation_df = "work_folder/data/GO/entrez_to_GO_{aspect}.csv"
    output:
        annotations_per_id_baits = "work_folder/data/GO/annotation_per_entrez_{aspect}_baits.csv",
        annotations_per_id_preys = "work_folder/data/GO/annotation_per_entrez_{aspect}_preys.csv"
    run:
        df_annot = pd.read_csv(input.annotation_df, sep="\t")
        df_studies_baits = pd.read_csv(input.b_count, sep="\t")
        df_studies_preys = pd.read_csv(input.p_count, sep="\t")

        '''
        Bait annotation counts
        '''
        df = pd.merge(df_studies_baits, df_annot, left_on="entrez_id_bait", right_on="entrez_id", how="right", suffixes=("_studies", "_annot"))
        df = df[df['entrez_id_bait'].notna()]
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id_baits, sep="\t", index=False)

        '''
        Prey annotation counts
        '''
        df = pd.merge(df_studies_preys, df_annot, left_on="entrez_id_prey", right_on="entrez_id", how="right", suffixes=("_studies", "_annot"))
        df = df[df['entrez_id_prey'].notna()]
        df.fillna(0, inplace=True)
        df.to_csv(output.annotations_per_id_preys, sep="\t", index=False)        
