from format_mitab import filter_mitab, reform_to_bait_prey 
import networkx as nx
import pandas as pd


rule get_intact:
    output:
        intact="work_folder/data/intact/human.txt"
    shell:
        """
        wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.txt -O {output.intact}
        """

rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab = "work_folder/data/intact/human.txt"
    output:
        formated = "work_folder/data/intact/bait_prey_publications.pq"
    run:
        mitab_df     = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        bait_prey_df.to_parquet(
            output.formated,
            index = None
        )

def detect_and_join_isoform(x):
    split = x.split("-")
    if len(split) > 1:
        if not split[1].isdigit():
            return "-".join(split)
    return split[0]
        
    
rule get_bait_prey_usage:
    input:
        bp_df = "work_folder/data/intact/bait_prey_publications.pq"
    output:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequecies.pq"
    run:
        bp_df = pd.read_parquet(input.bp_df)
        bp_df["uniprot_id_bait"] = bp_df["uniprot_id_bait"].apply(detect_and_join_isoform)
        bp_df["uniprot_id_prey"] = bp_df["uniprot_id_prey"].apply(detect_and_join_isoform)
        G_bp = nx.from_pandas_edgelist(bp_df,"uniprot_id_bait", "uniprot_id_prey",create_using=nx.DiGraph)
        degree_prey = pd.DataFrame(G_bp.in_degree())
        degree_prey["type"] = "prey"
        degree_bait = pd.DataFrame(G_bp.out_degree())
        degree_bait["type"] = "bait"

        bp_frequencies = pd.concat([degree_bait, degree_prey], ignore_index=True)
        bp_frequencies.columns = ["uniprot_id", "count", "type"]
        bp_frequencies.to_parquet(
            output.bp_frequencies
        )

rule get_bait_count_per_study:
    input:
        bp_df = "work_folder/data/intact/bait_prey_publications.pq"
    output:
        b_frequencies = "work_folder/data/intact/bait_count.csv"
    run:
        bp_df = pd.read_parquet(input.bp_df)
        bp_df["uniprot_id_bait"] = bp_df["uniprot_id_bait"].apply(detect_and_join_isoform)
        bp_df["uniprot_id_prey"] = bp_df["uniprot_id_prey"].apply(detect_and_join_isoform)
        
        b_count_per_study = bp_df.groupby(["uniprot_id_bait", "pubmed_id", "detection_method"], as_index=False).size()
        number_of_studies_per_bait_prey = b_count_per_study.groupby(["uniprot_id_bait"], as_index=False).size()
        number_of_studies_per_bait_prey.columns = ["uniprot_id_bait", "count"]

        number_of_studies_per_bait_prey.to_csv(output.b_frequencies, sep="\t", index=False)