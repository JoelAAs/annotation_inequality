import urllib.request
import re
from format_mitab_complete_dates import filter_mitab, reform_to_bait_prey

TARGET_GO_TERMS = "" 

# --- Configuration ---
SPECIES = "HUMAN" 

def get_all_GO_snapshots_downloads(wildcards):
    # This print will appear ONLY after the checkpoint finishes
    # print("\n[STATUS] Checkpoint triggered! Reading the newly generated file list...")
    
    with open(checkpoints.get_GO_ftp_list.get().output[0]) as f:
        files = [line.strip() for line in f if line.strip()]
    
    # Print how many files it found before it starts downloading
    # print(f"[STATUS] Found {len(files)} historical snapshots. Queuing download jobs...\n")
    
    return expand("work_folder/data/dates/GO/snapshots/{filename}", filename=files)

# Checkpoint: Scrape the HTTPS server
checkpoint get_GO_ftp_list:
    output:
        "work_folder/data/dates/GO/snapshot_list.txt"
    message: 
        "--> [JOB] Scraping EBI HTTPS server for available snapshots..."
    run:
        url = f"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/{SPECIES}/"
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            html = response.read().decode('utf-8')
        
        # 1. Grab EVERY link on the page
        all_links = re.findall(r'href="([^"]+)"', html)
        
        valid_files = set()
        for link in all_links:
            filename = link.split("/")[-1]
            
            # 2. Check for modern files AND legacy files that are compressed
            if ("gaf" in filename or "gene_association" in filename) and filename.endswith(".gz"):
                
                # 3. NEW: Explicitly exclude isoforms, RNAs, and complexes
                if not any(x in filename for x in ["isoform", "rna", "complex", "plus", "ref"]):
                    valid_files.add(filename)
        
        # 4. Write them to the text file
        with open(output[0], "w") as f:
            for file_name in sorted(valid_files):
                f.write(file_name + "\n")

# Rule: Download using the HTTPS link
rule download_GO_snapshot:
    output:
        "work_folder/data/dates/GO/snapshots/{filename}"
    params:
        # FIXED: Single curly braces here so the wildcard evaluates correctly!
        url=lambda wildcards: f"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/{SPECIES}/{wildcards.filename}"
    message: 
        "--> [JOB] Downloading snapshot: {wildcards.filename}"
    shell:
        # Removed -q so wget will print helpful error logs if it ever fails
        "wget -O {output} {params.url}"

rule mark_GO_downloads_complete:
    input:
        get_all_GO_snapshots_downloads
    message: 
        "--> [JOB] All downloads complete! Touching final flag file."
    output:
        touch("work_folder/data/dates/GO/downloads_complete.txt")

def get_all_parsed_snapshots(wildcards):
    # print("\n[STATUS] Checkpoint triggered! Reading the newly generated file list...")
    with open(checkpoints.get_GO_ftp_list.get().output[0]) as f:
        files = [line.strip() for line in f if line.strip()]
    # print(f"[STATUS] Found {len(files)} clean historical snapshots. Queuing jobs...\n")
    return expand("work_folder/data/dates/GO/parsed/{filename}.csv", filename=files)

rule extract_GO_snapshot_dates:
    input:
        snapshot = "work_folder/data/dates/GO/snapshots/{filename}",
        genes = "work_folder/data/network/raw_networks/network_genes_translated.txt"
    output:
        "work_folder/data/dates/GO/parsed/{filename}.csv"
    params:
        target_gos = TARGET_GO_TERMS
    script:
        "../pyScripts/dates/GO/extract_go_snapshot_dates.py"

rule compute_GO_annotations_first_dates:
    input:
        get_all_parsed_snapshots
    output:
        "work_folder/data/dates/GO/GO_first_annotation_dates.csv"
    script:
        "../pyScripts/dates/GO/aggregate_GO_annotations_dates.py"

rule inject_dates_to_GO_network:
    input:
        network = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
        dates = "work_folder/data/dates/GO/GO_first_annotation_dates.csv",
        mapping = "work_folder/data/network/raw_networks/hgnc_complete_set.txt"
    output:
        network_with_dates = "work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates.pkl"
    script:
        "../pyScripts/dates/GO/inject_GO_annotations_dates.py"

rule assign_dates_to_GO_missing_ones:
    input:
        network_with_dates = "work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates.pkl",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        network_with_all_dates = "work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates_complete.pkl"
    script:
        "../pyScripts/dates/GO/assign_dates_to_GO_missing_ones.py"

rule format_miTab_with_complete_dates:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab = "work_folder/data/intact/human.txt"
    output:
        formated = "work_folder/data/intact/bait_prey_publications_no_gene_complete_dates.pq"
    run:
        mg = mygene.MyGeneInfo()
        mitab_df = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        #bait_prey_df = bait_prey_df.iloc[:100000,:]
        bait_prey_df.to_parquet(
            output.formated        
            )

rule add_entrez_gene_complete_dates:
    input:
        formated_in = "work_folder/data/intact/bait_prey_publications_no_gene_complete_dates.pq"
    output:
        formated_out = "work_folder/data/intact/bait_prey_publications_complete_dates.pq"
    script:
        "../get_entrez_id.R"        

rule add_edge_dates_to_GO_network:
    input:
        network_with_all_dates = "work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates_complete.pkl", 
        bp_publications = "work_folder/data/intact/bait_prey_publications_complete_dates.pq"
    output:
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    script:
        "../pyScripts/dates/GO/add_edge_dates_to_GO_network.py"

rule find_GO_nodes_with_top_5_annotations:
    input:
        network_with_all_dates = "work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates_complete.pkl", 
        coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        nodes_with_top_5_annotations_df = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_{depth}_cutoff_{cutoff}.csv",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_{depth}_cutoff_{cutoff}.pkl" 
    script: 
        "../pyScripts/dates/GO/find_GO_nodes_with_top_5_annotations.py"

rule compare_past_present_and_future_GO_networks:
    input: 
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        networks_statistics = "work_folder/data/dates/GO/network_statistics/{aspect}_depth_{depth}_cutoff_{cutoff}_networks_statistics.csv"
    script: 
        "../pyScripts/dates/GO/compare_past_present_and_future_GO_networks.py"

rule visualize_past_present_and_future_GO_networks_statistics:
    input: 
        networks_statistics = "work_folder/data/dates/GO/network_statistics/{aspect}_depth_{depth}_cutoff_{cutoff}_networks_statistics.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        network_statistics_total_plots = "work_folder/data/dates/GO/plots/{aspect}_depth_{depth}_cutoff_{cutoff}_total_plots.png",
        network_statistics_annotated_plots = "work_folder/data/dates/GO/plots/{aspect}_depth_{depth}_cutoff_{cutoff}_annotated_plots.png"
    script:
        "../pyScripts/plotting/visualize_past_present_and_future_GO_networks_statistics.py"

## TODO refine
rule calculate_GO_temporal_correlation:
    input:
        networks_statistics = "work_folder/data/dates/GO/network_statistics/{aspect}_depth_{depth}_cutoff_{cutoff}_networks_statistics.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        correlation_results = "work_folder/data/dates/GO/correlation/{aspect}_depth_{depth}_cutoff_{cutoff}_temporal_correlation_results.csv"
    script:
        "../pyScripts/dates/GO/calculate_GO_temporal_correlation.py"

## TODO
rule visualize_GO_temporal_correlation:
    input:
        correlation_results = "work_folder/data/dates/GO/correlation/{aspect}_depth_{depth}_cutoff_{cutoff}_temporal_correlation_results.csv"
    output:
        correlation_plot = "work_folder/data/dates/GO/plots/{aspect}_depth_{depth}_cutoff_{cutoff}_temporal_correlation_plot.png"
    script:
        "../pyScripts/dates/GO/visualize_GO_temporal_correlation.py"