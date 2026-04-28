checkpoint create_bait_prey_publications_network:
    input: 
        bp_publications = "work_folder/data/intact/bait_prey_publications.pq",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output: 
        degree_frequencies = "work_folder/data/network/raw_networks/raw_network_degree_frequencies.csv",
        bp_network = "work_folder/data/network/raw_networks/bait_prey_publications_network.pkl"
    script: 
        "../pyScripts/network/create_bait_prey_publications_network.py"

rule plot_bait_prey_publications_network_degree_frequencies:
    input: 
        degree_frequencies = "work_folder/data/network/raw_networks/raw_network_degree_frequencies.csv"
    output: 
        degree_frequencies_plot = "work_folder/data/network/raw_networks/raw_network_degree_frequencies.png"
    script: 
        "../pyScripts/plotting/plot_bait_prey_publications_network_degree_frequencies.py"   

rule extract_raw_network_genes:
    input: 
        bp_network = "work_folder/data/network/raw_networks/bait_prey_publications_network.pkl"        
    output: 
        network_genes = "work_folder/data/network/raw_networks/network_genes.txt"
    run: 
        import pickle

        with open(input.bp_network, 'rb') as f:
            network_data = pickle.load(f)

        with open(output.network_genes, 'w') as f:
            for gene_id in network_data.nodes():
                f.write(f'{gene_id}\n')

        print(f"[STATUS] Extracted {len(network_data)} genes from the network!")

rule download_hgnc_mapping:
    output:
        "work_folder/data/network/raw_networks/hgnc_complete_set.txt"
    shell:
        "wget -O {output} https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

rule translate_network_genes:
    input:
        network_genes = "work_folder/data/network/raw_networks/network_genes.txt",
        mapping = "work_folder/data/network/raw_networks/hgnc_complete_set.txt"
    output:
        "work_folder/data/network/raw_networks/network_genes_translated.txt"
    script:
        "../pyScripts/network/translate_network_gene_ids.py"