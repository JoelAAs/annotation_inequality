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