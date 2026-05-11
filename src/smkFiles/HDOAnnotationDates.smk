## TODO retrieve the gene-doid association dates, inject them to the network and add also edge dates to it (there should already be the function) look on Gemini
checkpoint get_HDO_disease_gene_annotation_dates:
    output:
        europepmc_dir = directory("work_folder/data/dates/HDO/opentargets/evidence")
    shell:
        """
        mkdir -p {output.europepmc_dir}

        wget -e robots=off -r -np -nd -R "index.html*" https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_europepmc/ -P {output.europepmc_dir}
        """ 

rule get_opentargets_targets:
    output:
        targets_dir = directory("work_folder/data/dates/HDO/opentargets/targets")
    shell:
        """
        mkdir -p {output.targets_dir}
        wget -e robots=off -r -np -nd -R "index.html*" https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/target/ -P {output.targets_dir}
        """

rule get_opentargets_diseases:
    output:
        diseases_dir = directory("work_folder/data/dates/HDO/opentargets/diseases")
    shell:
        """
        mkdir -p {output.diseases_dir}
        wget -e robots=off -r -np -nd -R "index.html*" https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/disease/ -P {output.diseases_dir}
        """

rule download_clinvar:
    output:
        "work_folder/data/clinvar/submission_summary.txt.gz"
    shell:
        "wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"

rule compute_HDO_annotations_first_dates:
    input:
        evidence_folder = "work_folder/data/dates/HDO/opentargets/evidence",
        target_folder = "work_folder/data/dates/HDO/opentargets/targets",
        disease_folder = "work_folder/data/dates/HDO/opentargets/diseases"
    output:
        first_annotation_dates = "work_folder/data/dates/HDO/HDO_first_annotation_dates.csv"
    script:
        "../pyScripts/dates/HDO/compute_HDO_annotations_dates.py"

rule inject_dates_to_HDO_network:
    input:
        network = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        dates = "work_folder/data/dates/HDO/HDO_first_annotation_dates.csv"
    output:
        network_with_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates.pkl"
    script:
        "../pyScripts/dates/HDO/inject_HDO_annotations_dates.py"


rule assign_dates_to_HDO_missing_ones:
    input:
        network_with_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates.pkl",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        network_with_all_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates_complete.pkl"
    script:
        "../pyScripts/dates/HDO/assign_dates_to_HDO_missing_ones.py"

rule add_edge_dates_to_HDO_network:
    input:
        network_with_all_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates_complete.pkl", 
        bp_publications = "work_folder/data/intact/bait_prey_publications_complete_dates.pq"
    output:
        final_network = "work_folder/data/dates/HDO/networks_with_dates/HDO_final_network.pkl"
    script:
        "../pyScripts/dates/HDO/add_edge_dates_to_HDO_network.py"

rule find_HDO_nodes_with_top_5_annotations:
    input:
        network_with_all_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates_complete.pkl", 
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_5_elastic_net_coefficients_cutoff_20.csv"
    output: 
        nodes_with_top_5_annotations_df = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.csv",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.pkl" 
    script: 
        "../pyScripts/dates/HDO/find_HDO_nodes_with_top_5_annotations.py"

rule compare_past_present_and_future_HDO_networks:
    input: 
        final_network = "work_folder/data/dates/HDO/networks_with_dates/HDO_final_network.pkl",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        networks_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_networks_statistics.csv"
    script: 
        "../pyScripts/dates/HDO/compare_past_present_and_future_HDO_networks.py"

rule visualize_past_present_and_future_HDO_networks_statistics:
    input: 
        networks_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_networks_statistics.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output: 
        network_statistics_total_plots = "work_folder/data/dates/HDO/plots/depth_{depth}_cutoff_{cutoff}_total_plots.png",
        network_statistics_annotated_plots = "work_folder/data/dates/HDO/plots/depth_{depth}_cutoff_{cutoff}_annotated_plots.png"
    script:
        "../pyScripts/plotting/visualize_past_present_and_future_HDO_networks_statistics.py"

rule visualize_HDO_fractions_of_annotated_genes_during_the_years:
    input:
        networks_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_networks_statistics.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        fractions_plot = "work_folder/data/dates/HDO/plots/depth_{depth}_cutoff_{cutoff}_fractions_plot.png",
        fractions_plot_general = "work_folder/data/dates/HDO/plots/depth_{depth}_cutoff_{cutoff}_fractions_plot_general.png"
    script:
        "../pyScripts/plotting/visualize_HDO_fractions_of_annotated_genes_during_the_years.py"

rule compare_edges_evolution_in_HDO_networks_time_traveler:
    input: 
        final_network = "work_folder/data/dates/HDO/networks_with_dates/HDO_final_network.pkl",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics_time_traveler.csv"
    script:
        "../pyScripts/dates/HDO/compare_edge_evolution_in_HDO_networks_time_traveler.py"

rule visualize_edges_evolution_in_HDO_networks_time_traveler:
    input:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics_time_traveler.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        fractions_plot = "work_folder/data/dates/HDO/plots/edges_evolution/depth_{depth}_cutoff_{cutoff}_edges_evolution_fraction_plot_time_traveler.png"
    script:
        "../pyScripts/plotting/visualize_edges_evolution_in_HDO_networks_time_traveler.py"

rule compare_edges_evolution_in_HDO_networks_time_traveler_one_year_span:
    input: 
        final_network = "work_folder/data/dates/HDO/networks_with_dates/HDO_final_network.pkl",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics_time_traveler_one_year_span.csv"
    script:
        "../pyScripts/dates/HDO/compare_edge_evolution_in_HDO_networks_time_traveler_one_year_span.py"

rule visualize_edges_evolution_in_HDO_networks_time_traveler_one_year_span:
    input:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics_time_traveler_one_year_span.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        fractions_plot = "work_folder/data/dates/HDO/plots/edges_evolution/depth_{depth}_cutoff_{cutoff}_edges_evolution_fraction_plot_time_traveler_one_year_span.png"
    script:
        "../pyScripts/plotting/visualize_edges_evolution_in_HDO_networks_time_traveler_one_year_span.py"

rule compare_edges_evolution_in_HDO_networks:
    input: 
        final_network = "work_folder/data/dates/HDO/networks_with_dates/HDO_final_network.pkl",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics.csv"
    script:
        "../pyScripts/dates/HDO/compare_edge_evolution_in_HDO_networks.py"

rule visualize_edges_evolution_in_HDO_networks:
    input:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        fractions_plot = "work_folder/data/dates/HDO/plots/edges_evolution/depth_{depth}_cutoff_{cutoff}_edges_evolution_fraction_plot.png"
    script:
        "../pyScripts/plotting/visualize_edges_evolution_in_HDO_networks.py"

rule compare_edges_evolution_in_HDO_networks_one_year_span:
    input: 
        final_network = "work_folder/data/dates/HDO/networks_with_dates/HDO_final_network.pkl",
        nodes_with_top_5_annotations_pickle = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics_one_year_span.csv"
    script:
        "../pyScripts/dates/HDO/compare_edge_evolution_in_HDO_networks_one_year_span.py"

rule visualize_edges_evolution_in_HDO_networks_one_year_span:
    input:
        edges_evolution_statistics = "work_folder/data/dates/HDO/network_statistics/depth_{depth}_cutoff_{cutoff}_edges_evolution_statistics_one_year_span.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        fractions_plot = "work_folder/data/dates/HDO/plots/edges_evolution/depth_{depth}_cutoff_{cutoff}_edges_evolution_fraction_plot_one_year_span.png"
    script:
        "../pyScripts/plotting/visualize_edges_evolution_in_HDO_networks_one_year_span.py"

# --- RE-QUERYING SECTION ---

rule re_query_for_HDO_annotations_per_gene:
    input: 
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output: 
        annot_df = "work_folder/data/HDO/new_annotations_per_gene_with_ancestors.csv"
    script:
        "../pyScripts/dates/HDO/re_query_for_HDO_annotations_per_gene.R" 

# rule compute_HDO_annotations_first_dates_second_version:
#     input:
#         evidence_folder = "work_folder/data/dates/HDO/opentargets/evidence",
#         disease_folder = "work_folder/data/dates/HDO/opentargets/diseases"
#     output:
#         first_annotation_dates = "work_folder/data/dates/HDO/HDO_first_annotation_dates_second_version.csv"
#     script:
#         "../pyScripts/dates/HDO/compute_HDO_annotations_dates_second_version.py"

rule re_query_for_disgenet:
    input: 
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output: 
        annot_df = "work_folder/data/HDO/new_annotations_per_gene_with_ancestors_disgenet.csv"
    script:
        "../pyScripts/dates/HDO/re_query_for_disgenet.R" 