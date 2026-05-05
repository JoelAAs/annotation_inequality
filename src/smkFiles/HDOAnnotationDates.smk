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

## TODO
rule assign_dates_to_HDO_missing_ones:
    input:
        network_with_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates.pkl",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        network_with_all_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates_complete.pkl"
    script:
        "../pyScripts/dates/HDO/assign_dates_to_HDO_missing_ones.py"

## TODO
rule add_edge_dates_to_HDO_network:
    input:
        network_with_all_dates = "work_folder/data/dates/HDO/networks_with_dates/HDO_network_with_dates_complete.pkl", 
        bp_publications = "work_folder/data/intact/bait_prey_publications_complete_dates.pq"
    output:
        final_network = "work_folder/data/dates/HDOO/networks_with_dates/HDO_final_network.pkl"
    script:
        "../pyScripts/dates/HDO/add_edge_dates_to_HDO_network.py"