rule download_ctd:
    output:
        "work_folder/data/dates/HDO_new/CTD_genes_diseases.tsv.gz"
    shell:
        "wget -O {output} https://ctdbase.org/reports/CTD_genes_diseases.tsv.gz"

rule decompress_ctd:
    input:
        "work_folder/data/dates/HDO_new/CTD_genes_diseases.tsv.gz"
    output:
        "work_folder/data/dates/HDO_new/CTD_genes_diseases.tsv"
    shell:
        "gunzip -k {input}"

# Parse the DO ontology once -> ancestors map + DOID-to-MeSH/OMIM xref map
rule parse_do_ontology:
    input:
        obo = "work_folder/data/HDO/doid.obo"
    output:
        ancestors = "work_folder/data/dates/HDO_new/results/ancestors_map.pkl",
        xrefs = "work_folder/data/dates/HDO_new/results/doid_to_mesh_omim.json"
    script:
        "../pyScripts/dates/HDO/parse_do_ontology.py"
    
# From the network, keep only the leaf DOIDs per node (drop propagated ancestors)
# Also flatten into a (entrez_id, doid) pairs table
rule extract_leaf_pairs:
    input:
        network = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        ancestors = "work_folder/data/dates/HDO_new/results/ancestors_map.pkl"
    output:
        pairs = "work_folder/data/dates/HDO_new/results/leaf_pairs.csv"
    script:
        "../pyScripts/dates/HDO/extract_leaf_pairs.py"

# Join (entrez_id, doid) pairs against CTD to get supporting PMIDs
rule query_ctd:
    input:
        pairs = "work_folder/data/dates/HDO_new/results/leaf_pairs.csv",
        xrefs = "work_folder/data/dates/HDO_new/results/doid_to_mesh_omim.json",
        ctd = "work_folder/data/dates/HDO_new/CTD_genes_diseases.tsv"
    output:
        pair_pmids = "work_folder/data/dates/HDO_new/pair_to_pmids.csv"
    script:
        "../pyScripts/dates/HDO/query_ctd.py"

# Fetch the most precise available publication date for every unique PMID colected in the previous rule
rule fetch_pubmed_dates:
    input:
        pair_pmids = "work_folder/data/dates/HDO_new/pair_to_pmids.csv"
    params:
        email = "Gabriele.Oberti@ieo.it",
        batch_size = 200
    output:
        pmid_dates = "work_folder/data/dates/HDO_new/pmid_dates.csv"
    script:
        "../pyScripts/dates/HDO/fetch_pubmed_dates.py"

# For each (entrez_id, doid) pair, pick the earliest dated PMID
# Final table --> entrez_id, doid, date, granularity, pmid
rule assemble_final_dates_table:
    input:
        pair_pmids = "work_folder/data/dates/HDO_new/pair_to_pmids.csv",
        pmid_dates = "work_folder/data/dates/HDO_new/pmid_dates.csv"
    output:
        final = "work_folder/data/dates/HDO_new/entrez_doid_dates.csv"
    script:
        "../pyScripts/dates/HDO/assemble_final_dates_table.py"