rule fit_elastic_net_HDO:
    input:
        annotation_df = "work_folder/data/HDO/feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_elastic_net_coefficients.csv"
    script:
        "compute_HDO_elastic_net_coefficients.py"