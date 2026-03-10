rule join_dfs:
    input:
        annotation_counts = "work_folder/data/disgenet/anotation_per_uniprot.csv",
        bait_prey_interactions = "work_folder/data/intact/bait_prey_frequecies.pq"
    output:
        joined_baits = "work_folder/data/joined/interactions_annotations_baits.pq",
        joined_preys = "work_folder/data/joined/interactions_annotations_preys.pq"
    run:
        import pandas as pd

        ac_df = pd.read_csv(input.annotation_counts)
        bpi_df = pd.read_parquet(input.bait_prey_interactions)
        print("Left DF columns:", ac_df.columns.tolist())
        print("Right DF columns:", bpi_df.columns.tolist())

        # First filtering for baits
        bpi_df_baits = bpi_df[bpi_df['type'] == 'bait']
        joined_baits = bpi_df_baits.merge(ac_df, left_on = 'uniprot_id', right_on = 'uniprotid', how = 'outer')
        # Keeping just one column for the id
        joined_baits = joined_baits.drop('uniprotid', axis = 1)
        joined_baits.to_parquet(output.joined_baits)

        # Now filtering for preys
        bpi_df_preys = bpi_df[bpi_df['type'] == 'prey']
        joined_preys = bpi_df_preys.merge(ac_df, left_on = 'uniprot_id', right_on = 'uniprotid', how = 'outer')
        # Also here keeping just one of the two columns
        joined_preys = joined_preys.drop('uniprotid', axis = 1)
        joined_preys.to_parquet(output.joined_preys)
