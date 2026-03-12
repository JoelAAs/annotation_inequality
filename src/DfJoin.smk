rule join_dfs:
    input:
        annotation_counts = "work_folder/data/{database}/annotation_per_entrez.csv",
        bait_prey_interactions = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        joined_baits = "work_folder/data/{database}/complete_baits.csv",
        joined_preys = "work_folder/data/{database}/complete_preys.csv"
    run:
        import pandas as pd

        ac_df = pd.read_csv(input.annotation_counts, sep="\t")
        bpi_df = pd.read_parquet(input.bait_prey_interactions)
        ac_df['entrez_id'] = ac_df['entrez_id'].astype(str)
        bpi_df['entrez_id'] = bpi_df['entrez_id'].astype(str)

        # First filtering for baits
        #bpi_df_baits = bpi_df[bpi_df['type'] == 'bait']
        joined_baits = bpi_df.merge(ac_df, on = 'entrez_id', how = 'outer')
        joined_baits = joined_baits.rename(columns = {"count": "count_interactions"})
        joined_baits = joined_baits[joined_baits['type'] == 'bait']
        print(f'Baits number = {len(joined_baits)}')
        # Keeping just one column for the id
        joined_baits = joined_baits.drop('entrez_id_bait', axis = 1)
        joined_baits.to_csv(output.joined_baits, sep = "\t", index = False)

        # Now filtering for preys
        #bpi_df_preys = bpi_df[bpi_df['type'] == 'prey']
        joined_preys = bpi_df.merge(ac_df, on = 'entrez_id', how = 'outer')
        joined_preys = joined_preys.rename(columns = {"count": "count_interactions"})
        joined_preys = joined_preys[joined_preys['type'] == 'prey']
        print(f'Preys number = {len(joined_preys)}')
        # Also here keeping just one of the two columns
        joined_preys = joined_preys.drop('entrez_id_bait', axis = 1)
        joined_preys.to_csv(output.joined_preys, sep = "\t", index = False)
