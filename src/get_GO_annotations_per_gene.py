import pandas as pd
import mygene

ASPECTS = ['BP', 'MF', 'CC']
input_df = snakemake.input.bp_frequencies

mg = mygene.MyGeneInfo()
df = pd.read_parquet(input_df)
df = df.drop_duplicates(subset = 'entrez_id')

for aspect in ASPECTS:
    output_df = [p for p in snakemake.output.annotation_df if f"{aspect}_" in p][0]
    output_list = [p for p in snakemake.output.annotation_list if f"{aspect}_" in p][0]

    results = mg.querymany(
        df['entrez_id'].astype(str),
        scopes = 'entrezgene',
        fields = f'entrezgene,go.{aspect}',
        species = 'human',
        as_dataframe = True
    )

    go_df = results
    rows = []
    
    for entrez_id, data in results.iterrows():
    
        go_entries = data.get(f'go.{aspect}')
        
        if isinstance(go_entries, list):
            for entry in go_entries:
                rows.append({
                    'entrez_id': entrez_id,
                    'go_id': entry.get('id'),
                    'go_term_name': entry.get('term')
                })
        elif isinstance(go_entries, dict):
            rows.append({
                'entrez_id': entrez_id,
                'go_id': go_entries.get('id'),
                'go_term_name': go_entries.get('term')
            })

        
    print(f'GO {aspect} annotations obtained!')

    final_df = pd.DataFrame(rows).drop_duplicates(subset=['entrez_id', 'go_id', 'go_term_name'])    
    annot_list = final_df['go_term_name']
    annot_list = annot_list.drop_duplicates()

    final_df.to_csv(output_df, sep='\t', index=False)
    annot_list.to_csv(output_list, sep='\t', index=False)