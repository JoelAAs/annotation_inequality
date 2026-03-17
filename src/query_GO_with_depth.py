import pandas as pd 
import mygene
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
import os

# Args
input_pod = snakemake.input.bp_frequencies
ontologies = {
    'BP': snakemake.output.annotation_df_BP,
    'MF': snakemake.output.annotation_df_MF,
    'CC': snakemake.output.annotation_df_CC
}
depth = snakemake.wildcards.depth

# Input
pod_df = pd.read_parquet(input_pod)

# Processing
mg = mygene.MyGeneInfo()

cache_dir = os.path.expanduser('~/.cache/goatools')
os.makedirs(cache_dir, exist_ok = True)
obo_path = os.path.join(cache_dir, 'go-basic.obo')

if not os.path.exists(obo_path):
    obo_path = download_go_basic_obo(obo = obo_path)

go_dag = GODag(obo_path)

for ontology, outputfile in ontologies.items():
    # Query mygene
    results = mg.querymany(
        pod_df['entrez_id'].astype(str),
        scopes = 'entrezgene',
        fields = f'entrezgene,go.{ontology}.id',
        species = 'human',
        as_dataframe = True
    )

    # Save results, drop na values, and clean df
    go_df = results
    go_df['go_id'] = go_df[f'go.{ontology}'].apply(
        lambda x: [d['id'] for d in x] if isinstance(x, list) else []
    )
    go_df = go_df[['entrezgene', 'go_id']].explode('go_id')
    go_df = go_df.dropna()
    go_df = go_df.rename(columns = {
        'entrezgene': 'entrez_id',
    })

    go_df = go_df.drop_duplicates(subset = ['entrez_id', 'go_id'])

    data = []

    for _, go_id in go_df.iterrows():
        term = go_dag.get(go_id['go_id'])
        if term:
            data.append({
                'entrez_id': go_id['entrez_id'],
                "go_id": go_id['go_id'],
                'namespace': term.namespace,
                "depth": term.depth
            })
        else:
            data.append({
                'entrez_id': go_id['entrez_id'],
                "go_id": go_id['go_id'],
                'namespace': term.namespace,
                "depth": None
            })

    df = pd.DataFrame(data)
    df = df[df['depth'] <= 4]

    pod_df.drop_duplicates(subset = ['entrez_id'])
    go_annot_count = df.groupby('entrez_id').size().reset_index(name = 'count')

    final_df = (
        pod_df[['entrez_id']]
        .astype(str)
        .merge(go_annot_count, on = 'entrez_id', how = 'outer')
        .fillna({'count':0}).astype(int)
    )

    print(f'GO {ontology} annotation counts with depth {depth} done!')

    final_df.to_csv(outputfile, sep = '\t', index = False)