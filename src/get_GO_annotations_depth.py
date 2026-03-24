import pandas as pd
import numpy as np
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
import os

input_pod = snakemake.input.annot_df
output_pod = snakemake.output.annot_df_depth

cache_dir = os.path.expanduser('~/.cache/goatools')
os.makedirs(cache_dir, exist_ok = True)
obo_path = os.path.join(cache_dir, 'go-basic.obo')

if not os.path.exists(obo_path):
    obo_path = download_go_basic_obo(obo = obo_path)

go_dag = GODag(obo_path)

def get_depth(go_id):
    if pd.isna(go_id):
        return np.nan
    
    term =  go_dag.get(go_id)
    
    if term is None:
        return np.nan
    
    return term.depth

for inputfile, outputfile in zip(input_pod, output_pod):
    annot_df = pd.read_csv(inputfile, sep = '\t')
    aspect = inputfile.split('/')[-1].split('_')[0]

    print(f'Adding depth to GO {aspect} terms...')

    annot_df['depth'] = annot_df['go_id'].apply(get_depth)

    annot_df.to_csv(outputfile, sep = '\t', index = False)

    print(f'Depth to GO {aspect} terms added!')