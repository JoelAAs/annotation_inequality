import pandas as pd 

def get_go_dfs(go_dfs_list):
    go_dfs = {
        f.split('_')[-2]: pd.read_csv(f, sep = "\t")
        for f in go_dfs_list
    }

    return go_dfs

go_dfs_list = snakemake.input.go_dfs
hdo_df_input = snakemake.input.hdo_df
outputfile = snakemake.output.merged_df

go_dfs = get_go_dfs(go_dfs_list)
bp_df = go_dfs.get('BP')
mf_df = go_dfs.get('MF')
cc_df = go_dfs.get('CC')
hdo_df = pd.read_csv(hdo_df_input, sep = '\t')

dfs = [bp_df, mf_df, cc_df]
for df in dfs:
    df = df.drop(columns = ['entrez_id_bait', 'count_studies'], inplace = True)

bp_df = bp_df.rename(columns = {'count_annot': 'count_annot_bp'}, inplace = True)
mf_df = mf_df.rename(columns = {'count_annot': 'count_annot_mf'}, inplace = True)
cc_df = cc_df.rename(columns = {'count_annot': 'count_annot_cc'}, inplace = True)
hdo_df = hdo_df.rename(columns = {'count_annot': 'count_annot_hdo'})

merged_df = hdo_df
for df in dfs:
    merged_df = pd.merge(merged_df, df, on = 'entrez_id', how = 'left')

merged_df = merged_df.drop(columns = ['entrez_id_bait'])

merged_df.to_csv(outputfile, sep = '\t')
