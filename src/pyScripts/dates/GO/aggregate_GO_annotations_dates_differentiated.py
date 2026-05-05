import pandas as pd
import os

input_files = snakemake.input
out_iea = snakemake.output.dates_IEA
out_manual = snakemake.output.dates_MANUAL

print(f"[STATUS] Aggregating data from {len(input_files)} parsed snapshots...")

df_list = []
for file in input_files:
    if os.path.getsize(file) > 50: # Skip basically empty files
        df_list.append(pd.read_csv(file, sep="\t"))

def write_empty(outfile):
    # Write empty file with headers if no data is found
    with open(outfile, "w") as f:
        f.write("Gene\tGO_ID\tFirst_Date_Annotated\tEvidence\tQualifier\n")

if not df_list:
    write_empty(out_iea)
    write_empty(out_manual)
    print("[WARNING] No matching Gene/GO annotations found across any files.")
else:
    master_df = pd.concat(df_list, ignore_index=True)
    
    # Remove 'NOT' qualifiers
    initial_len = len(master_df)
    master_df = master_df[~master_df['Qualifier'].fillna('').str.contains('NOT', case=False)]
    print(f"[STATUS] Removed {initial_len - len(master_df)} annotations with 'NOT' qualifier.")
    
    # Sort chronologically, then drop duplicates keeping the oldest date
    master_df = master_df.sort_values(by="Date")
    master_df = master_df.drop_duplicates(subset=["Gene", "GO_ID"], keep="first")
    
    # Split into IEA and MANUAL DataFrames
    # Assuming your evidence code column is named 'Evidence'
    iea_df = master_df[master_df['Evidence'] == 'IEA'].copy()
    manual_df = master_df[master_df['Evidence'] != 'IEA'].copy()
    
    # Define a helper function to process and save each subset
    def process_and_save(df, outfile, label):
        if df.empty:
            write_empty(outfile)
            print(f"[WARNING] No records found for {label}.")
            return
            
        # Sort chronologically, then drop duplicates keeping the oldest date
        df = df.sort_values(by="Date")
        final_df = df.drop_duplicates(subset=["Gene", "GO_ID"], keep="first")
        
        # Drop missing dates, rename, cast to int, and save
        final_df = final_df.dropna(subset=["Date"])
        final_df = final_df.rename(columns={"Date": "First_Date_Annotated"})
        final_df["First_Date_Annotated"] = final_df["First_Date_Annotated"].astype(int)
        
        final_df.to_csv(outfile, sep="\t", index=False)
        print(f"[STATUS] Saved historical inception dates for {len(final_df)} {label} Gene-GO pairs.")

    # 4. Execute for both files
    print("\n--- Processing IEA ---")
    process_and_save(iea_df, out_iea, "IEA")
    
    print("\n--- Processing MANUAL ---")
    process_and_save(manual_df, out_manual, "MANUAL")

print("[STATUS] Processing complete!")