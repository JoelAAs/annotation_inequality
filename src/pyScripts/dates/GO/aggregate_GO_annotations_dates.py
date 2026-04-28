import pandas as pd
import os

input_files = snakemake.input
output_file = snakemake.output[0]

print(f"[STATUS] Aggregating data from {len(input_files)} parsed snapshots...")

df_list = []
for file in input_files:
    if os.path.getsize(file) > 50: # Skip basically empty files
        df_list.append(pd.read_csv(file, sep="\t"))

if not df_list:
    with open(output_file, "w") as f:
        f.write("Gene\tGO_ID\tFirst_Date_Annotated\n")
    print("[WARNING] No matching Gene/GO annotations found.")
else:
    master_df = pd.concat(df_list, ignore_index=True)
    
    # Sort chronologically, then drop duplicates keeping the oldest date
    master_df = master_df.sort_values(by="Date")
    final_df = master_df.drop_duplicates(subset=["Gene", "GO_ID"], keep="first")
    
    print(f"[STATUS] Found inception dates for {len(final_df)} before removing NaN dates.")

    final_df = final_df.dropna(subset=["Date"])
    final_df = final_df.rename(columns={"Date": "First_Date_Annotated"})
    final_df["First_Date_Annotated"] = final_df["First_Date_Annotated"].astype(int)
    final_df.to_csv(output_file, sep="\t", index=False)
    
    print(f"[STATUS] Found historical inception dates for {len(final_df)} Gene-GO pairs.")