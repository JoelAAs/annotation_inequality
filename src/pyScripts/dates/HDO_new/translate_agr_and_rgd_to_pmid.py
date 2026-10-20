import pandas as pd
import requests
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

file_input = snakemake.input.associations
file_output = snakemake.output.translated_refs

# Extract unique references
print("--- [INFO] Reading Parquet dataset... ---")
df = pd.read_parquet(file_input)

# Isolate everything that is NOT a PMID and NOT 'no_reference'
mask = (df['Reference'] != 'no_reference') & (~df['Reference'].str.startswith('PMID:', na=False))
refs_to_translate = df.loc[mask, 'Reference'].unique().tolist()

print(f"--- [INFO] Found {len(refs_to_translate)} UNIQUE references (AGRKB, RGD, etc.) to translate. ---")

# API connection setup
# Create a robust session that retries automatically on 5xx errors or rate limits
session = requests.Session()
retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
session.mount('https://', HTTPAdapter(max_retries=retries))

base_url = "https://www.alliancegenome.org/api/reference/"
results = []

# Querying
print("--- [INFO] Starting Alliance of Genome Resources API querying... ---")

for i, ref_curie in enumerate(refs_to_translate, 1):
    if i % 100 == 0:
        print(f"    > Processed {i}/{len(refs_to_translate)} references...")
        
    found_pmid = None
    fallback_date = None
    
    try:
        # API call (10-second timeout to avoid getting stuck)
        response = session.get(f"{base_url}{ref_curie}", timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # AGR returns JSON data. Look for cross_references.
            # Handle both 'crossReferences' and 'cross_references' keys for safety
            xrefs = data.get('crossReferences', []) or data.get('cross_references', [])
            
            for xref in xrefs:
                curie = xref.get('curie', '')
                if curie.startswith('PMID:'):
                    found_pmid = curie
                    break
            
            # If no PMID is found, save the AGR fallback date
            if not found_pmid:
                fallback_date = data.get('datePublished')
                
        else:
            # 404 means the reference cannot be resolved via this API
            pass
            
    except Exception as e:
        print(f"    [ERROR] Query failed for {ref_curie}: {e}")
        
    # Save the result for this iteration
    results.append({
        'Original_Reference': ref_curie,
        'Translated_PMID': found_pmid,
        'Fallback_Date': fallback_date
    })
    
    # 0.1-second pause: courtesy to the server (approx. 10 requests per second)
    time.sleep(0.1)

# Final stats and Save
df_results = pd.DataFrame(results)
recovered_pmids = df_results['Translated_PMID'].notna().sum()
recovered_dates = df_results['Fallback_Date'].notna().sum()
unresolved = len(df_results) - recovered_pmids - recovered_dates

print("\n--- [INFO] Translation completed! ---")
print(f"      > Recovered PMIDs: {recovered_pmids}")
print(f"      > Recovered fallback dates (no PMID): {recovered_dates}")
print(f"      > Unresolved references: {unresolved}")

# Ensure columns are string format before saving to Parquet
df_results = df_results.astype(str)
df_results.to_parquet(file_output, index=False)
print(f"--- [INFO] Translation dictionary saved to: {file_output} ---")