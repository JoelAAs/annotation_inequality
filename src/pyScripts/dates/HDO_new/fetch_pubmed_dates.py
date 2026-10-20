import pandas as pd
import time
from Bio import Entrez

file_input = snakemake.input.associations
file_output = snakemake.output.pmid_dates

# NCBI requires an email address for API usage
Entrez.email = "Gabriele.Oberti@ieo.it"

# Extract unique PMIDs
print("--- [INFO] Reading Parquet dataset... ---")
df = pd.read_parquet(file_input)

# Filter only rows with a PMID and remove the prefix to get the raw numbers
pmid_series = df.loc[df['Reference'].str.startswith('PMID:', na=False), 'Reference']
unique_pmids = pmid_series.str.replace('PMID:', '').unique().tolist()

print(f"--- [INFO] Found {len(unique_pmids)} UNIQUE PMIDs to fetch from NCBI. ---")

# Batching function
def chunker(seq, size):
    """Splits a list into smaller chunks of a specified size."""
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

# Query PubMed
date_results = []
batch_size = 200  # NCBI recommends batches of 200-500 to avoid timeouts
total_batches = (len(unique_pmids) // batch_size) + 1

print(f"--- [INFO] Starting PubMed data retrieval in {total_batches} batches... ---")

for i, batch in enumerate(chunker(unique_pmids, batch_size), 1):
    print(f"    > Processing batch {i}/{total_batches}...")
    
    try:
        # efetch retrieves full bibliographic data in XML format
        handle = Entrez.efetch(db="pubmed", id=",".join(batch), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Parse each returned article
        for article in records['PubmedArticle']:
            pmid_str = str(article['MedlineCitation']['PMID'])
            
            # Look for the exact date the paper entered the PubMed database
            pub_date_node = None
            for date_node in article['PubmedData']['History']:
                if date_node.attributes.get('PubStatus') == 'pubmed':
                    pub_date_node = date_node
                    break
            
            if pub_date_node:
                year = pub_date_node.get('Year', '')
                month = pub_date_node.get('Month', '01') 
                day = pub_date_node.get('Day', '01')     
                
                # Format as YYYY-MM-DD
                exact_date = f"{year}-{month.zfill(2)}-{day.zfill(2)}"
            else:
                exact_date = "Date_not_found"
            
            date_results.append({
                'Reference': f"PMID:{pmid_str}", 
                'Publication_Date': exact_date
            })
            
    except Exception as e:
        print(f"    [ERROR] Batch {i} failed: {e}")
    
    # Polite delay to prevent NCBI server bans (approx 2-3 requests/sec)
    time.sleep(0.4) 

# Save results
df_dates = pd.DataFrame(date_results)
print(f"\n--- [INFO] Download complete! Recovered dates for {len(df_dates)} PMIDs. ---")

# Ensure no duplicates exist
df_dates = df_dates.drop_duplicates(subset=['Reference'])

# Save translation dictionary to Parquet
df_dates.to_parquet(file_output, index=False)
print(f"--- [INFO] PubMed dates saved to: {file_output} ---")