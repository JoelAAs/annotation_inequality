import urllib.request
import re

# --- Configuration ---
SPECIES = "HUMAN" 

def get_all_GO_snapshots_downloads(wildcards):
    # This print will appear ONLY after the checkpoint finishes
    print("\n[STATUS] Checkpoint triggered! Reading the newly generated file list...")
    
    with open(checkpoints.get_GO_ftp_list.get().output[0]) as f:
        files = [line.strip() for line in f if line.strip()]
    
    # Print how many files it found before it starts downloading
    print(f"[STATUS] Found {len(files)} historical snapshots. Queuing download jobs...\n")
    
    return expand("work_folder/data/dates/GO/snapshots/{filename}", filename=files)

# Checkpoint: Scrape the HTTPS server instead of FTP
checkpoint get_GO_ftp_list:
    output:
        "work_folder/data/dates/GO/snapshot_list.txt"
    message: 
        "--> [JOB] Scraping EBI HTTPS server for available snapshots..."
    run:
        url = f"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/{SPECIES}/"
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            html = response.read().decode('utf-8')
        
        # 1. Grab EVERY link on the page
        all_links = re.findall(r'href="([^"]+)"', html)
        
        valid_files = set()
        for link in all_links:
            filename = link.split("/")[-1]
            
            # 2. Check for modern files AND legacy files that are compressed
            if ("gaf" in filename or "gene_association" in filename) and filename.endswith(".gz"):
                
                # 3. NEW: Explicitly exclude isoforms, RNAs, and complexes
                if not any(x in filename for x in ["isoform", "rna", "complex", "plus", "ref"]):
                    valid_files.add(filename)
        
        # 4. Write them to the text file
        with open(output[0], "w") as f:
            for file_name in sorted(valid_files):
                f.write(file_name + "\n")

# Rule: Download using the HTTPS link
rule download_GO_snapshot:
    output:
        "work_folder/data/dates/GO/snapshots/{filename}"
    params:
        # FIXED: Single curly braces here so the wildcard evaluates correctly!
        url=lambda wildcards: f"https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/{SPECIES}/{wildcards.filename}"
    message: 
        "--> [JOB] Downloading snapshot: {wildcards.filename}"
    shell:
        # Removed -q so wget will print helpful error logs if it ever fails
        "wget -O {output} {params.url}"

rule mark_GO_downloads_complete:
    input:
        get_all_GO_snapshots_downloads
    message: 
        "--> [JOB] All downloads complete! Touching final flag file."
    output:
        touch("work_folder/data/dates/GO/downloads_complete.txt")