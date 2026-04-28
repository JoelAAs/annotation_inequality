import pandas as pd
import pickle

## This function is not being used, it is here in case it will be required in the future!

# Load your new comprehensive network
path = "work_folder/data/dates/GO/networks_with_dates/BP_network_with_dates.pkl" 
with open(path, "rb") as f:
    G = pickle.load(f)

# The magic list comprehension
flat_data = [
    {
        "Entrez_ID": node, 
        "Bait_Count": attr.get("bait_count", 0), # Assuming this is still an attribute!
        "GO_ID": item["go_id"], 
        "Depth": item["depth"], 
        "Date": item["first_annotation_date"]
    }
    for node, attr in G.nodes(data=True) 
    for item in attr.get("go_annotations", [])
]

df = pd.DataFrame(flat_data)

# Show off your clean data
print(df.head())