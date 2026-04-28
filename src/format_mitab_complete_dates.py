import pandas as pd
import re

def _find_pattern(cell, pattern, single=True):
    # Added str(cell) just to be safe against nulls/floats
    match = re.search(pattern, str(cell))
    if match:
        matches = match.groups()
        if single:
            return matches[0]
        return match.groups()
    return None

def filter_mitab(filename):
    """
    Remove non-human proteins.
    Remove any non bait-prey interaction.
    Format publications ids and uniprot identifiers
    :param filename: miTab interaction file
    :return: filtered and formated pandas dataframe
    """
    mitab_df = pd.read_csv(filename, sep="\t")
    mitab_df = mitab_df[
        (
                mitab_df['Taxid interactor A'].str.contains("9606") &
                mitab_df['Taxid interactor B'].str.contains("9606")
        )
    ]
    mitab_df['RoleA'] = mitab_df['Experimental role(s) interactor A'].str.contains("bait")
    mitab_df['RoleB'] = mitab_df['Experimental role(s) interactor B'].str.contains("bait")

    # No bait-bait prey-prey pairs
    mitab_df = mitab_df[mitab_df['RoleA'] ^ mitab_df['RoleB']]

    mitab_df["pubmed_id"] = mitab_df["Publication Identifier(s)"].apply(
        _find_pattern, args=(r"pubmed:(\d+)",))
    mitab_df["IDA"] = mitab_df["#ID(s) interactor A"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    mitab_df["IDB"] = mitab_df["ID(s) interactor B"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    mitab_df["detection_method"] = mitab_df["Interaction detection method(s)"].apply(
        _find_pattern, args=(r"psi-mi:\"(MI:\d+)\"",))
    
    # NEW CODE: Capture full YYYY/MM/DD and remove the slashes to get YYYYMMDD
    mitab_df["date_created"] = mitab_df["Creation date"].apply(
        _find_pattern, args=(r"(\d{4}/\d{2}/\d{2})",)
    ).str.replace("/", "")
    
    return mitab_df


def reform_to_bait_prey(mitab_df):
    """
    Reform PPI pandas dataframe into bait-prey format
    :param mitab_df: input filtered
    :return: data frame with columns:
        uniprot_bait
        uniprot_prey
        pubmed_id
        detection_method
        date_created
    """
    # Updated column list to include 'date_created'
    mitab_df = mitab_df[["IDA", "IDB", "RoleA", "RoleB", "detection_method", "pubmed_id", "date_created"]].drop_duplicates(keep="first")
    n_ppis = mitab_df.shape[0]
    mitab_df = mitab_df.dropna()
    print(f"dropped {n_ppis - mitab_df.shape[0]} out of {n_ppis} due to missing values")
    mitab_df = mitab_df[mitab_df["RoleA"] ^ mitab_df["RoleB"]]
    
    reform_list = []
    for _, row in mitab_df.iterrows():
        if row["RoleA"]:
            bait = row["IDA"]
            prey = row["IDB"]
        else:
            bait = row["IDB"]
            prey = row["IDA"]

        reform_list.append({
            "uniprot_id_bait": bait,
            "uniprot_id_prey": prey,
            "pubmed_id": row["pubmed_id"],
            "detection_method": row["detection_method"].replace(":", "-"),
            "date_created": row["date_created"] # Storing the full YYYYMMDD date
        })

    bait_prey_df = pd.DataFrame(reform_list)
    bait_prey_df = bait_prey_df.drop_duplicates()
    
    return bait_prey_df