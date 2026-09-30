import json
import pandas as pd

def load_ctd(ctd_path):
    """Load CTD_genes_diseases.tsv and keep only rows with direct
    evidence and at least one PubMed reference. Rows without
    DirectEvidence are chemical-inferred associations, not direct
    gene-disease evidence, and are dropped."""
    ctd = pd.read_csv(
        ctd_path, sep="\t", comment="#",
        names=[
            "GeneSymbol", "GeneID", "DiseaseName", "DiseaseID",
            "DirectEvidence", "InferenceChemicalName", "InferenceScore",
            "OmimIDs", "PubMedIDs"
        ],
        dtype=str,
        low_memory=False,
    )
    ctd = ctd[ctd["DirectEvidence"].notna()]
    ctd = ctd.dropna(subset=["PubMedIDs"])
    ctd["GeneID"] = ctd["GeneID"].astype(str)
    return ctd


def merge_pmid_strings(series):
    """Union of all pipe-separated PMID lists in a pandas Series
    into a single sorted list of unique PMIDs."""
    all_pmids = set()
    for pmid_str in series:
        all_pmids.update(pmid_str.split("|"))
    return sorted(all_pmids)


def main():
    # --- Load inputs (paths injected by Snakemake) ---
    pairs_path = snakemake.input.pairs   # noqa: F821
    xrefs_path = snakemake.input.xrefs   # noqa: F821
    ctd_path = snakemake.input.ctd       # noqa: F821
    output_path = snakemake.output.pair_pmids  # noqa: F821

    pairs_df = pd.read_csv(pairs_path, dtype=str)

    with open(xrefs_path) as f:
        doid_to_mesh_omim = json.load(f)

    ctd = load_ctd(ctd_path)

    # --- Expand each DOID into its MeSH/OMIM xrefs, one row per xref ---
    pairs_df["mesh_omim_ids"] = pairs_df["doid"].map(doid_to_mesh_omim)

    n_before = len(pairs_df)
    exploded = pairs_df.explode("mesh_omim_ids").dropna(subset=["mesh_omim_ids"])
    n_no_xref = n_before - exploded["doid"].nunique()
    print(f"DOIDs with no MeSH/OMIM xref (cannot be matched in CTD): {n_no_xref}")

    # CTD's xref IDs may need normalization: our obo xrefs look like
    # "MESH:D001943" or "OMIM:114480". CTD's DiseaseID column follows
    # the same "MESH:xxx" / "OMIM:xxx" convention -- verify this matches
    # your actual obo xref format before relying on it as-is.
    merged = exploded.merge(
        ctd,
        left_on=["entrez_id", "mesh_omim_ids"],
        right_on=["GeneID", "DiseaseID"],
        how="inner",
    )

    print(f"Pairs with at least one CTD match: {merged[['entrez_id', 'doid']].drop_duplicates().shape[0]} "
          f"out of {n_before} total leaf pairs")

    # --- Aggregate: union of PMIDs per (entrez_id, doid) pair ---
    pair_to_pmids = (
        merged.groupby(["entrez_id", "doid"])["PubMedIDs"]
        .apply(merge_pmid_strings)
        .reset_index()
        .rename(columns={"PubMedIDs": "pmids"})
    )

    # Store pmids as a pipe-separated string for a clean flat CSV;
    # downstream steps will split it back into a list.
    pair_to_pmids["pmids"] = pair_to_pmids["pmids"].apply(lambda l: "|".join(l))

    pair_to_pmids.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()