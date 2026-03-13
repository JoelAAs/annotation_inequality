import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path

def get_bait_or_prey(input_path):
    input_path = Path(input_path)
    filename = input_path.stem
    bait_or_prey = filename.split("_")
    return bait_or_prey[-1]

def get_database(input_path):
    path = Path(input_path)
    database = path.parent.name
    return database

def get_aspect(input_path):
    name = input_path.split("/")[-1]
    if get_bait_or_prey(input_path) == 'baits':
        return name.replace("annotation_per_entrez_", "").replace("_baits.csv", "")
    else:
        return name.replace("annotation_per_entrez_", "").replace("_preys.csv", "")

def plot_distributions(df, name, database): 
    study_counts = df['count_studies'].dropna()
    annotations = df['count_annot'].dropna()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 5))

    ax1.hist(study_counts[study_counts > 0], bins = 50, alpha = 0.7, color = 'skyblue', edgecolor = 'navy')
    ax1.set_xlabel('Study counts')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'Study counts distribution ({name})')
    ax1.set_yscale('log')

    ax2.hist(annotations[annotations > 0], bins = 50, alpha = 0.7, color = 'skyblue', edgecolor = 'navy')
    ax2.set_xlabel('Annotation counts')
    ax2.set_ylabel('Frequency')
    ax2.set_title(f'Annotation counts distribution ({name})')
    ax2.set_yscale('log')

    plt.tight_layout()
    plt.savefig(f"work_folder/data/plots/{database}_plots/{name}distributions.png")

def get_pearson_correlation(df):
    x = df['count_annot']
    y = df['count_studies']
    
    r, pval = pearsonr(x, y)
    pval = pval if float(pval) > 1e-300 else 1e-300
    print(f'r = {r}\tp-value = {pval}')

    return r, pval

def make_scatterplot():
    inputs = snakemake.input
    outputs = snakemake.output

    for inputfile, outputfile in zip(inputs, outputs):
        df = pd.read_csv(inputfile, sep = "\t")
        print(f'Opening file {inputfile}')
        print(f'Saving image on {outputfile}')

        r, pval = get_pearson_correlation(df)
        bait_or_prey = get_bait_or_prey(inputfile)
        database = get_database(inputfile)
        if database == 'GO':
            aspect = get_aspect(inputfile)
        else:
            aspect = ''
        plot_distributions(df, bait_or_prey, database)

        # Plot with log counts
        plt.figure(figsize = (10, 8))
        plt.scatter(df['count_annot'], df['count_studies'], alpha = 0.6, s = 50)
        plt.text(0.05, 0.95, f'Pearson r = {r:.5f}\np-value = {pval:.2e}',
                 transform = plt.gca().transAxes, fontsize = 12,
                 bbox = dict(boxstyle = "round,pad=0.3", facecolor = "white"))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Annotation Counts (log scale)')
        plt.ylabel('Study Counts (log scale)')
        plt.grid(True, alpha = 0.3)
        plt.title(f'{database} {aspect} Pearson Correlation Annotation vs Study counts ({bait_or_prey})\nr = {r:.3f}')

        plt.savefig(outputfile, dpi = 300, bbox_inches = 'tight')
        plt.close()

        # Plot with linear counts
        plt.figure(figsize = (10, 8))
        plt.scatter(df['count_annot'], df['count_studies'], alpha = 0.6, s = 50)
        plt.text(0.05, 0.95, f'Pearson r = {r:.5f}\np-value = {pval:.2e}',
                 transform = plt.gca().transAxes, fontsize = 12,
                 bbox = dict(boxstyle = "round,pad=0.3", facecolor = "white"))
        plt.xlabel('Annotation Counts')
        plt.ylabel('Study Counts')
        plt.grid(True, alpha = 0.3)
        plt.title(f'Pearson Correlation Annotation vs Study counts ({bait_or_prey})\nr = {r:.3f}')
        plt.savefig(f"work_folder/data/plots/{database}_plots/{aspect}pearson_linear_counts_{bait_or_prey}.png", dpi = 300, bbox_inches = 'tight')
        plt.close()

make_scatterplot()