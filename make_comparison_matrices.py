import os
import shutil
import argparse
import pandas as pd
import numpy as np
import re
from Bio import Phylo
from io import StringIO

def copy_gff_files(main_folder, gff_output_folder):
    """
    Copy all .gff files from nested genome subdirectories into an output directory.
    This ensures that all .gff annotation files are available in a single location for parsing.
    """
    os.makedirs(gff_output_folder, exist_ok=True)

    for root, _, files in os.walk(main_folder):
        for file in files:
            if file.endswith(".gff"):
                src = os.path.join(root, file)
                dst = os.path.join(gff_output_folder, file)
                shutil.copy2(src, dst)

    print(f"All .gff files copied to {gff_output_folder}")

def parse_gff(file_path, target_genes):
    """
    Parse GFF file and detect presence of specific target genes.

    Args:
        file_path (str): Path to the .gff annotation file.
        target_genes (list): List of genes to check for.

    Returns:
        tuple: Genome name and a dictionary of gene presence/absence (0 or 1).
    """
    genome_name = os.path.basename(file_path).replace(".gff", "")
    gene_presence = {gene: 0 for gene in target_genes}

    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                if len(columns) > 8:
                    attributes = columns[8]
                    match = re.search(r"Name=([^;]+)", attributes)
                    if match:
                        gene_name = match.group(1)
                        if gene_name in target_genes:
                            gene_presence[gene_name] = 1

    return genome_name, gene_presence

def reorder_matrix_by_tree(df, tree_file_path):
    """
    Reorder a square matrix based on the tip order of a phylogenetic tree.

    Args:
        df (DataFrame): Square matrix of presence/absence.
        tree_file_path (str): Path to a Newick tree file.

    Returns:
        DataFrame: Matrix reordered to match the tree tip labels.
    """
    with open(tree_file_path, "r") as f:
        newick_text = f.read()
    tree = Phylo.read(StringIO(newick_text), "newick")
    tip_labels_ordered = [leaf.name.strip("'") for leaf in tree.get_terminals()]

    # Keep only genomes that exist in both the tree and the matrix
    valid_genomes = [g for g in tip_labels_ordered if g in df.index]
    df_ordered = df.loc[valid_genomes, valid_genomes]
    return df_ordered

def process_gff_folder(folder_path, target_genes, output_folder, tree_file, label):
    """
    Generate a gene presence/absence matrix and ordered gene-gene comparison matrices.

    Args:
        folder_path (str): Folder with GFF files.
        target_genes (list): List of genes to evaluate.
        output_folder (str): Where to save output CSVs.
        tree_file (str): Path to Newick tree for ordering.
        label (str): Identifier prefix (e.g., "HGT" or "VGT").
    """
    os.makedirs(output_folder, exist_ok=True)
    results = []

    # Collect and sort GFF files by creation time
    file_list = sorted([
        f for f in os.listdir(folder_path) if f.endswith(".gff")
    ], key=lambda x: os.path.getctime(os.path.join(folder_path, x)))

    # Parse each file and collect presence data
    for filename in file_list:
        file_path = os.path.join(folder_path, filename)
        genome_name, gene_presence = parse_gff(file_path, target_genes)
        results.append({"Genome": genome_name, **gene_presence})

    df = pd.DataFrame(results)
    df.set_index("Genome", inplace=True)

    # Save the presence/absence matrix
    matrix_path = os.path.join(output_folder, f"{label}_gene_presence_matrix.csv")
    df.to_csv(matrix_path)
    print(f"[{label}] Presence/absence matrix saved: {matrix_path}")

    # Create pairwise comparison matrices for each gene
    for gene in target_genes:
        # 1 if gene is present in both genomes, 0 otherwise
        gene_matrix = (df[gene].values[:, None] & df[gene].values[None, :]).astype(int)
        gene_df = pd.DataFrame(gene_matrix, index=df.index, columns=df.index)

        # Reorder the matrix based on the tree
        gene_df_ordered = reorder_matrix_by_tree(gene_df, tree_file)

        # Save reordered matrix
        output_filename = f"{gene}_comparison_matrix_ordered_by_tree.csv"
        output_path = os.path.join(output_folder, output_filename)
        gene_df_ordered.to_csv(output_path)
        print(f"[{label}] {gene} ordered comparison matrix saved: {output_path}")

def main(input_folder, output_folder, tree_file):
    """
    Main function: Copies files, processes HGT and VGT genes, outputs matrices.
    """
    gff_flat_folder = os.path.join(output_folder, "gff")
    copy_gff_files(input_folder, gff_flat_folder)

    # Define and process horizontally transferred genes
    hgt_genes = ['mecA', 'cmlA', 'tetA', 'tetM', 'vanA']
    hgt_output = os.path.join(output_folder, "hgt")
    process_gff_folder(gff_flat_folder, hgt_genes, hgt_output, tree_file, label="HGT")

    # Define and process vertically transferred genes
    vgt_genes = ['gapA', 'gyrA', 'gyrB', 'recA', 'rpoB']
    vgt_output = os.path.join(output_folder, "vgt")
    process_gff_folder(gff_flat_folder, vgt_genes, vgt_output, tree_file, label="VGT")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full pipeline: extract, process, and reorder gene matrices by tree.")
    parser.add_argument("--input_folder", required=True, help="Folder containing annotated subfolders with .gff/.ffn files")
    parser.add_argument("--output_folder", required=True, help="Destination folder for output")
    parser.add_argument("--tree_file", required=True, help="Path to the Newick tree file (.txt)")

    args = parser.parse_args()
    main(args.input_folder, args.output_folder, args.tree_file)
