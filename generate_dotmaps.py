import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_dotmap(matrix_folder):
    """
    Generates a combined dotmap image by overlaying all
    ordered gene comparison matrices found in the specified folder.

    Each gene is represented by a different color on the dotmap.

    Args:
        matrix_folder (str): Path to folder containing *_comparison_matrix_ordered_by_tree.csv files.
    """
    # Collect all relevant matrix files
    matrix_files = [
        os.path.join(matrix_folder, f)
        for f in os.listdir(matrix_folder)
        if f.endswith("_comparison_matrix_ordered_by_tree.csv")
    ]

    if not matrix_files:
        print(f"No ordered matrix files found in {matrix_folder}")
        return

    # Sort files for consistent plotting order
    matrix_files.sort()
    num_files = len(matrix_files)

    # Generate a colormap with a unique color for each gene
    colormap = plt.colormaps.get_cmap('nipy_spectral')
    colors = [colormap(i / num_files) for i in range(num_files)]

    plt.figure(figsize=(25, 25))  # Large figure for high-density matrices

    for i, path in enumerate(matrix_files):
        gene_name = os.path.basename(path).split("_comparison_matrix")[0]
        color = colors[i]

        try:
            df = pd.read_csv(path, index_col=0)
        except Exception as e:
            print(f"Error reading {path}: {e}")
            continue

        # Ensure binary matrix and get coordinates of 1s (gene presence in both strains)
        binary_df = (df > 0).astype(int)
        y_coords, x_coords = binary_df.values.nonzero()

        # Plot as square scatter points
        plt.scatter(x_coords, y_coords, color=color, s=6, marker='s', label=gene_name, alpha=0.9)

    # Label axes with strain names
    plt.xticks(ticks=range(len(binary_df.columns)), labels=binary_df.columns, rotation=90, fontsize=3)
    plt.yticks(ticks=range(len(binary_df.index)), labels=binary_df.index, fontsize=3)

    plt.title("Binární mapa genomových podobností")
    plt.legend(markerscale=2, fontsize=8, loc='upper right')
    plt.tight_layout()

    # Save as a PNG image
    output_path = os.path.join(matrix_folder, "combined_dotmap.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Dotmap saved to: {output_path}")

def main(hgt_folder, vgt_folder):
    """
    Main function: Generate dotmaps for HGT and VGT folders.
    """
    if hgt_folder:
        print(f"\nProcessing HGT folder: {hgt_folder}")
        create_dotmap(hgt_folder)

    if vgt_folder:
        print(f"\nProcessing VGT folder: {vgt_folder}")
        create_dotmap(vgt_folder)

if __name__ == "__main__":
    # Argument parsing for command line execution
    parser = argparse.ArgumentParser(description="Generate combined dotmaps for HGT and VGT comparison matrices.")
    parser.add_argument("--hgt_folder", required=True, help="Folder with ordered HGT matrices")
    parser.add_argument("--vgt_folder", required=True, help="Folder with ordered VGT matrices")

    args = parser.parse_args()
    main(args.hgt_folder, args.vgt_folder)
