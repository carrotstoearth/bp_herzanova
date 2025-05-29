# === Imports ===
import os
import sys
import argparse
import pandas as pd
import numpy as np
import tempfile
from PyQt5.QtWidgets import (
    QApplication, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QListWidgetItem, QAbstractItemView, QPushButton, QLabel
)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl, QTimer
import plotly.graph_objects as go
import plotly.express as px

class DotmapViewer(QWidget):
    """
    Main GUI class for the Gene Dotmap Viewer.
    Allows selection of HGT/VGT genes and plots their presence/absence across strains.
    Optionally overlays phylum information using color-coded heatmap stripes.
    """
    def __init__(self, hgt_folder, vgt_folder, phylum_file=None):
        super().__init__()
        self.setWindowTitle("Gene Dotmap Viewer")
        self.resize(1800, 1000)

        # Load HGT and VGT matrix files from their respective folders
        self.hgt_files = self.load_matrix_files(hgt_folder)
        self.vgt_files = self.load_matrix_files(vgt_folder)

        # Initialize phylum-related attributes
        self.strain_to_phylum = {}
        self.phylum_colors = {}
        self.has_phylum_info = False

        # Load phylum info if provided
        if phylum_file:
            try:
                phylum_df = pd.read_excel(phylum_file, header=1)
                # Extract phylum name using regex from field
                phylum_df["Phylum"] = phylum_df["GTDB-Tk"].str.extract(r"p__([^;]+)")
                phylum_df["Phylum"] = phylum_df["Phylum"].str.replace("_+", "_", regex=True).str.strip("_")
                phylum_df["Phylum"] = phylum_df["Phylum"].replace(regex=r"^Firmicutes.*", value="Firmicutes")

                # Create mapping from strain to phylum
                self.strain_to_phylum = dict(zip(phylum_df["Strain"], phylum_df["Phylum"]))

                # Assign colors to each unique phylum
                unique_phyla = sorted(set(self.strain_to_phylum.values()))
                palette = px.colors.qualitative.Set3
                self.phylum_colors = dict(zip(unique_phyla, palette[:len(unique_phyla)]))
                self.has_phylum_info = True
            except Exception as e:
                print(f"[Warning] Could not load phylum mapping file: {e}")

        # === GUI Layout ===
        layout = QHBoxLayout(self)

        # --- Left Panel: Controls ---
        left_layout = QVBoxLayout()

        # HGT gene list
        hgt_label = QLabel("HGT Genes:")
        hgt_label.setStyleSheet("font-size: 13pt; font-weight: bold;")
        left_layout.addWidget(hgt_label)

        self.hgt_list = QListWidget()
        self.hgt_list.setSelectionMode(QAbstractItemView.MultiSelection)
        for gene in sorted(self.hgt_files.keys()):
            item = QListWidgetItem(gene)
            item.setCheckState(False)
            self.hgt_list.addItem(item)
        left_layout.addWidget(self.hgt_list)

        # VGT gene list
        vgt_label = QLabel("VGT Genes:")
        vgt_label.setStyleSheet("font-size: 13pt; font-weight: bold;")
        left_layout.addWidget(vgt_label)

        self.vgt_list = QListWidget()
        self.vgt_list.setSelectionMode(QAbstractItemView.MultiSelection)
        for gene in sorted(self.vgt_files.keys()):
            item = QListWidgetItem(gene)
            item.setCheckState(False)
            self.vgt_list.addItem(item)
        left_layout.addWidget(self.vgt_list)

        # Plot button
        self.plot_button = QPushButton("Plot Selected Genes")
        self.plot_button.setStyleSheet("font-size: 10pt;")
        self.plot_button.clicked.connect(self.plot_selected)
        left_layout.addWidget(self.plot_button)

        # Loading message label
        self.loading_label = QLabel("")
        self.loading_label.setStyleSheet("color: black; font-size: 10pt;")
        left_layout.addWidget(self.loading_label)

        left_layout.addStretch()
        layout.addLayout(left_layout, 1)

        # --- Right Panel: Plot Area ---
        self.web_view = QWebEngineView()
        layout.addWidget(self.web_view, 4)
        self.setLayout(layout)

    def load_matrix_files(self, folder):
        """
        Load comparison matrix CSV files from a folder.
        Only includes files ending with '_comparison_matrix_ordered_by_tree.csv'.
        Returns a dictionary mapping gene name to file path.
        """
        return {
            os.path.basename(f).split("_comparison_matrix")[0]: os.path.join(folder, f)
            for f in os.listdir(folder)
            if f.endswith("_comparison_matrix_ordered_by_tree.csv")
        }

    def plot_selected(self):
        """
        Main function to generate the dotmap from selected genes.
        Reads matrices, builds a binary map, adds Plotly traces, overlays optional phylum bars.
        """
        self.loading_label.setText("Generating plot, please wait...")
        self.plot_button.setEnabled(False)
        QApplication.processEvents()

        # Get selected gene names from both lists
        selected_hgt = [self.hgt_list.item(i).text() for i in range(self.hgt_list.count()) if self.hgt_list.item(i).checkState()]
        selected_vgt = [self.vgt_list.item(i).text() for i in range(self.vgt_list.count()) if self.vgt_list.item(i).checkState()]
        selected_genes = [(g, self.hgt_files[g], 'HGT') for g in selected_hgt] + [(g, self.vgt_files[g], 'VGT') for g in selected_vgt]

        if not selected_genes:
            self.loading_label.setText("")
            self.plot_button.setEnabled(True)
            return

        fig = go.Figure()
        x_labels = y_labels = None

        # Gene color legend
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=0, color='rgba(0,0,0,0)'),
            legendgroup="genes",
            showlegend=True,
            name="<b>Gene Colors</b>"
        ))

        # Color palettes for HGT and VGT
        hgt_colors = px.colors.qualitative.Plotly
        vgt_colors = px.colors.qualitative.Dark24
        hgt_index = vgt_index = 0

        # Loop over selected genes and create a scatter trace for each
        for gene, path, category in selected_genes:
            df = pd.read_csv(path, index_col=0)
            binary_df = (df > 0).astype(int)  #  binary presence/absence
            y, x = binary_df.values.nonzero()

            # Save axis labels only once
            if x_labels is None:
                x_labels = list(df.columns)
            if y_labels is None:
                y_labels = list(df.index)

            color = hgt_colors[hgt_index % len(hgt_colors)] if category == "HGT" else vgt_colors[vgt_index % len(vgt_colors)]
            hgt_index += category == "HGT"
            vgt_index += category == "VGT"

            # Add the gene presence trace
            fig.add_trace(go.Scattergl(
                x=x,
                y=y,
                mode='markers',
                marker=dict(symbol='square', size=6, color=color, opacity=0.8),
                name=gene,
                legendgroup="genes",
                showlegend=True,
                text=[  # Tooltip text
                    (
                        f"<b>Gene:</b> {gene}"
                        f"<br><b>X:</b> {x_labels[x[j]]}" +
                        (f" ({self.strain_to_phylum[x_labels[x[j]]]}" if self.has_phylum_info and x_labels[x[j]] in self.strain_to_phylum else "") +
                        (")" if self.has_phylum_info and x_labels[x[j]] in self.strain_to_phylum else "") +
                        f"<br><b>Y:</b> {y_labels[y[j]]}" +
                        (f" ({self.strain_to_phylum[y_labels[y[j]]]}" if self.has_phylum_info and y_labels[y[j]] in self.strain_to_phylum else "") +
                        (")" if self.has_phylum_info and y_labels[y[j]] in self.strain_to_phylum else "")
                    )
                    for j in range(len(x))
                ],
                hoverinfo='text'
            ))

        # Optional: Add phylum color heatmaps on axes
        if self.has_phylum_info:
            phyla_keys = list(self.phylum_colors.keys())
            x_phyla = [self.strain_to_phylum.get(label, 'Unknown') for label in x_labels]
            y_phyla = [self.strain_to_phylum.get(label, 'Unknown') for label in y_labels]
            x_colors = [1 + phyla_keys.index(p) if p in phyla_keys else 0 for p in x_phyla]
            y_colors = [[1 + phyla_keys.index(p) if p in phyla_keys else 0] for p in y_phyla]

            # Define color scale
            colorscale = [[0.0, "#DDDDDD"]] + [
                [(i + 1) / len(self.phylum_colors), color]
                for i, color in enumerate(self.phylum_colors.values())
            ]

            # Horizontal and vertical color bars
            fig.add_trace(go.Heatmap(
                z=[x_colors] * 6,
                x=list(range(len(x_labels))),
                y=[-i for i in range(1, 7)],
                showscale=False,
                colorscale=colorscale,
                hoverinfo='skip'
            ))
            fig.add_trace(go.Heatmap(
                z=np.column_stack([y_colors] * 6).tolist(),
                x=[-i for i in range(1, 7)],
                y=list(range(len(y_labels))),
                showscale=False,
                colorscale=colorscale,
                hoverinfo='skip'
            ))

            # Add phylum color legend
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=0, color='rgba(0,0,0,0)'),
                legendgroup="phyla",
                showlegend=True,
                name="<b>Phylum Colors</b>"
            ))
            for phylum, color in self.phylum_colors.items():
                fig.add_trace(go.Scatter(
                    x=[None], y=[None],
                    mode='markers',
                    marker=dict(size=10, color=color),
                    legendgroup="phyla",
                    name=phylum,
                    showlegend=True
                ))

        # Final layout updates
        fig.update_layout(
            title=dict(
                text="HGT & VGT Gene Dotmap" + (" with Phylum Color Bars" if self.has_phylum_info else ""),
                font=dict(size=24, family="Arial", color="black"),
                x=0.5
            ),
            width=1800,
            height=1400,
            margin=dict(b=300, t=80, l=300, r=300),
            xaxis=dict(
                tickmode='array',
                tickvals=list(range(len(x_labels))),
                ticktext=x_labels,
                tickangle=90,
                tickfont=dict(size=2),
                range=[-7, len(x_labels) + 2]
            ),
            yaxis=dict(
                tickmode='array',
                tickvals=list(range(len(y_labels))),
                ticktext=y_labels,
                tickangle=0,
                tickfont=dict(size=2),
                range=[-7, len(y_labels) + 2]
            ),
            hovermode='closest',
            legend=dict(
                tracegroupgap=20,
                font=dict(size=18),
                itemsizing='constant',
                orientation='v',
                x=1.02,
                y=1
            ),
            hoverlabel=dict(
                font=dict(size=22, color='black'),
                bgcolor="white",
                bordercolor="black"
            )
        )

        # Write to temporary HTML and render in the web view
        with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as f:
            fig.write_html(f.name)
            self.web_view.load(QUrl.fromLocalFile(f.name))

        QTimer.singleShot(3000, lambda: self.loading_label.setText(""))
        self.plot_button.setEnabled(True)

# === Main ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gene Dotmap Viewer")
    parser.add_argument("--hgt_folder", required=True)
    parser.add_argument("--vgt_folder", required=True)
    parser.add_argument("--phylum_file", required=False)
    args = parser.parse_args()

    app = QApplication(sys.argv)
    viewer = DotmapViewer(args.hgt_folder, args.vgt_folder, args.phylum_file)
    viewer.show()
    sys.exit(app.exec_())
