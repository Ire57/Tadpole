import streamlit as st
import streamlit.components.v1 as components
import traceback
import matplotlib.pyplot as plt
import graphviz
import glob
import sbol3
import os
import RNA
# Import using absolute rout so it can be used as a package and launch the streamlit app
from tadpole import *
from tadpole.rna_cluster import matriz_rmsd, cluster_structures, organise_per_clusters
from tadpole.rna_cluster import compute_metrics, visualise_metrics
from tadpole.input_utils import parse_fasta_msa
from tadpole.structure import predict_secundary_structure
from tadpole.visualizacion import generate_pre_command_per_sequence
from tadpole.search import run_linker_search
from tadpole.search import count_rna1_rna3_pairings, constrained_mfe
from tadpole.structure import predict_secundary_structure, align_secondary_structure
from tadpole.rna_structures import get_base_pairs
from tadpole.genetic_algorithm import run_genetic_algorithm_search
from tadpole.conservation import calculate_conservation, clasify_conservation, conservation_category, calculate_complementarity_conservation
from tadpole.io_tools import create_zip_archive
from tadpole.SBOL import export_single_linker
from scipy.spatial.distance import pdist, squareform
import numpy as np
import shutil
import subprocess
import time
from io import StringIO
from contextlib import redirect_stdout
# from weasyprint import HTML # Weasyprint is typically for PDF generation and might require specific server setup not always available in Streamlit Cloud
import tempfile
import datetime
from collections import defaultdict
import base64 
import pandas as pd # Added pandas import as it's used later for dataframes

def plot_delta_mfe(results, output_dir="report_images", bins=20):
    """
    Generates a histogram of the ΔMFE distribution (mfe_constrained – mfe_unconstrained)
    and saves it as an image in output_dir/delta_mfe.png.

    This function calculates the difference in minimum free energy (ΔMFE) between constrained
    and unconstrained structures for a list of results, and visualizes the distribution
    using a histogram. The resulting plot is saved as a PNG file.

    :param results: A list of dictionaries, each containing the keys 'mfe_1' and 'mfe_2',
                    representing the unconstrained and constrained MFE values, respectively. (list of dict)
    :param output_dir: Directory where the plot image will be saved. Created if it doesn't exist. (str)
    :param bins: Number of bins to use in the histogram. (int)

    :returns: The full file path to the saved plot image, or None if there is no data to plot. (str or None)
    """
    # Extract ΔMFE values
    deltas = [res['mfe_2'] - res['mfe_1'] for res in results]
    if not deltas:
        return None  # nothing to plot

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "delta_mfe.png")

    # Create histogram figure
    plt.figure(figsize=(6, 4))
    plt.hist(deltas, bins=bins, edgecolor='black')
    plt.xlabel("ΔMFE = MFE_constrained − MFE_unconstrained (kcal/mol)")
    plt.ylabel("Number of designs")
    plt.title("ΔMFE Distribution")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    return out_path

def png_to_data_uri(png_path: str) -> str:
    """
    Reads a PNG file and returns a data URI string for embedding in HTML.

    This function takes the path to a PNG image file, encodes its contents
    in base64, and returns a properly formatted data URI. This is useful
    for embedding images directly into HTML without needing external file references.

    :param png_path: Path to the PNG image file. (str)

    :returns: A data URI string representing the PNG image, suitable for HTML embedding. (str)
    """
    with open(png_path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"

def plot_pairings_histogram(results, rna1, rna3, output_dir="report_images", bins=None):
    """
    Generates a histogram of the number of base pairings between RNA1 and RNA3
    for each design in `results`, and saves it to output_dir/pairings_count.png.

    This function calculates how many base pairs form between RNA1 and RNA3 in each
    RNA structure from the provided results. It then visualizes the distribution
    of these counts using a histogram, saved as a PNG image.

    :param results: A list of dictionaries, each containing RNA structure data and linker. (list of dict)
    :param rna1: The sequence of the first RNA fragment. (str)
    :param rna3: The sequence of the third RNA fragment. (str)
    :param output_dir: Directory where the histogram image will be saved. (str)
    :param bins: Number of histogram bins. If None, matplotlib chooses automatically. (int or None)

    :returns: The full path to the saved histogram image, or None if no data is available. (str or None)
    """
    # Recompute the number of RNA1–RNA3 pairings for each design
    counts = [
        count_rna1_rna3_pairings(
            res['structure_unconstrained'],  # or structure_constrained, depending on your use case
            rna1, rna3, res['linker']
        )
        for res in results
    ]
    if not counts:
        return None

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "pairings_count.png")

    # Create histogram figure
    plt.figure(figsize=(6, 4))
    plt.hist(counts, bins=bins, edgecolor='black')
    plt.xlabel("Number of RNA1–RNA3 base pairings")
    plt.ylabel("Number of designs")
    plt.title("Distribution of RNA1–RNA3 Base Pairings")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    return out_path

def build_full_html_report(
    results,
    report,
    rna1,
    rna3,
    struct1,
    constraint,
    mutable_rna1,
    watched_positions,
    use_mutations,
    mfe_delta,
    max_pairings,
    max_changes,
    num_mut,
    linker_min,
    linker_max,
    verbose,
    cluster_labels,
    search_method,
    image_output_dir="estructura_img",
    representative_img_bases=None
):
    """
    Generates a complete HTML report for RNA switch design results.

    This function assembles a structured, styled HTML report summarizing
    the outcomes of a computational RNA switch design pipeline. It includes
    global statistics, search parameters, result clustering, representative images,
    and automated analysis and conclusions.

    :param results: List of design results, each as a dictionary with relevant metrics. (list of dict)
    :param report: String of the report, not used in this function. (str)
    :param rna1: RNA1 sequence string. (str)
    :param rna3: RNA3 sequence string. (str)
    :param struct1: Secondary structure of RNA1 in dot-bracket notation. (str)
    :param constraint: Constraint string used for ON state simulation. (str)
    :param mutable_rna1: List of mutable positions in RNA1. (list of int)
    :param watched_positions: List of positions to track for functional conservation. (list of int)
    :param use_mutations: Whether mutations were allowed in RNA1. (bool)
    :param mfe_delta: Minimum required ΔMFE (ON–OFF energy difference). (float)
    :param max_pairings: Maximum allowed RNA1–RNA3 base pairings. (int)
    :param max_changes: Maximum allowed structural changes in RNA1. (int)
    :param num_mut: Maximum number of mutations allowed in RNA1. (int)
    :param linker_min: Minimum linker length. (int)
    :param linker_max: Maximum linker length. (int)
    :param verbose: Whether to include verbose information. (bool)
    :param cluster_labels: List of clustering labels for each design. (list of int)
    :param search_method: Label for the method used to perform the design search. (str)
    :param image_output_dir: Path to directory containing result image files. (str)
    :param representative_img_bases: List of base filenames for representative structure images. (list of str or None)

    :returns: A string containing the complete HTML content. (str)
    """

    # Paths to your PNGs for the global plots (if they exist)
    delta_uri = png_to_data_uri(f"{image_output_dir}/delta_mfe.png")
    pairings_uri = png_to_data_uri(f"{image_output_dir}/pairings_count.png")

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    total_found = len(results)

    # Ensure total_clusters correctly counts non-noise clusters
    total_clusters = len(set(c for c in cluster_labels if c != -1)) if len(cluster_labels) > 0 else 0

    # Calculate global metrics
    diffs = [res['mfe_2'] - res['mfe_1'] for res in results]
    avg_delta = sum(diffs) / len(diffs) if diffs else 0
    avg_pairs = sum(
        res.get('pairings_count', max_pairings) for res in results
    ) / total_found if total_found else 0
    
    # Group results by cluster
    cluster_dict = defaultdict(list)
    for idx, lbl in enumerate(cluster_labels):
        if lbl != -1:
            cluster_dict[lbl].append((idx, results[idx]))

    cluster_representatives = {}
    for cluster_id, items in sorted(cluster_dict.items()):
        # min() devuelve la tupla (idx, res) del enlazador representativo
        representative_tuple = min(items, key=lambda r: abs((r[1]['mfe_2'] - r[1]['mfe_1']) - mfe_delta))
        cluster_representatives[cluster_id] = representative_tuple

    # Create a sorted list of clusters based on the delta_mfe of their representative
    sorted_clusters_energy = sorted(
        cluster_representatives.items(),
       
        key=lambda item: abs((item[1][1]['mfe_2'] - item[1][1]['mfe_1']) - mfe_delta)
    )
    sorted_clusters_pairings = sorted(
        cluster_representatives.items(),
        
        key=lambda item: abs(count_rna1_rna3_pairings(
            item[1][1]['structure_unconstrained'],
            item[1][1].get('rna1_mutated_seq', item[1][1].get('rna1_sequence', 'N/A')),
            item[1][1]['rna3'],
            item[1][1]['linker']
        ) - max_pairings)
    )
    # Start building the HTML
    html = f"""
    <html>
    <head>
      <meta charset="utf-8"/>
      <style>/* Estilos para el modo oscuro */
    
        body {{
          margin: 40px; font-family: sans-serif; color: #1f2937;
          white-space: pre-wrap; word-wrap: break-word; overflow-wrap: break-word;
        }}
        h1 {{ color: #1e40af; }}
        h2 {{ color: #2563eb; border-bottom: 2px solid #2563eb; padding-bottom: 4px; }}
        h3 {{ color: #3b82f6; }}
        .section {{ margin-bottom: 30px; }}
        .linker-box {{ border:1px solid #d1d5db; border-radius:8px; padding:15px; background:#f9fafb; }}
        code {{ background:#f3f4f6; padding:3px 5px; border-radius:4px; font-family: monospace; }}
        .small {{ font-size:0.85rem; color:#6b7280; }}
        .metric-table td, .metric-table th, .pros-cons-table td, .pros-cons-table th {{ padding:8px 12px; border:1px solid #e5e7eb; }}
        .metric-table, .pros-cons-table {{ border-collapse: collapse; width:100%; }}
        .pros-cons-table th {{ background-color:#eff6ff; }}
        .img-preview {{ max-width:300px; margin:8px; vertical-align:top; }}
        ul {{ list-style-type: disc; margin-left: 20px; }}
        ol {{ list-style-type: decimal; margin-left: 20px; }}
      </style>
    </head>
    <body>
      <h1> SRE–Aptamer Switch Design Report</h1>
      <p class="small">Generated: {now}</p>

      <div class="section">
        <h2> Introduction and Objectives</h2>
        <p>In this study, we aim to convert a functional element of RNA1 (e.g., a frameshift or SCR) into a molecular "off→on" switch regulated by an aptamer (RNA2) and its ligand. This is achieved by designing a "linker" that connects RNA1 with RNA2, allowing RNA1 to toggle between an "off" (inactive) and an "on" (active) state in response to ligand binding.</p>
    </div>

      <div class="section">
        <h2> Search Parameters</h2>
        <ul>
          <li><strong>Search Method:</strong> {search_method}</li>
          <li>Mutations in SRE: {"Yes" if use_mutations else "No"}</li>
          <li>Mutable positions: {', '.join(map(str, mutable_rna1)) or '—'}</li>
          <li>Watched positions: {', '.join(map(str, watched_positions)) or '—'}</li>
          <li>Maximum allowed mutations: {num_mut or '—'}</li>
          <li>Linker length: {linker_min}–{linker_max} nt</li>
          <li>Minimum required ΔMFE: {mfe_delta:.1f} kcal/mol</li>
          <li>Maximum SRE-Aptamer pairings: {max_pairings}</li>
          <li>Maximum structural changes in SRE: {max_changes}</li>
          <li>Verbose: {"Yes" if verbose else "No"}</li>
        </ul>
      </div>
    """

    html += f"""
      <div class="section">
        <h2> Summary of Results</h2>
        <p>A total of <strong>{total_found}</strong> linkers meeting the search criteria were found, grouped into <strong>{total_clusters}</strong> distinct structural clusters. The average global metrics for the solutions found are:</p>
        <table class="metric-table">
          <tr><th>Total linkers found</th><td>{total_found}</td></tr>
          <tr><th>Total distinct clusters</th><td>{total_clusters}</td></tr>
          <tr><th>Average ΔMFE (on–off)</th><td>{avg_delta:.2f} kcal/mol</td></tr>
          <tr><th>Average SRE-Aptamers pairings</th><td>{avg_pairs:.1f}</td></tr>
        </table>
      </div>

            <div class="section">
            <h2> Detailed Analysis by Cluster</h2>
            <p>Below is a list of the representative sequences of all the clusters, ordered by the proximity of the representative linker's ΔMFE to the desired value ({mfe_delta:.1f} kcal/mol).

    """

    # Dynamically add sorted cluster analysis rows
    for cluster_id, representative_result in sorted_clusters_energy:
        rep_delta = representative_result[1]['mfe_2'] - representative_result[1]['mfe_1']
        html += f"""
                Cluster {cluster_id}:{rep_delta:.2f} kcal/mol
        """

    html += f"""
    <p>Below is a list of the representative sequences of all the clusters, ordered by the proximity of the number of pairings between the SRE and the Aptamer to the desired value ({max_pairings:.1f} ).    
    """
    for cluster_id, representative_result in sorted_clusters_pairings:
        html += f"""
                Cluster {cluster_id}:{count_rna1_rna3_pairings(representative_result[1]['structure_unconstrained'], representative_result[1].get('rna1_mutated_seq', representative_result[1].get('rna1_sequence', 'N/A')), representative_result[1]['rna3'], representative_result[1]['linker'])} 
        """
    html += f"""
    <p>How to choose the best option? </p>
    <div class="section">
        <p>The key aspects to consider for a successful ON/OFF system design are:</p>
        <ol>
          <li><strong>ΔMFE:</strong> You have above a list of the representatives of each cluster sorted by how close the ΔMFE values are to the desired value (input). The best option for this feature will be the first element of the list.</li>
          <li><strong>SRE-Aptamer Interactions:</strong> This interactions should be enough to disrupt the structure of the SRE o the OFF state, 
          but at the same time, too many pairings could make it impossible to separate the SRE from the Aptamer. If the input value for 'Maximum number of SRE-aptamer pairings'
           was balanced, the best option for this parameter are the top elements of the list above, as they are ordered by proximity to the input value. </li>
          <li><strong>Aptamer Entry and Binding Sites:</strong> The key binding sites of the aptamer should be mostly free (on the OFF state) to enable efficient ligand binding.</li>
          <li><strong>Disruption of the OFF State:</strong> Key nucleotides for the SRE's function must be identified and their functional structure disrupted in the unbound (OFF) state to ensure inactivity.</li>
          <li><strong>Maintenance of the ON State:</strong> Structures responsible for the SRE's function must be preserved in the bound (ON) state to ensure activity.</li>
        </ol>
        
      </div>
    """
    aptamer_uri = png_to_data_uri(f"tadpole_package/images/aptamer_report.jpg")
    html += f"""
                    <div style="display:inline-block; margin-right:20px; vertical-align:top;">
                    <p>If you are using the default aptamer (theophillyne) the binding sites are coloured in yellow in the following image: </p>
                      <img class="img-preview" src="{aptamer_uri}" alt="Unconstrained structure - Cluster {cluster_id}" />
                    </div>
                    """

    if not results:
        html += """
      <div class="section">
        <h2>⚠️ No Valid Results Found</h2>
        <p>No linker satisfying all criteria was identified. Suggestions:</p>
        <ul>
          <li>Relax the minimum ΔMFE to allow a smaller energy difference.</li>
          <li>Reduce RNA1–RNA3 pairing restrictions to permit more unwanted interactions but get more results.</li>
          <li>Expand the linker length range to explore a broader search space.</li>
          <li>Allow more mutations in RNA1 if greater design flexibility is desired.</li>
        </ul>
      </div>
        """
    else:
        cluster_dict = defaultdict(list)
        for idx, lbl in enumerate(cluster_labels):
            if lbl != -1:
                cluster_dict[lbl].append((idx, results[idx]))
        html_sections_img = ""
        if representative_img_bases:
            
            
            cluster_dict = defaultdict(list)
            for idx, lbl in enumerate(cluster_labels):
                if lbl != -1:
                    cluster_dict[lbl].append((idx, results[idx]))

            # Prepare a list of all clusters with their representative data
            cluster_data = []
            for cluster_id, items in sorted(cluster_dict.items()):
                # Find the representative linker using your logic
                idx, res = min(items, key=lambda r: abs((r[1]['mfe_2'] - r[1]['mfe_1']) - mfe_delta))

                # Store all the necessary information in a single tuple or dictionary
                # so we can iterate over it just once
                cluster_data.append({
                    'cluster_id': cluster_id,
                    'variants_count': len(items),
                    'representative_result': res,
                    'image_base': res.get('filename')  # Use the filename from the representative result
                })
            for i, img_base in enumerate(representative_img_bases):
                cluster_data[i]['image_base'] = img_base  # Use the filename from the representative result
                
            # Now, loop through the prepared data once to display everything
            for data in cluster_data:
                cluster_id = data['cluster_id']
                res = data['representative_result']
                img_base = data['image_base']
                variants_count = data['variants_count']

                html += f"""
                    <div class="section">
                        <h2>🧬 Cluster Details {cluster_id} — {variants_count} variants</h2>
                        <div class="linker-box">
                            <h3>Example Linker in Cluster {cluster_id}: <code>{res['linker']}</code></h3>
                            <p><strong>Full Sequence (RNA1+Linker+RNA3):</strong> <code>{res['sequence']}</code></p>
                            <p><strong>Mutations in RNA1:</strong> {res.get('mut1_info', 'N/A')}</p>
                            <p><strong>ΔMFE (Off→On):</strong> {(res['mfe_2']-res['mfe_1']):.2f} kcal/mol —
                              <strong>RNA1–RNA3 Pairings:</strong> {count_rna1_rna3_pairings(res['structure_unconstrained'], res.get('rna1_mutated_seq', res.get('rna1_sequence', 'N/A')), res['rna3'], res['linker'])}</p>
                        </div>
                    </div>
                    """

                # Check if a filename exists before trying to generate the image paths
                if img_base:
                    off_uri = png_to_data_uri(f"{image_output_dir}/linker_{linker_min}/{img_base}_unconstrained_plot.png")
                    on_uri = png_to_data_uri(f"{image_output_dir}/linker_{linker_min}/{img_base}_constrained_plot.png")

                    html += f"""
                    <div style="display:inline-block; margin-right:20px; vertical-align:top;">
                      <h3>Cluster {cluster_id}</h3>
                      <p><em>Off (unconstrained)</em></p>
                      <img class="img-preview" src="{off_uri}" alt="Unconstrained structure - Cluster {cluster_id}" />
                      <p><em>On (constrained)</em></p>
                      <img class="img-preview" src="{on_uri}" alt="Constrained structure - Cluster {cluster_id}" />
                    </div>
                    """
                else:
                    html += f"""
                    <div style="display:inline-block; margin-right:20px; vertical-align:top;">
                        <h3>Cluster {cluster_id}</h3>
                        <p>No representative image was found for this cluster.</p>
                    </div>
                    """
                    html_sections_img += "</div>"

        html += html_sections_img

        
        
        
            

        html += f"""
      <div class="section">
        <p>🔍 Full images and data for each linker are available at: <code>{image_output_dir}/</code></p>
      </div>
        """

    
    html += f"""
      

      <div class="section">
        <h2>💡 Possible Improvements and Next Steps</h2>
        <p>To further optimize the RNA switch design, consider the following recommendations (for more information on the parameters, visit the 'Help' section of Tadpole):</p>
        <ul>
          <li><strong>Adjust Search Parameters:</strong>
            <ul>
              <li>If the number of results is low or metrics are suboptimal, try relaxing the minimum ΔMFE, reducing SRE-Aptamer pairing restrictions, expanding the linker length range or adding mutations to the SRE. For the Geentic algorithm, try incrementing the population or generations</li>
              <li>For more complex systems or greater diversity, experiment by allowing more mutations in SRE or tuning genetic algorithm parameters: greater mutation rate leads to more diversity. </li>
            </ul>
          </li>
          <li><strong>Additional Cluster Analysis:</strong>
            <ul>
              <li>The representative linker for each cluster is the one that has the ΔMFE closest to the input. </li>
              <li>For the Genetic Algorithm, if there is low variety among clusters, consider incrementing the mutation rate and/or the population.</li>
            </ul>
          </li>
          <li><strong>Exploration of Alternative Aptamers:</strong>
            <ul>
              <li>If performance remains suboptimal, consider using or designing aptamers with different binding or stability characteristics that may better integrate with your SRE.</li>
            </ul>
          </li>
        </ul>
      </div>

    </body>
    </html>
    """

    return html

def run_command(cmd):
    """
    Executes a shell command and returns its standard output.

    This function runs a shell command using Python's subprocess module. If the command
    is successful, the standard output is returned. If the command fails, an error message
    is displayed using Streamlit and the function returns None.

    :param cmd: Shell command to execute. (str)

    :returns: The standard output of the command if successful, otherwise None. (str or None)
    """
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as e:
        st.error(f"Command failed:\n`{cmd}`\n{e.stderr}")
        return None

def stop_requested():
    """
    Checks whether a stop request has been triggered in the Streamlit session.

    This function looks for the "stop_requested" flag in Streamlit's session state
    to determine if the user has requested to stop a running process.

    :returns: True if a stop has been requested, False otherwise. (bool)
    """
    return st.session_state.get("stop_requested", False)

log_area = st.empty()

def update_log(msg):
    """
    Appends a message to the session log and updates the displayed log area.

    This function adds a new line to the "log_text" entry in Streamlit's session state,
    maintaining a running log of messages. It also updates the visible log area in the app.

    :param msg: The message string to append to the log. (str)
    """
    current_text = st.session_state.get("log_text", "")
    st.session_state["log_text"] = current_text + msg + "\n"
    log_area.code(st.session_state["log_text"], language="text")
