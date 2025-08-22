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

st.set_page_config(page_title="Tadpole", layout="wide", page_icon="logo.png") # Updated page title

# ---------------- CUSTOM STYLES ----------------
st.markdown("""
    <style>
            
    /* Global Styles */
    html, body, [class*="css"] {
        font-family: 'Helvetica Neue', sans-serif;
        background-color: #f8f8f8; /* Light background for the entire app */
        color: #333333; /* Default text color */
    }

    /* Primary Color Palette based on #e52920 */
   :root {
    --bg-light-red: #fdeeee;      /* Very light red/pink for backgrounds */
    --soft-red: #fbdad9;          /* Light red/pink for inactive elements, light card borders */
    --medium-light-red: #f7bcbb;  /* Slightly darker pink */
    --light-accent-red: #f28f8c;  /* Medium pink/light red for hover states */
    --primary-red: #ea534e;       /* Adjusted to #ea534e as per new palette, was #e52920 */
    --dark-accent-red: #e62e25;   /* Slightly darker primary red */
    --dark-red: #e41513;          /* Darker red for button hover/active */
    --very-dark-red: #be1818;     /* Even darker red */
    --deep-red: #9d1915;          /* Dark red/maroon for strong accents/borders */
    --text-darker-red: #53110e;   /* Very dark red for text/headers */
    
    --soft-gray: #e0e0e0;         /* Soft gray for general backgrounds/borders (kept as is) */
    --medium-gray: #cccccc;       /* Medium gray for borders/dividers (kept as is) */
    --dark-gray: #1a1a1a;         /* General dark text/header color (kept as is) */
    --text-color-light: #333333;  /* Dark text on light backgrounds (kept as is) */
    --text-color-dark: #ffffff;   /* White text on dark primary colors (kept as is) */
    --card-bg: #ffffff;           /* White background for cards/panels (kept as is) */
    --border-color: var(--soft-gray); /* Default border color (kept as is) */
    --shadow-color: rgba(0, 0, 0, 0.1); /* Subtle shadow (kept as is) */
}
    }

    /* Headers - Ensure good contrast */
    h1, h2, h3, h4, h5, h6 {
        color: var(--text-darker-red); /* Using the darkest red for headers */
        font-weight: 600;
    }
    h1 { font-size: 2.5em; color: var(--primary-red); margin-bottom: 0.5em; } /* Keep H1 primary red */
    h2 { font-size: 2em; margin-bottom: 0.5em; }
    h3 { font-size: 1.5em; margin-bottom: 0.5em; }
    

    /* Main Container Padding */
    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
        padding-left: 3rem;
        padding-right: 3rem;
    }

    /* Tabs Styling - Professional and clear */
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
        margin-bottom: 1.5rem;
        border-bottom: 2px solid var(--deep-red); /* Stronger separator for main tabs */
    }
    .stTabs [data-baseweb="tab"] {
        background-color: var(--soft-red); /* Light red for inactive tabs */
        border-radius: 8px 8px 0 0;
        padding: 10px 20px;
        color: var(--text-darker-red); /* Dark text on light red */
        font-weight: 500;
        transition: background-color 0.3s, color 0.3s;
        border: 1px solid var(--soft-red); /* Subtle border */
        border-bottom: none; /* No bottom border for inactive */
    }
    .stTabs [data-baseweb="tab"]:hover {
        background-color: var(--light-accent-red); /* Lighter red on hover */
        color: var(--text-color-dark); /* White text on hover */
    }
    .stTabs [aria-selected="true"] {
        background-color: var(--primary-red); /* Primary red for active tab */
        color: var(--text-color-dark); /* White text on active tab */
        border: 1px solid var(--primary-red); /* Border matches active color */
        border-bottom: none; /* No bottom border for active */
        font-weight: 600; /* Bolder for active tab */
    }

    /* Buttons - Clear action, professional look */
    .stButton>button {
        background-color: var(--primary-red);
        color: var(--text-color-dark); /* White text */
        border-radius: 8px;
        border: none;
        padding: 0.7em 1.2em;
        font-weight: 600;
        transition: background-color 0.3s, transform 0.2s;
        box-shadow: 0 2px 5px var(--shadow-color);
        cursor: pointer;
    }
    .stButton>button:hover {
        background-color: var(--dark-red); /* Darker red on hover */
        color: #FFF;
        transform: translateY(-2px);
    }
    .stButton>button:active {
        background-color: var(--dark-red);
        color: #FFF;
        transform: translateY(0);
        box-shadow: none;
    }

    /* File Uploader - Clear visual cue */
    .stFileUploader {
        border: 2px dashed var(--medium-gray); /* Dashed border for upload area */
        border-radius: 8px;
        padding: 1rem;
        text-align: center;
        background-color: var(--card-bg);
        margin-top: 1rem;
        margin-bottom: 1rem;
    }

    /* Expander/Collapsible Panels - Clean and distinguishable */
    
    .expander_div {
        background-color: var(--soft-gray); /* Soft gray header */
        border-radius: 8px;
        padding: 10px 15px;
        font-weight: 500;
        color: var(--text-color-light);
        transition: background-color 0.3s;
        border: 1px solid var(--soft-gray);
        font-size: 1.2rem;
    }
    .expander_div:hover {
        background-color: #d0d0d0; /* Slightly darker gray on hover */
    }
    .expander_subdiv {
        background-color: var(--card-bg); /* White content area */
        border: 1px solid var(--soft-gray);
        border-top: none; /* No top border to blend with header */
        border-radius: 0 0 8px 8px;
        padding: 15px;
    }

    /* Metrics - Highlighted key information */
    [data-testid="stMetric"] {
        background-color: var(--card-bg);
        border-radius: 8px;
        padding: 1rem;
        box-shadow: 0 2px 5px var(--shadow-color);
        border: 1px solid var(--border-color);
        margin-bottom: 1rem;
    }
    [data-testid="stMetricLabel"] {
        color: var(--primary-red); /* Red label for importance */
        font-weight: 600;
        font-size: 0.9em;
    }
    [data-testid="stMetricValue"] {
        color: var(--text-darker-red); /* Darker red for metric values */
        font-size: 1.8em;
        font-weight: 700;
    }
    [data-testid="stMetricDelta"] {
        color: #4CAF50; /* Green for positive delta, adjust for negative if needed */
    }

    /* Text Input & Text Area - Clean borders, clear focus */
    .stTextInput>div>div>input, .stTextArea>div>div>textarea {
        border: 1px solid var(--medium-gray); /* Medium gray border */
        border-radius: 8px;
        padding: 8px 12px;
        transition: border-color 0.3s, box-shadow 0.3s;
        color: var(--text-color-light);
        background-color: var(--card-bg); /* Explicitly set background to white */
    }
    .stTextInput>div>div>input:focus, .stTextArea>div>div>textarea:focus {
        border-color: var(--primary-red); /* Red border on focus */
        box-shadow: 0 0 0 1px var(--primary-red); /* Red shadow on focus */
        outline: none; /* Remove default outline */
    }

    /* Info/Warning/Error Messages - Clear visual cues */
    .stAlert {
        border-radius: 8px;
        margin-bottom: 1rem;
    }
    .stAlert > div > div {
        /* Default background and text color for info */
        background-color: #e6f7ff; 
        color: #004085;
        border-left: 5px solid #007bff;
    }
    .stAlert.st-warning > div > div { /* Specific for warning */
        background-color: #fff3cd; 
        color: #664d03;
        border-left: 5px solid #ffc107;
    }
    .stAlert.st-success > div > div { /* Specific for success */
        background-color: #d4edda; 
        color: #0f5132;
        border-left: 5px solid #28a745;
    }
    .stAlert.st-error > div > div { /* Specific for error */
        background-color: #f8d7da; 
        color: #842029;
        border-left: 5px solid #dc3545;
    }

    

    /* General containers for visual separation - "Panels" */
    .stContainer {
        background-color: var(--card-bg);
        border-radius: 8px;
        padding: 1.5rem;
        margin-bottom: 1.5rem;
        box-shadow: 0 2px 5px var(--shadow-color);
        border: 1px solid var(--border-color);
    }

    /* Dataframes - Clean and readable tables */
    .stDataFrame {
        border: 1px solid var(--medium-gray);
        border-radius: 8px;
        overflow: hidden; /* Ensures border-radius applies to content */
    }
    .stDataFrame > div > div {
        background-color: var(--card-bg);
    }
    .stDataFrame th {
        background-color: var(--soft-gray);
        color: var(--text-color-light);
        font-weight: 600;
        border-bottom: 1px solid var(--medium-gray);
    }
    .stDataFrame td {
        color: var(--text-color-light);
        border-bottom: 1px solid var(--soft-gray);
    }
    .stDataFrame tr:nth-child(even) {
        background-color: #f5f5f5; /* Subtle stripe for readability */
    }

    
     /* sidebar style*/
    [data-testid="stSidebar"] {
        background-color: #0e1117;
        border-right: 1px solid var(--medium-gray);
        padding-top: 1rem;
    }
    
    /* Sidebar header*/
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h1,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h2,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h3 {
        color: var(--deep-red);
        font-weight: 600;
        padding-bottom: 0.5rem;
        border-bottom: 1px solid var(--light-red);
        margin-bottom: 1.5rem;
    }
    
    /* Radio buttons on the sidebar container */
    [data-testid="stSidebar"] .stRadio {
        margin-bottom: 2rem;
    }
    
    /* Style for options buttons */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label > div {
        background-color: var(--sidebar-item-bg);
        border-radius: 8px;
        padding: 12px 15px;
        margin-bottom: 0.75rem;
        transition: all 0.2s ease;
        color: #FFF;
        font-weight: 500;
        border-left: 4px solid transparent;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    
    /* Hover state for options for radio button */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label:hover > div {
        background-color: var(--sidebar-item-hover);
        border-left: 4px solid var(--light-red);
        transform: translateX(3px);
    }
    
    /* Active state (selected) for radio buttom options */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] input[type="radio"]:checked + div {
        background-color: var(--sidebar-item-active);
        color: var(--text-color-dark);
        border-left: 4px solid var(--deep-red);
        box-shadow: 0 2px 8px rgba(229, 41, 32, 0.3);
        transform: translateX(5px);
    }
    
    /* Sidebar sections */
    [data-testid="stSidebar"] .sidebar-section {
        background-color: var(--card-bg);
        border-radius: 8px;
        padding: 1rem;
        margin-bottom: 1.5rem;
        border: 1px solid var(--medium-gray);
        box-shadow: 0 2px 5px rgba(0,0,0,0.05);
    }
    
    /* Section title sidebar */
    [data-testid="stSidebar"] .sidebar-section-title {
        color: var(--deep-red);
        font-weight: 600;
        font-size: 1rem;
        margin-bottom: 0.75rem;
        padding-bottom: 0.5rem;
        border-bottom: 1px solid var(--soft-red);
    }
    
    /* Sidebar buttons */
    [data-testid="stSidebar"] .stButton > button {
        background-color: var(--soft-red);
        color: #FFF;
        border: 1px solid var(--light-red);
        border-radius: 6px;
        font-weight: 500;
        transition: all 0.2s ease;
        width: 100%;
        margin-bottom: 0.5rem;
    }
    
    /* Hover state for botones on sidebar */
    [data-testid="stSidebar"] .stButton > button:hover {
        background-color: var(--light-red);
        color: #FFF;
        border-color: var(--primary-red);
        transform: translateY(-1px);
        box-shadow: 0 2px 5px rgba(229, 41, 32, 0.2);
    }
    
    
    
    
    /* Sidebar Metrix*/
    [data-testid="stSidebar"] [data-testid="stMetric"] {
        background-color: var(--soft-red);
        padding: 0.75rem;
        border-radius: 6px;
        margin-bottom: 0.75rem;
    }
    
    [data-testid="stSidebar"] [data-testid="stMetric"] label {
        color: var(--deep-red) !important;
        font-weight: 600 !important;
    }
    
    [data-testid="stSidebar"] [data-testid="stMetric"] [data-testid="metric-container"] {
        color: var(--text-color-light) !important;
    }
    
    /* Dataframes - Clean and readable tables */
    .stDataFrame {
        border: 1px solid var(--medium-gray);
        border-radius: 8px;
        overflow: hidden; /* Ensures border-radius applies to content */
    }
    .stDataFrame > div > div {
        background-color: var(--card-bg);
    }
    .stDataFrame th {
        background-color: var(--soft-gray);
        color: var(--text-color-light);
        font-weight: 600;
        border-bottom: 1px solid var(--medium-gray);
    }
    .stDataFrame td {
        color: var(--text-color-light);
        border-bottom: 1px solid var(--soft-gray);
    }
    .stDataFrame tr:nth-child(even) {
        background-color: #f5f5f5; /* Subtle stripe for readability */
    }
    /* Main header style */
    .main-header {
        background: linear-gradient(135deg, var(--header-gradient-start), var(--header-gradient-end));
        border-radius: 12px;
        box-shadow: 0 8px 24px var(--header-shadow);
        margin-bottom: 2rem;
        padding: 2.5rem 2rem;
        border: 1px solid var(--header-border);
        position: relative;
        overflow: hidden;
    }
    
    /* Header decorations */
    .header-decoration {
        position: absolute;
        top: 0;
        right: 0;
        width: 100%;
        height: 100%;
        background-image: 
            linear-gradient(45deg, var(--text-darker-red) 25%, transparent 25%),
            linear-gradient(-45deg, var(--text-darker-red) 25%, transparent 25%);
        background-size: 12px 12px;
        background-position: 0 0, 6px 0;
        opacity: 0.03;
        z-index: 0;
        clip-path: polygon(70% 0, 100% 0, 100% 100%, 90% 100%);
    }

    .header-decoration-2 {
        position: absolute;
        bottom: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-image: 
            linear-gradient(45deg, var(--soft-red) 25%, transparent 25%),
            linear-gradient(-45deg, var(--soft-red) 25%, transparent 25%);
        background-size: 10px 10px;
        background-position: 0 0, 5px 0;
        opacity: 0.02;
        z-index: 0;
        clip-path: polygon(0 30%, 0 100%, 30% 100%);
    }
    
    /* Header content */
    .header-content {
        position: relative;
        z-index: 1;
        text-align: center;
    }
    
    /* MAin title */
    .app-title {
        font-size: 3rem;
        font-weight: 700;
        color: var(--deep-red);
        margin-bottom: 0.75rem;
        letter-spacing: -0.5px;
        text-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    
    /* subtitle */
    .app-description {
        font-size: 1.15rem;
        color: var(--text-color-light);
        max-width: 800px;
        margin: 0 auto;
        line-height: 1.6;
    }
</style>
""", unsafe_allow_html=True)


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
    aptamer_uri = png_to_data_uri(f"images/aptamer_report.jpg")
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

# --- Content for "structural RNA element" sub-tab ---


# ===== Estilos auxiliares =====
def example_box(text):
    st.markdown(f"""
    <div style="
        border-left: 5px solid #17a2b8; 
        background-color: #e8f7fa; 
        padding: 10px; 
        border-radius: 8px; 
        margin: 1rem 0;
        color: #000;">
        <b>Example:</b><br>{text}
    </div>
    """, unsafe_allow_html=True)

def example_cont_box(text):
    st.markdown(f"""
    <div style="border-left: 5px solid #17a2b8; 
        background-color: #e8f7fa; 
            padding: 8px; 
            border-radius: 8px; 
            margin: 1rem;
            color: #000">
    {text}
    </div>
    """, unsafe_allow_html=True)



def transparency_note(text):
    st.markdown(f"""
    <div style="
        border-left: 8px solid #ffc107; 
        background-color: #fff3cd; 
        padding: 12px; 
        border-radius: 8px; 
        margin: 1.5rem 0;
        color: #000;">
        <i>Transparency Note:</i><br>{text}
    </div>
    """, unsafe_allow_html=True)

def input_block(title):
    st.markdown(f"""
    <h4 style="
        color: #0056b3; 
        margin-top: 1.5rem; 
        margin-bottom: 0.5rem;">
        {title}
    </h4>
    
    """, unsafe_allow_html=True)


# ===== SRE sub-tab =====
def structural_rna_element_tab(): #(With the MSA Tool)
    
        
    st.markdown("""<div style="font-size: 2.7rem; color: #ea534e; font-weight: bold; margin-top: 1rem;">
        MSA Tool:  
    </div>
    """, unsafe_allow_html=True)


    st.markdown("""**Load MSA file**""")
    uploaded_file = st.file_uploader("Formats: FASTA (.fasta, .fa) or Clustal (.aln)", type=["fasta", "fa", "aln"])
    
    st.markdown("""
    This section helps characterize the structure of an RNA element.
    If the element occurs naturally, an evolutionary analysis can be applied.
    While this is a widely used approach, no existing tool allows researchers to perform the analysis and visualize the results directly on the RNA structure, 
    which can be highly insightful for identifying conserved (and thus functionally important) structural regions.
    """)
    st.markdown("---")
    msa= []
    # Processing upload file
    if uploaded_file:
        msa = parse_fasta_msa(uploaded_file) 
        msa_name = uploaded_file.name.split('.')[0] 
        if len(msa) < 1: 
            st.error("The MSA file must contain at least one sequence to run the analysis.")
        else:
            st.success(f"'{uploaded_file.name}' file uploaded correctly with {len(msa)} sequences.")
    else:
        st.info("Please, upload an MSA file to proceed with the analysis.")

    # SUB-SUB-TABS para "structural RNA element"
    structural_rna_sub_sub_tabs = st.tabs([
        "General information",
        "Structural Clustering Analysis",
        "Visualise All Structures"
    ])

    # Pass msa and msa_name to sub-tab functions
    with structural_rna_sub_sub_tabs[0]:
        if msa and len(msa) >= 1:
            # --- START analysis content ---
            import shutil, os
            from matplotlib import pyplot as plt

            st.header(" First Sequence Analysis & Conservation")
            sequence = msa[0].replace("-", "")
            structure, energy = predict_secundary_structure(sequence)
            st.code(f"{structure} ({energy:.2f} kcal/mol)")
            st.markdown("""
                This is the predicted structure and energy of the first sequence of your MSA. 
            """)

            conservation_scores = calculate_conservation(msa)
            msa_len = len(msa[0])
            col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
            colors = {
                "always": (255, 0, 0),
                "except2": (139, 0, 0),
                "except4": (0, 0, 255),
                "others": (255, 255, 0)
            }

            output_dir_first = f"rna_outputs_{msa_name}_general"
            shutil.rmtree(output_dir_first, ignore_errors=True)
            os.makedirs(output_dir_first)

            pre_command = generate_pre_command_per_sequence(msa[0], col_categories, colors, lw=10)
            fold_file = os.path.join(output_dir_first, "rna_structure_0.fold")
            with open(fold_file, "w") as f:
                f.write(sequence + "\n" + structure + "\n")

            rnaplot_cmd = f'RNAplot --pre "{pre_command}" < {fold_file}'
            if run_command(rnaplot_cmd):
                try:
                    ps_file = os.path.join(output_dir_first, "rna_structure_0.ps")
                    png_file = os.path.join(output_dir_first, "rna_structure_0.png")
                    os.rename("rna.ps", ps_file)
                    gs_cmd = (
                        f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
                        f'-dEPSCrop -sOutputFile={png_file} {ps_file}'
                    )
                    run_command(gs_cmd)
                    if os.path.exists(png_file):
                        st.markdown("<p style='text-align:center; font-style: italic;'>Structure Cluster 0</p>", unsafe_allow_html=True)
                        st.image(png_file, use_column_width=True)
                except FileNotFoundError:
                    st.warning("RNAplot did not generate the expected output.")

            st.subheader(" Conservation Map")
            fig, ax = plt.subplots(figsize=(10, 2))
            ax.plot(conservation_scores, color="#2c7a7b", linewidth=2)
            ax.set_title("Per-Base Conservation Score")
            ax.set_xlabel("Alignment Position")
            ax.set_ylabel("Score")
            ax.grid(True)
            ax.set_ylim(0, 1)
            st.pyplot(fig) 
            st.markdown("""
                This is the conservation map. This indicates you the exact level of conservation of each nucleotide. 
            """)





            structure = (predict_secundary_structure(msa[0].replace('-','')))[0]
            pairs = get_base_pairs(align_secondary_structure(msa[0], structure))
            col_categories = {idx: calculate_complementarity_conservation(msa, idx, pairs) for idx in range(len(msa[0]))}
            st.markdown("This visualization reflects how well base pair complementarity from the reference (first) sequence is maintained across the MSA. A high level"
             "of conservation implies structural stability and the potential for conserved RNA secondary structures " \
             "(If these pairings are frequently conserved, it suggests the RNA structure may also be preserved).")
            colors = {
                    "always": (255, 0, 0),
                    "notpaired": (255,255,255),
                    "except2": (0, 0, 255),
                    "except4": (255, 255, 0),                    
                    "others": (255, 255, 0)
                }
            fold_path = os.path.join(output_dir_first, "pairing_conservation.fold")
            with open(fold_path, "w") as f:
                f.write(sequence + "\n" + structure + "\n")
            pre_command = generate_pre_command_per_sequence(msa[0], col_categories, colors, lw=10)
            rnaplot_cmd = f'RNAplot --pre "{pre_command}" < {fold_path}'
            run_command(rnaplot_cmd)

            ps_path = os.path.join(output_dir_first, "pairing_conservation.ps")
            png_path = os.path.join(output_dir_first, "pairing_conservation.png")
            if os.path.exists("rna.ps"):
                os.rename("rna.ps", ps_path)
                run_command(f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
                            f'-dEPSCrop -sOutputFile={png_path} {ps_path}')
                if os.path.exists(png_path):
                    st.image(png_path, caption=f".")

            # legend
            st.markdown("### Legend:")
            st.markdown(
                """
                <div style='display: flex; flex-direction: column; gap: 8px;'>
                    <div><span style='background-color: rgb(255, 0, 0); 
                                     display: inline-block; width: 20px; height: 20px;
                                     margin-right: 10px; border: 1px solid black;'></span>
                        <b>always</b> — present in all sequences </div>
                    <div><span style='background-color: rgb(0, 0, 255); 
                                     display: inline-block; width: 20px; height: 20px;
                                     margin-right: 10px; border: 1px solid black;'></span>
                        <b>almost always conserved</b> — present in all sequences except for up to two</div>
                    <div><span style='background-color: rgb(255, 255, 0); 
                                     display: inline-block; width: 20px; height: 20px;
                                     margin-right: 10px; border: 1px solid black;'></span>
                        <b>others</b> — less consens</div>
                    <div><span style='background-color: rgb(255, 255, 255); 
                                         display: inline-block; width: 20px; height: 20px;
                                         margin-right: 10px; border: 1px solid black;'></span>
                            <b>not paired</b></div>
                </div>
                """,
                unsafe_allow_html=True
            )
                        
            # --- END analysis content ---
        else:
            st.info("Once the file is uploaded, here you will find the predicted structure and energy of the first seuqence of the MSA (presumably the one you are interesetd on characterising)"
            "and a conservation map, with the conservation score per position.")
    
    with structural_rna_sub_sub_tabs[1]: # structural clustering
        if msa and len(msa) >= 2: # Clustering requires at least 2 sequences
            # --- START clustering content ---
            import shutil, os, time
            st.header("Structural Clustering")

            if st.button("Run clustering and analysis"):
                output_dir = f"rna_outputs_{msa_name}_clustering"
                shutil.rmtree(output_dir, ignore_errors=True)
                os.makedirs(output_dir)

                structures = [predict_secundary_structure(seq.replace("-", ""))[0] for seq in msa] # msa es lista de strings
                start = time.time()
                rmsd_matrix = matriz_rmsd(structures)
                labels = cluster_structures(rmsd_matrix, eps=3.0, min_samples=1)
                cluster_summary = organise_per_clusters(msa, structures, labels, rmsd_matrix, output_base=output_dir) # msa es lista de strings
                end = time.time()

                times = [(end - start) / len(structures)] * len(structures)
                metrics = compute_metrics(structures, times)
                visualise_metrics(metrics,output_dir)

                st.subheader("📈 Evaluation Metrics")
                st.markdown(f"""
                    <div class="metric-box">
                        <b>Total Structures:</b> {metrics["total_structures"]}<br>
                        <b>Unique Structures:</b> {metrics["unique_structures"]}<br>
                        <b>Total Time:</b> {metrics["total_time"]:.3f} sec
                    </div>
                """, unsafe_allow_html=True)

                st.subheader("📊 Cluster Summary")
                os.makedirs(output_dir, exist_ok=True)

                # To use sequence, col_categories y colors, we compute locally:
                sequence = msa[0]
                msa_len = len(msa[0]) 
                col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
                colors = {
                    "always": (255, 0, 0),
                    "except2": (0, 0, 255),
                    "except4": (255, 255, 0),
                    "others": (255, 255, 0)
                }

                # Modifica la sección dentro de la pestaña structural clustering

                for cluster_id, cluster_data in cluster_summary.items():
                    with st.expander(f"Cluster {cluster_id} — {cluster_data['num_variants']} structures"):
                        rep_structure = cluster_data['representative_structure']
                        st.markdown(f"Representative Structure: `{rep_structure}`")
                        for var in cluster_data.get("members", []):
                            st.markdown(f"Specific sequence:  `{var['sequence']}`")
                        
                        # ------------------------------------------------------------------------
                        # ESTO ES LO NUEVO: Cargar la imagen existente en vez de generarla
                        
                        # Obtener el índice del representante buscando su estructura en la lista original de estructuras.
                        # Asegúrate de que 'structures' (la lista de todas las estructuras) esté en el ámbito de este bucle.
                        try:
                            # Obtener el índice del representante
                            rep_index = structures.index(rep_structure)
                            
                            # Definir la ruta a la imagen pre-generada en la carpeta 'all_structures'
                            all_structures_dir = f'rna_outputs_{msa_name}_all_structures'
                            # La imagen se llama 'cimages_{index}.png' en la pestaña 'all_structures'
                            rep_image_path = os.path.join(all_structures_dir, f"cimages_{rep_index}.png")
                
                            # Verificar si la imagen existe antes de mostrarla
                            if os.path.exists(rep_image_path):
                                st.image(rep_image_path, caption=f"Representative Structure Cluster {cluster_id}")
                            else:
                                st.warning(f"No se encontró la imagen para el representante del clúster {cluster_id} en {rep_image_path}.")
                        except ValueError:
                            st.warning(f"No se pudo encontrar el índice del representante del clúster {cluster_id}. Es posible que la lista de estructuras haya cambiado.")
                            
                        # ------------------------------------------------------------------------
                
                        # legend
                        st.markdown("### Legend:")
                        st.markdown(
                            """
                            <div style='display: flex; flex-direction: column; gap: 8px;'>
                                <div><span style='background-color: rgb(255, 0, 0); 
                                                    display: inline-block; width: 20px; height: 20px;
                                                    margin-right: 10px; border: 1px solid black;'></span>
                                    <b>always</b> — present in all sequences </div>
                                <div><span style='background-color: rgb(0, 0, 255); 
                                                    display: inline-block; width: 20px; height: 20px;
                                                    margin-right: 10px; border: 1px solid black;'></span>
                                    <b>almost always conserved</b> — present in all sequences except for up to two</div>
                                <div><span style='background-color: rgb(255, 255, 0); 
                                                    display: inline-block; width: 20px; height: 20px;
                                                    margin-right: 10px; border: 1px solid black;'></span>
                                    <b>others</b> — less consens</div>
                            </div>
                            """,
                            unsafe_allow_html=True
                        )      
                        

                if os.path.exists(f"{output_dir}/diversity.png"):
                    st.image(f"{output_dir}/diversity.png", caption="Structural Diversity")


            # --- END clustering content ---
        else:
            st.info("This section allows you to see the structure of the different unique structures of your MSA, with the nucleotides colored by conservation.")

    with structural_rna_sub_sub_tabs[2]: # Visualise all structures
        if msa and len(msa) >= 1:
            
            # --- START visualise all structures content  ---
            output_dir = f'rna_outputs_{msa_name}_all_structures'
            msa_len = len(msa[0]) 
            os.makedirs(output_dir, exist_ok=True)
            col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
            colors = {
                    "always": (255, 0, 0),
                    "except2": (0, 0, 255),
                    "except4": (255, 255, 0),
                    "others": (255, 255, 0)
                }
            
            for i, msa_seq in enumerate(msa): 
                with st.expander(f"### Sequence {i}"):
                    sequence = msa_seq.replace("-", "")
                    structure = predict_secundary_structure(sequence)[0]
                    st.markdown(f"Structure: `{structure}`")
                    st.markdown(f"Energy: `{predict_secundary_structure(sequence)[1]}`")
                    st.markdown(f"Sequence: `{sequence}`")
                    pre_command = generate_pre_command_per_sequence(msa_seq, col_categories, colors, lw=10)

                    fold_file = os.path.join(output_dir, f"images_{i}.fold")
                    with open(fold_file, "w") as f:
                        f.write(sequence + "\n" + structure + "\n")

                    run_command(f'RNAplot --pre "{pre_command}" < {fold_file}')

                    ps_file = os.path.join(output_dir, f"images_{i}.ps")
                    png_file = os.path.join(output_dir, f"cimages_{i}.png")

                    if os.path.exists("rna.ps"):
                        os.rename("rna.ps", ps_file)
                        run_command(
                            f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
                            f'-dEPSCrop -sOutputFile={png_file} {ps_file}'
                        )

                        if os.path.exists(png_file):
                            st.image(png_file, caption=f"Sequence {i}")
                        else:
                            st.warning(f"Could not generate image for sequence {i}.")
                    else:
                        st.warning(f"No output from RNAplot for sequence {i}.")
                     # legend
                    st.markdown("## Legend:")
                    st.markdown(
                        """
                        <div style='display: flex; flex-direction: column; gap: 8px;'>
                            <div><span style='background-color: rgb(255, 0, 0); 
                                             display: inline-block; width: 20px; height: 20px;
                                             margin-right: 10px; border: 1px solid black;'></span>
                                <b>always</b> — present in all sequences </div>
                            <div><span style='background-color: rgb(0, 0, 255); 
                                             display: inline-block; width: 20px; height: 20px;
                                             margin-right: 10px; border: 1px solid black;'></span>
                                <b>almost always conserved</b> — present in all sequences except for up to two</div>
                            <div><span style='background-color: rgb(255, 255, 0); 
                                             display: inline-block; width: 20px; height: 20px;
                                             margin-right: 10px; border: 1px solid black;'></span>
                                <b>others</b> — less consens</div>
                            
                        </div>
                        """,
                        unsafe_allow_html=True
                    )
                # --- END visualise all structures content ---
        else:
            st.info("Please, upload an MSA file to proceed with the analysis.")




# --- Content for "Linker Finder" sub-tab ---
def linker_finder_tab():
    st.markdown("""<div style="font-size: 2.7rem; color: #ea534e; font-weight: bold; margin-top: 1rem;">
        Linker Finder  
    </div>
    """, unsafe_allow_html=True)
    # Main container
    with st.container():
        # Inputs
        with st.expander("Sequences and targeted structure", expanded=True):
            rna1 = st.text_area(
                label="Structural RNA element",
                value='ACCAGUGUGCGGAUGAUAACUACUGACGAAAGAGUCAUCGACUCAGUUAGUGGUUGGAUGUAGUCACAUUAGU',
                height=110,
                placeholder="Introduce the SRE sequence  (A, C, G, U).",
                help="More infromation on 'About the Parts -> Structural RNA element'"
            )

            rna3 = st.text_area(
                label="Aptamer",
                value='AGUUGGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACCAAAA',
                height=110,
                placeholder="Introduce the aptamer sequence (A, C, G, U).",
                help="More infromation on 'About the Parts -> Aptamer'"
            )

            struct1 = st.text_input(
                label="Targeted structure of the SRE element (dot-bracket)",
                value='((.(((((((......(((((((((((....(((((...))))).)))))))))))......).)))))).))',
                placeholder="Example: ((.((...)))).",
                help="Put here the theoretical structure of the SRE element (determined experimentally or through " \
                "a prediction software, like RNAFold). This is the structure that the SRE sequence is going to be forced to have on the constrained (ON) state"
            )
            st.markdown('</div>', unsafe_allow_html=True)

        with st.expander("⚙️ Advanced Parameters", expanded=False):

            # --- Search Method Selection ---
            search_method = st.radio(
                "Select Search Method",
                ("Genetic Algorithm", "Brute Force"),
                help="More information on 'Help -> Documentation'"
            )

            constraint = st.text_input(
                label="Restriction chain for the Aptamer sequence",
                value='........................................................................................<<<<<<<...............................',
                placeholder="Example: ...........<<<<<............",
                help="Restrictions for the folding of the constrained (ON) state. It must have the same " \
                "length as the RNA2 sequence. '.' means no restriction, '(' and ')' for corresponding pairs, 'x' " \
                "to avoid pairing that nucleotide and '>' to avoid pairing with later nts and '<' to avoid pairing with previous nucleotides."
            )
            mutable_rna1 = st.text_input(
                label="Mutable positions for SRE",
                value="1, 2, 4, 5, 6, 7, 8, 9",
                placeholder="Ex: 1, 2, 4, 5, 6",
                help="This positions should be nucleotides known for NOT having an impact on functionallity. " \
                "Notice that, if paired on the input targeted structure, the pairs will be mutated aswell to mantain complementarity."
            )
            watched_positions_str = st.text_input(
                label="Watched positions for SRE",
                value="56",
                placeholder="Ex: 56, 57",
                help="Positions key to functionality for the SRE element whose pairing "
                "(from respect to the input targeted structure) will be forced to change to ensure the OFF state disrupts functionality."
            )
        
            use_mutations = st.checkbox("Allow mutations on SRE", value=False)
        
            # compute mutation rang
            mutable_set = [int(x.strip()) for x in mutable_rna1.split(",") if x.strip().isdigit()]
            max_mutations_allowed = len(mutable_set) if mutable_set else 0
        
            if use_mutations and max_mutations_allowed > 0:
                num_mut = st.slider(
                    label=f"Maximum number of mutations(1 - {max_mutations_allowed})",
                    min_value=1,
                    max_value=max_mutations_allowed,
                    value=1,
                    step=1,
                    help="Maximum number of mutations allowed."
                )
            else:
                num_mut = 0
        
            # Slider para longitud linker entre 4 y 10
            linker_min = st.slider(
                label="Linker's minimum length",
                min_value=4,
                max_value=10,
                value=7,
                step=1,
                help="More infromation on 'About the Parts -> System in Detail'. The linker is the key element that will be randomised untill one is found such that it makes the ON/OFF system possible"
            )
            linker_max = st.slider(
                label="Linker's maximum length",
                min_value=linker_min,
                max_value=10,
                value=7,
                step=1,
                help="More infromation on 'About the Parts -> System in Detail'. Linkers of length from the minimum length up to this length will be tested."
            )
        
            mfe_delta = st.number_input(
                label="Minimum Energy difference (kcal/mol)",
                min_value=0.0, max_value=20.0, value=4.0, step=0.1,
                help="More infromation on 'About the Parts -> System in Detail'. Difference between the energy of the constrained" \
                "(ON) state and the unconstrained (OFF) state. It should be approximately half " \
                "the binding energy of the aptamer–ligand interaction."
            )
            max_pairings = st.number_input(
                label="Maximum number of SRE-aptamer pairings",
                min_value=0, max_value=20, value=10, step=1,
                help="More infromation on 'About the Parts -> System in Detail'. Binding with the ligand can't brake many pairings. To ensure the switch's functionality, " \
                "one should choose a low number of pairings. However, if this conditions is very restrictive, " \
                "the ON/OFF system might not be possible to construct, and thus it is allowed to modulate this value."
            )
            max_changes = st.number_input(
                label="Maximum changes on the SRE structure",
                min_value=0, max_value=20, value=6, step=1,
                help="The input targeted structure might not need to be much restrictive, " \
                "as maybe not the whole structure is important for functionality. Therefore " \
                "the user can indicate how many changes are supported. This might allow more linkers to be found."
            )
            verbose = st.checkbox("Show execution details", value=True)

            # --- GA Specific Parameters (conditionally displayed) ---
            if search_method == "Genetic Algorithm":
                st.markdown("---")
                st.subheader("Genetic Algorithm Parameters")
                ga_population_size = st.number_input(
                    label="Population Size",
                    min_value=10, max_value=500, value=100, step=10,
                    help="Number of individuals in each generation."
                )
                ga_generations = st.number_input(
                    label="Number of Generations",
                    min_value=10, max_value=1000, value=20, step=10,
                    help="Number of generations to run the algorithm."
                )
                ga_elitism_rate = st.slider(
                    label="Elitism Rate",
                    min_value=0.0, max_value=0.5, value=0.1, step=0.01,
                    help="Fraction of the best individuals that are directly passed to the next generation."
                )
                ga_mutation_rate_rna1 = st.slider(
                    label="SRE Mutation Rate",
                    min_value=0.0, max_value=0.4, value=0.05, step=0.005,
                    help="Probability of a base mutation in SRE per position."
                )
                ga_mutation_rate_linker = st.slider(
                    label="Linker Mutation Rate",
                    min_value=0.0, max_value=0.5, value=0.1, step=0.01,
                    help="Probability of a base mutation in the linker per position."
                )
                ga_tournament_size = st.number_input(
                    label="Tournament Size",
                    min_value=2, max_value=10, value=3, step=1,
                    help="Number of individuals competing in each tournament selection."
                )
                # For GA, linker length must be fixed, so we'll use linker_min
                st.info(f"For Genetic Algorithm, only the 'Linker's minimum length' ({linker_min}) will be used as the fixed linker length.")
                linker_length_for_ga = linker_min # Use the min length as fixed length for GA

            st.markdown('</div>', unsafe_allow_html=True)
        st.markdown(
            """
            <p style='color:red; font-weight:bold;'>
            WARNING: If you do not change this name, the results will be mixed and some might be rewritten.
            </p>
            """,
            unsafe_allow_html=True
        )

        folder = st.text_area(
                label="Folder name",
                value='linker_search',
                help="Folder name in which the files will be placed. "
            )
        st.markdown("---")
        # Execute and STOP buttons
        run_col1, run_col2 = st.columns([4, 1])
        with run_col1:
            st.markdown("""
                ⚠️ **Notice – Responsible Use Agreement**

                By pressing this button, you confirm that you will use this tool **solely for research or educational purposes** and in compliance with all applicable **biosafety** and **biosecurity** regulations.  

                **Misuse of this software for harmful purposes is strictly prohibited.**
                """)

            run_button = st.button("Execute Linker search")
        with run_col2:
            st.markdown("ATENTION, images and .txt will be saved on the corresponding folder but the text bellow will restart.")
            if st.button("⛔ Stop Search"):
                st.session_state.stop_requested = True

        def check_stop():
            return st.session_state.get("stop_requested", False)


        # Área de log output
        output_container = st.empty()
        st.session_state.setdefault("log_text", "")
        st.session_state.setdefault("stop_requested", False)


        def update_log(msg):
            current_text = st.session_state.get("log_text", "")
            st.session_state["log_text"] = current_text + msg + "\n"
            output_container.code(st.session_state["log_text"], language="text")

        # Initial variables
        results = []
        report = ""
        labels = []
        if run_button:
            st.session_state.stop_requested = False  # Restart when execute
            st.session_state["log_text"] = ""        # Clean prev logs
            output_container.code("", language="text")  # Clean container

            if not rna1.strip() or not rna3.strip() or not struct1.strip():
                st.warning("Please, introduce SRE and aptamer sequences and targeted structure")
                st.stop()

            try:
                mutable_set = [int(x.strip()) for x in mutable_rna1.split(",") if x.strip().isdigit()]
                watched_positions = [int(x.strip()) for x in watched_positions_str.split(",") if x.strip().isdigit()]
                
                if search_method == "Brute Force":
                    linker_lengths = range(linker_min, linker_max + 1)
                    with st.spinner("Searching for linkers (Brute Force)..."):
                        results, report_lines, labels = run_linker_search(
                            rna1=rna1.strip(),
                            rna3=rna3.strip(),
                            struct1=struct1.strip(),
                            constraint=constraint.strip(),
                            mutable_rna1=mutable_set,
                            watched_positions=watched_positions,
                            output_dir=folder,
                            use_mutations=use_mutations,
                            mfe_delta=mfe_delta,
                            max_pairings_rna1_rna3=max_pairings,
                            max_structure_changes=max_changes,
                            num_mut=num_mut,
                            linker_lengths=linker_lengths,
                            verbose=verbose,
                            check_stop_func=check_stop,
                            log_func=update_log
                        )
                    report = report_lines # Join the report lines generated by the GA wrapper
                    if results: # Check if any valid solutions were found
                        # The `results` list is already populated in the desired format by run_genetic_algorithm_search
                        update_log(f"Brute force search completed. Found {len(results)} valid linker(s).")
                    else:
                        results = []
                        labels = []
                        update_log("Brute force search completed. No valid linker found. Please, consider changing some parameters")  
                elif search_method == "Genetic Algorithm":
                    # Call the genetic algorithm wrapper function
                    with st.spinner("Searching for linkers (Genetic Algorithm)..."):

                        # The wrapper function already handles internal logging if verbose=True
                        # and prepares the results, report, and cluster labels.
                        results, report_lines, labels = run_genetic_algorithm_search(
                            rna1=rna1.strip(),
                            rna3=rna3.strip(),
                            struct1=struct1.strip(),
                            constraint=constraint.strip(),
                            mutable_rna1=mutable_set,
                            watched_positions=watched_positions,
                            output_dir=folder,
                            use_mutations=use_mutations,
                            mfe_delta=mfe_delta,
                            max_pairings_rna1_rna3=max_pairings,
                            max_structure_changes=max_changes,
                            num_mut=0, # Not directly used by GA, but kept for compatibility of the wrapper
                            verbose=verbose,
                            check_stop_func=check_stop,
                            log_func=update_log,
                            # --- Pass the GA-specific parameters that are defined in app.py ---
                            population_size=ga_population_size,
                            generations=ga_generations,
                            linker_length_ga=linker_length_for_ga, # This maps app.py's linker_length_for_ga to ga's linker_length_ga
                            elitism_rate=ga_elitism_rate,
                            mutation_rate_rna1=ga_mutation_rate_rna1,
                            mutation_rate_linker=ga_mutation_rate_linker,
                            tournament_size=ga_tournament_size
                        )
                    report = report_lines # Join the report lines generated by the GA wrapper

                    if results: # Check if any valid solutions were found
                        # The `results` list is already populated in the desired format by run_genetic_algorithm_search
                        update_log(f"Genetic Algorithm search completed. Found {len(results)} valid linker(s).")
                    else:
                        results = []
                        labels = []
                        update_log("Genetic Algorithm search completed. No valid linker found. Please try again")                
                st.session_state["results"] = results
                st.session_state["report"] = report
                st.session_state["cluster_labels"] = labels

            except Exception as e:
                st.error(f"❌ Error during execution: {e}")
                # Log the error as well
                update_log(f"ERROR: {e}")
                st.code(traceback.format_exc()) # This displays the full traceback in Streamlit
                update_log(f"FULL TRACEBACK: {traceback.format_exc()}") # This logs it internally


            if not results:
                st.warning("No valid linkers found. If you are using the GA, please try again. If you are using Brute Force, consider changing the parameters. ")

        # Mostrar resultados si existen
        if "results" in st.session_state and st.session_state["results"]:
            results = st.session_state["results"]
            report = st.session_state["report"] # This is still the raw text log
            labels = st.session_state.get("cluster_labels", [])
            mutable_set = [int(x.strip()) for x in mutable_rna1.split(",") if x.strip().isdigit()]
            watched_positions = [int(x.strip()) for x in watched_positions_str.split(",") if x.strip().isdigit()]
            linker_lengths = range(linker_min, linker_max + 1) # This is mostly for report generation; actual linker length used depends on method

            # These plotting functions are designed for multiple results.
            # For GA, they will receive a list with one item.
            img_path = plot_delta_mfe(results, output_dir=folder)
            # Pass the correct RNA1 sequence (mutated for GA) to plot_pairings_histogram if it uses it.
            # For GA, rna1_mutated_seq will be in results[0]['rna1_mutated_seq']
            rna1_for_plotting = results[0].get('rna1_mutated_seq', rna1) if search_method == "Genetic Algorithm" and results else rna1
            img_pairings = plot_pairings_histogram(results, rna1_for_plotting, rna3, output_dir=folder)
            
            st.success(f"✅ {len(results)} valid linkers found.")

            # Clustering
            cluster_dict = defaultdict(list)
            for idx, label in enumerate(labels):
                cluster_dict[label].append((idx, results[idx]))
            representative_img_bases = []
            
            if (len(results) > 1):
                for cluster_id, cluster_items in sorted(cluster_dict.items()):
                    st.markdown(f"## Cluster {cluster_id} ({len(cluster_items)} sequences)")

                    # Choosing representative
                    representative = min(cluster_items, key=lambda x: abs(abs(x[1]['mfe_2']-x[1]['mfe_2']) - mfe_delta))

                    idx, res = representative
                    representative_img_bases.append(res['linker'])
                    with st.expander(f" Representative Linker (#{idx + 1}): {res['linker']}"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**SRE's mutations:** {res['mut1_info']}")

                        # RNAplot creates images in the folder, the file names are derived from linker
                        # We need to explicitly check for the generated files based on the 'best_individual'
                        # and the save_and_plot_structures function
                        linker_for_img = res['linker']
                        linker_len_for_img = len(linker_for_img) # Get the actual linker length


                        unconstrained_image_path = f"{folder}/linker_{len(linker_for_img)}/{linker_for_img}_unconstrained_plot.png"
                        constrained_image_path = f"{folder}/{linker_for_img}_constrained_plot.png"

                        if os.path.exists(unconstrained_image_path):
                            st.image(unconstrained_image_path, caption=f"Unconstrained Structure for Linker {linker_for_img}")
                        if os.path.exists(constrained_image_path):
                            st.image(constrained_image_path, caption=f"Constrained Structure for Linker {linker_for_img}")

                

            # Adjust parameters for HTML report based on search method
            current_rna1_for_report = rna1.strip()
            current_num_mut_for_report = num_mut
            current_linker_min_for_report = linker_min
            current_linker_max_for_report = linker_max
            
            if search_method == "Genetic Algorithm" and results:
                # Try to get the mutated RNA1 sequence from results.
                # If 'rna1_mutated_seq' is not present, it will fallback to the original 'rna1' sequence.
                current_rna1_for_report = results[0].get('rna1_mutated_seq', rna1)
                # For GA, num_mut is implicitly handled by GA's internal mutation logic,
                # but we can pass 0 or a placeholder as max mutations for RNA1 as per GA parameters
                current_num_mut_for_report = f"GA (Rate: {ga_mutation_rate_rna1})" 
                current_linker_min_for_report = linker_length_for_ga
                current_linker_max_for_report = linker_length_for_ga # Fixed length

            # Generate the rich HTML report
            html_report_content = build_full_html_report(
                results=st.session_state["results"],
                report=st.session_state["report"], # Passing the raw log to be included in the HTML report
                rna1=current_rna1_for_report,
                rna3=rna3,
                struct1=struct1,
                constraint=constraint,
                mutable_rna1=mutable_set,
                watched_positions=watched_positions,
                use_mutations=use_mutations,
                mfe_delta=mfe_delta,
                max_pairings=max_pairings,
                max_changes=max_changes,
                num_mut=current_num_mut_for_report,
                linker_min=current_linker_min_for_report,
                linker_max=current_linker_max_for_report,
                verbose=verbose,
                cluster_labels=st.session_state["cluster_labels"],
                search_method=search_method, # Now included
                image_output_dir=folder,
                representative_img_bases=representative_img_bases
            )
            # results in a zip file
            zip_data = create_zip_archive(folder)

            # allow download
            st.download_button(
                label="Download results",
                data=zip_data,
                file_name="results_rna_switch.zip",
                mime="application/zip"
            )

            # Download in SBOL
            st.markdown("""
                ### Download your design in SBOL format 

                **SBOL (Synthetic Biology Open Language)** is the global standard for describing synthetic biology designs. 
                        By downloading your results in this XML format, you can share your RNA design in a way that is easily 
                        reproducible, reusable, and compatible with other tools in the scientific community.

                """)
            merged_doc = sbol3.Document()

            if results:
                for index, result in enumerate(results):
                    # Llama a la función para cada resultado y añade el componente al documento
                    merged_doc = export_single_linker(
                        doc=merged_doc,
                        result=result,
                        index=index
                    )
                
                # Convierte el documento completo a una cadena JSON-LD
                sbol_bytes = merged_doc.write_string("json-ld").encode("utf-8")
            
                # Muestra el botón de descarga solo si hay resultados
                st.download_button(
                    label="Download All Linkers as SBOL",
                    data=sbol_bytes,
                    file_name="all_linkers.jsonld",
                    mime="application/ld+json"
                )
            else:
                st.info("No hay resultados de linkers para exportar.")



            # 4. Limpia la carpeta temporal después de la descarga (opcional, pero recomendado)
            # Esto asegura que no se acumulen archivos innecesarios en el servidor.
            #if os.path.exists(folder):
            #    shutil.rmtree(folder)
            
            

            st.subheader("Final Report")
            with st.expander("View Full HTML Report", expanded=False):
                # The 'height' and 'scrolling' parameters are important for long reports
                st.components.v1.html(html_report_content, height=800, scrolling=True)

            # 3. Add a download button for the HTML file
            st.download_button(
                label="Download Report as HTML",
                data=html_report_content,
                file_name=f"{folder}.html",
                mime="text/html"
            )

            show_all = st.checkbox("Show all linkers")
            folder = f"{folder}/linker_{len((results[0])['linker'])}"
            if show_all:
                for idx, res in enumerate(results):
                    with st.expander(f" Linker {idx + 1}: {res['linker']}"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Unconstrained structure:** `{res['structure_unconstrained']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**Mutations on SRE:** {res['mut1_info']}")
                        
                        linker_for_img = res['linker']
                        linker_len_for_img = len(linker_for_img)
                        
                        
                        unconstrained_image_path = f"{folder}/{linker_for_img}_unconstrained_plot.png"
                        constrained_image_path = f"{folder}/{linker_for_img}_constrained_plot.png"
                        print(constrained_image_path)
                        if os.path.exists(unconstrained_image_path):
                            st.image(unconstrained_image_path, caption=f"Unconstrained Structure for Linker {linker_for_img}")
                        if os.path.exists(constrained_image_path):
                            st.image(constrained_image_path, caption=f"Constrained Structure for Linker {linker_for_img}")



# --- Create diagrams for the description section ---
def create_ga_diagram():
    # Create new object
    dot = graphviz.Digraph(comment='Genetic Algorithm Workflow (Complete)')
    dot.attr('node', shape='box', style='rounded', fontname='Helvetica')
    dot.attr('edge', fontname='Helvetica')
    dot.attr('graph', rankdir='TB')

    # Main ndodes
    dot.node('start', 'Start')
    dot.node('init_pop', 'Initialize Population')

    with dot.subgraph(name='cluster_main_loop') as c:
        c.attr(label='Main Loop')
        c.attr(style='rounded')
        c.node('eval_fitness', 'Evaluate Fitness for each Individual')
        c.node('store_solutions', 'Store Valid Solutions\n(Check all criteria)')
        c.node('sel_parents', 'Select Parents (Elitism & Tournament)')
        c.node('gen_offspring', 'Crossover & Mutation\n(Generate Offspring)')
        c.node('next_gen', 'Form Next Generation')
    
    dot.node('max_gen_check', 'Maximum number of\ngenerations reached?')
    dot.node('cluster_solutions', 'Cluster Solutions\nby Structure Similarity')
    dot.node('generate_report', 'Generate Final Report\n& Output Files')
    dot.node('end', 'End')

    # Connexions
    dot.edge('start', 'init_pop')
    dot.edge('init_pop', 'eval_fitness')
    
    # Main loop
    dot.edge('eval_fitness', 'store_solutions')
    dot.edge('store_solutions', 'sel_parents')
    dot.edge('sel_parents', 'gen_offspring')
    dot.edge('gen_offspring', 'next_gen')
    dot.edge('next_gen', 'eval_fitness', label='Loop')
    
    # Final part of the loop
    dot.edge('eval_fitness', 'max_gen_check')
    dot.edge('max_gen_check', 'cluster_solutions', label='Yes')
    dot.edge('max_gen_check', 'sel_parents', label='No')

    # End process
    dot.edge('cluster_solutions', 'generate_report')
    dot.edge('generate_report', 'end')

    return dot


def create_brute_force_diagram():
    dot = graphviz.Digraph(comment='Brute Force Workflow')
    dot.attr('node', shape='box')

    dot.node('A', 'Start')
    dot.node('B', 'Iterate through all Linker Lengths')
    dot.node('C', 'Generate all possible Linker sequences')
    dot.node('D', 'Select a Linker')

    with dot.subgraph(name='cluster_1') as c:
        c.attr(label='Evaluation Loop')
        c.attr(style='rounded')
        
        c.node('E', 'Evaluate OFF-State (Unconstrained)')
        c.node('F', 'Evaluate ON-State (Constrained)')
        c.node('G', 'Calculate MFE Delta')
        
        c.edge('E', 'F')
        c.edge('F', 'G')

    dot.node('H', 'Meets all Criteria?')
    dot.node('I', 'Store as Valid Solution')
    dot.node('J', 'Mutations Enabled?')
    dot.node('K', 'Try all possible SRE Mutations')
    dot.node('L', 'Mfe Delta met with Mutation?')
    dot.node('M', 'Store as Valid Solution (with Mutations)')

    dot.node('N', 'All Linkers Checked?')
    dot.node('O', 'Cluster Solutions & Generate Report')
    dot.node('P', 'End')
    
    # Connexions
    dot.edge('A', 'B')
    dot.edge('B', 'C')
    dot.edge('C', 'D')
    dot.edge('D', 'E')
    dot.edge('G', 'H')
    
    dot.edge('H', 'I', label='Yes')
    dot.edge('H', 'J', label='No')

    dot.edge('J', 'K', label='Yes')
    dot.edge('J', 'D', label='No') # Iterate again if there where not mutations

    dot.edge('K', 'L')
    dot.edge('L', 'M', label='Yes')
    dot.edge('L', 'D', label='No') # Back to linker selection
    
    dot.edge('M', 'D')
    dot.edge('I', 'D') # Load and back to linker selection
    
    dot.edge('D', 'N')
    dot.edge('N', 'O', label='Yes')
    dot.edge('O', 'P')
    
    dot.edge('N', 'D', label='No')

    return dot


def create_simple_workflow_diagram():
    dot = graphviz.Digraph(comment='RNA Switch Designer Workflow')
    dot.attr('node', shape='box', style='rounded', fontname='Helvetica', fontsize='12')
    dot.attr('edge', fontname='Helvetica', fontsize='10')

    # Nodes
    dot.node('A', 'User Input\n(SRE, Aptamer, Constraints)')
    dot.node('B', 'Search Algorithm\n(Brute Force or GA)')
    dot.node('C', 'Evaluation\n(OFF/ON State, MFE Delta)')
    dot.node('D', 'Valid Linker Found?')
    dot.node('E', 'Final Report & Output\n(Linkers, Structures, Energies)')
    dot.node('F', 'Download Results')

    # Connexions
    dot.edge('A', 'B', label='User provides input')
    dot.edge('B', 'C', label='Software runs search')
    dot.edge('C', 'D', label='Evaluates each solution')
    dot.edge('D', 'C', label='No', arrowhead='none', style='dashed')
    dot.edge('D', 'E', label='Yes')
    dot.edge('E', 'F', label='User can download')
    
    return dot


def help_tab():
    st.header("Help and Documentation")
    st.write("Find user guides and tutorials")

    help_sub_tabs = st.tabs([" User Guide", "Why build an RNA-switch"," Documentation"])

    # === USER GUIDE TAB ===
    with help_sub_tabs[0]:
        st.markdown("""
        <h1 style="text-align:center; margin-bottom: 0.5em;">Introduction to the RNA ON/OFF Switch Design Toolkit</h1>
        

       
        
        <h2 style="margin-top: 1em; color:#cf211b;">Main Feature: Model System</h2>
        <h4 style="margin-top: 1em; color:#f7bcbb;">In short, you input your Structural RNA Element (SRE) and Aptamer, and the tool helps
                     you turn them into a functional RNA switch. This means you can turn the function of your SRE on and off:</h4>
     
        

        In this image, you can see the system that builds this tool, a switch with ON (right) and OFF (left) states. 
        """, unsafe_allow_html=True)
        # Centering image with collumns
        col1, col2, col3 = st.columns([1, 4, 1])
        with col2:
            st.image("images/system_img.png", width=900)
        st.markdown("""
        **OFF state:**  
           The aptamer and the SRE are paired. This means the structure of the SRE is disrupted, and it will not perform its function; therefore, this is the OFF state.
                    
        **On State:**  
           The aptamer is bound to the ligand and forms an independent structure, releasing the SRE. The SRE can then form the proper structure to perform its function, therefore leading to an ON state.

    
        ### Core Methodology
        The design process is a multi-step pipeline built on computational and biological principles:
        - **Constraint-Based Design:** You provide a desired 'ON' state structure as a constraint to guide the search algorithms.
        - **Thermodynamic Prediction:** The tool predicts the most stable secondary structures for your designs. The key evaluation metric is **ΔMFE**, the energy difference between the 'ON' and 'OFF' states, which indicates the switching efficiency.
        - **Dual Search Algorithms:** The tool uses both **Brute-Force Search** (for short linkers) and a **Genetic Algorithm** (for more complex designs) to find the best possible linker sequences.
        - **Evaluation:** Each candidate design is rigorously assessed based on its **ΔMFE**, structural preservation of the SRE, and minimization of unwanted interactions.
        - **Structural Analysis:** The final designs are clustered to help you identify unique and structurally robust solutions.
           
        Below you can find a diagram of the general workflow of the software. 

        """, unsafe_allow_html=True)
        col1, col2, col3 = st.columns([2.5, 4, 1])
        with col2:
            st.graphviz_chart(create_simple_workflow_diagram())
       
        st.markdown(
            "<h2 style='font-size:22px; margin:0 0 6px 0; color:#fdeeee;'>"
            "What do we refer to as a Structural RNA Element (SRE)?</h2>",
            unsafe_allow_html=True
        )
        with st.expander(" See details", expanded=False):


            st.markdown("""
            **Definition**  
            A Structural RNA Element (SRE) is an RNA sequence region that folds into a specific secondary structure crucial for its biological function.  
            These elements often rely on their folded conformation to regulate processes such as Stop Codon Readthrough or frameshifting.
            """)

            # SECIS Example
            example_box("""
            The SECIS element is an RNA structure that stimulates stop codon readthrough.  
            Normally, translation terminates at a stop codon, releasing the protein.  
            However, the SECIS element can cause the ribosome to skip the stop codon, resulting in a longer protein product.  
            Its function relies entirely on its secondary structure, so by controlling this structure, we can modulate its regulatory effect.
            """)

           
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/SCR.png", width=900)

            st.subheader("Secondary Structure Considerations")
            st.markdown("""
            The SRE must have a known or predictable secondary structure, typically represented in dot-bracket notation.  
            If an experimentally validated structure is not available, users are encouraged to predict one using established tools—**RNAfold** from the Vienna RNA Package is highly recommended.  
            """)

           
            st.markdown("""
            <div style="font-size: 1.1em; font-weight: bold; margin-top: 1rem; color: #e41513">
            In case the SRE’s structure you aim to study is not well characterised, this software includes an evolutionary analysis using a multiple sequence alignment (MSA Tool below) to help characterise conserved structural features.
            </div>""", unsafe_allow_html=True)

            transparency_note("""
            Note that RNAfold has certain limitations:<br>
            - It predicts only pseudoknot-free secondary structures.<br>
            - Complex tertiary interactions and non-canonical base pairs are not modelled.<br>
            - Predictions are minimum free energy (MFE) approximations and may not capture all biologically relevant conformations.<br>
            Despite these limitations, RNAfold remains one of the most comprehensive and accessible tools available for local execution.
            """)


            example_box("Stop Codon Readthrough (SCR)  \nWithout SCR, translation usually ends when the ribosome encounters a Stop Codon.")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/NoSCR.png", width=900)

            example_cont_box("However, when an SCR element is present, its structure interacts with the ribosome and makes it skip the Stop Codon, resulting in an elongated protein.")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/SCR.png", width=900)

            example_cont_box("To build a switch, we add an aptamer." \
            "- In the absence of ligand: the aptamer pairs with the SCR element, disrupting its structure. The ribosome stops at the stop codon, producing only protein A." \
            "")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/OFF.png", width=1200)

            example_cont_box(""
            "- In the presence of ligand: the aptamer binds to it and forms its own structure, releasing the SCR to fold normally, stimulating readthrough and producing protein AB." \
            "")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/ON.png", width=1200)

        
        
        
        st.markdown(
            "<h2 style='font-size:22px; margin:0 0 6px 0; color:#fdeeee;'>"
            "What is an Aptamer?</h2>",
            unsafe_allow_html=True
        )
        with st.expander(" See details", expanded=False):
            st.write("""
            An aptamer is a short, structured RNA (or DNA) sequence capable of binding specifically to a target molecule, usually a small ligand. Aptamers are often used as biosensors or regulatory elements in synthetic biology, particularly in riboswitches (RNA systems that regulate gene expression in response to ligand binding).
            In these systems, the aptamer forms part of a larger RNA structure that adopts one conformation in the absence of the ligand and switches to a different conformation when the ligand is present.
            """)

            # create 3 columns to make image central
            col1, col2, col3 = st.columns([2, 4, 1])

            with col2:
                st.image("images/aptamer.png", width=400)


            st.write("""
            Aptamers function through structure-based recognition, meaning their ability to bind the ligand depends on folding into a specific 3D shape. Therefore, preserving or engineering the aptamer’s functional conformation is critical in any synthetic RNA system.
            """)

           
            

            st.subheader("Aptamer Used in This Tool")
            st.write("""
            By default, the tool includes support for the theophylline aptamer, one of the most well-characterised synthetic RNA aptamers in the literature.
            It binds the ligand with high specificity—over 10,000-fold more selective than for caffeine, a closely related molecule.
            This aptamer has been widely studied and used in various synthetic biology applications, including riboswitches in both prokaryotic and eukaryotic systems.
            For example, the aptamer variant used by Anzalone et al. (2016) was shown to effectively regulate frameshifting in mammalian cells, making it a strong candidate for this system.
            The structure of the theophylline aptamer is well defined, and its binding pocket and surrounding regions are known from both experimental data (e.g., NMR, crystallography) and computational predictions (RNAfold, etc.).
            """)
            # create 3 columns to make image central
            col1, col2, col3 = st.columns([1, 4, 1])

            with col2:
                st.image("images/aptamer3D.png", width=750)

            st.subheader("Other Compatible Aptamers")
            st.write("""
            While the tool uses the theophylline aptamer by default, it is also compatible with other ligand-binding aptamers, such as:
            - Tetracycline aptamer
            - Warfarin aptamer
            - Ciprofloxacin aptamer
            """)

            st.info("""
            When using a different aptamer, be sure to:
            - Provide the correct aptamer sequence.
            - Define the “Restriction chain for the aptamer sequence” based on its ligand-binding region.
            - Change the Minimum ΔMFE (explanation bellow)
            """)
            st.subheader("Simulating Ligand Binding in Aptamers")
            st.write("""
            Simulating ligand binding directly is a highly complex task, often requiring techniques like molecular dynamics or docking simulations.
            These are not practical when screening hundreds or thousands of sequence variants in a high-throughput or exploratory context.
            To address this, the tool implements a simplified yet effective strategy for simulating ligand binding using secondary structure constraints.

            **Structural Restriction Approach**  
            Instead of modelling ligand–aptamer interactions directly, this parameter enforces structural constraints to mimic ligand binding.
            Specifically, the nucleotides known to interact with the ligand are restricted from forming base pairs with upstream regions of the RNA.
            This forces the aptamer to fold into its ligand-bound conformation, reproducing the structural effect of binding without needing to simulate the ligand itself.
            This approach has proven sufficient for many practical applications, including recreating known aptamer conformations upon ligand binding.
            
            ### Minimum ΔMFE (ON-OFF Energy Difference)
            In short, you need to know the energy of the binding with the ligand, and put the Minimum ΔMFE as the half of that energy. Here is why:
            <span style="color:#1f77b4; font-weight:bold;">What is the Minimum ΔMFE?:</span>  
            The minimum required difference in Minimum Free Energy (MFE) between the OFF and ON states. MFE is a proxy for how stable a structure is.

            <span style="color:#ff7f0e; font-weight:bold;">Why it matters:</span>  
            You want the system to have two distinct and energetically stable conformations:  
            - One in the absence of ligand (OFF), where the SRE is disrupted  
            - One in the presence of ligand (ON), where the SRE folds correctly  

            If both states have similar energies, the switch may not behave reliably, flipping randomly or being stuck in one state. By setting a minimum ΔMFE, you enforce a thermodynamic bias that makes switching behavior more robust and predictable.

            <span style="color:#2ca02c; font-weight:bold;">Goal:</span>  
            Ensure there's a meaningful, energetic difference between ON and OFF that reflects real structural change. This enhances the functional specificity of the switch.

            In this image, you can see the energy of the different configurations.  
            """)
            
            
            st.image("images/energies.png", width=1200)

            st.markdown("""
            The key concept is that the binding with the ligand changes the energy of the system (in the case of the theophylline aptamer, the binding energy is -9.5 kcal/mol). Therefore, the structure of the OFF state can have a different energies depending on whether the ligand is present or not.

            Without the ligand, we want to have the OFF state, so we want the OFF state without the ligand (central figure) to be the most stable one (and therefore have a lower energy than the ON state with the ligand, left figure). With the ligand, we want to have the ON state, so we want the ON state (left figure) to be more stable than the OFF state with the ligand (right figure) (and therefore have a lower energy).

                     """)
            

        st.markdown(
            "<h2 style='margin-top: 1em; color:#cf211b;'>Input Requirements</h2>",
            unsafe_allow_html=True
        ) 
        st.markdown("""
            
            #### General Key Parameters:
            """)
        with st.expander(" See details", expanded=False):
            st.markdown("""
                ### 1. Linker Length

                <span style="color:#1f77b4; font-weight:bold;">What it is:</span>  
                The number of nucleotides between the end of the SRE and the start of the aptamer.

                <span style="color:#ff7f0e; font-weight:bold;">Why it matters:</span>  
                The linker physically connects the aptamer to the regulatory element. Its length influences the folding potential of the full RNA molecule in both ON and OFF states:  
                - Too short, and it may prevent the aptamer and SRE from interacting sufficiently to form the OFF state.  
                - Too long, and it may introduce too much flexibility or unwanted structural configurations, reducing the precision of switching.

                <span style="color:#2ca02c; font-weight:bold;">Goal:</span>  
                Find a length that allows:  
                - Effective interaction between aptamer and SRE in the OFF state (leading to SRE disruption),  
                - But enables SRE recovery when the aptamer is ligand-bound (ON state).  

                The tool allows you to scan across a range of linker lengths to identify the most effective ones for switching behaviour.


                ### 2. Maximum Number of SRE–Aptamers Pairings

                <span style="color:#1f77b4; font-weight:bold;">What it is:</span>  
                A constraint on how many base pairs are allowed to form between the SRE and the aptamer in the OFF state.

                <span style="color:#ff7f0e; font-weight:bold;">Why it matters:</span>  
                In the OFF state, the aptamer must partially bind the SRE to prevent it from forming its active structure. However:  
                - Too few pairings might not be enough to effectively disrupt the SRE.  
                - Too many pairings could over-stabilise the OFF state, making the ON state hard to reach even when the ligand is present.

                <span style="color:#2ca02c; font-weight:bold;">Goal:</span>  
                Strike a balance: allow enough pairing for functional repression, but leave room for switching. This variable helps tune the strength of repression and influences whether ligand binding can successfully revert the system to the ON state.  
                This parameter is especially useful when testing different SRE or aptamer sequences, as their pairing potential may vary significantly.
            """, unsafe_allow_html=True)

            # Mutations Block
            st.markdown("""
                 ## SRE's key parameters:
                        """, unsafe_allow_html=True)
            input_block("Mutations") 
            st.markdown("""
                This software allows the user to introduce mutations in the SRE sequence to explore functional variants.  
                If a nucleotide involved in a base pair is mutated, its paired base is automatically mutated to a complementary one, preserving the secondary structure.""")
            # SECIS example
            example_box("""
                    Take the SECIS element as our SRE. The SECIS element has a certain structure (on the image bellow) that allows it to perform its function.
                        Without that certain structure, it does not work. Not all parts of the structure are equally important, though. Certain positions (green) are key for functions.  
                    As the bottom positions are not important for function, they can be mutated without disrupting activity.
                    """)
            col1, col2, col3 = st.columns([4, 4, 1])
            with col2:
                st.image("images/SECIS.png", width=100)
            st.markdown("""  
                It is recommended to limit mutations to non-critical nucleotides to avoid compromising function.
            """)

        # Key nucleotides Block
        st.markdown("""
     
            #### Key Nucleotides for Function:
            """)
        with st.expander(" See details", expanded=False):
            st.markdown("""
                The software allows specifying “watched positions”: nucleotides known to be critical for function.  
            """)
            example_box("""
                In the case of the SECIS element, the green positions are key for for function.  
            """)


            st.markdown("""
                Watched positions are used as functional markers:  
                - If their structure relative to the input structure is maintained, the SRE is considered functional.  
                - If disrupted, the SRE is considered non-functional.

                <b>Hint:</b> Avoid including all key nucleotides in ‘Watched positions’.  
                Even if one maintains pairing, the overall structure could still change.  
                Experiment with different sets to achieve accurate results.
                """, unsafe_allow_html=True)

            # Allow changes
            input_block("Allow Changes on SRE’s Structure") 
            st.markdown("""
                While certain nucleotides and structural features of an SRE are essential, others may be flexible or non-critical.  
                This tool includes a parameter: **Maximum changes on the SRE structure**, defining tolerated deviation from the original structure.

                - A strict value (0) means the structure must remain identical.  
                - A higher value allows more flexibility.
                """)

            # Final Block
            st.markdown("""
                <h3>Why These Inputs Matter</h3>
                RNA-based systems are complex; small structural changes can have a big impact, but not all regions matter equally.  
                This tool lets you define what matters most, preserving function where it counts while exploring a broad design space.
                """, unsafe_allow_html=True)


        st.markdown("""
     
            #### Search methods:
            """)
        with st.expander(" See details", expanded=False):
            st.markdown("""
                           
        
            #### Brute Force
            This method **tries every single possible combination** to see which ones work.  
            It is simple and guarantees you’ll find *all* possible solutions, but if there are too many combinations, it can take an extremely long time.

            **⚠️ Important:**  
            If you allow mutations in Brute Force, the number of possibilities grows enormously, and the search can become practically impossible to finish.  
            If you need mutations, or the number of possibilities is high (for a long linker for example) **use the Genetic Algorithm instead**.

            ---

            #### Genetic Algorithm (GA)
            A Genetic Algorithm works a bit like **natural selection in biology**:
            1. **Start** with a group of random possible solutions (called *population*).
            2. **Test** each one to see how well it works (*fitness*).
            3. **Keep** the best ones.
            4. **Mix** them together (reproduction) and make small random changes (*mutations*) to create a new generation.
            5. **Repeat** this process until a good solution is found.

            It doesn’t try *every* possibility, but it’s much faster to find solutions when the number of possibilities is huge.

            ---

            ###  Genetic Algorithm Parameters

            **Population Size – How many different candidates are tested at the same time.**  
            - 🔼 **Higher** → Explores more possibilities at once, increasing chances of finding good solutions but each generation takes longer to compute.  
            - 🔽 **Lower** → Faster per generation, but less diversity and higher risk of getting stuck in mediocre solutions.  

            **Number of Generations – How many times the process repeats.**  
            - 🔼 **Higher** → More time for the algorithm to refine solutions, but longer total runtime.  
            - 🔽 **Lower** → Quicker results, but solutions may be less optimised.  

            **Elitism Rate – Percentage of the best candidates kept unchanged for the next round.**  
            - 🔼 **Higher** → Keeps more top solutions, ensuring quality, but reduces exploration of new possibilities.  
            - 🔽 **Lower** → More exploration, but may lose very good solutions along the way.  

            **SRE Mutation Rate – Chance of changing each base in SRE during evolution.**  
            - 🔼 **Higher** → More exploration and diversity in SRE sequences, but can destroy promising structures.  
            - 🔽 **Lower** → More stability in SRE, but less chance to discover unexpected improvements.  

            **Linker Mutation Rate – Chance of changing each base in the linker.**  
            - 🔼 **Higher** → More variation in linker sequences and potential new folding patterns, but also more instability.  
            - 🔽 **Lower** → More structural stability, but fewer creative linker solutions.  

            **Tournament Size – How many candidates compete at once before picking the winner.**  
            - 🔼 **Higher** → Picks stronger candidates more often, speeding convergence but risking loss of diversity.  
            - 🔽 **Lower** → Keeps more diversity, but takes longer to reach high-quality solutions.  

            **Note:** For GA, the linker length is **fixed** and comes from the "Linker's minimum length" setting.
            """)
        
        st.markdown("""
            Below you can see a diagram of the two possible Search Methods (see explanation on the expander above):""")
        st.markdown("""
            Genetic Algorithm Workflow""")
        col1, col2, col3 = st.columns([1, 4, 0.5])
        with col2:
            st.graphviz_chart(create_ga_diagram())
        st.markdown("""
            Brute Force implementation""")
        col1, col2, col3 = st.columns([0.5, 4, 0.5])
        with col2:
            st.graphviz_chart(create_brute_force_diagram())

    with help_sub_tabs[1]:
        
        st.markdown("""
        

        <h2>Unlocking RNA control beyond the old playbook</h2>
        <p>
        Synthetic biology thrives on the ability to reprogram living systems with speed, precision, and predictability.<br>
        For years, gene regulation has relied on tools like transcription factor switches, CRISPR interference, RNAi, and inducible promoters, each powerful, but each with trade-offs: transcriptional delays, off-target effects, bulky components, or hard-to-deliver inducers.<br>
        What if you could bypass these limitations?<br>
        Acting at the translational level offers faster, more direct control over protein output, and RNA-based devices make this possible in an extremely compact form. Aptamers, short RNA sequences that bind specific small molecules, are already a staple in this space. But until now, they’ve been almost exclusively paired with RBS or IRES elements, limiting the range of control they could achieve.
        </p>

        <h3>A new design space: Structural RNA Elements</h3>
        <p>
        The real leap forward comes from coupling aptamers to SREs: Some RNA sequences carry out their function purely through structure, forming shapes that control ribosome behaviour, alter reading frames, or enable readthrough. These SREs are like “mechanical parts” in the RNA world, but traditionally they’ve been fixed in one mode by their sequence.<br>
        Our toolkit changes that.<br>
        We make any SRE conditional by coupling it to a ligand-responsive aptamer. The aptamer acts as a switch: in the absence of the ligand, it forces the SRE into an inactive fold; in the presence of the ligand, it releases the SRE to fold into its active form.<br>
        This turns static RNA motifs into dynamic, programmable ON/OFF switches.
        </p>

        
        
        """ , unsafe_allow_html=True)
        
        # Before / After 
        col_before, col_after = st.columns(2)

        with col_before:
            st.markdown("""
            <h3 style="text-align:center; color:#d9534f;">Before this toolkit</h3>
            <p style="text-align: justify;">
            Designing aptamer–SRE fusions meant exploring hundreds to thousands of linker variants, aptamer orientations, and sequence tweaks — each requiring structural prediction and ON/OFF scoring. This level of exhaustive modelling was beyond what most labs could realistically do by hand, making the approach impractical in day-to-day research.<br><br>
            It also required specialist expertise in RNA folding, free-energy modelling, and aptamer integration — time spent on design mechanics instead of developing applications.
            </p>
            """, unsafe_allow_html=True)

        with col_after:
            st.markdown("""
            <h3 style="text-align:center; color:#5cb85c;">After this toolkit</h3>
            <p style="text-align: justify;">
            Any team can design custom ligand-responsive RNA switches in minutes, without RNA expertise.<br><br>
            It works with any aptamer and any structural RNA element, from frameshift control to Stop Codon Readthrough and beyond.<br><br>
            The platform automates modelling, scoring, and visualisation, turning what was once a niche experimental challenge into a standard, accessible part of the synthetic biology toolbox.
            </p>
            """, unsafe_allow_html=True)

        st.markdown("---")

    with help_sub_tabs[2]:  # Documentation
        st.header("Modules and Documentation")
        st.markdown("""
           

            ### The TADPOLE Modules: What They Do

            While you won't need to interact with these directly, this is a high-level overview of the main components that make the tool run:

            * **`conservation.py`**: Analyzes sequence conservation to help ensure the stability and functionality of your designs.
            * **`genetic_algorithm.py`**: This is the "smart" engine for longer linker designs. It intelligently searches for optimal sequences.
            * **`input_utils.py`**: Ensures that all the files and sequences you upload are correctly formatted for the tool to use.
            * **`io_tools.py`**: Manages all the input and output, including saving your results as reports and images.
            * **`rna_cluster.py`**: Groups your final designs by structural similarity to help you identify unique solutions.
            * **`rna_mutation.py`**: Manages the sequence changes needed to explore new designs.
            * **`rna_structures.py`**: Provides the foundational logic for working with RNA structures and base pairs.
            * **`search.py`**: Runs the simple search for short linkers and evaluates how well your designs meet the criteria.
            * **`structure.py`**: Predicts the secondary shape of your RNA designs, which is crucial for their function.
            * **`visualisation.py`**: Generates all the helpful plots and images you receive in your final output folder.

            ---

            ### Source Code and Contribution

            TADPOLE is a project developed by Team Barcelona-UB 2025. You can view, fork, download the application and run in local, and contribute to the source code on our GitLab repository:

            * **GitLab Repository**: [https://gitlab.igem.org/2025/software-tools/barcelona-ub/](https://gitlab.igem.org/2025/software-tools/barcelona-ub/)

            ---

            ### Authors and Acknowledgements

            This software was developed by **Team Barcelona-UB 2025**.

            We extend special thanks to the developers of ViennaRNA, SBOL, and the open-source community for providing the foundational tools and libraries that made this project possible.
        
        """,
            unsafe_allow_html=True
        )
# ---------------- MAIN APPLICATION LAYOUT ----------------

#Header
html_content = """
<style>
        /* Primary Color Palette based on #e52920 */
        :root {
            --bg-light-red: #fdeeee;           /* Very light red/pink for backgrounds */
            --soft-red: #fbdad9;               /* Light red/pink for inactive elements, light card borders */
            --medium-light-red: #f7bcbb;       /* Slightly darker pink */
            --light-accent-red: #f28f8c;       /* Medium pink/light red for hover states */
            --primary-red: #ea534e;            /* Adjusted to #ea534e as per new palette, was #e52920 */
            --dark-accent-red: #e62e25;        /* Slightly darker primary red */
            --dark-red: #e41513;               /* Darker red for button hover/active */
            --very-dark-red: #be1818;          /* Even darker red */
            --deep-red: #9d1915;               /* Dark red/maroon for strong accents/borders */
            --text-darker-red: #53110e;        /* Very dark red for text/headers */
            
            --soft-gray: #e0e0e0;              /* Soft gray for general backgrounds/borders (kept as is) */
            --medium-gray: #cccccc;            /* Medium gray for borders/dividers (kept as is) */
            --dark-gray: #1a1a1a;              /* General dark text/header color (kept as is) */
            --text-color-light: #333333;       /* Dark text on light backgrounds (kept as is) */
            --text-color-dark: #ffffff;        /* White text on dark primary colors (kept as is) */
            --card-bg: #ffffff;                /* White background for cards/panels (kept as is) */
            --border-color: var(--soft-gray);  /* Default border color (kept as is) */
            --shadow-color: rgba(0, 0, 0, 0.1); /* Subtle shadow (kept as is) */
        }

        /* Headers - Ensure good contrast */
        h1, h2, h3, h4, h5, h6 {
            color: var(--text-darker-red); /* Using the darkest red for headers */
            font-weight: 600;
        }
        h1 { 
            font-size: 2.5em; 
            color: var(--primary-red); /* Keep H1 primary red */
            margin-bottom: 0.5em; 
        } 
        h2 { font-size: 2em; margin-bottom: 0.5em; }
        h3 { font-size: 1.5em; margin-bottom: 0.5em; }


        /* Existing styles adapted to new variables */
        body {
            font-family: 'Poppins', sans-serif;
            margin: 0;
        }
        
        .main-header {
            padding: 2rem 1rem;
            background-color: var(--primary-red);
            border-radius: 2rem;
            text-align: center;
        }
        .app-title {
            color: var(--medium-light-red);
            margin-bottom: 0.5rem;
        }
        .app-description, .secondary-line {
            font-size: 1.4rem;
            color: #000000;
            margin: 0.3rem 0;
        }
        .rotating-quote {
            font-size: 1.3rem;
            font-style: italic;
            margin-top: 1.2rem;
            color: var(--text-darker-red);
            min-height: 2rem;
            transition: opacity 0.3s ease-out, transform 0.3s ease-out;
            opacity: 1; 
            transform: translateY(0); 
        }

        /* --- RESPONSIVE STYLES --- */
        @media (max-width: 768px) {
            h1 {
                font-size: 2em;
            }
            .main-header {
                padding: 1.5rem 0.5rem;
            }
            .app-description, .secondary-line {
                font-size: 1.2rem;
            }
            .rotating-quote {
                font-size: 1rem;
                min-height: 1.5rem;
            }
        }
        
        @media (max-width: 480px) {
            h1 {
                font-size: 1.5em;
            }
            .main-header {
                border-radius: 1rem;
            }
            .app-description, .secondary-line {
                font-size: 1rem;
            }
        }
    </style>
</head>
<body>
    <div class="main-header">
        <div class="header-content">
            <h1 class="app-title">TADPOLE: Build your own RNA switch.</h1>
            <p class="app-description">Get ON/OFF systems that reprogram translation through shape, not sequence.</p>
            <p class="secondary-line">Input your sequence. TADPOLE finds the path to control.</p>
            <p class="rotating-quote" id="quote"></p>
        </div>
    </div>

    <script>
        const quotes = [
            "Forget promoters. Think structure.<br>The future of regulation folds differently.",
            "What if translation was the new frontier of control?<br>Welcome to the RNA switch revolution.",
            "Biology has always known how to self-regulate.<br>We’re just learning to speak its language.",
            "Fast, reversible, translation-level regulation?<br>We’re already skipping limits. Join us, beyond the finish line."
        ];
        let quoteIndex = 0;
        const quoteElement = document.getElementById("quote");

        function rotateQuote() {
            quoteElement.style.opacity = 0;
            quoteElement.style.transform = 'translateY(10px)'; 

            setTimeout(() => {
                quoteElement.innerHTML = quotes[quoteIndex];
                quoteIndex = (quoteIndex + 1) % quotes.length;
                
                quoteElement.style.transform = 'translateY(-10px)'; 

                setTimeout(() => {
                    quoteElement.style.opacity = 1;
                    quoteElement.style.transform = 'translateY(0)'; 
                }, 10); 

            }, 300); 
        }

        document.addEventListener('DOMContentLoaded', () => {
            rotateQuote();
            setInterval(rotateQuote, 6000); 
        });
    </script>

"""

st.components.v1.html(html_content, height=350)

# Sidebar for navegation
with st.sidebar:
    st.image("images/logo.png")
    st.header("TADPOLE")
    selected_main_tab = st.radio(
        "Navegation",
        ["Tools", "Help"],
        index=1 # default: "Model System"
    )

# Main content (acoding to sidebar selection)
if selected_main_tab == "Help":
    help_tab()
elif selected_main_tab == "Tools":

    tools_sub_tabs = st.tabs(["Model System", "Structural RNA element"])

    with tools_sub_tabs[0]: # "Model system"
        linker_finder_tab() 
    with tools_sub_tabs[1]: # "structural RNA element"
        structural_rna_element_tab() 

   

