import streamlit as st
import streamlit.components.v1 as components
import traceback
import matplotlib.pyplot as plt
import os
import RNA
from rna_cluster import matriz_rmsd, cluster_structures, organise_per_clusters
from rna_cluster import compute_metrics, visualise_metrics
from input_utils import parse_fasta_msa
from structure import predict_secundary_structure
from visualizacion import generate_pre_command_per_sequence
from search import run_linker_search
from search import count_rna1_rna3_pairings, constrained_mfe
from structure import predict_secundary_structure, align_secondary_structure 
from rna_structures import get_base_pairs
from genetic_algorithm import run_genetic_algorithm_search 
from conservation import calculate_conservation, clasify_conservation, conservation_category, calculate_complementarity_conservation
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
        transform: translateY(-2px);
    }
    .stButton>button:active {
        background-color: var(--dark-red);
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
        background-color: var(--sidebar-bg);
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
        color: var(--text-color-light);
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
        color: var(--deep-red);
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

def structure_distance(s1, s2):
    """
    Computes the Hamming distance between two RNA secondary structures.

    This function compares two RNA secondary structures (in dot-bracket notation)
    and returns the number of positions where they differ. It assumes both input
    structures are of equal length.

    :param s1: First RNA secondary structure string in dot-bracket notation. (str)
    :param s2: Second RNA secondary structure string in dot-bracket notation. (str)

    :returns: The Hamming distance (number of mismatched characters) between the two structures. (int)
    """
    # Simple Hamming distance between two secondary structures
    return sum(a != b for a, b in zip(s1, s2))


#Write Report
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
    Builds a complete HTML report for RNA switch design results.

    This function assembles a structured, styled HTML report summarizing 
    the outcomes of a computational RNA switch design pipeline. It includes 
    global statistics, search parameters, result clustering, representative images,
    and automated analysis and conclusions.

    :param results: List of design results, each as a dictionary with relevant metrics. (list of dict)
    :param report: Not used in this function but passed externally, possibly for future integration. (any)
    :param rna1: RNA1 sequence string. (str)
    :param rna3: RNA3 sequence string. (str)
    :param struct1: Secondary structure of RNA1 in dot-bracket notation. (str)
    :param constraint: Constraint string used for ON state simulation. (str)
    :param mutable_rna1: List of mutable positions in RNA1. (list of int)
    :param watched_positions: List of positions to track for functional conservation. (list of int)
    :param use_mutations: Whether mutations were allowed in RNA1 (should be renamed to `use_mutations`). (bool)
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

    # Paths to your PNGs for the global plots
    delta_uri       = png_to_data_uri(f"{image_output_dir}/delta_mfe.png")
    pairings_uri    = png_to_data_uri(f"{image_output_dir}/pairings_count.png")

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    total_found = len(results)
    # Ensure total_clusters correctly counts non-noise clusters
    total_clusters = len(set(cluster_labels)) - (1 if -1 in set(cluster_labels) else 0) if (cluster_labels is not None and len(cluster_labels) > 0) else 0

    # Calculate global metrics
    diffs = [res['mfe_2'] - res['mfe_1'] for res in results]
    avg_delta = sum(diffs)/len(diffs) if diffs else 0
    avg_pairs = sum(
        res.get('pairings_count', max_pairings) for res in results
    ) / total_found if total_found else 0

    

    # Group results by cluster
    cluster_groups = defaultdict(list)
    for r in results:
        cluster = r.get('cluster')
        if cluster is not None:
            cluster_groups[cluster].append(r)

    linker_analysis_data = {}

    for cluster, group in cluster_groups.items():
        avg_delta = sum(r['mfe_2'] - r['mfe_1'] for r in group) / len(group)
        avg_pairs = sum(r.get('pairings_count', 0) for r in group) / len(group)
        avg_mfe1 = sum(r['mfe_1'] for r in group) / len(group)

        pros = []
        cons = []

        # ΔMFE
        if avg_delta > 4.0:
            pros.append(f"Good energy difference ({avg_delta:.1f} kcal/mol)")
        elif avg_delta < 2.0:
            cons.append(f"Low energy difference ({avg_delta:.1f} kcal/mol)")
        else:
            pros.append(f"Moderate energy difference ({avg_delta:.1f} kcal/mol)")

        # Pairings
        if avg_pairs < 4:
            pros.append("Low number of pairings between the aptamer and the RNA1 element")
        elif avg_pairs > 6:
            cons.append("High number of pairings between the aptamer and the RNA1 element")
        else:
            pros.append("Moderate pairings between the aptamer and the RNA1 element")

        # OFF structure confidence
        if avg_mfe1 < -20:
            pros.append("High confidence on the OFF structure ensuring no readthrough")
        elif avg_mfe1 > -10:
            cons.append("Low confidence on the OFF structure ensuring no readthrough")

        # Optional: Aptamer accessibility
        if all('accessibility' in r for r in group):
            avg_access = sum(r['accessibility'] for r in group) / len(group)
            if avg_access > 0.7:
                pros.append("Aptamer entry and binding site mostly free")
            elif avg_access < 0.3:
                cons.append("Aptamer entry and binding site may be hard to access")

        linker_analysis_data[cluster] = {
            "Pros": pros,
            "Cons": cons
        }

    # Start building the HTML
    html = f"""
    <html>
    <head>
      <meta charset="utf-8"/>
      <style>
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
      <h1>🔬 RNA1–Aptamer Switch Design Report</h1>
      <p class="small">Generated: {now}</p>

      <div class="section">
        <h2>📚 Introduction and Objectives</h2>
        <p>In this study, we aim to convert a functional element of RNA1 (e.g., a frameshift or SCR) into a molecular "off→on" switch regulated by an aptamer (RNA2) and its ligand. This is achieved by designing a "linker" that connects RNA1 with RNA2, allowing RNA1 to toggle between an "off" (inactive) and an "on" (active) state in response to ligand binding.</p>
        <p>The key aspects to consider for a successful ON/OFF system design are:</p>
        <ol>
          <li><strong>OFF State Energy:</strong> The energy of the unbound (OFF) state should be lower than that of the bound (ON) state without the ligand. Ideally, this difference should be about half the ligand binding energy, ensuring that the bound state is the most stable only when the ligand is present.</li>
          <li><strong>RNA1-RNA3 Interactions:</strong> Pairings between the aptamer (RNA2) and the RNA1 structural element (RNA3) should not exceed a certain threshold to allow the necessary conformational change.</li>
          <li><strong>Aptamer Entry and Binding Sites:</strong> The key entry and binding sites of the aptamer should be mostly free to enable efficient ligand binding.</li>
          <li><strong>Disruption of the OFF State:</strong> Key nucleotides for the RNA1 structural element’s function must be identified and their functional structure disrupted in the unbound (OFF) state to ensure inactivity.</li>
          <li><strong>Maintenance of the ON State:</strong> Structures responsible for function must be preserved in the bound (ON) state to ensure activity.</li>
        </ol>
        <h3>Simulation of the Bound (ON) State with RNAFold</h3>
        <p>Since RNAFold cannot directly simulate interactions between RNA molecules or ligand binding, a constraint-based approach is used to simulate the bound (ON) state. This involves imposing known base-pairing constraints that the aptamer adopts when bound to its ligand.</p>
        <p>In our case, less restrictive constraints were imposed to increase confidence that the resulting structure would naturally occur. Specifically, nucleotides 9, 10, and 11 of the aptamer binding site were constrained to pair only with downstream nucleotides. This simulates how, upon theophylline binding to the aptamer, these nucleotides can no longer pair with the upstream structural element (RNA1), forcing the switch to the "on" state.</p>
      </div>

      <div class="section">
        <h2>⚙️ Search Parameters</h2>
        <ul>
          <li><strong>Search Method:</strong> {search_method}</li>
          <li>Mutations in RNA1: {"Yes" if use_mutations else "No"}</li>
          <li>Mutable positions: {', '.join(map(str, mutable_rna1)) or '—'}</li>
          <li>Watched positions: {', '.join(map(str, watched_positions)) or '—'}</li>
          <li>Maximum allowed mutations: {num_mut}</li>
          <li>Linker length: {linker_min}–{linker_max} nt</li>
          <li>Minimum required ΔMFE: {mfe_delta:.1f} kcal/mol</li>
          <li>Maximum RNA1–RNA3 pairings: {max_pairings}</li>
          <li>Maximum structural changes in RNA1: {max_changes}</li>
          <li>Verbose: {"Yes" if verbose else "No"}</li>
        </ul>
      </div>
    """
    #  <div class="section">
    #    <h2>📊 Summary of Results</h2>
    #    <p>A total of <strong>{total_found}</strong> linkers meeting the search criteria were found, grouped into <strong>{total_clusters}</strong> distinct structural clusters. The average global metrics for the solutions found are:</p>
    #    <table class="metric-table">
    #      <tr><th>Total linkers found</th><td>{total_found}</td></tr>
    #      <tr><th>Total distinct clusters</th><td>{total_clusters}</td></tr>
    #      <tr><th>Average ΔMFE (on–off)</th><td>{avg_delta:.2f} kcal/mol</td></tr>
    #      <tr><th>Average RNA1–RNA3 pairings</th><td>{avg_pairs:.1f}</td></tr>
    #    </table>
    #  </div>
#
    #  <div class="section">
    #    <h2>🔬 Detailed Analysis by Cluster</h2>
    #    <p>Below is an analysis of observed pros and cons for different clusters, based on examples and design criteria:</p>
    #    <table class="pros-cons-table">
    #      <thead>
    #        <tr>
    #          <th>Cluster</th>
    #          <th>Pros</th>
    #          <th>Cons</th>
    #        </tr>
    #      </thead>
    #      <tbody>
    #"""
    ## Dynamically add cluster analysis rows
    #for length in sorted(linker_analysis_data.keys()):
    #    data = linker_analysis_data[length]
    #    pros_list = "".join([f"<li>{p}</li>" for p in data['Pros']])
    #    cons_list = "".join([f"<li>{c}</li>" for c in data['Cons']])
    #    html += f"""
    #        <tr>
    #          <td>{length}</td>
    #          <td><ul>{pros_list}</ul></td>
    #          <td><ul>{cons_list}</ul></td>
    #        </tr>
    #    """

    html += """
          </tbody>
        </table>
      </div>
    """

    if not results:
        html += """
      <div class="section">
        <h2>⚠️ No Valid Results</h2>
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
            if lbl != -1:  # Exclude noise points from clusters
                cluster_dict[lbl].append((idx, results[idx]))

        html_sections_img = ""
        if representative_img_bases:
            html_sections_img += f"""
            <div class="section">
                <h2>🖼️ Representative Structures by Cluster</h2>
                <p>Below are visual examples of the 'off' (unconstrained) and 'on' (constrained) structures for each identified cluster.</p>
                """
            for i, img_base in enumerate(representative_img_bases):
                off_uri = png_to_data_uri(f"{image_output_dir}/linker_{linker_min}/{img_base}_unconstrained_plot.png")
                on_uri = png_to_data_uri(f"{image_output_dir}/linker_{linker_min}/{img_base}_constrained_plot.png")

                html_sections_img += f"""
                <div style="display:inline-block; margin-right:20px; vertical-align:top;">
                  <h3>Cluster {i}</h3>
                  <p><em>Off (unconstrained)</em></p>
                  <img class="img-preview" src="{off_uri}" alt="Unconstrained structure - Cluster {i}" />
                  <p><em>On (constrained)</em></p>
                  <img class="img-preview" src="{on_uri}" alt="Constrained structure - Cluster {i}" />
                </div>
                """
            html_sections_img += "</div>"

        html += html_sections_img

        for cluster_id, items in sorted(cluster_dict.items()):
            if cluster_id == -1:
                continue

            idx, res = items[0]  # Representative example for this cluster

            html += f"""
      <div class="section">
        <h2>🧬 Cluster Details {cluster_id} — {len(items)} variants</h2>
        <div class="linker-box">
          <h3>Example Linker in Cluster {cluster_id}: <code>{res['linker']}</code></h3>
          <p><strong>Full Sequence (RNA1+Linker+RNA3):</strong> <code>{res['sequence']}</code></p>
          <p><strong>Mutations in RNA1:</strong> {res.get('mut1_info', 'N/A')}</p>
          <p><strong>ΔMFE (Off→On):</strong> {(res['mfe_2']-res['mfe_1']):.2f} kcal/mol — 
              <strong>RNA1–RNA3 Pairings:</strong> {res.get('pairings_count','—')}</p>
          </div>
      </div>
            """

        html += f"""
      <div class="section">
        <p>🔍 Full images and data for each linker are available at: <code>{image_output_dir}/</code></p>
      </div>
        """

    stability_ok = avg_delta >= mfe_delta
    low_interaction = avg_pairs <= max_pairings
    html += f"""
      <div class="section">
        <h2>🔎 Automated Conclusions</h2>
        <ul>
          <li><strong>Relative “off”→“on” stability (ΔMFE):</strong> {'✅ Met' if stability_ok else '❌ Insufficient'} (Average ΔMFE = {avg_delta:.2f} kcal/mol, required = {mfe_delta:.1f} kcal/mol)</li>
          <li><strong>Limited RNA1–RNA3 Interactions:</strong> {'✅ Met' if low_interaction else '❌ Excess'} (Average pairings = {avg_pairs:.1f}, max allowed = {max_pairings})</li>
          <li><strong>Functional integrity of RNA1 (watched positions):</strong> Verifiable in {total_found} cases. Review individual results for each watched position is recommended.</li>
          <li><strong>Cluster diversity:</strong> {'High' if total_clusters>3 else 'Low'}. {'Multiple structural topologies identified.' if total_clusters>1 else 'A single dominant topology found.'}.</li>
        </ul>
      </div>

      <div class="section">
        <h2>💡 Possible Improvements and Next Steps</h2>
        <p>To further optimize the RNA switch design, consider the following recommendations:</p>
        <ul>
          <li><strong>Adjust Search Parameters:</strong>
            <ul>
              <li>If the number of results is low or metrics are suboptimal, try relaxing the minimum ΔMFE, reducing RNA1–RNA3 pairing restrictions, or expanding the linker length range.</li>
              <li>For more complex systems or greater diversity, experiment by allowing more mutations in RNA1 or tuning genetic algorithm parameters (population size, generations, mutation/crossover rates).</li>
            </ul>
          </li>
          <li><strong>Refine Watched Positions:</strong>
            <ul>
              <li>Ensure `watched_positions` cover all critical functional domains of RNA1 whose disruption/maintenance is essential for switch function. Consider adding or modifying these positions based on biological knowledge.</li>
            </ul>
          </li>
          <li><strong>Experimental Validation:</strong>
            <ul>
              <li>Predicted designs should be experimentally validated to confirm their "off→on" behavior and ligand binding efficiency.</li>
              <li>Consider targeted mutagenesis experiments on linkers or RNA1 to validate structural and functional predictions.</li>
            </ul>
          </li>
          <li><strong>Additional Cluster Analysis:</strong>
            <ul>
              <li>Examine representative structures and linker sequences in each cluster to identify common motifs or desirable features.</li>
              <li>Clusters with high internal diversity or those with better average metrics may be promising candidates for testing.</li>
            </ul>
          </li>
          <li><strong>Exploration of Alternative Aptamers:</strong>
            <ul>
              <li>If performance remains suboptimal, consider using or designing aptamers with different binding or stability characteristics that may better integrate with your RNA1.</li>
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


# ===== Pestaña principal =====
def structural_rna_element_tab():
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

        # Ejemplo SECIS
        example_box("""
        The SECIS element is an RNA structure that stimulates stop codon readthrough.  
        Normally, translation terminates at a stop codon, releasing the protein.  
        However, the SECIS element can cause the ribosome to skip the stop codon, resulting in a longer protein product.  
        Its function relies entirely on its secondary structure, so by controlling this structure, we can modulate its regulatory effect.
        """)

        # Imagen centrada
        col1, col2, col3 = st.columns([1, 4, 1])
        with col2:
            st.image("images/SCR.png", width=900)

        st.subheader("Secondary Structure Considerations")
        st.markdown("""
        The SRE must have a known or predictable secondary structure, typically represented in dot-bracket notation.  
        If an experimentally validated structure is not available, users are encouraged to predict one using established tools—**RNAfold** from the Vienna RNA Package is highly recommended.  
        """)

        # Nota especial
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

        st.subheader("Input Requirements")

        # Bloque Mutations
        input_block("Mutations") 
        st.markdown("""
        This software allows the user to introduce mutations in the SRE sequence to explore functional variants.  
        If a nucleotide involved in a base pair is mutated, its paired base is automatically mutated to a complementary one, preserving the secondary structure.""")
        # Ejemplo imagen SECIS
        example_box("""
        In the SECIS element, certain positions (green) are key for promoting readthrough.  
        As the bottom positions are not important for function, they can be mutated without disrupting activity.
        """)
        col1, col2, col3 = st.columns([4, 4, 1])
        with col2:
            st.image("images/SECIS.png", width=100)
        st.markdown("""  
        It is recommended to limit mutations to non-critical nucleotides to avoid compromising function.
        """)

        # Bloque Key nucleotides
        input_block("Key Nucleotides for Function") 
        st.markdown("""
        The software allows specifying “watched positions”: nucleotides known to be critical for function.  
        """)
        example_box("""
        In the case of the SECIS element, the green positions are key for promoting readthrough.  
        """)
        

        st.markdown("""
        Watched positions are used as functional markers:  
        - If their structure relative to the input structure is maintained, the SRE is considered functional.  
        - If disrupted, the SRE is considered non-functional.

        <b>Hint:</b> Avoid including all key nucleotides in ‘Watched positions’.  
        Even if one maintains pairing, the overall structure could still change.  
        Experiment with different sets to achieve accurate results.
        """, unsafe_allow_html=True)

        # Bloque Allow changes
        input_block("Allow Changes on SRE’s Structure") 
        st.markdown("""
        While certain nucleotides and structural features of an SRE are essential, others may be flexible or non-critical.  
        This tool includes a parameter: **Maximum changes on the RNA1 structure**, defining tolerated deviation from the original structure.

        - A strict value (0) means the structure must remain identical.  
        - A higher value allows more flexibility.
        """)

        # Bloque final
        st.markdown("""
        <h3>Why These Inputs Matter</h3>
        RNA-based systems are complex; small structural changes can have a big impact, but not all regions matter equally.  
        This tool lets you define what matters most, preserving function where it counts while exploring a broad design space.
        """, unsafe_allow_html=True)








    st.markdown("""<div style="font-size: 2.7rem; color: #ea534e; font-weight: bold; margin-top: 1rem;">
        MSA Tool:  
    </div>
    """, unsafe_allow_html=True)


    # --- Reubicación del cargador de archivos y la sección "About" del antiguo sidebar ---
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
            # --- START tab1_analysis content ---
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
                        
            # --- END tab1_analysis content ---
        else:
            st.info("Once the file is uploaded, here you will find the predicted structure and energy of the first seuqence of the MSA (presumably the one you are interesetd on characterising)"
            "and a conservation map, with the conservation score per position.")
    # --- START tab1_analysis content ---
    with structural_rna_sub_sub_tabs[1]: # structural clustering
        if msa and len(msa) >= 2: # Clustering requires at least 2 sequences
            # --- START tab2_clustering content ---
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
                sequence = msa[0].replace("-", "") 
                msa_len = len(msa[0]) 
                col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
                colors = {
                    "always": (255, 0, 0),
                    "except2": (0, 0, 255),
                    "except4": (255, 255, 0),
                    "others": (255, 255, 0)
                }

                for cluster_id, cluster_data in cluster_summary.items():
                    with st.expander(f"Cluster {cluster_id} — {cluster_data['num_variants']} structures"):
                        rep = cluster_data['representative_structure']
                        st.markdown(f"Representative Structure: `{rep}`")
                        for var in cluster_data.get("members", []):
                            st.markdown(f"Specific sequence:  `{var['sequence']}`")

                        fold_path = os.path.join(output_dir, f"cluster_{cluster_id}.fold")
                        with open(fold_path, "w") as f:
                            f.write(sequence + "\n" + rep + "\n")

                        pre_command = generate_pre_command_per_sequence(msa[0], col_categories, colors, lw=10)
                        rnaplot_cmd = f'RNAplot --pre "{pre_command}" < {fold_path}'
                        run_command(rnaplot_cmd)

                        ps_path = os.path.join(output_dir, f"cluster_{cluster_id}.ps")
                        png_path = os.path.join(output_dir, f"cluster_{cluster_id}.png")
                        if os.path.exists("rna.ps"):
                            os.rename("rna.ps", ps_path)
                            run_command(f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
                                        f'-dEPSCrop -sOutputFile={png_path} {ps_path}')
                            if os.path.exists(png_path):
                                st.image(png_path, caption=f"Structure Cluster {cluster_id}")
                        
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


            # --- END tab2_clustering content ---
        else:
            st.info("This section allows you to see the structure of the different unique structures of your MSA, with the nucleotides colored by conservation.")

    with structural_rna_sub_sub_tabs[2]: # Visualizar Todas las Estructuras
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

                    fold_file = os.path.join(output_dir, f"cluster_images_{i}.fold")
                    with open(fold_file, "w") as f:
                        f.write(sequence + "\n" + structure + "\n")

                    run_command(f'RNAplot --pre "{pre_command}" < {fold_file}')

                    ps_file = os.path.join(output_dir, f"cluster_images_{i}.ps")
                    png_file = os.path.join(output_dir, f"cluster_images_{i}.png")

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
                label="RNA Segment 1",
                value='ACCAGUGUGCGGAUGAUAACUACUGACGAAAGAGUCAUCGACUCAGUUAGUGGUUGGAUGUAGUCACAUUAGU',
                height=110,
                placeholder="Introduce the RNA1 sequence  (A, C, G, U).",
                help="This is the RNA element whose function (to be modulated) depends on structure."
            )

            rna3 = st.text_area(
                label="RNA Segment 2",
                value='AGUUGGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACCAAAA',
                height=110,
                placeholder="Introduce the RNA2 sequence (A, C, G, U).",
                help="This is the aptamer's sequence."
            )

            struct1 = st.text_input(
                label="Targeted structure of the RNA1 element (dot-bracket)",
                value='((.(((((((......(((((((((((....(((((...))))).)))))))))))......).)))))).))',
                placeholder="Example: ((.((...)))).",
                help="Put here the theoretical structure of the RNA1 element (determined experimentally or through " \
                "a prediction software, like RNAFold). This is the structure that the RNA1 sequence is going to be forced to have on the constrained (ON) state"
            )
            st.markdown('</div>', unsafe_allow_html=True)

        with st.expander("⚙️ Advanced Parameters", expanded=False):

            # --- Search Method Selection ---
            search_method = st.radio(
                "Select Search Method",
                ("Brute Force", "Genetic Algorithm"),
                help="Choose between exhaustive search (Brute Force) or an evolutionary algorithm (Genetic Algorithm)."
            )

            constraint = st.text_input(
                label="Restriction chain for the RNA2 sequence",
                value='........................................................................................<<<<<<<...............................',
                placeholder="Example: ...........<<<<<............",
                help="Restrictions for the folding of the constrained (ON) state. It must have the same " \
                "length as the RNA2 sequence. '.' means no restriction, '(' and ')' for corresponding pairs, 'x' " \
                "to avoid pairing that nucleotide and '>' to avoid pairing with later nts and '<' to avoid pairing with previous nucleotides."
            )
            mutable_rna1 = st.text_input(
                label="Mutable positions for RNA1",
                value="1, 2, 4, 5, 6, 7, 8, 9",
                placeholder="Ex: 1, 2, 4, 5, 6",
                help="This positions should be nucleotides known for NOT having an impact on functionallity. " \
                "Notice that, if paired on the input targeted structure, the pairs will be mutated aswell to mantain complementarity."
            )
            watched_positions_str = st.text_input(
                label="Watched positions for RNA1",
                value="56",
                placeholder="Ex: 56, 57",
                help="Positions key to functionality for the RNA1 element whose pairing "
                "(from respect to the input targeted structure) will be forced to change to ensure the OFF state disrupts functionality."
            )
        
            use_mutations = st.checkbox("Allow mutations on RNA1", value=True)
        
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
                help="The linker is the key element that will be randomised untill one is found such that it makes the ON/OFF system possible"
            )
            linker_max = st.slider(
                label="Linker's maximum length",
                min_value=linker_min,
                max_value=10,
                value=7,
                step=1,
                help="Linkers of length from the minimum length up to this length will be tested."
            )
        
            mfe_delta = st.number_input(
                label="Minimum Energy difference (kcal/mol)",
                min_value=0.0, max_value=20.0, value=0.0, step=0.1,
                help="Difference between the energy of the constrained" \
                "(ON) state and the unconstrained (OFF) state. It should be approximately half " \
                "the binding energy of the aptamer–ligand interaction."
            )
            max_pairings = st.number_input(
                label="Maximum number of RNA1-RNA3 pairings",
                min_value=0, max_value=20, value=15, step=1,
                help="Binding with the ligand can't brake many pairings. To ensure the switch's functionality, " \
                "one should choose a low number of pairings. However, if this conditions is very restrictive, " \
                "the ON/OFF system might not be possible to construct, and thus it is allowed to modulate this value."
            )
            max_changes = st.number_input(
                label="Maximum changes on the RNA1 structure",
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
                    min_value=10, max_value=500, value=50, step=10,
                    help="Number of individuals in each generation."
                )
                ga_generations = st.number_input(
                    label="Number of Generations",
                    min_value=10, max_value=1000, value=100, step=10,
                    help="Number of generations to run the algorithm."
                )
                ga_elitism_rate = st.slider(
                    label="Elitism Rate",
                    min_value=0.0, max_value=0.5, value=0.1, step=0.01,
                    help="Fraction of the best individuals that are directly passed to the next generation."
                )
                ga_mutation_rate_rna1 = st.slider(
                    label="RNA1 Mutation Rate",
                    min_value=0.0, max_value=0.2, value=0.02, step=0.005,
                    help="Probability of a base mutation in RNA1 per position."
                )
                ga_mutation_rate_linker = st.slider(
                    label="Linker Mutation Rate",
                    min_value=0.0, max_value=0.5, value=0.05, step=0.01,
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
        folder = st.text_area(
                label="Folder name",
                value='linker_search',
                placeholder="WARNING: If you do not change this name, " \
                "the results will be mixed and some might be rewriten.",
                help="Folder name in which the files will be placed. "
            )
        st.markdown("---")
        # Execute and STOP buttons
        run_col1, run_col2 = st.columns([4, 1])
        with run_col1:
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
                st.warning("Please, introduce RNA1 and RNA3 sequences and targeted structure")
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
                        update_log("Brute force search completed. No valid linker found.")  
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
                            linker_lengths=linker_length_for_ga, # Pass this as linker_length_ga to the function's arg
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
                        update_log("Genetic Algorithm search completed. No valid linker found.")                
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
                st.warning("No valid linkers found.")
                st.markdown("### 📝 Summarised search report:")
                st.text(report)

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
                    representative = min(cluster_items, key=lambda x: x[1]['mfe_2'])  

                    idx, res = representative
                    representative_img_bases.append(res['linker'])
                    with st.expander(f" Representative Linker (#{idx + 1}): {res['linker']}"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**RNA1 mutations:** {res['mut1_info']}")

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

            # This is where you generate your rich HTML report
            html = build_full_html_report(
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

            st.markdown("### 📝 Full Search Report:")
            

            #pdf_path = generar_pdf_desde_html(html)

            #with open(pdf_path, "rb") as f:
            #    st.download_button(
            #        label="📥 Download full PDF report",
            #        data=f.read(),
            #        file_name="informe_linkers.pdf",
            #        mime="application/pdf"
            #    )

            show_all = st.checkbox("Show all linkers")
            if show_all:
                for idx, res in enumerate(results):
                    with st.expander(f" Linker {idx + 1}: {res['linker']}"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Unconstrained structure:** `{res['structure_unconstrained']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**Mutations on RNA1:** {res['mut1_info']}")
                        
                        linker_for_img = res['linker']
                        linker_len_for_img = len(linker_for_img)
                        folder = f"{folder}/linker_{linker_len_for_img}"
                        
                        unconstrained_image_path = f"{folder}/{linker_for_img}_unconstrained_plot.png"
                        constrained_image_path = f"{folder}/{linker_for_img}_constrained_plot.png"
                        
                        if os.path.exists(unconstrained_image_path):
                            st.image(unconstrained_image_path, caption=f"Unconstrained Structure for Linker {linker_for_img}")
                        if os.path.exists(constrained_image_path):
                            st.image(constrained_image_path, caption=f"Constrained Structure for Linker {linker_for_img}")

# --- Content for "Aptamers" sub-tab ---
def aptamers_tab():
    st.subheader("What is an Aptamer?")
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
    """)

    st.info("""
    **Important Notes for Users**  
    - If you're using a different aptamer than the default one, you must adjust this parameter based on:
        - The length of the aptamer.
        - The position of ligand-binding nucleotides.
    - The user is responsible for verifying that, in the constrained structure of the design, the aptamer folds into the expected structure in the presence of the ligand.
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
    - Check if the secondary structure of the bound configuration matches the aptamer’s structure on the constrained structure of the system.
    - Define the “Restriction chain for the aptamer sequence” based on its ligand-binding region.
    """)

def system_detail():
    st.header("The system in detail")
    st.markdown("""
    Let’s break down how the three key parameters shape the system's functionality:

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

    ### 2. Minimum ΔMFE (ON-OFF Energy Difference)

    <span style="color:#1f77b4; font-weight:bold;">What it is:</span>  
    The minimum required difference in Minimum Free Energy (MFE) between the OFF and ON states. MFE is a proxy for how stable a structure is.

    <span style="color:#ff7f0e; font-weight:bold;">Why it matters:</span>  
    You want the system to have two distinct and energetically stable conformations:  
    - One in the absence of ligand (OFF), where the SRE is disrupted  
    - One in the presence of ligand (ON), where the SRE folds correctly  

    If both states have similar energies, the switch may not behave reliably, flipping randomly or being stuck in one state. By setting a minimum ΔMFE, you enforce a thermodynamic bias that makes switching behavior more robust and predictable.

    <span style="color:#2ca02c; font-weight:bold;">Goal:</span>  
    Ensure there's a meaningful, energetic difference between ON and OFF that reflects real structural change. This enhances the functional specificity of the switch.

    In this image, you can see the energy of the different configurations.  
    The key concept is that the binding with the ligand changes the energy of the system (in the case of the theophylline aptamer, the binding energy is -9.5 kcal/mol). Therefore, the same structure (in this case, the OFF state) can have a different energy depending on whether the ligand is present or not.

    Without the ligand, we want to have the OFF state, so we want the OFF state without the ligand (central figure) to be the most stable one (and therefore have a lower energy than the ON state with the ligand, left figure). With the ligand, we want to have the ON state, so we want the ON state (left figure) to be more stable than the OFF state with the ligand (right figure) (and therefore have a lower energy).

    ### 3. Maximum Number of RNA1–RNA3 Pairings

    <span style="color:#1f77b4; font-weight:bold;">What it is:</span>  
    A constraint on how many base pairs are allowed to form between RNA1 (the SRE) and RNA3 (the aptamer) in the OFF state.

    <span style="color:#ff7f0e; font-weight:bold;">Why it matters:</span>  
    In the OFF state, the aptamer must partially bind the SRE to prevent it from forming its active structure. However:  
    - Too few pairings might not be enough to effectively disrupt the SRE.  
    - Too many pairings could over-stabilise the OFF state, making the ON state hard to reach even when the ligand is present.

    <span style="color:#2ca02c; font-weight:bold;">Goal:</span>  
    Strike a balance: allow enough pairing for functional repression, but leave room for switching. This variable helps tune the strength of repression and influences whether ligand binding can successfully revert the system to the ON state.  
    This parameter is especially useful when testing different SRE or aptamer sequences, as their pairing potential may vary significantly.
    """, unsafe_allow_html=True)





def help_tab():
    st.header("Help and Documentation")
    st.write("Find user guides and tutorials")

    help_sub_tabs = st.tabs([" User Guide", " Documentation"])

    # === USER GUIDE TAB ===
    with help_sub_tabs[0]:
        st.markdown("""
        <h1 style="text-align:center; margin-bottom: 0.5em;">Introduction to the RNA ON/OFF Switch Design Toolkit</h1>
        """, unsafe_allow_html=True)

        st.markdown("""
        <h2 style="margin-top: 1em;">How does the toolkit work?</h2>
        <p>The web interface is divided into intuitive sections, each supporting a key part of the design process:</p>

        <h3 style="margin-top: 1em; color:#4a90e2;">Model System</h3>
        <p>The core module. It helps you explore candidate linker sequences that allow your aptamer and structural RNA to interact correctly in the OFF state, and separate correctly in the ON state (when the ligand is present).<br>
        It simulates folding behaviour, compares energy states, and allows for visualisation of the structure of the system (ON and OFF).</p>

        <h3 style="margin-top: 1em; color:#4a90e2;">About the parts</h3>
        <h4 style="margin-top: 1em; color:#f7bcbb;"> SRE (Structural RNA Element) Analysis</h4>
        <p>This section helps you understand SREs. It also has a tool to make an evolutionary analysis from a multi-sequence alignment (MSA), to help identify the functional regions of your structural element, relevant for the design of the system.</p>

        <h4 style="margin-top: 1em; color:#f7bcbb;">Aptamer</h4>
        <p>This section helps you understand what Aptamers are and the key concepts to have in mind.</p>

        <h3 style="margin-top: 1em; color:#4a90e2;">Documentation</h3>
        <p>All calculations are powered by <strong>RNAfold</strong>, and results are presented with plots, base-pairing histograms, energy comparisons, and automatically generated reports.</p>

        <hr style="margin-top:2em; margin-bottom:2em;">
        """ , unsafe_allow_html=True)

        # Imagen centralizada con columnas
        col1, col2, col3 = st.columns([1.5, 4, 1])
        with col2:
            st.image("images/system.png", width=700)

        st.markdown("""
        **ON state vs OFF state:**

        - **OFF state:**  
          The aptamer and the SRE are paired. This means the structure of the SRE is disrupted, and it will not perform its function — this is the OFF state.

        - **ON state:**  
          The aptamer is bound to the ligand and forms an independent structure, releasing the SRE. The SRE can then form the proper structure to perform its function, therefore leading to an ON state.

        ---
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
            st.image("images/OFF.png", width=900)

        example_cont_box(""
        "- In the presence of ligand: the aptamer binds to it and forms its own structure, releasing the SCR to fold normally, stimulating readthrough and producing protein AB." \
        "")
        col1, col2, col3 = st.columns([1, 4, 1])
        with col2:
            st.image("images/ON.png", width=900)

        st.markdown("---")

        st.markdown("""
        ### What the toolkit does:  
        This toolkit automates the design of the entire system:
        - Searches for optimal linker sequences.
        - Tests structural behaviour.
        - Simulates ligand binding.
        - Proposes concrete sequences that meet precise energetic and functional criteria.

        <p style="color:#f28f8c; font-weight:bold; font-size:1.2em; margin-top:1em;">
        It allows researchers to build working ON/OFF systems based on any SRE and aptamer combination, not just predefined templates.
        </p>

        ---
                    
        <h2>What is the value?</h2>

        <h3>Unlocking RNA control beyond the old playbook</h3>
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

        <h3>From concept to working switch... in minutes</h3>
        <p>
        Our toolkit automates every step of the design process:<br>
        - Smart linker selection – Finds spacer sequences that connect the aptamer and SRE, determining their folding.<br>
        - Two-state RNA folding prediction – Models the RNA’s shape in both ON and OFF states.<br>
        - Performance scoring – Ranks designs based on predicted stability and how well key structural features are preserved.<br>
        - Clear, interactive results – Delivers visualisations, scores, and ready-to-test sequences with no coding required.
        </p>
        
        """ , unsafe_allow_html=True)
        
        # Before / After esquema en columnas
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

    with help_sub_tabs[1]:  # Documentation
        st.header("Modules and Functions Documentation")

        st.subheader("conservation.py")
        st.markdown("""
        This module contains functions to analyze sequence conservation in alignments.

        **calculate_complementarity_conservation(msa, idx, pairs):**  
        Description: Calculates complementarity conservation at a given position (`idx`) of a multiple sequence alignment (MSA). It identifies the paired position and evaluates consistency of complementary base pairs across all MSA sequences.  
        Parameters:  
        - `msa`: A multiple sequence alignment as a list of strings.  
        - `idx`: The zero-based index of the nucleotide position to analyze.  
        - `pairs`: A list of tuples representing paired base positions.  
        Returns:  
        A string describing the level of complementarity conservation. Possible values include `"always_conserved"`, `"not_paired"`, `"-"`, or `"others"`.
        """)

        st.subheader("genetic_algorithm.py")
        st.markdown("""
        This module implements the genetic algorithm to design RNA switches.

        **get_pair_table(struct):**  
        Description: Generates a pairing table for a given RNA secondary structure in dot-bracket notation.  
        Parameters:  
        - `struct`: The RNA secondary structure string in dot-bracket notation.  
        Returns:  
        A 1-based list representing the pairing table.

        **initialize_population(rna1_orig, rna3_orig, linker_length, population_size, mutable_positions_rna1, struct1_orig):**  
        Description: Initializes a population of individuals for the genetic algorithm. Each individual is an RNA switch candidate with a full sequence, mutated rna1 segment, linker, and rna3 segment.  
        Parameters:  
        - `rna1_orig`: Original rna1 sequence.  
        - `rna3_orig`: Original rna3 sequence.  
        - `linker_length`: Desired length of linker sequences.  
        - `population_size`: Number of individuals to generate.  
        - `mutable_positions_rna1`: List of zero-based positions in rna1 where mutations are allowed.  
        - `struct1_orig`: Original rna1 structure in dot-bracket notation.  
        Returns:  
        A list of dictionaries, each representing an individual.

        **fitness(individual, rna1_orig, rna3_orig, struct1_orig, constraint_orig, watched_positions_orig, mfe_delta, max_pairings, max_structure_changes):**  
        Description: Calculates the fitness score of an individual. Evaluates the RNA switch based on criteria such as observed position changes, minimizing rna1-rna3 pairings, preserving rna1 structure, and minimum free energy (MFE) difference.  
        Parameters:  
        - `individual`: Dictionary representing the individual to evaluate.  
        - `rna1_orig`: Original rna1 sequence.  
        - `rna3_orig`: Original rna3 sequence.  
        - `struct1_orig`: Original rna1 structure.  
        - `constraint_orig`: Dot-bracket constraint string guiding folding toward the "ON" state.  
        - `watched_positions_orig`: List of positions in rna1 expected to show structural change in the "OFF" state.  
        - `mfe_delta`: Minimum acceptable energy difference between "ON" and "OFF" states.  
        - `max_pairings`: Maximum allowed base pairs between rna1 and rna3.  
        - `max_structure_changes`: Maximum allowed character differences between original rna1 structure and rna1 substructure within the full construct.
        """)

        st.subheader("input_utils.py")
        st.markdown("""
        This module manages loading and parsing alignment files.

        **parse_fasta_msa(uploaded_file):**  
        Description: Parses an uploaded file as a multiple sequence alignment (MSA) in FASTA or Clustal format. Performs validations to ensure the alignment is suitable for downstream analyses.  
        Parameters:  
        - `uploaded_file`: A file-like object expected to have a `read()` method.
        """)
        st.subheader("io_tools.py")
        st.markdown("""
        This module handles data input/output, including saving and plotting structures.

        **save_and_plot_structures(seq, structure_unconstr, structure_constr, rna1, linker, rna3, mut1_info, mfe_1, mfe_2, folder_prefix='propuestas'):**  
        Description: Saves and plots the unconstrained and constrained secondary structures of an RNA switch. Creates folder structure, saves detailed info, and generates files for ViennaRNA's RNAplot tool.  
        Parameters:  
        - `seq`: The full RNA sequence.  
        - `structure_unconstr`: The unconstrained structure in dot-bracket notation.  
        - `structure_constr`: The constrained structure in dot-bracket notation.  
        - `rna1`, `linker`, `rna3`: Sequences of the switch segments.  
        - `mut1_info`: Information about mutations in rna1.  
        - `mfe_1`, `mfe_2`: MFE values for the structures.  
        - `folder_prefix`: Output folder prefix, default is "propuestas".
        """)

        st.subheader("rna_cluster.py")
        st.markdown("""
        This module clusters RNA structures based on base-pair similarity.

        **parse_pairs(structure):**  
        Description: Parses a dot-bracket structure string and extracts all base pairs.  
        Parameters:  
        - `structure`: RNA structure in dot-bracket notation.  
        Returns:  
        A list of tuples, each tuple representing a base pair.

        **rmsd_pairs(pairs1, pairs2):**  
        Description: Calculates the root mean square deviation (RMSD) between two sets of base pairs, quantifying structural similarity between two RNA secondary structures.  
        Parameters:  
        - `pairs1`: List of base pairs from the first structure.  
        - `pairs2`: List of base pairs from the second structure.

        **rmsd_matrix(structures):**  
        Description: Computes an RMSD distance matrix for a set of structures, where each element [i, j] is the RMSD between structure i and j.  
        Parameters:  
        - `structures`: List of RNA structures in dot-bracket notation.  
        Returns:  
        A NumPy matrix of RMSD distances.

        **cluster_structures(rmsd_matrix, eps=1.0):**  
        Description: Performs clustering of structures using DBSCAN algorithm based on RMSD distance matrix.  
        Parameters:  
        - `rmsd_matrix`: RMSD distance matrix.  
        - `eps`: Maximum distance between two samples for them to be considered in the same neighborhood. Default is 1.0.  
        Returns:  
        A list of cluster labels where -1 indicates noise.

        **organise_per_clusters(structures, clusters):**  
        Description: Organizes structures into a dictionary where keys are cluster labels and values are lists of structures in that cluster.  
        Parameters:  
        - `structures`: List of structures.  
        - `clusters`: List of cluster labels.  
        Returns:  
        Dictionary of structures grouped by cluster.

        **compute_metrics(clusters_dict):**  
        Description: Computes various metrics for each cluster such as number of elements and sequences, storing them in a dictionary.  
        Parameters:  
        - `clusters_dict`: Dictionary of structures grouped by cluster.  
        Returns:  
        Dictionary of cluster metrics.

        **visualise_metrics(clusters_metrics):**  
        Description: Visualizes the computed metrics for each cluster by creating bar charts or other visualizations.  
        Parameters:  
        - `clusters_metrics`: Dictionary of cluster metrics.
        """)

        st.subheader("rna_mutation.py")
        st.markdown("""
        This module contains functions for mutating RNA sequences.

        **mutate_sequence(seq, mutable_positions, n_mut):**  
        Description: Generates all possible mutated sequences applying a specific number of mutations (`n_mut`) at allowed positions (`mutable_positions`).  
        Parameters:  
        - `seq`: Original nucleotide sequence.  
        - `mutable_positions`: List of zero-based positions where mutations are allowed.  
        - `n_mut`: Exact number of mutations to apply.  
        Returns:  
        A generator yielding tuples containing the new mutated sequence and a list of applied mutations.
        """)

        st.subheader("rna_structures.py")
        st.markdown("""
        This module provides utility functions for working with RNA structures.

        **get_pairing_dict(structure):**  
        Description: Parses an RNA secondary structure in dot-bracket notation and returns a dictionary mapping each paired base index to its partner's index.  
        Parameters:  
        - `structure`: RNA secondary structure in dot-bracket notation.  
        Returns:  
        A dictionary of base pairs.
        """)

        st.subheader("search.py")
        st.markdown("""
        This module includes functions to evaluate structures and pairings.

        **count_rna1_rna3_pairings(struct, rna1, rna3, linker):**  
        Description: Counts the number of base pairs formed between rna1 and rna3 segments within a full secondary structure. This helps evaluate unwanted interactions.  
        Parameters:  
        - `struct`: Dot-bracket string of the full structure.  
        - `rna1`, `rna3`, `linker`: Sequences of the segments.  
        Returns:  
        Total number of base pairs formed between rna1 and rna3.

        **check_rna1_structure_preserved(rna1, linker, rna3, structure_constrained, structure_rna1):**  
        Description: Calculates the number of differences between the original rna1 structure and the predicted rna1 substructure when part of the full constrained construct. Checks if rna1 segment maintains its desired structure in the "ON" state.  
        Parameters:  
        - `rna1`, `linker`, `rna3`: Sequences of the segments.  
        - `structure_constrained`: Dot-bracket string of the full structure with constraint.  
        - `structure_rna1`: Original dot-bracket structure of just the rna1 segment.  
        Returns:  
        Number of differing characters between the two structures.
        """)

        st.subheader("structure.py")
        st.markdown("""
        This module focuses on RNA secondary structure prediction.

        **predict_secundary_structure(sequence):**  
        Description: Predicts the minimum free energy (MFE) secondary structure of an RNA sequence using ViennaRNA's RNA.fold function.  
        Parameters:  
        - `sequence`: RNA nucleotide sequence as a string.  
        Returns:  
        A tuple containing the predicted secondary structure in dot-bracket notation and the MFE in kcal/mol.
        """)

        st.subheader("visualisation.py")
        st.markdown("""
        This module handles visualization of RNA structures.

        **generate_rnaplot_with_colours(sequence, secundary_structure, groups, tag='conservation', output_dir='outputs'):**  
        Description: Generates a plot of the RNA secondary structure with specific positions colored according to predefined conservation categories. Uses ViennaRNA's RNAplot tool.  
        Parameters:  
        - `sequence`: RNA nucleotide sequence.  
        - `secundary_structure`: Dot-bracket notation of the structure.  
        - `groups`: A dictionary where keys are category names and values are lists of 1-based nucleotide indices belonging to that category.  
        - `tag`: Prefix for output filenames.  
        - `output_dir`: Directory to save output files.  
        Returns:  
        A PIL image object of the generated PNG plot, or None if failed.
        """)
        

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
.main-header {
    padding: 2rem 1rem;
    background-color: var(--primary-red); /* Using new light red background */
    border-radius: 2rem;
    text-align: center;
    font-family: 'Poppins', sans-serif; /* Ensure 'Poppins' is imported if used, otherwise fallback */
}
.app-title {
    /* This will now be controlled by the h1 style above */
    /* font-size: 2.2rem; - Removed as h1 defines 2.5em */
    color: var(--medium-light-red);
    margin-bottom: 0.5rem; /* Kept for specific spacing */
}
.app-description, .secondary-line {
    font-size: 1.4rem;
    color:#000000; /* Using general dark text color */
    margin: 0.3rem 0;
}
.rotating-quote {
    font-size: 1.3rem;
    font-style: italic;
    margin-top: 1.2rem;
    color: var(--text-darker-red); /* Using a lighter accent red for the quotes */
    min-height: 2rem;
    transition: opacity 0.3s ease-out, transform 0.3s ease-out;
    opacity: 1; 
    transform: translateY(0); 
}
</style>

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

rotateQuote();
setInterval(rotateQuote, 6000); 
</script>
"""

st.components.v1.html(html_content, height=350)

# Sidebar for navegation
with st.sidebar:
    st.image("images/logo.png")
    st.header("TADPOLE")
    selected_main_tab = st.radio(
        "Navegation",
        ["Model System", "About the Parts", "Help"],
        index=2 # default: "Model System"
    )

# Main content (acoding to sidebar selection)
if selected_main_tab == "Help":
    help_tab()
elif selected_main_tab == "About the Parts":

    about_parts_sub_tabs = st.tabs(["Structural RNA element", "Aptamers", 'System in Detail'])

    with about_parts_sub_tabs[0]: # "structural RNA element"
        structural_rna_element_tab() 

    with about_parts_sub_tabs[1]: # "Aptamers"
        aptamers_tab()
    with about_parts_sub_tabs[2]: # "Aptamers"
        system_detail()
elif selected_main_tab == "Model System":
    linker_finder_tab()

