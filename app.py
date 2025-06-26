import streamlit as st
import traceback
import RNA
from rna_cluster import matriz_rmsd, cluster_structures, organise_per_clusters
from rna_cluster import compute_metrics, visualise_metrics
from input_utils import parse_fasta_msa
from structure import predict_secundary_structure
from visualizacion import generate_pre_command_per_sequence
from search import run_linker_search
from search import count_rna1_rna3_pairings, constrained_mfe
from genetic_algorithm import run_genetic_algorithm_search 
from conservacion import calculate_conservation, clasify_conservation, conservation_category
from scipy.spatial.distance import pdist, squareform
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import os
import shutil
import subprocess
import time
import matplotlib.pyplot as plt
from io import StringIO
from contextlib import redirect_stdout
import tempfile
import datetime
from collections import defaultdict
import base64 

st.set_page_config(page_title="RNA Conservation & Structure Analysis", layout="wide")

# ---------------- CUSTOM STYLES ----------------
st.markdown("""
    <style>
    html, body, [class*="css"] {
        font-family: 'Helvetica Neue', sans-serif;
        background-color: #f4f6f9;
        color: #1f2937;
    }
    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
        max-width: 1200px;
    }
    h1, h2, h3 {
        color: #2c5282;
    }
    .css-1d391kg, .css-1cypcdb {
        background-color: #e2e8f0;
        padding: 2rem 1rem;
        border-radius: 0 1rem 1rem 0;
    }
    .stButton>button {
        background-color: #2b6cb0;
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.5rem 1.5rem;
        font-weight: bold;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #2c5282;
        color: #fefefe;
    }
    .stTabs [role="tablist"] {
        gap: 0.5rem;
        border-bottom: 2px solid #e2e8f0;
        margin-bottom: 1rem;
    }
    .stTabs [role="tab"] {
        background: #cbd5e0;
        color: #1a202c;
        border-radius: 8px 8px 0 0;
        padding: 0.75rem 1.25rem;
        font-weight: 600;
    }
    .stTabs [aria-selected="true"] {
        background: white;
        color: #2c5282;
        border-bottom: 2px solid white;
    }
    .metric-box {
        background: #ff5733;
        padding: 1.2rem;
        border-radius: 12px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        margin-bottom: 1rem;
    }
    .st-expanderHeader {
        font-weight: 600;
        font-size: 1.1rem;
        color: #2c5282;
    }
    .log-area code {
        background: #f7fafc;
        border-left: 4px solid #2b6cb0;
        padding-left: 1rem;
        display: block;
        margin-bottom: 1rem;
    }
    input, textarea {
        background-color: ##000000 !important;
        border-radius: 0.5rem !important;
    }
    </style>
""", unsafe_allow_html=True)


# ---------------- HEADER ----------------
st.markdown("""
<div style="
    background: linear-gradient(135deg, #cbd5e0, #edf2f7);
    padding: 2rem;
    border-radius: 1rem;
    box-shadow: 0 8px 20px rgba(0,0,0,0.08);
    margin-bottom: 1rem;
    text-align: center;
">
    <h1 style='font-size: 2.7rem; color: #1a365d; margin-bottom: 0.5rem;'>🧬 RNA Conservation & Structure Studio</h1>
    <p style='font-size: 1.15rem; color: #2d3748; max-width: 800px; margin: auto;'>Visualize RNA secondary structures, cluster alignments, explore conservation, and discover functional linkers.</p>
</div>
""", unsafe_allow_html=True)

# ---------------- SIDEBAR ----------------
with st.sidebar:
    st.header("Upload MSA File")
    uploaded_file = st.file_uploader("Formats: FASTA (.fasta, .fa) or Clustal (.aln)", type=["fasta", "fa", "aln"])
    st.markdown("---")
    st.markdown("### About")
    st.markdown("""
Upload your multiple sequence alignment (MSA) file to enable structural clustering, visualization, and conservation scoring.
The linker analysis tool is available without uploading data.
    """)
    st.markdown("---")

#Create HTML and PDF
# def generar_pdf_desde_html(html_content):
#     with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmpfile:
#         HTML(string=html_content).write_pdf(tmpfile.name)
#         return tmpfile.name
#Images for report


def plot_delta_mfe(results, output_dir="report_images", bins=20):
    """
    Genera un histograma de la distribución de ΔMFE (mfe_constrained – mfe_unconstrained)
    y lo guarda en output_dir/delta_mfe.png.
    """
    # Extraer deltas
    deltas = [res['mfe_2'] - res['mfe_1'] for res in results]
    if not deltas:
        return None  # nada que plotear

    # Asegurar que exista la carpeta
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "delta_mfe.png")

    # Crear figura
    plt.figure(figsize=(6,4))
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
    Lee un fichero PNG y devuelve un data URI para incrustar en HTML.
    """
    # with open(png_path, "rb") as f:
    #     b64 = base64.b64encode(f.read()).decode("ascii")
    # return f"data:image/png;base64,{b64}"

import matplotlib.pyplot as plt
import os

def plot_pairings_histogram(results, rna1, rna3, output_dir="report_images", bins=None):
    """
    Genera un histograma del número de emparejamientos entre RNA1 y RNA3
    para cada diseño en `results`, y lo guarda en output_dir/pairings_count.png.
    
    `bins`: número de barras; si es None, matplotlib elige automáticamente.
    """
    # Recalcular conteo de emparejamientos para cada diseño
    counts = [
        count_rna1_rna3_pairings(
            res['structure_unconstrained'],  # ó structure_constrained, según prefieras
            rna1, rna3, res['linker']
        )
        for res in results
    ]
    if not counts:
        return None

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "pairings_count.png")

    plt.figure(figsize=(6,4))
    plt.hist(counts, bins=bins, edgecolor='black')
    plt.xlabel("Número de emparejamientos RNA1–RNA3")
    plt.ylabel("Número de diseños")
    plt.title("Distribución de emparejamientos entre RNA1 y RNA3")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    return out_path

def structure_distance(s1, s2):
    # Simple distancia Hamming entre dos estructuras secundarias
    return sum(a != b for a, b in zip(s1, s2))

def plot_structure_dendrogram(results, output_path="report_images/structure_dendrogram.png"):
    # Extraer estructuras constrained
    structures = [res['structure_constrained'] for res in results]

    # Comprobar que todas tengan la misma longitud (necesario para Hamming)
    assert all(len(s) == len(structures[0]) for s in structures), "Las estructuras deben tener la misma longitud"

    # Crear matriz de distancias por Hamming
    dist_matrix = squareform(pdist(np.array(structures)[:, None], 
                                   lambda u, v: structure_distance(u[0], v[0])))

    # Clustering jerárquico
    linked = linkage(dist_matrix, method='average')

    # Graficar dendrograma
    plt.figure(figsize=(10, 6))
    dendrogram(linked,
               labels=[res['linker'] for res in results],
               orientation='top',
               distance_sort='descending',
               show_leaf_counts=True)
    plt.title("🔬 Dendrograma de estructuras (con restricción)")
    plt.xlabel("Linker")
    plt.ylabel("Distancia estructural (Hamming)")
    plt.tight_layout()

    # Guardar imagen
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close()
#Write Report


# Assuming png_to_data_uri is defined elsewhere and accessible

def construir_html_reporte_completo(
    results, report,
    rna1, rna3, struct1, constraint,
    mutable_rna1, watched_positions, use_mutaciones,
    mfe_delta, max_pairings, max_changes, num_mut,
    linker_min, linker_max, verbose, cluster_labels,
    search_method, # Added search_method parameter
    image_output_dir="estructura_img", # General output directory for images
    representative_img_bases=None 
):
    # Paths to your PNGs for the global plots
    delta_uri       = png_to_data_uri("report_images/delta_mfe.png")
    pairings_uri    = png_to_data_uri("report_images/pairings_count.png")
    dendo_uri       = png_to_data_uri("report_images/structure_dendrogram.png")

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    total_found = len(results)
    # Ensure total_clusters correctly counts non-noise clusters
    total_clusters = len(set(cluster_labels)) - (1 if -1 in set(cluster_labels) else 0) if (cluster_labels is not None and cluster_labels.size > 0) else 0

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
          <li>Mutations in RNA1: {"Yes" if use_mutaciones else "No"}</li>
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

        html_secciones_img = ""
        if representative_img_bases:
            html_secciones_img += f"""
            <div class="section">
                <h2>🖼️ Representative Structures by Cluster</h2>
                <p>Below are visual examples of the 'off' (unconstrained) and 'on' (constrained) structures for each identified cluster.</p>
                """
            for i, img_base in enumerate(representative_img_bases):
                off_uri = png_to_data_uri(f"propuestas/linker_{linker_min}/{img_base}_unconstrained_plot.png")
                on_uri = png_to_data_uri(f"propuestas/linker_{linker_min}/{img_base}_constrained_plot.png")

                html_secciones_img += f"""
                <div style="display:inline-block; margin-right:20px; vertical-align:top;">
                  <h3>Cluster {i}</h3>
                  <p><em>Off (unconstrained)</em></p>
                  <img class="img-preview" src="{off_uri}" alt="Unconstrained structure - Cluster {i}" />
                  <p><em>On (constrained)</em></p>
                  <img class="img-preview" src="{on_uri}" alt="Constrained structure - Cluster {i}" />
                </div>
                """
            html_secciones_img += "</div>"

        html += html_secciones_img

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
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as e:
        st.error(f"Command failed:\n`{cmd}`\n{e.stderr}")
        return None


def stop_requested():
    return st.session_state.get("stop_requested", False)

log_area = st.empty()

def update_log(msg):
    current_text = st.session_state.get("log_text", "")
    st.session_state["log_text"] = current_text + msg + "\n"
    log_area.code(st.session_state["log_text"], language="text")

def tab1_analysis(msa, msa_name):
    import shutil, os
    from matplotlib import pyplot as plt

    st.header(" First Sequence Analysis & Conservation")
    sequence = msa[0].replace("-", "")
    structure, energy = predict_secundary_structure(sequence)
    st.code(f"{structure} ({energy:.2f} kcal/mol)")

    conservation_scores = calculate_conservation(msa)
    msa_len = len(msa[0])
    col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
    colors = {
        "always": (255, 0, 0),
        "except2": (139, 0, 0),
        "except4": (0, 0, 255),
        "others": (255, 255, 0)
    }

    output_dir_first = f"rna_outputs_{msa_name}_first"
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

def tab2_clustering(msa, msa_name):
    import shutil, os, time
    st.header(" Full Structural Clustering")

    if st.button("Run clustering and analysis"):
        output_dir = f"rna_outputs_{msa_name}"
        shutil.rmtree(output_dir, ignore_errors=True)
        os.makedirs(output_dir)

        structures = [predict_secundary_structure(seq.replace("-", ""))[0] for seq in msa]
        start = time.time()
        rmsd_matrix = matriz_rmsd(structures)
        labels = cluster_structures(rmsd_matrix, eps=3.0, min_samples=1)
        cluster_summary = organise_per_clusters(msa, structures, labels, rmsd_matrix, output_base="rna_clusters")
        end = time.time()

        times = [(end - start) / len(structures)] * len(structures)
        metrics = compute_metrics(structures, times)
        visualise_metrics(metrics)

        st.subheader("📈 Evaluation Metrics")
        st.markdown(f"""
            <div class="metric-box">
                <b>Total Structures:</b> {metrics["total_structures"]}<br>
                <b>Unique Structures:</b> {metrics["unique_structures"]}<br>
                <b>Total Time:</b> {metrics["total_time"]:.3f} sec
            </div>
        """, unsafe_allow_html=True)

        st.subheader("📊 Cluster Summary")
        cluster_img_dir = "cluster_images"
        os.makedirs(cluster_img_dir, exist_ok=True)

        # Para usar sequence, col_categories y colors, calculamos localmente:
        sequence = msa[0].replace("-", "")
        msa_len = len(msa[0])
        col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
        colors = {
            "always": (255, 0, 0),
            "except2": (139, 0, 0),
            "except4": (0, 0, 255),
            "others": (255, 255, 0)
        }

        for cluster_id, cluster_data in cluster_summary.items():
            with st.expander(f"Cluster {cluster_id} — {cluster_data['num_variants']} structures"):
                rep = cluster_data['representative_structure']
                st.markdown(f"** Representative Structure:** `{rep}`")

                fold_path = os.path.join(cluster_img_dir, f"cluster_{cluster_id}.fold")
                with open(fold_path, "w") as f:
                    f.write(sequence + "\n" + rep + "\n")

                pre_command = generate_pre_command_per_sequence(msa[0], col_categories, colors, lw=10)
                rnaplot_cmd = f'RNAplot --pre "{pre_command}" < {fold_path}'
                run_command(rnaplot_cmd)

                ps_path = os.path.join(cluster_img_dir, f"cluster_{cluster_id}.ps")
                png_path = os.path.join(cluster_img_dir, f"cluster_{cluster_id}.png")
                if os.path.exists("rna.ps"):
                    os.rename("rna.ps", ps_path)
                    run_command(f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
                                f'-dEPSCrop -sOutputFile={png_path} {ps_path}')
                    if os.path.exists(png_path):
                        st.image(png_path, caption=f"Structure Cluster {cluster_id}")

                for var in cluster_data.get("members", []):
                    if var["type"] == "variante":
                        st.code(f"[{var['index']}] {var['structure']}")

        if os.path.exists("rna_clusters/diversity.png"):
            st.image("rna_clusters/diversity.png", caption="Structural Diversity")


def generar_y_mostrar_estructuras(msa, output_dir="cluster_images"):
    msa_len = len(msa[0])
    os.makedirs(output_dir, exist_ok=True)
    col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}
    colors = {
        "always": (255, 0, 0),
        "except2": (139, 0, 0),
        "except4": (0, 0, 255),
        "others": (255, 255, 0)
    }
    for i, msa_seq in enumerate(msa):
        sequence = msa_seq.replace("-", "")
        structure = predict_secundary_structure(sequence)[0]

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

def linker_finder_tab():
    st.markdown("""
    <style>
    /* Fuente y fondo general */
    html, body, [class*="css"] {
        font-family: 'Inter', 'Helvetica Neue', sans-serif;
        background-color: #f9fafb;
        color: #1f2937;
    }

    /* Contenedor principal */
    .block-container {
        max-width: 900px;
        margin: auto;
        padding: 2rem 2rem 4rem 2rem;
    }

    /* Títulos */
    h1, h2, h3 {
        color: #2563eb; /* azul profesional */
        font-weight: 700;
        letter-spacing: 0.03em;
    }

    /* Secciones con borde y sombra suave */
    .section-box {
        background: white;
        border-radius: 1rem;
        padding: 1.8rem 2.4rem;
        margin-bottom: 1.8rem;
        box-shadow: 0 3px 12px rgba(100, 116, 139, 0.12);
        transition: box-shadow 0.3s ease;
    }
    .section-box:hover {
        box-shadow: 0 6px 20px rgba(59, 130, 246, 0.3);
    }

    /* Labels personalizados */
    label {
        font-weight: 600;
        font-size: 1rem;
        color: #374151;
        margin-bottom: 0.3rem;
        display: block;
    }

    /* Inputs y textareas */
    textarea, input[type="text"], input[type="number"] {
        background-color: #f3f4f6;
        border: 1.8px solid #d1d5db;
        border-radius: 0.6rem;
        padding: 0.6rem 1rem;
        font-size: 1rem;
        transition: border-color 0.3s ease;
        width: 100%;
        box-sizing: border-box;
        font-family: 'Inter', sans-serif;
        color: #111827;
        resize: vertical;
        min-height: 38px;
    }
    textarea:focus, input[type="text"]:focus, input[type="number"]:focus {
        outline: none;
        border-color: #2563eb;
        background-color: #eff6ff;
        box-shadow: 0 0 8px rgb(37 99 235 / 0.4);
    }

    /* Botones */
    .stButton>button {
        background-color: #2563eb;
        color: white;
        border-radius: 0.75rem;
        padding: 0.65rem 1.6rem;
        font-size: 1.1rem;
        font-weight: 700;
        box-shadow: 0 3px 10px rgba(37, 99, 235, 0.4);
        transition: background-color 0.3s ease, box-shadow 0.3s ease;
        width: 100%;
        cursor: pointer;
    }
    .stButton>button:hover {
        background-color: #1e40af;
        box-shadow: 0 6px 20px rgba(30, 64, 175, 0.6);
    }

    /* Checkbox */
    div[role="checkbox"] {
        margin-right: 0.5rem;
    }
    label[for^="checkbox"] {
        font-weight: 600;
        color: #374151;
    }

    /* Expander header */
    .st-expanderHeader {
        font-weight: 700 !important;
        font-size: 1.15rem !important;
        color: #2563eb !important;
    }

    /* Contenedor log output */
    .log-area code {
        background: #f1f5f9;
        border-left: 5px solid #2563eb;
        padding: 1rem 1.5rem;
        display: block;
        max-height: 220px;
        overflow-y: auto;
        font-size: 0.9rem;
        font-family: 'Source Code Pro', monospace;
        border-radius: 0 0.6rem 0.6rem 0;
        white-space: pre-wrap;
        color: #1e293b;
        margin-top: 1rem;
        box-shadow: inset 0 0 6px rgba(0,0,0,0.05);
    }

    /* Divisiones */
    hr {
        border: none;
        border-top: 1.5px solid #e5e7eb;
        margin: 2rem 0;
    }

    /* Columns para botones y checkbox */
    .stColumns > div {
        padding-right: 0.8rem;
    }

    /* Placeholder estilo */
    ::placeholder {
        color: #9ca3af;
        opacity: 1;
        font-style: italic;
    }
    </style>
    """, unsafe_allow_html=True)
    st.header(" Linker Finder")

    # Contenedor principal


def linker_finder_tab():
    st.header(" Linker Finder")

    # Contenedor principal
    with st.container():
        # Sección Inputs RNA
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
        
            use_mutaciones = st.checkbox("Allow mutations on RNA1", value=True)
        
            # Parsear las posiciones mutables para calcular rango mutaciones
            mutable_set = [int(x.strip()) for x in mutable_rna1.split(",") if x.strip().isdigit()]
            max_mutations_allowed = len(mutable_set) if mutable_set else 0
        
            if use_mutaciones and max_mutations_allowed > 0:
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

        st.markdown("---")
        # Botones Ejecutar y Detener en columnas
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

        # Variables iniciales
        results = []
        report = ""
        labels = []
        if run_button:
            
            st.session_state.stop_requested = False  # Reinicia al ejecutar
            st.session_state["log_text"] = ""        # Limpia logs anteriores
            output_container.code("", language="text")  # Limpia el contenedor

            if not rna1.strip() or not rna3.strip() or not struct1.strip():
                st.warning("Please, introduce RNA1 and RNA3 sequences and targeted structure")
                st.stop()

            try:
                mutable_set = [int(x.strip()) for x in mutable_rna1.split(",") if x.strip().isdigit()]
                watched_positions = [int(x.strip()) for x in watched_positions_str.split(",") if x.strip().isdigit()]
                
                if search_method == "Brute Force":
                    linker_lengths = range(linker_min, linker_max + 1)
                    with st.spinner("Searching for linkers (Brute Force)..."):
                        results, report_lines = run_linker_search(
                            rna1=rna1.strip(),
                            rna3=rna3.strip(),
                            struct1=struct1.strip(),
                            constraint=constraint.strip(),
                            mutable_rna1=mutable_set,
                            watched_positions=watched_positions,
                            use_mutaciones=use_mutaciones,
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
                        # Ensure run_genetic_algorithm_search is imported from genetic_algorithm.py

                        # The wrapper function already handles internal logging if verbose=True
                        # and prepares the results, report, and cluster labels.
                        results, report_lines, labels = run_genetic_algorithm_search(
                            rna1=rna1.strip(),
                            rna3=rna3.strip(),
                            struct1=struct1.strip(),
                            constraint=constraint.strip(),
                            mutable_rna1=mutable_set,
                            watched_positions=watched_positions,
                            use_mutaciones=use_mutaciones,
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
                # ADD OR ENSURE THESE TWO LINES ARE PRESENT:
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
            img_path = plot_delta_mfe(results)
            # Pass the correct RNA1 sequence (mutated for GA) to plot_pairings_histogram if it uses it.
            rna1_for_plotting = results[0].get('rna1_mutated_seq', rna1) if search_method == "Genetic Algorithm" and results else rna1
            img_pairings = plot_pairings_histogram(results, rna1_for_plotting, rna3)
    
            st.success(f"✅ {len(results)} valid linkers found.")

            # --- Display overall plots if they generated paths ---
            if img_path and os.path.exists(img_path):
                st.image(img_path, caption="Delta MFE Distribution")
            if img_pairings and os.path.exists(img_pairings):
                st.image(img_pairings, caption="RNA1-RNA3 Pairings Histogram")
            # ----------------------------------------------------

            st.markdown("---")
            st.header("Candidate Linker Solutions")
            labels = list(labels)
            if labels and len(set(labels)) > 1: # If clustering was performed and there's more than one cluster
                st.subheader("Representative Linkers (by Cluster)")
                cluster_dict = defaultdict(list)
                for idx, label in enumerate(labels):
                    cluster_dict[label].append((idx, results[idx]))

                for cluster_id, cluster_items in sorted(cluster_dict.items()):
                    if cluster_id == -1: # Noise cluster
                        st.markdown(f"## Noise Cluster ({len(cluster_items)} sequences)")
                    else:
                        st.markdown(f"## Cluster {cluster_id} ({len(cluster_items)} sequences)")

                    # Choose representative: e.g., the one with the best constrained MFE
                    representative = min(cluster_items, key=lambda x: x[1]['mfe_2']) # menor es mejor
                    idx, res = representative[0], representative[1] # Unpack idx and result dict

                    with st.expander(f" Representative Linker (#{idx + 1}): {res['linker']} (Cluster {cluster_id})"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Unconstrained structure:** `{res['structure_unconstrained']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**RNA1 mutations:** {res['mut1_info']}")

                        image_paths = res.get('image_paths', {}) # Safely get the dictionary of paths
                        if image_paths:
                            if 'unconstrained' in image_paths and os.path.exists(image_paths['unconstrained']):
                                st.image(image_paths['unconstrained'], caption=f"Unconstrained Structure for Linker {res['linker']}")
                            else:
                                st.warning(f"Unconstrained image not found for Rep. Linker {res['linker']}: {image_paths.get('unconstrained', 'N/A')}")
                    
                            if 'constrained' in image_paths and os.path.exists(image_paths['constrained']):
                                st.image(image_paths['constrained'], caption=f"Constrained Structure for Linker {res['linker']}")
                            else:
                                st.warning(f"Constrained image not found for Rep. Linker {res['linker']}: {image_paths.get('constrained', 'N/A')}")
                        else:
                            st.info(f"No image paths stored for Representative Linker {res['linker']}.")
            else: # No clustering or only one result/cluster (or all in noise)
                st.subheader("Best Candidate Linker (No Clustering)")
                if results:
                    res = results[0] # Assuming results are sorted by fitness, so first is best
                    with st.expander(f" Top Linker: {res['linker']}"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Unconstrained structure:** `{res['structure_unconstrained']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**RNA1 mutations:** {res['mut1_info']}")

                        # >>> Display images for the best linker <<<
                        image_paths = res.get('image_paths', {})
                        if image_paths:
                            if 'unconstrained' in image_paths and os.path.exists(image_paths['unconstrained']):
                                st.image(image_paths['unconstrained'], caption=f"Unconstrained Structure for Linker {res['linker']}")
                            else:
                                st.warning(f"Unconstrained image not found for Best Linker {res['linker']}: {image_paths.get('unconstrained', 'N/A')}")
                            if 'constrained' in image_paths and os.path.exists(image_paths['constrained']):
                                st.image(image_paths['constrained'], caption=f"Constrained Structure for Linker {res['linker']}")
                            else:
                                st.warning(f"Constrained image not found for Best Linker {res['linker']}: {image_paths.get('constrained', 'N/A')}")
                        else:
                            st.info(f"No image paths stored for Best Linker {res['linker']}.")

            # --- Dendrogram (moved here for better flow after main candidate display) ---
            if labels and len(set(labels)) > 1: # Only plot if clustering actually happened and there's more than one cluster
                plot_structure_dendrogram(results) # This might not be meaningful for a single GA result
            # -------------------------------------------------------------------------

            st.markdown("---")
            show_all = st.checkbox("Show all linkers")
            if show_all:
                st.subheader("All Valid Linkers Found")
                for idx, res in enumerate(results): # Iterate through ALL results
                    linker_name_for_expander = res.get('linker', f"Linker {idx + 1}")
                    with st.expander(f" Linker {idx + 1}: {linker_name_for_expander} (Fitness: {res.get('fitness', 'N/A'):.2f})"):
                        st.markdown(f"**Sequence:** `{res['sequence']}`")
                        st.markdown(f"**Unconstrained structure:** `{res['structure_unconstrained']}`")
                        st.markdown(f"**Constrained structure:** `{res['structure_constrained']}`")
                        st.markdown(f"**MFE (unconstrained):** {res['mfe_1']:.2f} kcal/mol")
                        st.markdown(f"**MFE (constrained):** {res['mfe_2']:.2f} kcal/mol")
                        if res.get("mut1_info"):
                            st.markdown(f"**Mutations on RNA1:** {res['mut1_info']}")

                        # >>> Display images using the paths stored in 'image_paths' for EACH linker <<<
                        image_paths = res.get('image_paths', {}) # Get the dictionary of paths for this linker
                
                        if image_paths: # Check if the dictionary is not empty
                            if 'unconstrained' in image_paths and os.path.exists(image_paths['unconstrained']):
                                st.image(image_paths['unconstrained'], caption=f"Unconstrained Structure for Linker {res['linker']}")
                            else:
                                st.warning(f"Unconstrained image not found for Linker {res['linker']}: {image_paths.get('unconstrained', 'N/A')}")
                    
                            if 'constrained' in image_paths and os.path.exists(image_paths['constrained']):
                                st.image(image_paths['constrained'], caption=f"Constrained Structure for Linker {res['linker']}")
                            else:
                                st.warning(f"Constrained image not found for Linker {res['linker']}: {image_paths.get('constrained', 'N/A')}")
                        else:
                            st.info(f"No image paths stored for Linker {res['linker']}. Images might not have been generated or stored correctly.")

            st.markdown("### 📝 Full Search Report:")
            st.text(report) # Display the full report text generated by the search function


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
            html = construir_html_reporte_completo(
                results=st.session_state["results"],
                report=st.session_state["report"], # Passing the raw log to be included in the HTML report
                rna1=current_rna1_for_report,
                rna3=rna3,
                struct1=struct1,
                constraint=constraint,
                mutable_rna1=mutable_set,
                watched_positions=watched_positions,
                use_mutaciones=use_mutaciones,
                mfe_delta=mfe_delta,
                max_pairings=max_pairings,
                max_changes=max_changes,
                num_mut=current_num_mut_for_report,
                linker_min=current_linker_min_for_report,
                linker_max=current_linker_max_for_report,
                verbose=verbose,
                cluster_labels=st.session_state["cluster_labels"],
                representative_img_bases=representative_img_bases, # This might need adjusting if your HTML report needs actual image paths
                search_method=search_method # Now included
            )
# ---------------- TABS PRINCIPALES ----------------
tabs = st.tabs([" Linker Finder", " MSA-Based Analyses"])

# ---------------- TAB LINKER: SIEMPRE DISPONIBLE ----------------
with tabs[0]:
    linker_finder_tab()

# ---------------- TAB MSA ----------------
with tabs[1]:
    if not uploaded_file:
        st.info("📂 Please upload a multiple sequence alignment (MSA) file to access structural and conservation analysis.")
        st.stop()
    if uploaded_file is not None:
        msa_name = uploaded_file.name.split('.')[0]  # nombre sin extensión
    else:
        msa_name = "default_name"  # o algún valor por defecto

    # PARSEAR MSA
    msa = parse_fasta_msa(uploaded_file)  # Usa tu función de parseo real aquí
    if len(msa) < 2:
        st.error("The file must contain at least two sequences.")
        st.stop()

    # SUBTABS MSA
    sub_tabs = st.tabs([
        " First Sequence Analysis", 
        " Full Structural Clustering", 
        " Visualize All Structures"
    ])

    with sub_tabs[0]:
        tab1_analysis(msa, msa_name)

    with sub_tabs[1]:
        tab2_clustering(msa, msa_name)


    with sub_tabs[2]:
        generar_y_mostrar_estructuras(msa, output_dir=msa_name)


