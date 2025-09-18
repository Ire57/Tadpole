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
            if app_functions.run_command(rnaplot_cmd):
                try:
                    ps_file = os.path.join(output_dir_first, "rna_structure_0.ps")
                    png_file = os.path.join(output_dir_first, "rna_structure_0.png")
                    os.rename("rna.ps", ps_file)
                    gs_cmd = (
                        f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
                        f'-dEPSCrop -sOutputFile={png_file} {ps_file}'
                    )
                    app_functions.run_command(gs_cmd)
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
            app_functions.run_command(rnaplot_cmd)

            ps_path = os.path.join(output_dir_first, "pairing_conservation.ps")
            png_path = os.path.join(output_dir_first, "pairing_conservation.png")
            if os.path.exists("rna.ps"):
                os.rename("rna.ps", ps_path)
                app_functions.run_command(f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
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

                    app_functions.run_command(f'RNAplot --pre "{pre_command}" < {fold_file}')

                    ps_file = os.path.join(output_dir, f"images_{i}.ps")
                    png_file = os.path.join(output_dir, f"cimages_{i}.png")

                    if os.path.exists("rna.ps"):
                        os.rename("rna.ps", ps_file)
                        app_functions.run_command(
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
