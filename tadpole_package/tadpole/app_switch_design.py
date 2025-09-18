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



# --- Content for "Switch Design" sub-tab ---
def switch_design():
    st.markdown("""<div style="font-size: 2.7rem; color: #ea534e; font-weight: bold; margin-top: 1rem;">
        Switch Design  
    </div>
    """, unsafe_allow_html=True)
    # Main container
    with st.container():
        # Inputs
        with st.expander("Essential Inputs", expanded=True):
            rna1 = st.text_area(
                label="Functional RNA Element (FRE) Sequence",
                value='ACCAGUGUGCGGAUGAUAACUACUGACGAAAGAGUCAUCGACUCAGUUAGUGGUUGGAUGUAGUCACAUUAGU',
                height=110,
                placeholder="Introduce the FRE sequence  (A, C, G, U).",
                help="More infromation on https://2025.igem.wiki/barcelona-ub/software"
            )

            rna3 = st.text_area(
                label="Conformational RNA Element (CRE) Sequence",
                value='AGUUGGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACCAAAA',
                height=110,
                placeholder="Introduce the CRE sequence (A, C, G, U).",
                help="More infromation on https://2025.igem.wiki/barcelona-ub/software"
            )

            struct1 = st.text_input(
                label="Targeted Structure of the FRE (dot-bracket)",
                value='((.(((((((......(((((((((((....(((((...))))).)))))))))))......).)))))).))',
                placeholder="Example: ((.((...)))).",
                help="Put here the theoretical structure of the FRE element (determined experimentally or through " \
                "a prediction software, like RNAFold). This is the structure that the FRE sequence is going to be forced to have on the constrained (ON) state"
            )
            watched_positions_str = st.text_input(
                label="Watched positions for FRE",
                value="56",
                placeholder="Ex: 56, 57",
                help="Positions key to functionality for the FRE element whose pairing "
                "(from respect to the input targeted structure) will be forced to change to ensure the OFF state disrupts functionality."
            )
            constraint = st.text_input(
                label="Restriction Chain for the CRE (constraint)",
                value='........................................................................................<<<<<<<...............................',
                placeholder="Example: ...........<<<<<............",
                help="Restrictions for the folding of the constrained (ON) state. It must have the same " \
                "length as the hole sequence. '.' means no restriction, '(' and ')' for corresponding pairs, 'x' " \
                "to avoid pairing that nucleotide and '>' to avoid pairing with later nts and '<' to avoid pairing with previous nucleotides."
            )
            st.markdown('</div>', unsafe_allow_html=True)

        with st.expander(" Advanced Parameters", expanded=False):

            # --- Search Method Selection ---
            search_method = st.radio(
                "Select Search Method",
                ("Genetic Algorithm", "Brute Force"),
                help="More infromation on https://2025.igem.wiki/barcelona-ub/software"
            )

            
            mutable_rna1 = st.text_input(
                label="Mutable positions for SRE",
                value="1, 2, 4, 5, 6, 7, 8, 9",
                placeholder="Ex: 1, 2, 4, 5, 6",
                help="This positions should be nucleotides known for NOT having an impact on functionallity. " \
                "Notice that, if paired on the input targeted structure, the pairs will be mutated aswell to mantain complementarity."
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
                help="The linker is the key element that will be randomised untill one is found such that it makes the ON/OFF system possible"
            )
            linker_max = st.slider(
                label="Linker's maximum length",
                min_value=linker_min,
                max_value=10,
                value=7,
                step=1,
                help=" Linkers of length from the minimum length up to this length will be tested."
            )
            max_changes = st.number_input(
                label="Maximum changes on the SRE structure",
                min_value=0, max_value=20, value=6, step=1,
                help="The input targeted structure might not need to be much restrictive, " \
                "as maybe not the whole structure is important for functionality. Therefore " \
                "the user can indicate how many changes are supported. This might allow more linkers to be found."
            )
            max_pairings = st.number_input(
                label="Maximum number of SRE-aptamer pairings",
                min_value=0, max_value=20, value=10, step=1,
                help="The external signal to the CRE can't brake many pairings. To ensure the switch's functionality, " \
                "one should choose a low number of pairings. However, if this conditions is very restrictive, " \
                "the ON/OFF system might not be possible to construct, and thus it is allowed to modulate this value."
            )
            mfe_delta = st.number_input(
                label="Minimum Energy difference (kcal/mol)",
                min_value=0.0, max_value=20.0, value=4.0, step=0.1,
                help="Difference between the energy of the constrained" \
                "(ON) state and the unconstrained (OFF) state. In the case of aptamers, it should be approximately half " \
                "the binding energy of the aptamer–ligand interaction."
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
        random_seed = st.number_input(
            label="Seed for the Genetic Algorithm",
            value=42,
            step=1,
            help="The Genetic Algorithm (GA) is deterministic. "
                 "Setting a seed ensures reproducibility of the results."
        )

        folder = 'results_switch_design'
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
            # Clean temporary folder 
            if os.path.exists(folder):
                shutil.rmtree(folder)
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
                            tournament_size=ga_tournament_size, random_seed=random_seed
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
                st.warning("No valid systems found. If you are using the GA, please try again with a different seed. If you are using Brute Force, consider changing the parameters. ")

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
            img_path = app_functions.plot_delta_mfe(results, output_dir=folder)
            # Pass the correct RNA1 sequence (mutated for GA) to plot_pairings_histogram if it uses it.
            # For GA, rna1_mutated_seq will be in results[0]['rna1_mutated_seq']
            rna1_for_plotting = results[0].get('rna1_mutated_seq', rna1) if search_method == "Genetic Algorithm" and results else rna1
            img_pairings = app_functions.plot_pairings_histogram(results, rna1_for_plotting, rna3, output_dir=folder)
            
            st.success(f"✅ {len(results)} valid systems found.")

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
            html_report_content = app_functions.build_full_html_report(
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


            igem_type_iis_sites = {
            "BbsI": ["GAAGAC", "GTCTTC"],
            "BsaI": ["GGTCTC", "GAGACC"],
            "BsmBI": ["CGTCTC", "GAGACG"],
            "BspQI": ["GCTCTTC", "GAAGAGC"],
            "BtgZI": ["GCGATG", "CATCGC"],
            "Esp3I": ["CGTCTC", "GAGACG"], # Es la misma que BsmBI
            "PaqCI": ["CACCTGC", "GCAGGTG"],
            "SapI": ["GCTCTTC", "GAAGAGC"] # Es la misma que BspQI
            }
        
            st.markdown("""
            ## Restriction sites
            More information about restriction sites: [iGEM scars database](https://parts.igem.org/partsdb/scars.cgi)  
            """)
            st.markdown("### Considered type IIS restriction enzymes:")
            for enzyme, sites in igem_type_iis_sites.items():
                st.markdown(f"- **{enzyme}** → {', '.join(sites)}")
            for idx, res in enumerate(results):
                full_sequence = res['sequence']
                for enzyme, sites in igem_type_iis_sites.items():
                    for site in sites:
                        if site.upper() in full_sequence.upper():
                            st.markdown(f"¡WARNING, A restriction site was found for the enzime {enzyme} ({site}) in the sequence {res['sequence']}.")
                            break  
            st.markdown("""Restriction sites cheking completed.""")
            
            
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
                    merged_doc = export_single_linker(
                        doc=merged_doc,
                        result=result,
                        index=index
                    )
                
                sbol_bytes = merged_doc.write_string("json-ld").encode("utf-8")
            
                # Show button if there are results
                st.download_button(
                    label="Download All Linkers as SBOL",
                    data=sbol_bytes,
                    file_name="all_linkers.jsonld",
                    mime="application/ld+json"
                )
            else:
                st.info("No results.")
  
            

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
