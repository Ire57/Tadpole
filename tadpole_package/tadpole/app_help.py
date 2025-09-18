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
            st.graphviz_chart(app_presentation.create_simple_workflow_diagram())
       
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
            app_presentation.example_box("""
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

            app_presentation.transparency_note("""
            Note that RNAfold has certain limitations:<br>
            - It predicts only pseudoknot-free secondary structures.<br>
            - Complex tertiary interactions and non-canonical base pairs are not modelled.<br>
            - Predictions are minimum free energy (MFE) approximations and may not capture all biologically relevant conformations.<br>
            Despite these limitations, RNAfold remains one of the most comprehensive and accessible tools available for local execution.
            """)


            app_presentation.example_box("Stop Codon Readthrough (SCR)  \nWithout SCR, translation usually ends when the ribosome encounters a Stop Codon.")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/NoSCR.png", width=900)

            app_presentation.example_cont_box("However, when an SCR element is present, its structure interacts with the ribosome and makes it skip the Stop Codon, resulting in an elongated protein.")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/SCR.png", width=900)

            app_presentation.example_cont_box("To build a switch, we add an aptamer." \
            "- In the absence of ligand: the aptamer pairs with the SCR element, disrupting its structure. The ribosome stops at the stop codon, producing only protein A." \
            "")
            col1, col2, col3 = st.columns([1, 4, 1])
            with col2:
                st.image("images/OFF.png", width=1200)

            app_presentation.example_cont_box(""
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
            app_presentation.input_block("Mutations") 
            st.markdown("""
                This software allows the user to introduce mutations in the SRE sequence to explore functional variants.  
                If a nucleotide involved in a base pair is mutated, its paired base is automatically mutated to a complementary one, preserving the secondary structure.""")
            # SECIS example
            app_presentation.example_box("""
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
            app_presentation.example_box("""
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
            app_presentation.input_block("Allow Changes on SRE’s Structure") 
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
            st.graphviz_chart(app_presentation.create_ga_diagram())
        st.markdown("""
            Brute Force implementation""")
        col1, col2, col3 = st.columns([0.5, 4, 0.5])
        with col2:
            st.graphviz_chart(app_presentation.create_brute_force_diagram())

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
