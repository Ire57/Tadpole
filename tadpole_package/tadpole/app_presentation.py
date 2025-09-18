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

