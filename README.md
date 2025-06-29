# RNA Conservation and Secondary Structure Analysis

## Project Objectives

This web application allows you to:

Analyze multiple sequence alignments (MSA) of RNA sequences to:

- Calculate positional conservation across the alignment.
- Predict and visualize RNA secondary structure.
- Generate colored graphical representations integrating conservation and structure.
- Facilitate biological interpretation through interactive visualizations and quantitative analysis.

Linker Finder: This tool integrates a genetic algorithm to optimize RNA linker sequences for functional RNA switches. It works by:

- Generating an initial population of random linker sequences.

- Evaluating each linker based on predicted RNA folding in OFF and ON states, structural constraints, and energy differences (switch efficiency).

- Selecting top-performing linkers for crossover and mutation to create new candidates.

- Iterating over multiple generations to improve solutions.

- Outputting the best linker sequences with their structural and fitness data.

---
---
## Web Application Access

While Streamlit offers a direct application deployment service, the inherent complexity of our application—specifically its reliance on libraries such as ViennaRNA and Ghostscript, which necessitate compilation beyond standard Python environments—precluded its direct deployment via Streamlit for this phase. Consequently, a demonstration version of the tool has been deployed on **Render** (https://render.com/) to ensure its accessibility for evaluation. It is important to note that this current deployment operates within Render's limited resource model, which may result in longer initial page load times. For optimal performance and to effectively test the application's core functionality, it is recommended to utilise the Genetic Algorithm (GA) option, configuring it with approximately 10 generations and a low energy difference (MFE Delta). It should be noted that a different deployment strategy may be employed for the application's final public release. Access the live demonstration here: \href{https://toolkit-m4g6.onrender.com/}{Access the demo here}

---
## Requirements and Installation
### Requirements
Python 3.8 or higher

ViennaRNA Package (includes RNAplot) installed and accessible from the command line

Ghostscript (gs) installed and accessible from the command line

System libraries required by WeasyPrint (for PDF/HTML rendering)

### Installation

1. Install system dependencies
Ubuntu/Debian:

sudo apt-get update
sudo apt-get install -y python3-pip python3-venv git ghostscript \
    libpango1.0-0 libgdk-pixbuf2.0-0 libcairo2

2. Install ViennaRNA Package
Download and install ViennaRNA following the official instructions:
https://www.tbi.univie.ac.at/RNA/
sudo apt-get install viennarna


3. Clone this repository:

git clone https://github.com/Ire57/Toolkit.git
cd Toolkit

4. Create and activate a Python virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate    # On Windows: venv\Scripts\actívate

5. Install Python dependencies
pip install --upgrade pip
pip install -r requirements.txt


Verify installation

Check ViennaRNA tools:
RNAplot --help


Check Ghostscript:
gs --version

Check Python packages:
python -c "import RNA, streamlit, weasyprint; print('All packages loaded successfully')"




## Quick Start Guide

1. Run the app:

streamlit run app.py

2. Upload an MSA file in the sidebar. Supported formats: FASTA (.fasta, .fa) or Clustal (.aln).

The app will validate the file and display:

- Predicted secondary structure for the first sequence.
- Conservation visualisation mapped onto the structure.
- Optionally, enable viewing structures for all sequences.

3. Explore generated images to analyze sequence conservation and structure.

## Scientific Background

Conservation in RNA Alignments
Conservation measures how consistent bases are at a given position across related RNA sequences. It is calculated using Shannon entropy: lower entropy means higher conservation.

RNA Secondary Structure and RNAplot
RNA folds into structures by base pairing (helices, loops, etc.). Secondary structure is computationally predicted (here we use ViennaRNA).

RNAplot is a tool that produces visual diagrams of RNA secondary structure in PostScript format, illustrating base pairings and folding.

## References and Documentation

ViennaRNA Package: https://www.tbi.univie.ac.at/RNA/

Streamlit: https://streamlit.io/

Conservation and Shannon Entropy:
Schneider TD, Stephens RM. "Sequence logos: a new way to display consensus sequences." Nucleic Acids Res. 1990.

Ghostscript: https://www.ghostscript.com/

Questions or suggestions? Contact: [ireneagudozamora@gmail.com]
