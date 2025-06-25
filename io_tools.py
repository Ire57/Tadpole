import os
from PIL import Image
from IPython.display import display

def save_and_plot_structures(seq, structure_unconstr, structure_constr,
                             rna1, linker, rna3, mut1_info, mfe_1, mfe_2,
                             folder_prefix="propuestas"):
    """
    Saves and plots the unconstrained and constrained secondary structures of an RNA switch. 

    This function takes a full RNA sequence, its two corresponding structures (unconstrained
    and constrained), and related information (RNA1, linker, RNA3 sequences, RNA1 mutation info,
    and MFE values). It then performs the following actions:
    1. Creates a dedicated folder structure (e.g., 'propuestas/linker_X').
    2. Saves detailed information (linker length, RNA1 mutations, MFE, structure, sequence)
       for both unconstrained and constrained states into separate text files.
    3. Generates `.fold` files, which are input files for ViennaRNA's RNAplot tool.
    4. Invokes `RNAplot` to generate PostScript (`.ps`) plots of the structures.
       These plots highlight the RNA1 (green), linker (red), and RNA3 (blue) segments
       using `omark` commands.
    5. Converts the generated `.ps` files to `.png` images using Ghostscript (`gs`).
    6. Displays the generated `.png` images directly in the output (e.g., in a Jupyter Notebook
       or IPython environment).

    :param seq: The full concatenated RNA sequence (rna1_mutated + linker + rna3) (str).
    :param structure_unconstr: The dot-bracket secondary structure of the RNA sequence
                               when folded without constraints (OFF state) (str).
    :param structure_constr: The dot-bracket secondary structure of the RNA sequence
                             when folded under specific constraints (ON state) (str).
    :param rna1: The RNA1 sequence (can be original or mutated) (str).
    :param linker: The linker sequence (str).
    :param rna3: The RNA3 sequence (str).
    :param mut1_info: Information about mutations applied to RNA1, typically a list of
                      (position, old_base, new_base) tuples (list).
    :param mfe_1: The Minimum Free Energy (MFE) of the unconstrained structure (float).
    :param mfe_2: The Minimum Free Energy (MFE) of the constrained structure (float).
    :param folder_prefix: The base directory for saving results. A subfolder based on
                          linker length will be created inside it. Defaults to "propuestas" (str).

    :returns: None. This function performs side effects (file saving and image display).
    """
    linker_len = len(linker)
    full_len = len(seq)
    rna1_len = len(rna1)
    
    folder = f"{folder_prefix}/linker_{linker_len}"
    os.makedirs(folder, exist_ok=True)

    for tag, structure, mfe in [
        ("unconstrained", structure_unconstr, mfe_1),
        ("constrained", structure_constr, mfe_2)
    ]:
        filename = f"{folder}/result_{linker}_{tag}.txt"
        with open(filename, "w") as f:
            f.write(f"Linker = {linker_len}\n")
            f.write(f"Mutaciones RNA1: {mut1_info}\n")
            f.write(f"mfe = {mfe:.2f}\n")
            f.write(structure + "\n")
            f.write(seq + "\n")

        # RNAplot
        fold_file = f"{folder}/{linker}_{tag}.fold"
        with open(fold_file, "w") as f:
            f.write(seq + "\n")
            f.write(structure + "\n")

        # Command for RNAplot to color different segments
        # 1-rna1_len (RNA1) -> GREEN
        # rna1_len+1 - rna1_len+linker_len (Linker) -> RED
        # rna1_len+linker_len+1 - full_len (RNA3) -> BLUE
        pre_command = f"1 {rna1_len} 10 GREEN omark " \
                      f"{rna1_len+1} {rna1_len+linker_len} 10 RED omark " \
                      f"{rna1_len+linker_len+1} {full_len} 10 BLUE omark"

        # Execute RNAplot command
        os.system(f'RNAplot --pre "{pre_command}" < {fold_file}')
        
        # Define paths for PostScript and PNG files
        ps_file = f"{folder}/{linker}_{tag}_plot.ps"
        png_file = f"{folder}/{linker}_{tag}_plot.png"
        
        # Move the output 'rna.ps' file to the designated ps_file path
        os.system(f"mv rna.ps {ps_file}")
        
        # Convert PostScript to PNG using Ghostscript
        # -dSAFER: safest operating mode
        # -dBATCH: exit after processing
        # -dNOPAUSE: do not pause between pages
        # -sDEVICE=pngalpha: output to PNG with alpha channel
        # -r150: resolution 150 DPI
        # -dEPSCrop: crop to bounding box
        # -sOutputFile: output file name
        os.system(f"gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 -dEPSCrop -sOutputFile={png_file} {ps_file}")
        
        # Display the generated PNG image in IPython environments (like Jupyter Notebook)
        display(Image.open(png_file))