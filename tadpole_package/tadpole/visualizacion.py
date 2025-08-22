import os
import subprocess
from PIL import Image

def generate_rnaplot_with_colours(sequence, secundary_structure, grups, tag="conservation", output_dir="outputs"):
    """
    Generates an RNA secondary structure plot with specific positions colored
    according to predefined conservation categories.

    This function utilizes ViennaRNA's `RNAplot` tool to visualize an RNA structure.
    It takes a sequence, its dot-bracket structure, and a dictionary of grouped
    positions (e.g., by conservation level). It then uses `RNAplot`'s `--pre`
    option to mark and color individual nucleotides based on their assigned group,
    converts the PostScript output to PNG, and returns the image.

    :param sequence: The RNA nucleotide sequence (str).
    :param secundary_structure: The dot-bracket notation of the RNA's secondary structure (str).
    :param grups: A dictionary where keys are category names (e.g., "siempre_conservado")
                   and values are lists of 1-based indices of nucleotides belonging to that category (dict).
    :param tag: A string prefix for output filenames (e.g., "conservation_plot.png"). Defaults to "conservation" (str).
    :param output_dir: The directory where output files (.fold, .ps, .png) will be saved. Defaults to "outputs" (str).

    :returns: A PIL Image object of the generated PNG plot if successful, otherwise None (PIL.Image.Image or None).
    """
    # Define RGB colors for each conservation category
    colours = {
        "always": (230, 0, 0),
        "except2": (139, 0, 0),
        "except4": (0, 0, 255),
        "others": (255, 255, 0)
    }
    category_map = {
        "always_conserved": "always",
        "always_conserved_except_two": "except2",
        "always_conserved_except_4": "except4",
        "others": "others"
    }
    lw = 10 # Line width for marking points on the plot

    # Build the `RNAplot --pre` command string for coloring individual positions
    pre_command = ""
    for category, indices in grups.items():
        color_key = category_map.get(category, "others")
        r, g, b = colours.get(color_key, (255, 255, 0)) # Get color, default to yellow if category not found
        for i in indices:
            # Format: 'position position line_width Red Green Blue omark'
            pre_command += f"{i} {i} {lw} {r} {g} {b} omark "

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define paths for input and output files
    fold_file = os.path.join(output_dir, f"{tag}.fold")
    ps_name = os.path.join(output_dir, f"{tag}_plot.ps")
    png_name = os.path.join(output_dir, f"{tag}_plot.png")

    # Write the sequence and structure to the .fold file, which RNAplot reads
    with open(fold_file, "w") as f:
        f.write(f"{sequence}\n{secundary_structure}\n")

    try:
        # Run RNAplot command using subprocess.
        # `shell=True` allows running the command string directly.
        # `check=True` raises CalledProcessError if the command returns a non-zero exit code.
        subprocess.run(f'RNAplot --pre "{pre_command.strip()}" < {fold_file}', shell=True, check=True)
        
        # RNAplot always outputs to "rna.ps", so rename it to the desired PS filename
        os.rename("rna.ps", ps_name)

        # Construct and run the Ghostscript command to convert PS to PNG
        gs_cmd = [
            "gs", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pngalpha",
            "-r150", "-dEPSCrop", f"-sOutputFile={png_name}", ps_name
        ]
        subprocess.run(gs_cmd, check=True) # Execute Ghostscript command
    except subprocess.CalledProcessError as e:
        print(f"Error during the generation of RNAplot or during conversion with Ghostscript: {e}")
        return None # Return None if there's an error during command execution

    try:
        # Open the generated PNG image using PIL and return it
        img = Image.open(png_name)
        return img
    except Exception as e:
        print(f"Could not open the generated image: {e}")
        return None # Return None if the PNG file cannot be opened


def generate_pre_command_per_sequence(msa_seq, col_categories, colors, lw):
    """
    Generates an `RNAplot --pre` command string to color and label positions
    of an RNA sequence based on their conservation categories.

    This function is designed to create a custom `--pre` command for `RNAplot`
    that takes into account gaps in the MSA. It assigns colors to non-gap positions
    based on their `col_categories` and adds numerical labels for every 10th
    non-gap position. The positions are adjusted to be 1-based, as required by `RNAplot`.

    :param msa_seq: A single sequence from an MSA, potentially containing gaps ('-') (str).
    :param col_categories: A dictionary mapping 0-based MSA indices to their
                           conservation category names (e.g., "siempre_conservado") (dict).
    :param colors: A dictionary mapping simplified category names (e.g., "siempre", "menos2")
                   to RGB color tuples (dict).
    :param lw: The line width for marking positions in the plot (int).

    :returns: A single string suitable for use with `RNAplot --pre` option,
              containing commands to color and label specific nucleotides (str).
    """
    comands = [] # List to build individual RNAplot pre-commands
    idx_sin_gaps = 0 # Counter for 0-based index in the gap-free sequence

    # Iterate through the MSA sequence, considering gaps
    for idx_msa, base in enumerate(msa_seq):
        if base == "-":
            continue # Skip gap positions

        # Get the conservation category for the current MSA position (ignoring gaps)
        # Default to "resto" if the category isn't found (e.g., for noise in clustering)
        category = col_categories.get(idx_msa, "others")

        # Map the full category name to a simplified name used in the `colors` dictionary
        if category == "always_conserved":
            name = "always"
        elif category == "always_conserved_except_two":
            name = "except2"
        elif category == "always_conserved_except_4":
            name = "except4"
        elif category == "not_paired":
            name = "notpaired"
        else:
            name = "others"

        r, g, b = colors[name] # Get RGB color for the current category
        pos = idx_sin_gaps + 1 # Convert to 1-based index for RNAplot

        # Add the 'omark' command to color the current position
        comands.append(f"{pos} {pos} {lw} {r} {g} {b} omark")
        idx_sin_gaps += 1 # Increment the gap-free sequence index

    # Add labels for every 10th position in the gap-free sequence
    for i in range(10, idx_sin_gaps + 1, 10):
        dx, dy = 0.2, 0.2 # Offset for label placement
        comands.append(f'{i} {dx} {dy} ({i}) Label') # RNAplot Label command

    return " ".join(comands) # Join all commands into a single string