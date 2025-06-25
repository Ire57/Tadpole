import numpy as np
from sklearn.cluster import DBSCAN
import os
import json
from PIL import Image
import matplotlib.pyplot as plt

def parse_pairs(structure):
    """
    Parses a dot-bracket RNA secondary structure string and extracts all base pairs.

    This function supports multiple types of pairing symbols (parentheses, square brackets,
    curly braces, and angle brackets) to identify base-paired nucleotides. It uses a stack-based
    approach to correctly match opening and closing symbols and identify the 0-indexed
    positions of the paired bases.

    :param structure: The RNA secondary structure in dot-bracket notation (str).

    :returns: A list of tuples, where each tuple `(i, j)` represents a base pair
              between the 0-indexed positions `i` and `j` in the sequence, with `i < j`.
              The list is sorted by the first element of the pairs (list of tuple).
    """
    stack = []
    pairs = []
    # Defines the mapping between opening and closing pairing symbols
    pair_symbols = {'(': ')', '[': ']', '{': '}', '<': '>'}
    # Separate stacks for each type of opening symbol to handle nested structures correctly
    stacks = {k: [] for k in pair_symbols.keys()}
    
    for i, char in enumerate(structure):
        if char in pair_symbols: # If it's an opening symbol
            stacks[char].append(i) # Push its index onto the corresponding stack
        else: # If it's a closing symbol or an unpaired dot
            for open_sym, close_sym in pair_symbols.items():
                if char == close_sym: # If the current character matches a closing symbol
                    if stacks[open_sym]: # Check if there's a corresponding opening symbol on the stack
                        j = stacks[open_sym].pop() # Pop the index of the opening symbol
                        pairs.append((j, i)) # Record the base pair (open_pos, close_pos)
                    # If the stack for this symbol is empty, it's an unmatched closing symbol,
                    # which indicates an invalid structure or an intended unpaired base.
                    break # Move to the next character after processing the closing symbol
    return sorted(pairs) # Return the sorted list of base pairs


def rmsd_pairs(pairs1, pairs2):
    """
    Calculates the Root Mean Square Deviation (RMSD) between two sets of base pairs.

    This function quantifies the structural similarity between two RNA secondary structures
    by comparing their identified base pairs. It assumes an ordered comparison of pairs
    and computes the squared Euclidean distance between corresponding pair indices.

    :param pairs1: A list of base pairs from the first structure,
                   e.g., `[(i1, j1), (i2, j2), ...]` (list of tuple).
    :param pairs2: A list of base pairs from the second structure,
                   e.g., `[(k1, l1), (k2, l2), ...]` (list of tuple).

    :returns: The RMSD value (float). Returns 0.0 if both lists are empty.
              Returns a large value (1e6) if one list is empty and the other isn't,
              or if the calculated RMSD is not finite (e.g., NaN due to empty diffs).
    """
    # If both lists are empty, RMSD is 0
    if not pairs1 and not pairs2:
        return 0.0

    # If one list is empty and the other isn't, they are completely different
    n = min(len(pairs1), len(pairs2))
    if n == 0:
        return 1e6  # Return a large value to indicate significant difference

    diffs = []
    # Compare pairs up to the length of the shorter list
    for k in range(n):
        p1 = pairs1[k]
        p2 = pairs2[k]
        # Calculate squared Euclidean distance between the pair coordinates
        diffs.append((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    # Calculate RMSD
    rmsd = np.sqrt(np.mean(diffs))
    
    # Handle cases where RMSD might be non-finite (e.g., if diffs is somehow empty, though 'n' prevents this)
    if not np.isfinite(rmsd):
        return 1e6
    return rmsd


def matriz_rmsd(estructuras):
    """
    Computes a symmetric RMSD distance matrix for a list of RNA secondary structures.

    This function calculates the RMSD between every unique pair of structures
    provided, forming an N x N distance matrix. It first parses the base pairs
    for each structure and then uses `rmsd_pairs` for pairwise comparisons.
    Includes print statements for real-time progress monitoring.

    :param estructuras: A list of RNA secondary structure strings in dot-bracket notation (list of str).

    :returns: A square NumPy array representing the symmetric RMSD distance matrix
              between all input structures (numpy.ndarray).
    """
    n = len(estructuras)
    # Parse base pairs for all structures once to avoid redundant computations
    pairs_list = [parse_pairs(s) for s in estructuras]
    matriz = np.zeros((n, n)) # Initialize an N x N matrix with zeros
    
    # Fill the upper triangle (and implicitly the lower due to symmetry)
    for i in range(n):
        for j in range(i + 1, n): # Start 'j' from 'i+1' for unique pairs
            print(f"Comparing structures {i} and {j}:")
            print(f"  Pairs {pairs_list[i]}")
            print(f"  Pairs {pairs_list[j]}")
            d = rmsd_pairs(pairs_list[i], pairs_list[j]) # Calculate RMSD
            print(f"  RMSD: {d}")
            matriz[i, j] = d # Assign to upper triangle
            matriz[j, i] = d # Assign to lower triangle for symmetry
    return matriz


def cluster_structures(rmsd_matrix, eps=100.0, min_samples=1):
    """
    Applies DBSCAN clustering to an RMSD distance matrix of RNA structures.

    DBSCAN (Density-Based Spatial Clustering of Applications with Noise) groups
    together structures that are closely packed together, marking as outliers
    (noise) those structures that lie alone in low-density regions. The clustering
    is performed using the precomputed RMSD distances.

    :param rmsd_matrix: A square NumPy array representing the RMSD distance matrix
                        between RNA structures (numpy.ndarray).
    :param eps: The maximum distance between two samples for one to be considered
                as in the neighborhood of the other (epsilon parameter for DBSCAN).
                Defaults to 3.0 (float).
    :param min_samples: The number of samples (or total weight) in a neighborhood
                        for a point to be considered as a core point. This includes
                        the point itself (minimum samples parameter for DBSCAN).
                        Defaults to 1, effectively treating each point as a core point
                        if `eps` is large enough, or creating many small clusters (int).

    :returns: A NumPy array where each element is the cluster label for the
              corresponding structure in the input matrix. -1 indicates noise points (numpy.ndarray).
    """
    db = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = db.fit_predict(rmsd_matrix)
    return labels


def organise_per_clusters(msa, estructuras, labels, rmsd_matrix, output_base="rna_clusters"):
    """
    Organizes RNA sequences and structures into folders based on their cluster assignments.

    This function processes the clustering results, identifies a representative structure
    for each cluster (the one with the lowest average RMSD to other members of its cluster),
    and saves all structures within their respective cluster directories. It also generates
    a summary JSON file and prepares data for Streamlit visualization.

    :param msa: A list of Multiple Sequence Alignment (MSA) sequences. These sequences
                are often gapped, so hyphens are removed before saving (list of str).
    :param estructuras: A list of RNA secondary structure strings corresponding to `msa` (list of str).
    :param labels: A NumPy array of cluster labels, where each label corresponds to an
                   individual structure in `estructuras`. -1 indicates noise (numpy.ndarray).
    :param rmsd_matrix: The precomputed RMSD distance matrix used for clustering (numpy.ndarray).
    :param output_base: The base directory where all cluster-related output will be saved.
                        Defaults to "rna_clusters" (str).

    :returns: A dictionary summarizing the clusters, including representative structures,
              number of variants, and detailed information about each member (dict).
    """
    # Create base output directory if it doesn't exist
    if not os.path.exists(output_base):
        os.makedirs(output_base)
    
    # Directory for representative structures (for easy access)
    rep_dir = os.path.join(output_base, "representantes")
    os.makedirs(rep_dir, exist_ok=True)

    summary = {} # Stores the final summary of clusters
    grups = {} # Temporarily groups indices by cluster label
    
    # Group indices based on their cluster label
    for idx, label in enumerate(labels):
        grups.setdefault(label, []).append(idx)

    # Process each cluster (or noise group, labeled -1)
    for group, indices in grups.items():
        grupo_dir = os.path.join(output_base, f"grupo_{group}")
        variantes_dir = os.path.join(grupo_dir, "variantes")
        os.makedirs(variantes_dir, exist_ok=True) # Create directories for the current cluster

        # Calculate representative for the current group
        # Create a sub-matrix of RMSD distances only for members of this group
        sub_matriz = rmsd_matrix[np.ix_(indices, indices)]
        # Calculate the average RMSD of each member to all other members in its group
        average = np.mean(sub_matriz, axis=1)
        # The representative is the one with the lowest average RMSD (most central)
        rep_idx_local = np.argmin(average)
        rep_idx = indices[rep_idx_local] # Get the global index of the representative

        # Save representative structure to a .fold file
        seq_rep = msa[rep_idx].replace("-", "") # Remove gaps from MSA sequence
        struct_rep = estructuras[rep_idx]
        rep_file = os.path.join(rep_dir, f"rep_grupo_{group}.fold")
        with open(rep_file, "w") as f:
            f.write(seq_rep + "\n" + struct_rep + "\n")

        variantes_info = [] # List to store info for all members of the current group
        # Save each variant's sequence and structure
        for i in indices:
            seq_var = msa[i].replace("-", "")
            struct_var = estructuras[i]
            var_file = os.path.join(variantes_dir, f"var_{i}.fold")
            with open(var_file, "w") as f:
                f.write(seq_var + "\n" + struct_var + "\n")

            # Store information about each member for the summary
            tipo = "representante" if i == rep_idx else "variante"
            variantes_info.append({
                "index": i,
                "sequence": seq_var,
                "structure": struct_var,
                "type": tipo
            })

        # Add current group's summary to the main summary dictionary
        summary[str(group)] = {
            "representative_index": rep_idx,
            "representative_structure": struct_rep,
            "num_variants": len(indices),
            "members": variantes_info
        }

    # Save the complete cluster summary to a JSON file
    resumen_file = os.path.join(output_base, "resumen_clusters.json")
    with open(resumen_file, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def summary_image_per_group(output_base="rna_clusters"):
    """
    Generates a summary image for each cluster, typically by copying the PNG plot
    of the representative structure into the cluster's directory.

    This function iterates through all identified clusters. For each cluster,
    it attempts to find the pre-generated PNG image of its representative structure
    (which would typically be created by an external plotting tool like RNAplot
    and named consistently). If found, it copies this PNG to the cluster's
    main directory as 'resumen.png' for easy overview. If the PNG is not found,
    it falls back to copying the representative's `.fold` file as a text summary.

    :param output_base: The base directory where cluster-related output (including
                        representative and group folders) is stored. Defaults to "rna_clusters" (str).

    :returns: None. This function performs side effects (file copying).
    """
    rep_dir = os.path.join(output_base, "representantes")
    # Get a list of all directories that represent a cluster group
    groups = [d for d in os.listdir(output_base) if d.startswith("grupo_")]

    for group in groups:
        grupo_path = os.path.join(output_base, group)

        # Extract only the numerical ID of the group (e.g., '15' from 'grupo_15')
        grupo_id = group.replace("group_", "")

        # Construct the expected path for the representative's .fold file and its corresponding .png plot
        rep_fold = os.path.join(rep_dir, f"rep_group_{grupo_id}.fold")
        rep_png = rep_fold.replace(".fold", "_plot.png") # Assumes _plot.png naming convention from RNAplot

        if os.path.exists(rep_png):
            # If the PNG exists, open it and save it as 'resumen.png' in the group's directory
            im = Image.open(rep_png)
            im.save(os.path.join(grupo_path, "summary.png"))
        else:
            # If the PNG is not found, as a fallback, copy the .fold file contents to a .txt file
            if os.path.exists(rep_fold):
                with open(rep_fold) as f:
                    txt = f.read()
                with open(os.path.join(grupo_path, "summary.txt"), "w") as f:
                    f.write(txt)


def compute_metrics(structures, tiempos):
    """
    Calculates key performance metrics for the RNA structure analysis and clustering.

    This function quantifies aspects like the total number of structures processed,
    the number of unique structures found (before clustering), and the total
    computational time spent on the analysis.

    :param structures: A list of RNA secondary structure strings. Used to determine
                        the number of unique structures (list of str).
    :param tiempos: A list of time values (e.g., seconds per folding/evaluation step).
                    Used to calculate total time (list of float).

    :returns: A dictionary containing the calculated metrics:
              - "total_structures": Total number of structures analyzed.
              - "unique_structures": Number of distinct secondary structures found.
              - "total_time": Sum of all provided time values (float).
    """
    total = len(structures) # Total number of structures provided
    unique = len(set(structures)) # Number of unique structure strings
    
    total_time = sum(tiempos) # Sum all individual time measurements

    metrics = {
        "total_structures": total,
        "unique_structures": unique,
        "total_time": total_time
    }

    return metrics

def visualise_metrics(metrics):
    """
    Generates and saves a bar plot visualizing key diversity metrics of the RNA structures.

    This function creates a simple bar chart comparing the total number of structures
    with the number of unique structures found, providing an insight into the diversity
    of the generated or analyzed RNA structures. The plot is saved as a PNG file.

    :param metrics: A dictionary containing calculated metrics, expected to have
                    "total_structures" and "unique_structures" keys (dict).

    :returns: None. This function generates and saves a plot as a side effect.
    """
    # Create a bar plot for diversity
    plt.figure(figsize=(5, 3)) # Set figure size for the plot
    plt.bar(["Total", "Unique"], [metrics["total_structures"], metrics["unique_structures"]], color=["skyblue", "blue"])
    plt.title("Diversity of Structures") # Set plot title
    plt.ylabel("Quantity") # Set y-axis label
    plt.tight_layout() # Adjust plot parameters for a tight layout
    plt.savefig("rna_clusters/diversity.png") # Save the plot to a file