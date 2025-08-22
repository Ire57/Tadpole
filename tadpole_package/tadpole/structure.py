import RNA


def predict_secundary_structure(sequence):
    """
    Predicts the minimum free energy (MFE) secondary structure of an RNA sequence.

    This function uses the ViennaRNA package's `RNA.fold` function to compute
    the most thermodynamically stable secondary structure for a given RNA sequence
    and its corresponding free energy.

    :param sequence: The RNA nucleotide sequence as a string (str).

    :returns: A tuple containing:
              - The predicted secondary structure in dot-bracket notation (str).
              - The minimum free energy (MFE) of the predicted structure in kcal/mol (float).
    """
    secundary_structure, energy = RNA.fold(sequence)
    return secundary_structure, energy


def align_secondary_structure(sequence: str, structure: str) -> str:
    """
    Aligns a given secondary structure string to a sequence that may contain gaps.

    This function is crucial when you have a secondary structure (which typically
    doesn't contain gaps) and a corresponding sequence that might have gaps
    (e.g., from a Multiple Sequence Alignment). It inserts gap characters ('-')
    into the structure string at positions corresponding to gaps in the sequence,
    ensuring that the final aligned structure has the same length as the gapped sequence.

    :param sequence: The RNA nucleotide sequence, potentially containing gaps ('-') (str).
    :param structure: The dot-bracket secondary structure string, which should
                      correspond to the `sequence` *without* its gaps (str).

    :raises ValueError: If the length of the `structure` does not match the
                        length of the `sequence` after removing all gaps from it.

    :returns: The secondary structure string aligned to the gapped sequence,
              with '-' characters inserted at gap positions (str).
    """
    seq_nogap = sequence.replace("-", "")
    # Validate that the structure length matches the gap-free sequence length
    if len(seq_nogap) != len(structure):
        raise ValueError(f"The length of the struture ({len(structure)}) does not match the length "
                         f"of the sequence without gaps ({len(seq_nogap)}).")

    aligned_structure = []
    struct_idx = 0 # Index for iterating through the original structure string
    for nt in sequence:
        if nt == "-":
            aligned_structure.append("-") # Add a gap to the structure if the sequence has one
        else:
            aligned_structure.append(structure[struct_idx]) # Add the structural character
            struct_idx += 1 # Move to the next character in the original structure
    return "".join(aligned_structure) # Join the list back into a string


def adjust_range_structure(structure: str) -> str:
    """
    Placeholder function intended for adjusting ranges or correcting RNA secondary structures.

    This function provides a template for operations that might involve iterating
    through a secondary structure string and potentially modifying it based on
    positional logic. It includes a basic example of checking index bounds to
    prevent `IndexError` when accessing adjacent characters. The current implementation
    does not perform any actual structural modifications.

    :param structure: The RNA secondary structure in dot-bracket notation (str).

    :returns: The (potentially adjusted) secondary structure string (str).
    """
    structure_ready = list(structure) # Convert to list for mutability
    length = len(structure_ready)

    for i in range(length):
        # Example: if you need to access i+1, first check that it's within range
        if i + 1 < length:
            # Add your logic here to modify structure_lista[i] or structure_lista[i+1]
            # For instance, you might change '(' to '.' under certain conditions.
            pass
        # More logic here with similar boundary checks if needed for other indices (e.g., i-1, i+2)

    return "".join(structure_ready) # Join the list back into a string


