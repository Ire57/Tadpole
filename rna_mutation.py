import itertools
from rna_structures import get_pair_table

def mutate_sequence(seq, mutable_positions, n_mut):
    """
    Generates all possible mutated sequences by applying 'n_mut' changes
    at specified mutable positions.

    This function explores all combinations of 'n_mut' positions from the
    'mutable_positions' list and all possible nucleotide replacements at these
    chosen positions. It yields each unique mutated sequence along with the
    details of the mutations applied.

    :param seq: The original nucleotide sequence (str).
    :param mutable_positions: A list of 0-indexed positions within the sequence
                              where mutations are allowed (list of int).
    :param n_mut: The exact number of mutations to apply to the sequence (int).

    :yields: A tuple containing:
             - A new string representing the mutated sequence.
             - A list of tuples, where each inner tuple describes a mutation as
               (position, new_base). Note: The original base is not directly
               captured here, as it's a mutation *from* the original sequence.
               Consider if you need (position, old_base, new_base) here.
               (tuple: (str, list of tuple))
    """
    bases = ['A', 'U', 'G', 'C']
    # Iterate through all combinations of 'n_mut' positions from the mutable_positions list
    for positions in itertools.combinations(mutable_positions, n_mut):
        # For each combination of positions, iterate through all possible nucleotide replacements
        # (A, U, G, C) for those 'n_mut' positions.
        for replacements in itertools.product(bases, repeat=n_mut):
            seq_list = list(seq) # Create a mutable list from the original sequence
            
            # Flag to check if any chosen replacement is the same as the original base
            mutation_is_identical = False
            for pos, base in zip(positions, replacements):
                if seq_list[pos] == base:
                    # If a selected new base is identical to the current base at that position,
                    # this combination of mutations isn't truly unique/different, so skip it.
                    mutation_is_identical = True
                    break
                seq_list[pos] = base # Apply the mutation to the sequence list
            
            # If no mutation was identical to the original base at its position, yield the result
            if not mutation_is_identical:
                # The 'else' clause for the inner 'for' loop means it executes if 'break' was NOT called.
                # This ensures we only yield if all mutations in the current combination are actual changes.
                yield ''.join(seq_list), list(zip(positions, [seq[p] for p in positions], replacements)) # Modified to include old base

def mutate_with_complement(seq, struct, mutations_list):
    """
    Applies a given set of mutations to a nucleotide sequence and, if applicable,
    mutates their complementary partners to preserve base pairing potential.

    This function is designed for targeted mutations where maintaining the
    secondary structure's integrity (specifically, specific base pairs) is
    important. It leverages a pairing table derived from the RNA structure
    to identify and modify complementary bases.

    :param seq: The original or current nucleotide sequence as a string (str).
    :param struct: The dot-bracket string representing the secondary structure
                   of the 'seq'. This structure is used to determine base pairs
                   for complementary mutations (str).
    :param mutations_list: A list of tuples, where each tuple specifies a mutation.
                           The format of each tuple is (position, new_base), or
                           if coming from `mutate_sequence` modified above, it could be
                           (position, old_base, new_base). For this function's logic,
                           only (position, new_base) is strictly needed (list of tuple).

    :returns: The mutated sequence as a string, with complementary bases potentially
              also modified to preserve base pairing (str).
    """
    seq_list = list(seq) # Convert sequence to a mutable list for modifications
    pt = get_pair_table(struct) # Get the pairing table from the structure
    complement_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'} # Map for canonical base pairing

    # Iterate through each mutation provided in the list
    for mutation_item in mutations_list:
        # The mutation_item can be (pos, new_base) or (pos, old_base, new_base)
        # We need the position and the new base to apply.
        pos = mutation_item[0]
        new_base = mutation_item[-1] # Safely get the new base (last element)
        
        seq_list[pos] = new_base # Apply the primary mutation
        
        # Check if the mutated position is paired in the structure
        # ptable returns 1-based index of partner, 0 if unpaired.
        paired_pos_1_based = pt[pos + 1] 
        
        if paired_pos_1_based != 0: # If the position is paired
            paired_pos_0_based = paired_pos_1_based - 1 # Convert to 0-based index for Python list
            
            # Ensure the paired position is valid (within sequence bounds and not self-pairing)
            # The condition `paired_pos_0_based > pos` means we only consider the "forward" pair
            # to avoid double-mutating the same pair if it's already in mutations_list.
            # `pair > i` (or `paired_pos_0_based > pos` here) might be intended to avoid
            # redundant processing or certain types of pairs, but a more robust check
            # for general complementary mutation would just be `paired_pos_0_based != pos`.
            # I will keep your original `pair > i` logic equivalent for consistency.
            if paired_pos_0_based > pos: # Ensure it's a valid pair to a *different* downstream position
                complementary_base = complement_map.get(new_base) # Get the required complement for the new base
                
                # If a complement exists and the base at the paired position is not already the correct complement
                if complementary_base and seq_list[paired_pos_0_based] != complementary_base:
                    seq_list[paired_pos_0_based] = complementary_base # Mutate the complementary base

    return ''.join(seq_list) # Return the final mutated sequence as a string