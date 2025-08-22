import RNA

def get_pairing_dict(structure):
    """
    Parses a dot-bracket RNA secondary structure string and returns a dictionary
    mapping each paired base's index to its partner's index.

    This function uses a stack-based approach to identify matching parentheses
    in the dot-bracket notation, effectively creating a lookup table for base pairs.
    It's useful for quickly finding the partner of any paired base.

    :param structure: The RNA secondary structure in dot-bracket notation (str).

    :returns: A dictionary where keys are 0-indexed positions of paired bases,
              and values are the 0-indexed positions of their respective partners.
              For example, `{0: 9, 9: 0, 1: 8, 8: 1}` (dict).
    """
    stack = []  # Used to store the indices of opening parentheses
    pairs = {}  # Dictionary to store pairing information
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)  # Push opening parenthesis index onto the stack
        elif c == ')':
            # If a closing parenthesis is found, pop the last opening parenthesis's index
            j = stack.pop()
            pairs[i] = j  # Store the pair: closing_idx -> opening_idx
            pairs[j] = i  # Store the pair: opening_idx -> closing_idx
    return pairs

def get_base_pairs(structure):
    """
    Parses a dot-bracket RNA secondary structure string and returns a list of
    tuples, where each tuple represents a base pair.

    This function identifies all canonical base pairs (represented by matching
    parentheses) in the dot-bracket notation. It's useful for iterating through
    pairs or comparing sets of pairs between structures.

    :param structure: The RNA secondary structure in dot-bracket notation (str).

    :returns: A list of tuples, where each tuple `(i, j)` represents a base pair
              between the 0-indexed positions `i` and `j` in the sequence, with `i < j`.
              For example, `[(0, 9), (1, 8)]` (list of tuple).
    """
    stack = []  # Used to store the indices of opening parentheses
    pairs = []  # List to store base pair tuples
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)  # Push opening parenthesis index onto the stack
        elif c == ')':
            # If a closing parenthesis is found, pop the last opening parenthesis's index
            j = stack.pop()
            pairs.append((j, i))  # Store the pair as (opening_idx, closing_idx)
        else:
            continue
    return pairs

def get_pair_table(struct):
    """
    Generates a pairing table from a dot-bracket secondary structure string.

    A pairing table is an array where `ptable[i]` contains the 1-based index of the
    base paired with `i` (if `i` is 1-based) or 0 if `i` is unpaired. This is a
    utility function from the ViennaRNA package, providing a standard representation
    for quick structural lookups.

    :param struct: The RNA secondary structure in dot-bracket notation (str).

    :returns: A NumPy array representing the ViennaRNA pairing table. The array is
              1-indexed, meaning `ptable[i]` gives the partner of base `i` (1-based).
              Unpaired bases are indicated by `0` (numpy.ndarray).
    """
    return RNA.ptable(struct)

def constrained_mfe(seq, constraint):
    """
    Predicts the minimum free energy (MFE) secondary structure of an RNA sequence
    under specific structural constraints.

    This function utilizes the ViennaRNA package's `RNA.fold_compound` to apply
    custom constraints to the folding process. This is essential for modeling
    RNA behavior where certain regions are forced to be paired or unpaired,
    mimicking environmental conditions or interactions.

    :param seq: The RNA nucleotide sequence as a string (str).
    :param constraint: A string representing the structural constraint in dot-bracket
                       notation, where positions to be constrained are specified
                       (e.g., `'((...)).'` or `'.xxx.'` for forced unpaired).
                       The length of the constraint string must match the sequence length (str).

    :returns: The minimum free energy (MFE) of the predicted structure
              under the given constraints, in kcal/mol (float).
    """
    fc = RNA.fold_compound(seq)  # Create a fold_compound object for the sequence
    # Apply the structural constraint
    # RNA.CONSTRAINT_DB_DEFAULT uses default constraint enforcement (e.g., '(' means force paired)
    fc.constraints_add(constraint, RNA.CONSTRAINT_DB_DEFAULT)
    return fc.mfe()  # Calculate and return the MFE under the applied constraints