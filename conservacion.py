from collections import Counter
import math



def calculate_conservation(msa):
    """
    For each column in the MSA, computes the probability of the most frequent character.

    :param msa: A list of strings (equal-length sequences).
    :returns: A list of float values, each representing the probability of the most
              typical character in the corresponding column.
    """
    num_seqs = len(msa)
    probs = []

    for i in range(len(msa[0])):  # For each column index
        column = [seq[i] for seq in msa]
        counts = Counter(column)
        max_count = max(counts.values())
        probs.append(max_count / num_seqs)

    return probs

def conservation_category(msa, idx):
    """
    Categorizes the conservation level of a specific column in an MSA.

    This function determines a qualitative category for the conservation of a column
    based on the frequency of its most common character, excluding gaps.
    It defines categories such as "always conserved",
    "always_conserved_except_two" (always conserved minus up to 2 differences), etc.

    :param msa: A list of strings representing the Multiple Sequence Alignment (list of str).
    :param idx: The 0-indexed position of the column to categorize (int).

    :returns: A string indicating the conservation category:
              - "always_conserved": If all non-gap characters in the column are identical.
              - "always_conserved_except_two": If there are 1 or 2 differences from the most common character.
              - "always_conserved_except_4": If there are 3 or 4 differences from the most common character.
              - "others": For all other cases (more than 4 differences).
              (str)
    """
    # Filter out gaps ('-') from the column before calculating conservation
    column = [seq[idx] for seq in msa if seq[idx] != '-']

    if not column: # Handle case of column with only gaps
        return "others" # Or a specific category for all-gap columns if needed

    most_common_char, _ = Counter(column).most_common(1)[0] # Get the most common character and its count
    
    # Calculate the number of characters that are NOT the most common
    diff = len(column) - column.count(most_common_char)
    
    if diff == 0:
        return "always_conserved"
    elif diff <= 2:
        return "always_conserved_except_two"
    elif diff <= 4:
        return "always_conserved_except_4"
    else:
        return "others"


def clasify_conservation(msa):
    """
    Classifies all non-gap positions in a Multiple Sequence Alignment into conservation categories.

    This function iterates through each position of the MSA that is not a gap in the
    first sequence (assumed to be a representative for gap presence). For each such
    position, it determines its conservation category using `conservation_category`
    and groups the 1-based indices of these positions accordingly.

    :param msa: A list of strings representing the Multiple Sequence Alignment (list of str).

    :returns: A dictionary where keys are conservation categories ("always_conserved",
              "always_conserved_except_two", "always_conserved_except_4", "others") and
              values are lists of 1-based indices of the columns belonging to that category (dict).
    """
    # Identify positions that are not gaps in the first sequence (assuming it's representative)
    non_gap_indices = [j for j in range(len(msa[0])) if msa[0][j] != "-"]
    
    # Calculate the category for each non-gap position
    col_categories = {idx: conservation_category(msa, idx)
                      for idx in non_gap_indices}
    
    # Initialize dictionary to hold grouped positions
    grups = {
        "always_conserved": [],
        "always_conserved_except_two": [],
        "always_conserved_except_4": [],
        "others": []
    }

    # Populate the groups with 1-based indices
    # 'pos' is the 0-based index of the item in 'gap_free_positions', 'j' is the actual 0-based MSA index
    for j in non_gap_indices: # Iterating directly over non_gap_indices is cleaner here
        cat = col_categories[j]
        grups[cat].append(j + 1)  # VARNA uses 1-based indexing, so add 1

    return grups