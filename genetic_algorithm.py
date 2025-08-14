import random
import RNA
import copy
import os

# Import functions from your existing modules
from rna_structures import get_pairing_dict, constrained_mfe, get_base_pairs
from io_tools import save_and_plot_structures
from rna_cluster import matriz_rmsd, cluster_structures, organise_per_clusters, compute_metrics, visualise_metrics
from search import count_rna1_rna3_pairings, check_rna1_structure_preserved

import search # Import the whole module to access its functions

def get_pair_table(struct):
    """
    Generates a pairing table for a given RNA secondary structure.

    This function takes an RNA secondary structure represented in dot-bracket notation
    and returns a pairing table. The pairing table is an array where `ptable[i]`
    stores the 1-based index of the nucleotide that position `i` (also 1-based) is paired with.
    If position `i` is unpaired, `ptable[i]` will be 0. This table is useful for
    quickly querying base-pairing partners within a structure.

    :param struct: An RNA secondary structure string in dot-bracket notation (e.g., "((.))") (str).

    :returns: A list (or array) representing the pairing table. The list is 1-indexed,
              meaning `ptable[0]` is unused, and `ptable[i]` corresponds to the pairing
              partner of the i-th base (1-based index) in the input structure. (list of int)
    """
    return RNA.ptable(struct)



def initialize_population(rna1_orig, rna3_orig, linker_length, population_size, mutable_positions_rna1, struct1_orig):
    """
    Initializes a population of individuals for the Genetic Algorithm.

    Each individual in the population represents a candidate RNA switch design,
    structured as a dictionary containing its full sequence, mutated SRE segment,
    linker sequence, Aptamer segment, and a record of SRE mutations. This function
    introduces initial diversity by applying a small rate of random mutations
    to the SRE segment (only at allowed positions) and generating random linker sequences.
    If an SRE base is mutated and it's part of a pair in the original SRE structure,
    its complementary base is also mutated to attempt to preserve pairing potential.

    :param rna1_orig: The original, unmutated SRE sequence (str).
    :param rna3_orig: The original Aptamer sequence (str).
    :param linker_length: The desired fixed length for the linker sequences
                          generated for each individual (int).
    :param population_size: The total number of individuals to generate
                            for the initial population (int).
    :param mutable_positions_rna1: A list of 0-indexed positions in SRE where
                                   mutations are allowed during initialization and evolution (list of int).
    :param struct1_orig: The dot-bracket structure of the original SRE. Used to identify
                         complementary base pairs for coordinated mutations during initialization (str).

    :returns: A list of dictionaries, where each dictionary represents a newly
              created individual (candidate RNA switch) for the starting population (list of dict).
    """
    population = []
    bases = ['A', 'U', 'G', 'C']
    
    # Get pair table for original SRE structure to guide complementary mutations
    # Assuming get_pair_table is a helper function that calls RNA.ptable
    pt_rna1_orig = RNA.ptable(struct1_orig) # Using RNA.ptable directly as per previous examples

    # Define an initial mutation rate for population diversity
    initial_rna1_mutation_rate = 0.05 

    for _ in range(population_size):
        current_rna1_list = list(rna1_orig)
        mutations_info = [] # To store information about mutations applied in SRE
        
        # Apply initial random mutations to SRE only at allowed mutable positions
        for pos_in_rna1 in mutable_positions_rna1:
            if random.random() < initial_rna1_mutation_rate:
                old_base = current_rna1_list[pos_in_rna1]
                new_base = random.choice([b for b in bases if b != old_base])
                current_rna1_list[pos_in_rna1] = new_base
                mutations_info.append((pos_in_rna1, old_base, new_base)) # Store pos, old_base, new_base
                
                # --- !!!!! Apply complementary mutation if paired in original structure !!!!!---
                # ViennaRNA's ptable returns 1-based indexing for pairs, 0 for unpaired
                paired_pos_1_based = pt_rna1_orig[pos_in_rna1 + 1] 
                if paired_pos_1_based != 0: # If the base is paired in the original structure
                    paired_pos_0_based = paired_pos_1_based - 1 # Convert to 0-based for Python list

                    # Ensure the paired position is valid (within RNA1 bounds and not self-pairing)
                    if 0 <= paired_pos_0_based < len(rna1_orig) and paired_pos_0_based != pos_in_rna1:
                        complement_base_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
                        complement_base = complement_base_map.get(new_base)
                        
                        # Only mutate if the complement exists and the current base is not already the complement
                        if complement_base and current_rna1_list[paired_pos_0_based] != complement_base:
                            old_complement_base = current_rna1_list[paired_pos_0_based]
                            current_rna1_list[paired_pos_0_based] = complement_base
                            mutations_info.append((paired_pos_0_based, old_complement_base, complement_base))


        mutated_rna1_seq = ''.join(current_rna1_list)
        
        # Generate a random linker of the specified length
        linker_seq = ''.join(random.choice(bases) for _ in range(linker_length))
        
        full_seq = mutated_rna1_seq + linker_seq + rna3_orig

        individual = {
            'seq': full_seq,
            'rna1_mutated': mutated_rna1_seq, # The SRE sequence with applied mutations
            'linker': linker_seq,
            'rna3': rna3_orig, # rna3 original, not mutated
            'rna1_mutations_info': mutations_info # Information about the mutations in SRE
        }
        population.append(individual)
    return population

def fitness(individual, rna1_orig, rna3_orig, struct1_orig, constraint_orig, watched_positions_orig,
            mfe_delta, max_pairings, max_structure_changes):
    """
    Calculates the fitness score of an RNA switch individual within the genetic algorithm.

    This function evaluates an individual based on several key criteria to determine its
    suitability as a functional RNA switch. It applies an "all-or-nothing" approach for
    primary viability conditions, returning a very low (negative infinity) fitness if
    any critical criterion is not met. Otherwise, it calculates a score primarily based
    on the thermodynamic efficiency of the switch.

    The criteria include:
    1.  **Watched Positions Change:** Ensures that specified positions in SRE change their
        pairing status (paired vs. unpaired) OR their **pairing partner** when folded into the
        unconstrained (OFF) state, compared to the original SRE structure. This is crucial
        for the switch mechanism. All watched positions must demonstrate such a change.
    2.  **Minimizing SRE-Aptamer Pairings:** Limits undesired base pairings between the SRE
        and Aptamer segments in the unconstrained (OFF) state, which could interfere with the switch.
    3.  **SRE Structure Preservation (ON State):** Assesses how well the SRE segment maintains
        its desired structure when the full construct is folded under a specific constraint
        (representing the "ON" state). Significant deviations are penalized.
    4.  **MFE Delta Condition:** Requires a sufficient thermodynamic energy difference
        (MFE_ON - MFE_OFF) to ensure a clear energetic distinction between the ON and OFF states,
        which is vital for switch efficiency.

    :param individual: A dictionary representing the individual (candidate RNA switch) to be evaluated.
                       Must contain 'seq' (full RNA sequence), 'linker', and 'rna1_mutated' keys (dict).
    :param rna1_orig: The original, unmutated SRE sequence. Used for context in evaluations (str).
    :param rna3_orig: The original RNA3 sequence. Used for constructing the full sequence and
                      evaluating SRE-Aptamer pairings (str).
    :param struct1_orig: The dot-bracket structure of the original SRE sequence. Used as a
                         reference for evaluating changes in 'watched_positions_orig' and
                         'max_structure_changes' in the ON state (str).
    :param constraint_orig: The dot-bracket constraint string used to guide the folding of the
                            full RNA sequence towards the desired "ON" state (str).
    :param watched_positions_orig: A list of 0-indexed positions within SRE that are specifically
                                   monitored to ensure they undergo a structural change
                                   (paired to unpaired, or vice-versa, or change partner) in the OFF state (list of int).
    :param mfe_delta: The minimum acceptable energy difference (MFE_ON - MFE_OFF) in kcal/mol.
                      A higher value indicates a more thermodynamically stable switch (float).
    :param max_pairings: The maximum number of allowed base pairs between the SRE and Aptamer
                         segments in the unconstrained (OFF) state. Exceeding this leads to a
                         disqualified individual (int).
    :param max_structure_changes: The maximum allowed number of character differences between
                                  the SRE sub-structure within the constrained (ON) state
                                  and the original SRE structure. Exceeding this leads to a
                                  disqualified individual (int).

    :returns: The fitness value of the individual (float). Returns `-float('inf')` if the
              individual fails any of the critical "all-or-nothing" conditions. Otherwise,
              it returns the difference `(MFE_ON - MFE_OFF)`, where a higher value indicates
              better fitness.
    """
    seq_full = individual['seq']
    linker = individual['linker']
    rna1_mutated = individual['rna1_mutated']

    # Step 1: Unconstrained folding (OFF state)
    struct_full, mfe_1 = RNA.fold(seq_full) 

    # Check: change in watched_positions in mutated rna1
    pairing1_orig = get_pairing_dict(struct1_orig) 
    pairing_full_structure = get_pairing_dict(struct_full) 
    
    num_watched_pos_changed = 0
    for pos in watched_positions_orig:
        # Ensure position is within rna1 bounds
        if pos < len(rna1_orig): 
            # Check if pairing status (paired/unpaired) changed
            is_paired_orig = pos in pairing1_orig
            is_paired_current = pos in pairing_full_structure

            if is_paired_orig != is_paired_current:  # Changed paired/unpaired status
                num_watched_pos_changed += 1
            elif is_paired_orig and is_paired_current:  # Both paired, check if partner changed
                if pairing1_orig[pos] != pairing_full_structure[pos]:
                    num_watched_pos_changed += 1

    # "All-or-nothing" condition for watched_positions
    # All watched positions must have changed.
    if len(watched_positions_orig) > 0 and num_watched_pos_changed < len(watched_positions_orig):
        return -float('inf')  # Maximum penalty

    # Check: undesired pairings between SRE and Aptamer
    pairings_rna1_rna3_count = count_rna1_rna3_pairings(struct_full, rna1_mutated, rna3_orig, linker)
    if pairings_rna1_rna3_count > max_pairings:
        return -float('inf')  # Maximum penalty

    # Step 2: Constrained folding (ON state)
    struct_constr, mfe_2 = constrained_mfe(seq_full, constraint_orig)

    # Check: preservation of original SRE structure under constraint
    changes_rna1_in_constrained = search.check_rna1_structure_preserved(rna1_mutated, linker, rna3_orig, struct_constr, struct1_orig)
    if changes_rna1_in_constrained > max_structure_changes:
        return -float('inf')  # Maximum penalty

    # Check: MFE energy difference condition
    if mfe_1 >= mfe_2 - mfe_delta:  # Means mfe_2 - mfe_1 < mfe_delta
        return -float('inf')  # Maximum penalty

    # If all "all-or-nothing" conditions are met, calculate actual fitness
    # Higher fitness is better. We favor a larger energy difference (mfe_2 - mfe_1)
    score = (mfe_2 - mfe_1)  # Energy difference is main factor if conditions met

    return score


def select_tournament(population, fitness_scores, tournament_size=3):
    """
    Selects individuals from a population using tournament selection.

    In tournament selection, a small group of individuals (the tournament) is randomly
    chosen from the population. The individual with the highest fitness within this
    tournament is then selected to be a parent for the next generation. This process
    is repeated until enough individuals are selected to form the new breeding pool.

    :param population: A list of individuals (dictionaries) representing the current generation (list of dict).
    :param fitness_scores: A list of numerical fitness scores, where each score corresponds
                           to an individual at the same index in the `population` list (list of float).
    :param tournament_size: The number of individuals that compete in each tournament.
                            A larger size increases selection pressure. Defaults to 3 (int).

    :returns: A list of selected individuals (dictionaries) that will serve as parents
              for the next generation. The size of this list will typically match
              the original population size. (list of dict)
    """
    selected = []
    # Ensure population is large enough for tournament size
    if len(population) < tournament_size:
        tournament_size = len(population)
    if tournament_size == 0 and len(population) > 0:  # If tournament size is 0 but population exists
        return random.sample(population, len(population))  # Select all randomly

    for _ in range(len(population)):  # Select as many individuals as original population for next gen
        if len(population) < tournament_size:
            # Fallback if population too small for tournament size
            selected.extend(random.sample(population, min(len(population), 2)))
            continue

        # Randomly choose individuals for the tournament
        tournament_competitors = random.sample(list(zip(population, fitness_scores)), tournament_size)
        
        # Find the best in the tournament
        best_competitor = max(tournament_competitors, key=lambda x: x[1])
        selected.append(best_competitor[0])  # Add only the individual, not the fitness

    return selected


def mutate(individual, rna1_orig_seq, struct1_orig, mutable_positions_rna1, linker_length, mutation_rate_rna1=0.02, mutation_rate_linker=0.05):
    """
    Applies random mutations to an RNA switch individual's SRE and linker segments.

    Mutations in SRE are applied **only at specified mutable positions**. If a mutated base
    is part of a canonical pair in the original SRE structure, an attempt is made
    to mutate its complementary partner as well, aiming to **preserve pairing potential**.
    Linker mutations are applied randomly across its entire length.

    :param individual: A dictionary representing the individual to be mutated. It must contain
                       'rna1_mutated' (current SRE sequence), 'linker' (current linker sequence),
                       'rna3' (RNA3 sequence), and 'SRE_mutations_info' (list of past SRE mutations) keys (dict).
    :param rna1_orig_seq: The original, unmutated SRE sequence. Used for checking bounds
                          and ensuring consistent SRE length (str).
    :param struct1_orig: The dot-bracket structure of the original SRE. This is crucial for
                         identifying complementary base pairs to perform coordinated mutations (str).
    :param mutable_positions_rna1: A list of 0-indexed positions within SRE where mutations
                                   are permitted. Mutations will only occur at these specified indices (list of int).
    :param linker_length: The fixed length of the linker segment. Used to iterate through
                          linker positions for mutation (int).
    :param mutation_rate_rna1: The probability (per mutable base) that a mutation occurs
                               at a given position in SRE. Defaults to 0.02 (float).
    :param mutation_rate_linker: The probability (per base) that a mutation occurs at a
                                 given position in the linker. Defaults to 0.05 (float).

    :returns: A new dictionary representing the mutated individual. This includes the
              updated full sequence ('seq'), mutated SRE ('rna1_mutated'), mutated linker ('linker'),
              unchanged RNA3 ('rna3'), and an updated record of SRE mutations ('rna1_mutations_info') (dict).
    """
    mutated_rna1_list = list(individual['rna1_mutated'])
    current_linker_list = list(individual['linker'])
    bases = ['A', 'U', 'G', 'C']
    
    # Get pair table for the original SRE structure to guide complementary mutations
    # This uses ViennaRNA's 1-based indexing for pairs (0 if unpaired).
    pt_rna1_orig = get_pair_table(struct1_orig) 
    
    # Initialize mutation info with a copy of the individual's existing mutations.
    # This ensures that all mutations are tracked, not just those from the current step.
    new_rna1_mutations_info = list(individual['rna1_mutations_info']) 

    # --- Mutate SRE ---
    # Iterate only through the explicitly allowed mutable positions in SRE
    for pos_in_rna1 in mutable_positions_rna1:
        if random.random() < mutation_rate_rna1:
            old_base = mutated_rna1_list[pos_in_rna1]
            # Select a new base that is different from the old one
            new_base = random.choice([b for b in bases if b != old_base])
            mutated_rna1_list[pos_in_rna1] = new_base
            new_rna1_mutations_info.append((pos_in_rna1, old_base, new_base))
            
            # --- Mutate the complementary base if it exists and is part of a pair in the original SRE structure ---
            # Retrieve the 1-based index of the paired position from the original SRE structure's pair table
            paired_pos_1_based = pt_rna1_orig[pos_in_rna1 + 1] 
            
            # Check if 'pos_in_rna1' is indeed paired in the original structure (paired_pos_1_based != 0)
            if paired_pos_1_based != 0:
                # Convert the 1-based paired position to a 0-based index for Python list
                paired_pos_0_based = paired_pos_1_based - 1 
                
                # Ensure the paired position is valid (within RNA1 bounds and not self-pairing)
                if (paired_pos_0_based != pos_in_rna1 and 
                    paired_pos_0_based >= 0 and 
                    paired_pos_0_based < len(rna1_orig_seq)): 
                    
                    # Define a mapping for complementary bases
                    complement_base_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
                    complement_base = complement_base_map.get(new_base) # Get complement for the NEWLY mutated base
                    
                    # If a complement exists and the current base at the paired position is not already the complement
                    if complement_base and mutated_rna1_list[paired_pos_0_based] != complement_base:
                        old_complement_base = mutated_rna1_list[paired_pos_0_based]
                        mutated_rna1_list[paired_pos_0_based] = complement_base
                        new_rna1_mutations_info.append((paired_pos_0_based, old_complement_base, complement_base))

    # --- Mutate linker ---
    # Iterate through all positions in the linker
    for i in range(linker_length):
        if random.random() < mutation_rate_linker:
            old_base = current_linker_list[i]
            new_base = random.choice([b for b in bases if b != old_base])
            current_linker_list[i] = new_base

    # Join the lists back into strings
    mutated_rna1_seq = ''.join(mutated_rna1_list)
    mutated_linker_seq = ''.join(current_linker_list)
    
    # Create a new individual dictionary to avoid modifying the original in-place
    new_individual = individual.copy()
    new_individual['rna1_mutated'] = mutated_rna1_seq
    new_individual['linker'] = mutated_linker_seq
    # Reconstruct the full sequence with the newly mutated SRE and linker
    new_individual['seq'] = new_individual['rna1_mutated'] + new_individual['linker'] + new_individual['rna3']
    new_individual['rna1_mutations_info'] = new_rna1_mutations_info # Update the SRE mutations record

    return new_individual


def crossover(parent1, parent2, rna1_length, linker_length):
    """
    Performs a single-point crossover operation between two parent RNA switch individuals.

    This function generates two new child individuals by exchanging segments of their linker
    sequences. The SRE (mutated) and RNA3 sequences of each child are inherited directly
    from their respective parents without modification during this operation.

    :param parent1: A dictionary representing the first parent individual. It should contain at
                    least 'linker', 'rna1_mutated', 'rna3', and 'rna1_mutations_info' keys (dict).
    :param parent2: A dictionary representing the second parent individual, with the same structure as parent1 (dict).
    :param rna1_length: The length of the SRE segment. This parameter is used to correctly
                        handle sequence construction, though not directly for crossover logic itself (int).
    :param linker_length: The length of the linker segment. This is crucial for determining
                          valid crossover points (int).

    :returns: A tuple containing two dictionaries, each representing a newly generated child individual.
              Each child dictionary includes its full sequence ('seq'), inherited mutated SRE ('rna1_mutated'),
              the new linker ('linker'), inherited RNA3 ('rna3'), and SRE mutation information ('rna1_mutations_info'). (tuple of dict)
    """
    
    linker1 = list(parent1['linker'])
    linker2 = list(parent2['linker'])

    # Crossover point in the linker (ensuring at least one element in each segment)
    crossover_point = random.randint(1, linker_length - 1) if linker_length > 1 else 0

    # Create new linkers by swapping segments
    child1_linker = linker1[:crossover_point] + linker2[crossover_point:]
    child2_linker = linker2[:crossover_point] + linker1[crossover_point:]

    # Build new individuals by combining components
    child1_seq = parent1['rna1_mutated'] + ''.join(child1_linker) + parent1['rna3']
    child2_seq = parent2['rna1_mutated'] + ''.join(child2_linker) + parent2['rna3']

    child1 = {
        'seq': child1_seq,
        'rna1_mutated': parent1['rna1_mutated'],
        'linker': ''.join(child1_linker),
        'rna3': parent1['rna3'],
        'rna1_mutations_info': parent1['rna1_mutations_info']  # SRE mutation info comes from the parent
    }
    child2 = {
        'seq': child2_seq,
        'rna1_mutated': parent2['rna1_mutated'],
        'linker': ''.join(child2_linker),
        'rna3': parent2['rna3'],
        'rna1_mutations_info': parent2['rna1_mutations_info']
    }

    return child1, child2


def genetic_algorithm(rna1_orig, rna3_orig, struct1_orig, constraint_orig, mutable_positions_rna1, watched_positions_orig,
                      population_size=50, generations=100, linker_length=7, elitism_rate=0.1,
                      mfe_delta=0, max_pairings=5, max_structure_changes=6,
                      mutation_rate_rna1=0.02, mutation_rate_linker=0.05,
                      tournament_size=3, verbose=True, log_func=None, check_stop_func=None):
    """
    Implements the core Genetic Algorithm (GA) for designing optimal RNA switches.

    This function evolves a population of RNA switch candidates over multiple generations
    by applying genetic operators (selection, crossover, mutation) based on a
    defined fitness function. It aims to find RNA linker sequences and RNA1 mutations
    that satisfy specific structural and thermodynamic criteria for a functional switch.

    :param rna1_orig: The original SRE sequence (str).
    :param rna3_orig: The RNA3 sequence (str).
    :param struct1_orig: The original dot-bracket structure of SRE (representing the OFF state) (str).
    :param constraint_orig: The dot-bracket constraint for the SRE + linker + RNA3 sequence
                            to guide the folding towards the desired ON state (str).
    :param mutable_positions_rna1: A list of 0-indexed positions in SRE where mutations are allowed during GA evolution (list of int).
    :param watched_positions_orig: A list of 0-indexed positions in SRE that must undergo
                                   a change in their pairing status between the OFF and ON states (list of int).
    :param population_size: The number of individuals (candidate designs) in each generation of the GA. Defaults to 50 (int).
    :param generations: The total number of generations the genetic algorithm will run. Defaults to 100 (int).
    :param linker_length: The fixed length of the linker sequences that the genetic algorithm will optimize. Defaults to 7 (int).
    :param elitism_rate: The proportion of the best-performing individuals (elite) from the current generation
                          that are directly carried over to the next generation without modification. Defaults to 0.1 (float).
    :param mfe_delta: The minimum required difference in MFE (MFE_ON - MFE_OFF) for a
                      valid switch, used within the fitness function. Defaults to 0 (float).
    :param max_pairings: The maximum allowed number of base pairings between SRE and RNA3 in the
                          unconstrained (OFF) state, used within the fitness function. Defaults to 5 (int).
    :param max_structure_changes: The maximum allowed number of changes in the SRE sub-structure within the
                                   constrained (ON) state compared to its original structure, used within the fitness function. Defaults to 6 (int).
    :param mutation_rate_rna1: The probability of a single nucleotide mutation occurring in the SRE segment
                               of an individual during the mutation phase. Defaults to 0.02 (float).
    :param mutation_rate_linker: The probability of a single nucleotide mutation occurring in the linker segment
                                  of an individual during the mutation phase. Defaults to 0.05 (float).
    :param tournament_size: The number of individuals randomly selected from the population for
                            tournament selection. The individual with the highest fitness among them is chosen as a parent. Defaults to 3 (int).
    :param verbose: If True, prints detailed progress, diagnostics for the best individual found, and generation summaries. Defaults to True (bool).
    :param log_func: An optional callable function to use for logging messages throughout the GA process.
                      This allows redirecting logs to a Streamlit app or a file. Defaults to print if None (callable, optional).
    :param check_stop_func: An optional callable function that, when invoked, should return True
                             if the genetic algorithm process needs to be halted prematurely (e.g., user interruption). Defaults to None (callable, optional).

    :returns: A list of dictionaries, where each dictionary represents a valid RNA switch solution
              found within the final population. Each dictionary includes the 'individual' data
              (linker, SRE mutations), its 'fitness', the 'seq_full', 'struct_unconstr', 'mfe_1',
              'struct_constr', and 'mfe_2'.
    """
    if log_func is None:
        log_func = print

    population = initialize_population(rna1_orig, rna3_orig, linker_length, population_size, mutable_positions_rna1, struct1_orig)
    rna1_len = len(rna1_orig)

    best_overall_individual = None
    best_overall_fitness = -float('inf')

    # --- NEW LIST TO STORE ALL UNIQUE VALID SOLUTIONS FOUND ---
    # We use a dictionary to store the best version (by fitness) of each unique complete sequence
    all_unique_valid_solutions = {}

    for generation in range(generations):
        if check_stop_func and check_stop_func():
            log_func("Genetic Algorithm stopped by user request.")
            break

        fitness_scores = []
        # Evaluate the current population and also record valid solutions
        for individual_index, individual in enumerate(population):
            fit = fitness(individual, rna1_orig, rna3_orig, struct1_orig, constraint_orig, watched_positions_orig,
                          mfe_delta, max_pairings, max_structure_changes)
            fitness_scores.append(fit)

            # --- RECORDING VALID SOLUTIONS DURING GENERATIONS ---
            if fit > -float('inf'):  # If it is a valid solution
                current_full_sequence = individual['seq']

                # Calculate MFEs and structures for this individual (needed for storage)
                unconstrained_struct, mfe_1 = RNA.fold(current_full_sequence)
                constrained_struct, mfe_2 = constrained_mfe(current_full_sequence, constraint_orig)

                # Create a dict with the data of this valid solution
                current_solution_data = {
                    'individual': copy.deepcopy(individual),  # IMPORTANT: deep copy
                    'fitness': fit,
                    'full_sequence': current_full_sequence,
                    'unconstrained_structure': unconstrained_struct,
                    'mfe_1': mfe_1,
                    'constrained_structure': constrained_struct,
                    'mfe_2': mfe_2
                }

                # Logic to save unique solutions or the best of duplicates (by full_sequence):
                if current_full_sequence not in all_unique_valid_solutions or \
                   fit > all_unique_valid_solutions[current_full_sequence]['fitness']:
                    all_unique_valid_solutions[current_full_sequence] = current_solution_data
                    if verbose:
                        log_func(f"New valid solution found in Gen {generation + 1}: Linker={individual['linker']}, RNA1_mut={individual['rna1_mutated']}, Fitness={fit:.2f}")

        scored_population = list(zip(population, fitness_scores))
        scored_population.sort(key=lambda x: x[1], reverse=True)  # Sort by fitness descending

        current_best_individual, current_best_fitness = scored_population[0]
        if current_best_fitness > best_overall_fitness:
            best_overall_fitness = current_best_fitness
            best_overall_individual = copy.deepcopy(current_best_individual) 

            if verbose:
                log_func(f"\n--- New best solution found in Gen {generation+1} ---")
                log_func(f"Best fitness: {best_overall_fitness:.2f}")
                log_func(f"Linker: {best_overall_individual['linker']}")
                log_func(f"SRE mutations: {best_overall_individual['rna1_mutations_info']}")

                # --- DIAGNOSTIC OF BEST SOLUTION BEFORE PLOTTING ---
                log_func("\n--- Diagnostic of the BEST SOLUTION before plotting ---")
                seq_diag = best_overall_individual['seq']
                struct_full_diag, mfe_1_diag = RNA.fold(seq_diag)
                struct_constr_diag, mfe_2_diag = constrained_mfe(seq_diag, constraint_orig)

                pairing1_orig_for_diag = get_pairing_dict(struct1_orig)
                pairing_full_structure_for_diag = get_pairing_dict(struct_full_diag)

                num_watched_pos_changed_for_diag = 0
                for pos in watched_positions_orig:
                    if pos < len(rna1_orig):
                        is_paired_orig = pos in pairing1_orig_for_diag
                        is_paired_current = pos in pairing_full_structure_for_diag
                        if is_paired_orig != is_paired_current:
                            num_watched_pos_changed_for_diag += 1
                        elif is_paired_orig and is_paired_current:
                            if pairing1_orig_for_diag[pos] != pairing_full_structure_for_diag[pos]:
                                num_watched_pos_changed_for_diag += 1

                log_func(f"✅ Watched positions changed: {num_watched_pos_changed_for_diag}/{len(watched_positions_orig)}")

                pairings_rna1_rna3_count_for_diag = count_rna1_rna3_pairings(
                    struct_full_diag,
                    best_overall_individual['rna1_mutated'],
                    best_overall_individual['linker'], 
                    best_overall_individual['rna3'] 
                )
                log_func(f"✅ SRE-Aptamer pairings: {pairings_rna1_rna3_count_for_diag} (Max allowed: {max_pairings})")

                changes_rna1_in_constrained_for_diag = check_rna1_structure_preserved(
                    best_overall_individual['rna1_mutated'],
                    best_overall_individual['linker'],
                    best_overall_individual['rna3'],
                    struct_constr_diag,
                    struct1_orig
                )
                log_func(f"✅ SRE changes (ON state): {changes_rna1_in_constrained_for_diag} (Max allowed: {max_structure_changes})")

                mfe_diff_for_diag = mfe_2_diag - mfe_1_diag
                log_func(f"✅ MFE difference (MFE_ON - MFE_OFF): {mfe_diff_for_diag:.2f} (Required: > {mfe_delta})")
                log_func("--------------------------------------------------\n")

                
                log_func("--------------------------------------------------")

        # Elitism selection
        num_elite = int(population_size * elitism_rate)
        elite = [ind for ind, _ in scored_population[:num_elite]]

        # Tournament selection for breeding
        selected_for_breeding = select_tournament([ind for ind, _ in scored_population], fitness_scores, tournament_size)

        # Create new population
        new_population = elite  # Elites are passed directly

        while len(new_population) < population_size:
            if len(selected_for_breeding) < 2:  # Ensure enough parents for crossover
                if len(elite) > 0:  # If not enough for tournament, duplicate elite randomly
                    new_population.append(copy.deepcopy(random.choice(elite))) # Copia profunda aquí también
                else:  # If no elite, generate random individuals
                    new_population.append(initialize_population(rna1_orig, rna3_orig, linker_length, 1, mutable_positions_rna1, struct1_orig)[0])
                continue

            parent1 = random.choice(selected_for_breeding)
            parent2 = random.choice(selected_for_breeding)

            # Crossover
            child1, child2 = crossover(parent1, parent2, rna1_len, linker_length)

            # Mutation
            child1 = mutate(child1, rna1_orig, struct1_orig, mutable_positions_rna1, linker_length, mutation_rate_rna1, mutation_rate_linker)
            child2 = mutate(child2, rna1_orig, struct1_orig, mutable_positions_rna1, linker_length, mutation_rate_rna1, mutation_rate_linker)

            new_population.extend([child1, child2])

        # Ensure population does not exceed max size
        population = new_population[:population_size]

        if verbose and generation % 5 == 0:  # Print progress every 5 generations
            avg_fitness = sum(fitness_scores) / len(fitness_scores) if fitness_scores else 0
            log_func(f"Generation {generation+1}/{generations}. Best fitness this generation: {current_best_fitness:.2f}. "
                     f"Average fitness: {avg_fitness:.2f}")

    # --- Collect all unique valid individuals found throughout the entire GA run ---
    # Convertir el diccionario de soluciones únicas a una lista de diccionarios
    final_solutions_to_return = list(all_unique_valid_solutions.values())

    # Sort final solutions by fitness
    final_solutions_to_return.sort(key=lambda x: x['fitness'], reverse=True)


    

    if verbose:
        log_func("\nGenetic Algorithm Finished.")
        log_func(f"Total unique valid solutions found across all generations: {len(final_solutions_to_return)}")

    return final_solutions_to_return



def run_genetic_algorithm_search(rna1, rna3, struct1, constraint,
                                 mutable_rna1, watched_positions, output_dir,
                                 use_mutations=True, mfe_delta=0,
                                 max_pairings_rna1_rna3=5,
                                 max_structure_changes=6, num_mut=0,
                                 linker_lengths=None, 
                                 verbose=True, check_stop_func=None, log_func=print,
                                 population_size=50, generations=100, linker_length_ga=7, 
                                 elitism_rate=0.1, mutation_rate_rna1=0.02,
                                 mutation_rate_linker=0.05, tournament_size=3):
    """
    Orchestrates the search for RNA switch designs using a Genetic Algorithm (GA).

    This function acts as a wrapper around the core genetic algorithm, setting up
    parameters, initiating the GA run, processing its results, saving all found
    valid designs, and generating a comprehensive report with optional structural clustering.

    :param rna1: The original SRE sequence (str).
    :param rna3: The RNA3 sequence (str).
    :param struct1: The original dot-bracket structure of SRE (representing the OFF state) (str).
    :param constraint: The dot-bracket constraint for the SRE + linker + RNA3 sequence
                       to guide the folding towards the desired ON state (str).
    :param mutable_rna1: A list of 0-indexed positions in SRE that are allowed to mutate (list of int).
    :param watched_positions: A list of 0-indexed positions in SRE to check for structural changes.
                              These positions should ideally change their pairing status between the OFF and ON states (list of int).
    :param use_mutations: If True, the GA will explore mutations in SRE. Defaults to True (bool).
    :param mfe_delta: The minimum required difference in MFE (MFE_ON - MFE_OFF) for a valid switch.
                      Higher values indicate a stronger switch. Defaults to 0 (float).
    :param max_pairings_rna1_rna3: Maximum allowed base pairings between SRE and RNA3 in the
                                   unconstrained (OFF) state. Defaults to 5 (int).
    :param max_structure_changes: Maximum allowed changes in the SRE sub-structure within the
                                  constrained (ON) state compared to its original structure. Defaults to 6 (int).
    :param num_mut: (Currently unused by GA's internal mutation rates, kept for compatibility if needed).
                    Defaults to 0 (int).
    :param linker_lengths: This parameter is primarily for compatibility with other search methods.
                           For the GA, it typically uses the first value or assumes a single linker length.
                           Defaults to None, which will be interpreted from `linker_length_ga`. (range or list of int, optional).
    :param verbose: If True, prints detailed progress and log messages during the search. Defaults to True (bool).
    :param check_stop_func: An optional callable function that, when invoked, returns True if the
                            GA process needs to be halted prematurely. Useful for UI integration (callable, optional).
    :param log_func: An optional callable function to use for logging messages. This allows
                     redirecting logs to a Streamlit app (e.g., st.write) or a file.
                     Defaults to print (callable, optional).
    :param population_size: The number of individuals in each generation of the genetic algorithm. Defaults to 50 (int).
    :param generations: The total number of generations the genetic algorithm will run. Defaults to 100 (int).
    :param linker_length_ga: The specific length of the linker sequences that the genetic algorithm will target. Defaults to 7 (int).
    :param elitism_rate: The proportion of the best individuals from the current generation that
                         are directly passed to the next generation without mutation or crossover. Defaults to 0.1 (float).
    :param mutation_rate_rna1: The probability of a mutation occurring in the SRE segment of an individual. Defaults to 0.02 (float).
    :param mutation_rate_linker: The probability of a mutation occurring in the linker segment of an individual. Defaults to 0.05 (float).
    :param tournament_size: The number of individuals chosen randomly from the population for
                            tournament selection. The individual with the highest fitness among them is selected. Defaults to 3 (int).

    :returns: A tuple containing:
              - all_final_results (list): A list of dictionaries, where each dictionary represents a
                                          valid linker solution found by the GA, including its sequence,
                                          structures, MFEs, and mutation information.
              - final_report_lines (str): A comprehensive textual report summarizing the GA search,
                                          including parameters and a summary of found linkers.
              - cluster_labels (list): A list of cluster labels assigned to each valid linker based
                                       on structural similarity, derived from post-search clustering.
    """
    
    all_final_results = []
    final_report_lines = []
    cluster_labels = [] 

    final_report_lines.append("### Genetic Algorithm Linker Search Report\n")
    final_report_lines.append(f"**Original SRE:** {rna1}\n")
    final_report_lines.append(f"**Aptamer:** {rna3}\n")
    final_report_lines.append(f"**Directed SRE Structure (Dot-Bracket):** {struct1}\n")
    final_report_lines.append(f"**Constraint for ON state:** {constraint}\n")
    final_report_lines.append(f"**Mutable positions in SRE:** {mutable_rna1}\n")
    final_report_lines.append(f"**'Watched' positions in SRE:** {watched_positions}\n")
    final_report_lines.append(f"**Allow mutations in SRE:** {use_mutations}\n")
    final_report_lines.append(f"**Minimum MFE Delta:** {mfe_delta:.2f} kcal/mol\n")
    final_report_lines.append(f"**Max SRE-Aptamer pairings:** {max_pairings_rna1_rna3}\n")
    final_report_lines.append(f"**Max SRE structure changes (ON state):** {max_structure_changes}\n")
    # Add GA specific parameters to the report
    final_report_lines.append(f"**GA Population Size:** {population_size}\n")
    final_report_lines.append(f"**GA Generations:** {generations}\n")
    final_report_lines.append(f"**GA Elitism Rate:** {elitism_rate}\n")
    final_report_lines.append(f"**GA SRE Mutation Rate:** {mutation_rate_rna1}\n")
    final_report_lines.append(f"**GA Linker Mutation Rate:** {mutation_rate_linker}\n")
    final_report_lines.append(f"**GA Tournament Size:** {tournament_size}\n")
    
    final_report_lines.append("\n---\n")

    # For now only one linker_length is supported (it gets the first value of the input string)
    actual_linker_length_for_ga = linker_lengths[0] if isinstance(linker_lengths, range) else linker_lengths
    
    log_func(f"\n🚀 Starting Genetic Algorithm for linker_length = {actual_linker_length_for_ga}...")
    
    # Receive the list of valid solutions
    valid_solutions_from_ga = genetic_algorithm(
        rna1_orig=rna1,
        rna3_orig=rna3,
        struct1_orig=struct1,
        constraint_orig=constraint,
        mutable_positions_rna1=mutable_rna1,
        watched_positions_orig=watched_positions,
        population_size=population_size, 
        generations=generations,         
        linker_length=linker_length_ga,  
        elitism_rate=elitism_rate,       
        mfe_delta=mfe_delta,
        max_pairings=max_pairings_rna1_rna3,
        max_structure_changes=max_structure_changes,
        mutation_rate_rna1=mutation_rate_rna1, 
        mutation_rate_linker=mutation_rate_linker, 
        tournament_size=tournament_size, 
        verbose=verbose,
        log_func=log_func,
        check_stop_func=check_stop_func
    )

    if valid_solutions_from_ga:
        log_func(f"✅ Found {len(valid_solutions_from_ga)} valid linker(s) by GA.")
        
        # Transform the list of valid solutions into the expected results structure
        for sol in valid_solutions_from_ga:
            ind = sol['individual'] # The individual object
            rna1_part_len = len(rna1) # Length of original rna1
            all_final_results.append({
                'linker': ind['linker'],
                'mut1_info': ind['rna1_mutations_info'],
                'sequence': sol['full_sequence'],
                'structure_unconstrained': sol['unconstrained_structure'],
                'structure_constrained': sol['constrained_structure'],
                'mfe_1': float(sol['mfe_1']),
                'mfe_2': float(sol['mfe_2']),
                'fitness': float(sol['fitness']),
                # Add rna1_mutated_seq for plotting/reporting if needed by outer functions
                'rna1_mutated_seq': sol['full_sequence'][0:rna1_part_len],
                'rna3': ind['rna3'] # Ensure rna3 is in the result for plotting functions
                #'pairings_count' : count_rna1_rna3_pairings(sol['unconstrained_structure'], sol['full_sequence'][0:rna1_part_len], ind['rna3'], ind['linker'])
            })
        # Sort results by fitness (higher is better) if desired
        all_final_results.sort(key=lambda x: x['fitness'], reverse=True)

        # ---: Save images for EVERY valid linker found ---
        log_func("\nSaving images for all valid linkers found by GA...")
        for i, res in enumerate(all_final_results):

            rna1_mutated_segment = res['sequence'][0:len(rna1)] 

            save_and_plot_structures(
                seq=res['sequence'],
                structure_unconstr=res['structure_unconstrained'],
                structure_constr=res['structure_constrained'],
                rna1=rna1_mutated_segment, 
                linker=res['linker'],
                rna3=res['rna3'], 
                mut1_info=res['mut1_info'],
                mfe_1=res['mfe_1'],
                mfe_2=res['mfe_2'],
                folder_prefix=output_dir
            )
            log_func(f"Saved images for linker: {res['linker']}")

        final_report_lines.append(f"\n### Genetic Algorithm Results ({len(all_final_results)} Valid Linker(s) Found)\n")
        
        # Add summary of all results to the report
        for i, res in enumerate(all_final_results):
            final_report_lines.append(f"**Result {i+1}:**\n")
            final_report_lines.append(f"   Linker: {res['linker']}\n")
            final_report_lines.append(f"   Fitness (MFE_ON - MFE_OFF): {res['fitness']:.2f}\n")
            final_report_lines.append(f"   SRE Mutations: {res['mut1_info']}\n")
            final_report_lines.append(f"   MFE OFF: {res['mfe_1']:.2f} kcal/mol, MFE ON: {res['mfe_2']:.2f} kcal/mol\n")
            final_report_lines.append("---\n")

        # Call clustering and plotting functions from rna_cluster and io_tools
        if len(all_final_results) > 1:
            msa_for_cluster = [res['sequence'] for res in all_final_results]
            estructuras_for_cluster = [res['structure_unconstrained'] for res in all_final_results]

            rmsd_matrix = matriz_rmsd(estructuras_for_cluster)
            labels = cluster_structures(rmsd_matrix, eps=3.0, min_samples=1)
            unique_labels = set(labels)
            num_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0) if labels.size > 0 else 0
            
            log_func(f"Clustering complete. Found {num_clusters} clusters.")
            
            cluster_labels = labels
            
            resumen_cluster = organise_per_clusters(msa_for_cluster, estructuras_for_cluster, labels, rmsd_matrix, output_base=f"{output_dir}/clusters")

            metricas = compute_metrics(estructuras_for_cluster, times=[])
            visualise_metrics(metricas, output_dir)
        else:
            log_func("Only one valid solution found, clustering skipped.")
            cluster_labels = [0] if all_final_results else []
            
    else:
        log_func(f"❌ No valid linkers found by Genetic Algorithm for linker_length = {actual_linker_length_for_ga}.")
        final_report_lines.append("\nNo valid linkers were found during the Genetic Algorithm search.\n")

    return all_final_results, "\n".join(final_report_lines), cluster_labels
