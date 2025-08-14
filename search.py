import itertools
import RNA
from rna_structures import get_pairing_dict, constrained_mfe, get_base_pairs
from rna_mutation import mutate_sequence, mutate_with_complement
from io_tools import save_and_plot_structures
from rna_cluster import matriz_rmsd, cluster_structures, organise_per_clusters, compute_metrics, visualise_metrics

def count_rna1_rna3_pairings(struct, rna1, rna3, linker):
    """
    Counts the number of base pairings that occur between the RNA1 and RNA3 segments
    within a given full RNA secondary structure.

    This function helps assess unwanted interactions between RNA1 and RNA3, which
    could interfere with the desired RNA switch mechanism.

    :param struct: The dot-bracket string representing the full secondary structure
                   of the RNA1-linker-RNA3 construct (str).
    :param rna1: The sequence of the RNA1 segment (str).
    :param rna3: The sequence of the RNA3 segment (str).
    :param linker: The sequence of the linker connecting RNA1 and RNA3 (str).
    :returns: The total number of base pairs formed between RNA1 and RNA3 (int).
    """
    pairs = get_base_pairs(struct)
    pairings_count = 0
    len1 = len(rna1)
    len_linker = len(linker)
    for pos1, pos2 in pairs:
        if (pos1 < len1 and pos2 >= len1 + len_linker) or \
           (pos2 < len1 and pos1 >= len1 + len_linker):
            pairings_count += 1
    return pairings_count

def check_rna1_structure_preserved(rna1, linker, rna3, structure_constrained, structure_rna1):
    """
    Calculates the number of differences between the original RNA1 structure
    and its predicted substructure when part of the full constrained construct.

    This helps verify if the RNA1 segment maintains its intended structure
    (or deviates minimally) when the switch is in its "ON" state under
    the given structural constraints.

    :param rna1: The sequence of the RNA1 segment (str).
    :param linker: The sequence of the linker connecting RNA1 and RNA3 (str).
    :param rna3: The sequence of the RNA3 segment (str).
    :param structure_constrained: The dot-bracket string of the full
                                  RNA1-linker-RNA3 construct when folded
                                  with a structural constraint (representing the ON state) (str).
    :param structure_rna1: The original dot-bracket structure of just the
                           RNA1 segment (representing its ideal OFF state structure or target ON structure) (str).
    :returns: The number of differing characters between the original RNA1 structure
              and its substructure within the constrained full construct (int).
    """
    len1 = len(rna1)
    substructure = structure_constrained[:len1]
    num_changes = sum(1 for a, b in zip(structure_rna1, substructure) if a != b)
    return num_changes

def run_linker_search(rna1, rna3, struct1, constraint,
                      mutable_rna1, watched_positions, output_dir,
                      use_mutations=True, mfe_delta=0,
                      max_pairings_rna1_rna3=5,
                      max_structure_changes=6, num_mut=0,
                      linker_lengths=range(7, 10),
                      verbose=True, check_stop_func=None, log_func=print):
    """
    Performs a comprehensive search for optimal RNA linkers and, optionally, RNA1
    mutations to design functional RNA switches.

    This function iterates through a defined range of linker lengths and all possible
    sequences for each length (brute-force approach). For each combination, it
    evaluates critical thermodynamic and structural criteria to identify valid
    candidates for RNA switches.

    :param rna1: The original RNA1 sequence (str).
    :param rna3: The RNA3 sequence (str).
    :param struct1: The original dot-bracket structure of RNA1 (representing the OFF state) (str).
    :param constraint: The dot-bracket constraint for the RNA1 + linker + RNA3 sequence
                       to guide the folding towards the desired ON state (str).
    :param mutable_rna1: A list of 0-indexed positions in RNA1 where mutations are allowed (list of int).
    :param watched_positions: A list of 0-indexed positions in RNA1 that must undergo
                              a change in their pairing status (paired/unpaired) from
                              the OFF state (struct1) to the OFF state of the full construct
                              (RNA1-linker-RNA3) (list of int).
    :param use_mutations: If True, the function will also explore mutations in RNA1
                           up to 'num_mut' for each linker. Defaults to True (bool).
    :param mfe_delta: The minimum required difference in MFE (MFE_ON - MFE_OFF) for a
                      linker to be considered valid. Higher values indicate a stronger switch (float).
                      Solutions where MFE_ON - MFE_OFF < mfe_delta are discarded. Defaults to 0.
    :param max_pairings_rna1_rna3: The maximum allowed number of base pairings between
                                   RNA1 and RNA3 in the unconstrained (OFF) state of the
                                   full RNA1-linker-RNA3 construct. Solutions exceeding
                                   this count are discarded. Defaults to 5 (int).
    :param max_structure_changes: The maximum allowed number of changes in the RNA1 sub-structure
                                  within the constrained (ON) state compared to its original
                                  structure (struct1). Solutions exceeding this are discarded.
                                  Defaults to 6 (int).
    :param num_mut: The maximum number of mutations to try for RNA1 if 'use_mutations' is True.
                    Defaults to 0 (int).
    :param linker_lengths: A range or list of integers specifying the lengths of linkers
                           to generate and test. Defaults to range(7, 10).
    :param verbose: If True, prints detailed progress and log messages during the search (bool).
                    Defaults to True.
    :param check_stop_func: An optional callable function that, when invoked, should return
                            True if the search process needs to be halted prematurely. Useful
                            for integrating with user interfaces (callable, optional).
    :param log_func: An optional callable function to use for logging messages. This allows
                     redirecting logs to a Streamlit app (e.g., st.write) or a file.
                     Defaults to print (callable, optional).

    :returns: A tuple containing:
              - results (list): A list of dictionaries, where each dictionary represents a
                                valid linker solution with its sequence, structures, MFEs,
                                and mutation information.
              - report (str): A comprehensive textual report summarizing the search process,
                              including statistics on discarded linkers and search parameters.
              - labels (list): A list of cluster labels assigned to each valid linker based
                               on structural similarity, derived from post-search clustering.
    """
    discarded_counts = {
        "pairings_exceed": 0,
        "watched_positions_no_change": 0,
        "structure_changes_exceed": 0,
        "mfe_delta_not_met": 0,
        "mutations_filtered": 0,
    }

    # Store discarted linkers for the report
    discarded_examples = {
        "pairings_exceed": [],
        "watched_positions_no_change": [],
        "structure_changes_exceed": [],
        "mfe_delta_not_met": [],
        "mutations_filtered": [],
    }


    results = []
    linkers_tried = set()
    count_found = 0
    count = 0
    stop_requested = False
    
    for linker_length in linker_lengths:
        all_linkers = [''.join(p) for p in itertools.product('AUGC', repeat=linker_length)]
        if verbose:
            log_func(f"🧪 Testing {len(all_linkers)} linkers of length {linker_length}")
        count=0
        for linker in all_linkers:
            
            if check_stop_func:
                val = check_stop_func()
                if val:
                    log_func("🛑 Search stoped by user.")
                    stop_requested = True
                    break
            if (count%1000==0):
                log_func(f'Tryed {count} linkers')
            count = count+1
            if linker in linkers_tried:
                continue
            linkers_tried.add(linker)

            seq_full = rna1 + linker + rna3
            struct_full, mfe_1 = RNA.fold(seq_full)

            pairing1 = get_pairing_dict(struct1)
            pairing_test = get_pairing_dict(struct_full)

            changed = all(pairing1.get(pos) != pairing_test.get(pos) for pos in watched_positions)
            pairings_count = count_rna1_rna3_pairings(struct_full, rna1, rna3, linker)
            if pairings_count > max_pairings_rna1_rna3:
                discarded_counts["pairings_exceed"] += 1
                if len(discarded_examples["pairings_exceed"]) < 3:
                    discarded_examples["pairings_exceed"].append(linker)
                continue
            if not changed:
                discarded_counts["watched_positions_no_change"] += 1
                if len(discarded_examples["watched_positions_no_change"]) < 3:
                    discarded_examples["watched_positions_no_change"].append(linker)
                continue
           

            struct_constr, mfe_2 = constrained_mfe(seq_full, constraint)
            changes_rna1 = check_rna1_structure_preserved(rna1, linker, rna3, struct_constr, struct1)
            
            if changes_rna1 > max_structure_changes:
                discarded_counts["structure_changes_exceed"] += 1
                if len(discarded_examples["structure_changes_exceed"]) < 3:
                    discarded_examples["structure_changes_exceed"].append(linker)
                continue
            
            if not use_mutations:
                if mfe_1 >= mfe_2 - mfe_delta:
                    discarded_counts["mfe_delta_not_met"] += 1
                    if len(discarded_examples["mfe_delta_not_met"]) < 3:
                        discarded_examples["mfe_delta_not_met"].append(linker)
                    continue
                elif mfe_1 < mfe_2 - mfe_delta:
                    log_func(f"✅ Valid linker: {linker}, ΔMFE={mfe_2 - mfe_1:.2f}")
                    results.append({
                        'linker': linker,
                        'mut1_info': [],
                        'sequence': seq_full,
                        'rna1' : rna1,
                        'rna3' : rna3,
                        'structure_unconstrained': struct_full,
                        'structure_constrained': struct_constr,
                        'mfe_1': mfe_1,
                        'mfe_2': mfe_2
                    })
                    save_and_plot_structures(seq_full, struct_full, struct_constr, rna1, linker, rna3, [], mfe_1, mfe_2, output_dir)
                    count_found += 1
                continue

            for n_mut in range(0, num_mut):
                for rna1_mut, mut1_info in mutate_sequence(rna1, mutable_rna1, n_mut):
                    mutated_rna1 = mutate_with_complement(rna1_mut, struct1, mut1_info)
                    seq_full = mutated_rna1 + linker + rna3
                    struct_full, mfe_1 = RNA.fold(seq_full)
                    pairing_full = get_pairing_dict(struct_full)

                    changed = all(pairing1.get(pos) == pairing_full.get(pos) for pos in watched_positions)
                    pairings_count = count_rna1_rna3_pairings(struct_full, rna1, rna3, linker)

                    if not changed or pairings_count > max_pairings_rna1_rna3:
                        continue

                    struct_constr, mfe_2 = constrained_mfe(seq_full, constraint)
                    changes_rna1 = check_rna1_structure_preserved(rna1, linker, rna3, struct_constr, struct1)
                    if changes_rna1 > max_structure_changes:
                        continue

                    if mfe_1 < mfe_2 - mfe_delta:
                        log_func(f"🧬 Mutated valid linker: {linker} with mutations {mut1_info}, ΔMFE={mfe_2 - mfe_1:.2f}")
                        results.append({
                            'linker': linker,
                            'mut1_info': mut1_info,
                            'sequence': seq_full,
                            'rna1_mutated_seq' : mutated_rna1,
                            'rna3' : rna3,
                            'structure_unconstrained': struct_full,
                            'structure_constrained': struct_constr,
                            'mfe_1': mfe_1,
                            'mfe_2': mfe_2
                        })
                        save_and_plot_structures(seq_full, struct_full, struct_constr, mutated_rna1, linker, rna3, mut1_info, mfe_1, mfe_2)
                        count_found += 1
                        break
        if stop_requested:
            break 
        print('strart clustering')

            # Call clustering and plotting functions from rna_cluster and io_tools
        if len(results) > 1:
            print('entered clustering')
            msa_for_cluster = [res['sequence'] for res in results]
            estructuras_for_cluster = [res['structure_unconstrained'] for res in results]
            print('entered clustering 1')
            rmsd_matrix = matriz_rmsd(estructuras_for_cluster)
            labels = cluster_structures(rmsd_matrix, eps=3.0, min_samples=1)
            print('entered clustering 2')
            unique_labels = set(labels)
            num_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0) if labels.size > 0 else 0
            print('entered clustering 3')
            log_func(f"Clustering complete. Found {num_clusters} clusters.")
            print('entered clustering 4')
            cluster_labels = labels
            
            resumen_cluster = organise_per_clusters(msa_for_cluster, estructuras_for_cluster, labels, rmsd_matrix, output_base=f"{output_dir}/clusters")

            metricas = compute_metrics(estructuras_for_cluster, times=[])
            visualise_metrics(metricas, output_dir)
        else:
            log_func("Only one valid solution found, clustering skipped.")
            cluster_labels = [0] if results else []
            

    report = f"--- Linker Finder Report ---\n"
    report += f"Total evaluated linkers: {len(linkers_tried)}\n"
    report += f"Valid found linkers: {len(results)}\n\n"

    for key in discarded_counts:
        report += f"Discarded Linkers for {key.replace('_',' ')}: {discarded_counts[key]}\n"
        if discarded_examples[key]:
            report += f"Examples: {discarded_examples[key]}\n\n"

    if discarded_counts['mfe_delta_not_met'] > len(linkers_tried)*0.8:
        report += "\n⚠️ MAny linkers do not met the ΔMFE criteria. Consider reducing that parameter.\n"
    

    return results, report, cluster_labels
