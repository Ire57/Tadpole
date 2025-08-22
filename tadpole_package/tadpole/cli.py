import argparse
import os
import sys
import traceback

# Imports all necessary modules from the tadpole package.
from tadpole.genetic_algorithm import run_genetic_algorithm_search
from tadpole.search import run_linker_search
from tadpole.structure import predict_secundary_structure
from tadpole.io_tools import save_and_plot_structures
from tadpole.input_utils import parse_fasta_msa
from tadpole.conservation import clasify_conservation
from tadpole.visualizacion import generate_rnaplot_with_colours
from tadpole.rna_cluster import matriz_rmsd, cluster_structures, organise_per_clusters, compute_metrics, visualise_metrics

# A simple function to log messages to the console.
def log_func(message):
    print(message, file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='TADPOLE: An RNA switch design toolkit.')
    
    # Configures all available subcommands.
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    ## Subcommand: `ga_search`
    ga_parser = subparsers.add_parser('ga_search', help='Runs the genetic algorithm.')
    ga_parser.add_argument('--rna1_seq', required=True, help='The RNA1 (SRE) sequence.')
    ga_parser.add_argument('--rna3_seq', required=True, help='The RNA3 (aptamer) sequence.')
    ga_parser.add_argument('--struct1', required=True, help='The RNA1 (SRE) structure in dot-bracket notation.')
    ga_parser.add_argument('--constraint', required=True, help='The aptamer structural constraint for the ON state.')
    ga_parser.add_argument('--linker_length_ga', type=int, default=7, help='The length of the linker to use.')
    ga_parser.add_argument('--num_mut', type=int, default=0, help='The number of mutations to introduce into the SRE.')
    ga_parser.add_argument('--output_dir', default='outputs', help='The directory to save results.')
    ga_parser.add_argument(
        '--mutable_rna1', 
        nargs='*', # Collects a list of arguments
        type=int,  # Ensures each item in the list is an integer
        default=[], 
        help='A list of mutable positions in RNA1.'
    )
    ga_parser.add_argument(
        '--watched_positions', 
        nargs='*', # Collects a list of arguments
        type=int,  # Ensures each item in the list is an integer
        default=[], 
        help='The watched positions in the SRE.'
    )
    ga_parser.add_argument(
        '--no_mutations',
        action='store_true',
        help='Disable mutations (default: mutations enabled).'
    )
    ga_parser.add_argument('--mfe_delta', type=float, default=0, help='The delta MFE value.')
    ga_parser.add_argument('--max_pairings_rna1_rna3', type=int, default=5, help='Maximum number of pairings between RNA1 and RNA3.')
    ga_parser.add_argument('--max_structure_changes', type=int, default=6, help='Maximum number of changes in structure.')
    ga_parser.add_argument('--verbose', action='store_true', default=True, help='Enable verbose output.')
    ga_parser.add_argument('--population_size', type=int, default=50, help='Population size for the genetic algorithm.')
    ga_parser.add_argument('--generations', type=int, default=100, help='Number of generations for the genetic algorithm.')
    ga_parser.add_argument('--elitism_rate', type=float, default=0.1, help='Elitism rate for the genetic algorithm.')
    ga_parser.add_argument('--mutation_rate_rna1', type=float, default=0.02, help='Mutation rate for RNA1.')
    ga_parser.add_argument('--mutation_rate_linker', type=float, default=0.05, help='Mutation rate for the linker.')
    ga_parser.add_argument('--tournament_size', type=int, default=3, help='Tournament size for selection.')
    
    
    # Subcommand: `linker_search`
    linker_parser = subparsers.add_parser('linker_search', help='Runs the exhaustive linker search.')
    linker_parser.add_argument('--rna1_seq', required=True, help='The RNA1 (SRE) sequence.')
    linker_parser.add_argument('--rna3_seq', required=True, help='The RNA3 (aptamer) sequence.')
    linker_parser.add_argument('--struct1', required=True, help='The RNA1 (SRE) structure in dot-bracket notation.')
    linker_parser.add_argument('--constraint', required=True, help='The aptamer structural constraint for the ON state.')
    linker_parser.add_argument('--linker_min', type=int, default=7, help='Minimum linker length to search for.')
    linker_parser.add_argument('--linker_max', type=int, default=10, help='Maximum linker length to search for.')
    
    linker_parser.add_argument(
        '--mutable_rna1', 
        nargs='*',
        type=int,  
        default=[], 
        help='A list of mutable positions in RNA1.'
    )
    linker_parser.add_argument(
        '--watched_positions', 
        nargs='*',
        type=int, 
        default=[], 
        help='The watched positions in the SRE.'
    )
    linker_parser.add_argument('--use_mutations', action='store_true', default=False, help='Allows mutations in the search.')
    linker_parser.add_argument('--mfe_delta', type=float, default=2.0, help='The MFE delta to filter structures.')
    linker_parser.add_argument('--max_pairings_rna1_rna3', type=int, default=10, help='The maximum number of pairings between RNA1 and RNA3.')
    linker_parser.add_argument('--max_structure_changes', type=int, default=6, help='Maximum number of changes in structure.')
    linker_parser.add_argument('--num_mut', type=int, default=0, help='The number of mutations to introduce into the SRE.')
    linker_parser.add_argument('--verbose', action='store_true', default=True, help='Enable verbose output.')
    linker_parser.add_argument('--output_dir', default='outputs', help='The directory to save results.')

    ## Subcommand: `predict_structure`
    predict_parser = subparsers.add_parser('predict_structure', help='Predicts the MFE secondary structure.')
    predict_parser.add_argument('--sequence', required=True, help='The RNA sequence.')
    predict_parser.add_argument('--output_dir', default='outputs', help='The directory to save results.')

    ## Subcommand: `plot_structure`
    plot_parser = subparsers.add_parser('plot_structure', help='Plots an RNA secondary structure from a sequence and a dot-bracket notation.')
    plot_parser.add_argument('--seq', required=True, help='The RNA sequence.')
    plot_parser.add_argument('--struct_unconstrained', required=True, help='The unconstrained RNA structure in dot-bracket notation.')
    plot_parser.add_argument('--struct_constrained', default=None, help='The constrained RNA structure in dot-bracket notation.')
    plot_parser.add_argument('--rna1', default=None, help='The RNA1 sequence segment.')
    plot_parser.add_argument('--linker', default=None, help='The linker sequence segment.')
    plot_parser.add_argument('--rna3', default=None, help='The RNA3 sequence segment.')
    plot_parser.add_argument('--mut1_info', default=None, help='Mutation information for RNA1.')
    plot_parser.add_argument('--mfe_1', type=float, default=None, help='MFE of the unconstrained structure.')
    plot_parser.add_argument('--mfe_2', type=float, default=None, help='MFE of the constrained structure.')
    plot_parser.add_argument('--output_dir', default='outputs', help='The directory to save the plot.')
    

    ## Subcommand: `analyze_conservation`
    conservation_parser = subparsers.add_parser('analyze_conservation', help='Analyzes and plots MSA conservation.')
    conservation_parser.add_argument('--msa_file', required=True, help='The path to the MSA file in FASTA or Clustal format.')
    conservation_parser.add_argument('--output_dir', default='outputs', help='The directory where the results will be saved.')

    ## Subcommand: `cluster_structures`
    cluster_parser = subparsers.add_parser('cluster_structures', help='Clusters similar structures from a file.')
    cluster_parser.add_argument('--sequences', required=True, help='A file containing sequences (one per line).')
    cluster_parser.add_argument('--structures', required=True, help='A file containing dot-bracket structures (one per line).')
    cluster_parser.add_argument('--output_dir', default='outputs', help='The directory where the results will be saved.')
    cluster_parser.add_argument('--eps', type=float, default=100.0, help='The maximum distance between two samples for one to be considered as in the neighborhood of the other. (DBSCAN epsilon)')
    cluster_parser.add_argument('--min_samples', type=int, default=1, help='The number of samples in a neighborhood for a point to be considered as a core point. (DBSCAN min_samples)')

    # --- Main logic for executing the commands ---
    args = parser.parse_args()
    
    if hasattr(args, 'output_dir'):
        os.makedirs(args.output_dir, exist_ok=True)

    if args.command == 'ga_search':
        try:
            run_genetic_algorithm_search(
                rna1=args.rna1_seq, 
                rna3=args.rna3_seq, 
                struct1=args.struct1,
                constraint=args.constraint,
                mutable_rna1=args.mutable_rna1,
                watched_positions=args.watched_positions,
                output_dir=args.output_dir,
                num_mut=args.num_mut,
                linker_length_ga=args.linker_length_ga,
                log_func=log_func, 
                use_mutations = not args.no_mutations,
                mfe_delta=args.mfe_delta,
                max_pairings_rna1_rna3=args.max_pairings_rna1_rna3,
                max_structure_changes=args.max_structure_changes,
                verbose=args.verbose,
                population_size=args.population_size,
                generations=args.generations,
                elitism_rate=args.elitism_rate,
                mutation_rate_rna1=args.mutation_rate_rna1,
                mutation_rate_linker=args.mutation_rate_linker,
                tournament_size=args.tournament_size
            )
            print(f"✅ Búsqueda GA completada. Los resultados se encuentran en la carpeta '{args.output_dir}'.")
        except Exception as e:
            print(f"❌ Ocurrió un error en el comando 'ga_search': {e}", file=sys.stderr)
            traceback.print_exc()

    elif args.command == 'linker_search':
        try:
            # Generate the linker_lengths range to match app.py's behavior
            linker_lengths_range = range(args.linker_min, args.linker_max + 1)
            
            # Pass this range to the function
            run_linker_search(
                rna1=args.rna1_seq,
                rna3=args.rna3_seq,
                struct1=args.struct1,
                constraint=args.constraint,
                # Pass the generated range here
                linker_lengths=linker_lengths_range,
                use_mutations=args.use_mutations,
                mfe_delta=args.mfe_delta,
                max_pairings_rna1_rna3=args.max_pairings_rna1_rna3,
                output_dir=args.output_dir,
                log_func=log_func,
                mutable_rna1=args.mutable_rna1, 
                watched_positions=args.watched_positions,
                max_structure_changes=args.max_structure_changes,
                num_mut=args.num_mut,
                verbose=args.verbose
            )
        except Exception as e:
            print(f"❌ Ocurrió un error en el comando 'linker_search': {e}", file=sys.stderr)
            traceback.print_exc()

    elif args.command == 'predict_structure':
        try:
            structure, energy = predict_secundary_structure(args.sequence)
            print(f"✅ Estructura MFE: {structure}")
            print(f"✅ Energía MFE: {energy} kcal/mol")
        except Exception as e:
            print(f"❌ Ocurrió un error en el comando 'predict_structure': {e}", file=sys.stderr)
            traceback.print_exc()

    elif args.command == 'plot_structure':
        try:
            save_and_plot_structures(
                seq=args.seq,
                structure_unconstr=args.struct_unconstrained,
                structure_constr=args.struct_constrained,
                rna1=args.rna1,
                linker=args.linker,
                rna3=args.rna3,
                mut1_info=args.mut1_info,
                mfe_1=args.mfe_1,
                mfe_2=args.mfe_2,
                folder_prefix=args.output_dir, 
            )
            print(f"✅ El gráfico se ha guardado en '{args.output_dir}'.")
        except Exception as e:
            print(f"❌ Ocurrió un error en el comando 'plot_structure': {e}", file=sys.stderr)
            traceback.print_exc()

    elif args.command == 'analyze_conservation':
        try:
            # 1. Open the file in binary read mode ('rb').
            #    This creates a file object that has the .read() and .decode() methods.
            with open(args.msa_file, 'rb') as f:
                msa = parse_fasta_msa(f)
            
            # 2. The sequence for the plot must be without hyphens
            sequence_for_plot = msa[0].replace('-', '')
            secundary_structure_for_plot, _ = predict_secundary_structure(sequence_for_plot)
            
            # 3. Classify the conservation in the gapped MSA
            conservation_groups_gapped = clasify_conservation(msa)
            
            # 4. Map the indices from the gapped MSA to the ungapped sequence
            conservation_groups_for_plot = {}
            # Create a map to translate gapped indices to ungapped indices
            gap_free_map = {}
            gap_free_index = 1
            for i, base in enumerate(msa[0]):
                if base != '-':
                    gap_free_map[i + 1] = gap_free_index
                    gap_free_index += 1

            for category, indices in conservation_groups_gapped.items():
                
                # Map the original indices to the new indices
                new_indices = []
                for idx in indices:
                    # Only add the index if it exists in the ungapped sequence
                    if idx in gap_free_map:
                        new_indices.append(gap_free_map[idx])

                if new_indices:
                    conservation_groups_for_plot[category] = new_indices
            
            # 5. Generate the plot with the corrected indices
            generate_rnaplot_with_colours(
                sequence=sequence_for_plot, 
                secundary_structure=secundary_structure_for_plot, 
                grups=conservation_groups_for_plot, 
                tag=os.path.basename(args.msa_file).split('.')[0], 
                output_dir=args.output_dir
            )
            print(f"✅ Conservation analysis complete. The plot has been saved to '{args.output_dir}'.")
        except Exception as e:
            print(f"❌ An error occurred in the 'analyze_conservation' command: {e}", file=sys.stderr)
            traceback.print_exc()

    elif args.command == 'cluster_structures':
        try:
            with open(args.sequences, 'r') as f:
                sequences = [line.strip() for line in f]
            with open(args.structures, 'r') as f:
                structures = [line.strip() for line in f]
            
            rmsd_matrix = matriz_rmsd(structures)
            # Use the values from the command line arguments
            labels = cluster_structures(rmsd_matrix, eps=args.eps, min_samples=args.min_samples)
            unique_labels = set(labels)
            num_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
            
            # The results are organized and saved.
            organise_per_clusters(sequences, structures, labels, rmsd_matrix, output_base=f"{args.output_dir}/clusters")
            metrics = compute_metrics(structures, times=[])
            visualise_metrics(metrics, args.output_dir)
            
            print(f"✅ Clustering completado. Se encontraron {num_clusters} clústeres. Los resultados se encuentran en '{args.output_dir}'.")
        except Exception as e:
            print(f"❌ Ocurrió un error en el comando 'cluster_structures': {e}", file=sys.stderr)
            traceback.print_exc()

if __name__ == '__main__':
    main()