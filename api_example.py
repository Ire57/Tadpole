from flask import Flask, request, jsonify
import sys
import os

import numpy as np

from tadpole import predict_secundary_structure, run_genetic_algorithm_search
#from structure import predict_secundary_structure
#from genetic_algorithm import run_genetic_algorithm_search

app = Flask(__name__)

@app.route('/predict_secondary_structure', methods=['POST'])
def predict_structure():
    """
    API endpoint to predict the secondary structure of an RNA sequence.
    This endpoint expects a JSON payload with a 'sequence' key.
    
    Example JSON payload:
    {
      "sequence": "AACGUGACCUAUCCCAUUACG"
    }
    """

    if not request.json or 'sequence' not in request.json:
        return jsonify({"error": "Invalid request. Please provide a JSON object with a 'sequence' key."}), 400

    rna_sequence = request.json['sequence']

    if not isinstance(rna_sequence, str) or not rna_sequence:
        return jsonify({"error": "The 'sequence' must be a non-empty string."}), 400

    try:
        predicted_structure, mfe = predict_secundary_structure(rna_sequence)
        
        response_data = {
            "sequence": rna_sequence,
            "structure": predicted_structure,
            "mfe": mfe,
        }

        return jsonify(response_data), 200

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return jsonify({"error": "An internal server error occurred."}), 500


@app.route('/run_ga_search', methods=['POST'])
def run_ga_search():
    """
    API endpoint to run the genetic algorithm search for RNA switch designs.

    Example JSON payload:
    {
      "rna1": "ACCAGUGUGCGGAUGAUAACUACUGACGAAAGAGUCAUCGACUCAGUUAGUGGUUGGAUGUAGUCACAUUAGU",
      "rna3": "AGUUGGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACCAAAA",
      "struct1": "((.(((((((......(((((((((((....(((((...))))).)))))))))))......).)))))).))",
      "constraint": "........................................................................................<<<<<<<...............................",
      "mutable_rna1": [1, 2, 4, 5, 6, 7, 8, 9],
      "watched_positions": [56],
      "output_dir": "ga_results"
    }
    """
    if not request.json:
        return jsonify({"error": "Invalid request. Please provide a JSON payload."}), 400

    # Required parameters check
    required_params = ['rna1', 'rna3', 'struct1', 'constraint', 'mutable_rna1', 'watched_positions', 'output_dir']
    if not all(param in request.json for param in required_params):
        return jsonify({"error": f"Missing one or more required parameters. Required: {', '.join(required_params)}"}), 400

    try:
        # Extract parameters with defaults
        rna1 = request.json['rna1']
        rna3 = request.json['rna3']
        struct1 = request.json['struct1']
        constraint = request.json['constraint']
        mutable_rna1 = request.json['mutable_rna1']
        watched_positions = request.json['watched_positions']
        output_dir = request.json['output_dir']

        # Optional parameters with default values from function signature
        use_mutations = request.json.get('use_mutations', True)
        mfe_delta = request.json.get('mfe_delta', 0)
        max_pairings_rna1_rna3 = request.json.get('max_pairings_rna1_rna3', 5)
        max_structure_changes = request.json.get('max_structure_changes', 6)
        num_mut = request.json.get('num_mut', 0)
        
        population_size = request.json.get('population_size', 50)
        generations = request.json.get('generations', 50)
        linker_length_ga = request.json.get('linker_length_ga', 7)
        elitism_rate = request.json.get('elitism_rate', 0.1)
        mutation_rate_rna1 = request.json.get('mutation_rate_rna1', 0.02)
        mutation_rate_linker = request.json.get('mutation_rate_linker', 0.05)
        tournament_size = request.json.get('tournament_size', 3)

        # Call the genetic algorithm function
        all_final_results, final_report_lines, cluster_labels = run_genetic_algorithm_search(
            rna1=rna1,
            rna3=rna3,
            struct1=struct1,
            constraint=constraint,
            mutable_rna1=mutable_rna1,
            watched_positions=watched_positions,
            output_dir=output_dir,
            use_mutations=use_mutations,
            mfe_delta=mfe_delta,
            max_pairings_rna1_rna3=max_pairings_rna1_rna3,
            max_structure_changes=max_structure_changes,
            num_mut=num_mut,
            population_size=population_size,
            generations=generations,
            linker_length_ga=linker_length_ga,
            elitism_rate=elitism_rate,
            mutation_rate_rna1=mutation_rate_rna1,
            mutation_rate_linker=mutation_rate_linker,
            tournament_size=tournament_size,
        )

        # Prepare and return the response
        response_data = {
            "search_method": "genetic_algorithm",
            "rna1": rna1,
            "rna3": rna3,
            "struct1": struct1,
            "constraint": constraint,
            "mutable_rna1": mutable_rna1,
            "watched_positions": watched_positions,
            "all_final_results": all_final_results,
            "kwargs": {
                "use_mutations": use_mutations,
                "mfe_delta": mfe_delta,
                "max_pairings_rna1_rna3": max_pairings_rna1_rna3,
                "max_structure_changes": max_structure_changes,
                "num_mut": num_mut,
                "population_size": population_size,
                "generations": generations,
                "linker_length_ga": linker_length_ga,
                "elitism_rate": elitism_rate,
                "mutation_rate_rna1": mutation_rate_rna1,
                "mutation_rate_linker": mutation_rate_linker,
                "tournament_size": tournament_size,
            },
            # 'cluster_labels' can return [] or Numpy array.
            # As an array is not JSON serializable, convert to list if is an array.
            # If an empty list is supplied, converting to array first ensures it has 'tolist' method.
            "cluster_labels": np.array(cluster_labels).tolist(),
        }
        
        return jsonify(response_data), 200

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return jsonify({"error": "An internal server error occurred."}), 500



if __name__ == '__main__':
    # Get the port from the environment, defaulting to 5000.
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)