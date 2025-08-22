import requests
import json

url = "http://127.0.0.1:5000/"

sequence = "AACGUGACCUAUCCCAUUACG"

def test_api_secondary_structure():
    """Sends a POST request and prints the response."""
    payload = {"sequence": sequence}

    try:
        response = requests.post(url+"predict_secondary_structure", json=payload)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        print(json.dumps(response.json(), indent=2))
        
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")

def test_api_ga():
    """Sends a POST request and prints the response."""
    payload =     {
        "rna1": "ACCAGUGUGCGGAUGAUAACUACUGACGAAAGAGUCAUCGACUCAGUUAGUGGUUGGAUGUAGUCACAUUAGU",
        "rna3": "AGUUGGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACCAAAA",
        "struct1": "((.(((((((......(((((((((((....(((((...))))).)))))))))))......).)))))).))",
        "constraint": "........................................................................................<<<<<<<...............................",
        "mutable_rna1": [1, 2, 4, 5, 6, 7, 8, 9],
        "watched_positions": [56],
        "output_dir": "ga_results"
    }


    try:
        response = requests.post(url+"run_ga_search", json=payload)
        # Raise an exception for bad status codes (4xx or 5xx)
        response.raise_for_status()
        print(json.dumps(response.json(), indent=2))
        
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    test_api_secondary_structure()
    test_api_ga()
