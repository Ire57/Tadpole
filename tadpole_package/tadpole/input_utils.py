from Bio import AlignIO
from io import StringIO

def parse_fasta_msa(uploaded_file):
    """
    Parses an uploaded file, attempting to read it as a FASTA or Clustal Multiple Sequence Alignment (MSA).
    
    This function first decodes the binary content of the uploaded file into a string.
    It then tries to parse this content as a FASTA alignment using Biopython's AlignIO. If
    that fails, it attempts to parse it as a Clustal alignment. After successful parsing,
    it performs several validation checks to ensure the alignment is valid for downstream
    MSA-based analyses.
    
    :param uploaded_file: An object representing the uploaded file, expected to have a
                          `read()` method that returns binary content (e.g., `streamlit.uploaded_file` object).
    
    :raises ValueError:
        - If the uploaded file is empty.
        - If the file content is neither a valid FASTA nor a valid Clustal format.
        - If the parsed alignment contains fewer than two sequences.
        - If not all sequences in the alignment have the same length.
        - If any of the sequences in the alignment are empty.
    
    :returns: A list of strings, where each string is a sequence from the validated MSA (list of str).
    """
    # Decode binary file content to text
    content = uploaded_file.read().decode("utf-8").strip()
    if not content:
        raise ValueError("The uploaded file is empty.")
    
    handle = StringIO(content)
    
    try:
        alignment = AlignIO.read(handle, "fasta")  # Try parsing as FASTA
    except Exception:
        # If FASTA parsing fails, try Clustal
        handle.seek(0)  # Reset buffer position for second read attempt
        try:
            alignment = AlignIO.read(handle, "clustal")
        except Exception:
            raise ValueError("The file is neither a valid FASTA nor a valid Clustal format.")
    
    # Validate number of sequences
    if len(alignment) < 2:
        raise ValueError("The file must contain at least two sequences for an MSA.")
    
    # Validate all sequences have the same length
    lengths = {len(record.seq) for record in alignment}
    if len(lengths) != 1:
        raise ValueError("Sequences are not of equal length. This is not a valid MSA.")
    
    # Validate sequences are not empty
    if 0 in lengths:
        raise ValueError("Sequences cannot be empty.")
    
    # Return list of sequences as strings
    return [str(record.seq) for record in alignment]


def get_msa_input():
    """
    This function is Not used on the Tadpole Software, it is to be used by developers that want to run 
    the code by terminal instead of the Tadpole interface, when building upon it, for example.
    Prompts the user to enter Multiple Sequence Alignment (MSA) sequences in FASTA format
    via the console and parses them.

    This function provides instructions for the expected input format and then
    captures a single line of user input. It expects a simplified FASTA-like
    format where sequence headers start with '>' followed by a name, and the
    sequence itself is on the same line, space-separated from the name.

    :returns: A list of strings, where each string is a sequence entered by the user (list of str).
    """
    print(
        "Enter your MSA sequences in FASTA format"
        "(e.g., '>dm6 AGGUC >dm7 UUCGA ...'):"
        )
    msa_input = input("Enter FASTA sequence: ").strip()
    # Splits the input string by '>' to separate individual sequences, then slices to remove empty first element
    sequences = msa_input.split('>')[1:] 
    # For each sequence block, splits by space and takes the second element (the sequence string itself)
    msa = [seq.split()[1] for seq in sequences]
    return msa



