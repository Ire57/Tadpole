import sbol3
import uuid
import json


def export_single_linker(doc, result, index=0):
    """
    Exports a single linker design to an SBOL3 Document.

    This function takes a dictionary of linker design results and converts
    them into SBOL3 objects (Sequence and Component). It adds these objects
    to a provided SBOL3 Document, ensuring each has a unique identifier.
    The function also includes additional design information as a JSON string
    in the Component's description property.

    :param doc: The SBOL3 Document to which the new component will be added. (sbol3.Document)
    :param result: A dictionary containing all design information for a single linker. (dict)
    :param index: A numerical index used to help generate unique identifiers. (int)

    :returns: The modified SBOL3 Document with the new linker component added. (sbol3.Document)
    """
    # Extract fields from result
    sequence_str = result.get('rna1_mutated_seq', 'rna1')
    structure_unconstrained = result.get('structure_unconstrained', '')
    structure_constrained = result.get('structure_constrained', '')
    mfe_unconstrained = result.get('mfe_1', None)
    mfe_constrained = result.get('mfe_2', None)
    mut_info = result.get("mut1_info", "")

    # Generate unique IDs
    comp_id = f"linker_{index}_{uuid.uuid4().hex}"
    seq_id = f"{comp_id}_sequence_{uuid.uuid4().hex}"

    # Create Sequence
    sequence = sbol3.Sequence(
        identity=seq_id,
        elements=sequence_str,
        encoding="https://identifiers.org/edam:format_1207"
    )
    doc.add(sequence)

    # Create Component
    component = sbol3.Component(
        identity=comp_id,
        types=["http://www.biopax.org/release/biopax-level3.owl#RNA"]
    )
    component.roles = ["http://identifiers.org/so/SO:0000356"]
    component.sequences = [sequence]

    # Store extra features in description as JSON
    extra_info = {
        "structure_unconstrained": structure_unconstrained,
        "structure_constrained": structure_constrained,
        "mfe_unconstrained": mfe_unconstrained,
        "mfe_constrained": mfe_constrained,
        "mutations_sre": mut_info
    }
    component.description = json.dumps(extra_info)

    doc.add(component)
    return doc


