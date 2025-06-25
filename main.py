from input_utils import get_msa_input, get_msa_name
from conservacion import calculate_conservation, clasify_conservation, conservation_category
from visualizacion import generate_rnaplot_with_colours, generate_pre_command_per_sequence
from structure import predict_secundary_structure
import os

# 1. Entrada
msa = get_msa_input()
msa_name = get_msa_name()

# 2. Conservación
conservation_scores = calculate_conservation(msa)

# 3. Estructura
sequence = msa[0].replace("-", "")
estructura_secundaria, energia = predict_secundary_structure(sequence)

# 4. Agrupación por conservación
grupos = clasify_conservation(msa)
# 5. Visualización coloreada
img = generate_rnaplot_with_colours(sequence, estructura_secundaria, grupos, tag=msa_name)

#FOR ALL SEQUENCES

output_dir = f"rna_outputs_{msa_name}"
os.makedirs(output_dir, exist_ok=True)

lw = 10
msa_len = len(msa[0])

# Categorizar columnas
col_categories = {idx: conservation_category(msa, idx) for idx in range(msa_len)}

colors = {
    "siempre": (255, 0, 0),
    "menos2": (139, 0, 0),
    "menos4": (0, 0, 255),
    "resto": (255, 255, 0)
}
# Recorrer cada secuencia del MSA
for i, msa_seq in enumerate(msa):
    sequence = msa_seq.replace("-", "")
    estructura_secundaria, energia = predict_secundary_structure(sequence)
    print(f"Secuencia {i}: {estructura_secundaria} ({energia:.2f} kcal/mol)")

    pre_command = generate_pre_command_per_sequence(msa_seq, col_categories, colors, lw)

    fold_filename = os.path.join(output_dir, f"rna_structure_{i}.fold")
    with open(fold_filename, "w") as f:
        f.write(sequence + "\n" + estructura_secundaria + "\n")

    os.system(f'RNAplot --pre "{pre_command}" < {fold_filename}')

    ps_name = os.path.join(output_dir, f"rna_structure_{i}.ps")
    png_name = os.path.join(output_dir, f"rna_structure_{i}.png")
    os.system(f"mv rna.ps {ps_name}")
    os.system(
        f'gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 '
        f'-dEPSCrop -sOutputFile={png_name} {ps_name}'
    )