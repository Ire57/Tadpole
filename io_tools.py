import os
import shutil
import tempfile
from PIL import Image

def save_and_plot_structures(seq, structure_unconstr, structure_constr,
                             rna1, linker, rna3, mut1_info, mfe_1, mfe_2):
    """
    Igual que antes, pero:
    - Usa carpeta temporal
    - Devuelve lista de rutas a imágenes generadas
    """
    linker_len = len(linker)
    full_len = len(seq)
    rna1_len = len(rna1)

    folder = tempfile.mkdtemp(prefix=f"linker_{linker_len}_")

    generated_images = []

    for tag, structure, mfe in [
        ("unconstrained", structure_unconstr, mfe_1),
        ("constrained", structure_constr, mfe_2)
    ]:
        # Guardar info
        filename = os.path.join(folder, f"result_{linker}_{tag}.txt")
        with open(filename, "w") as f:
            f.write(f"Linker = {linker_len}\n")
            f.write(f"Mutaciones RNA1: {mut1_info}\n")
            f.write(f"mfe = {mfe:.2f}\n")
            f.write(structure + "\n")
            f.write(seq + "\n")

        # Input para RNAplot
        fold_file = os.path.join(folder, f"{linker}_{tag}.fold")
        with open(fold_file, "w") as f:
            f.write(seq + "\n")
            f.write(structure + "\n")

        pre_command = f"1 {rna1_len} 10 GREEN omark " \
                      f"{rna1_len+1} {rna1_len+linker_len} 10 RED omark " \
                      f"{rna1_len+linker_len+1} {full_len} 10 BLUE omark"

        # Ejecutar RNAplot
        os.system(f'RNAplot --pre "{pre_command}" < {fold_file}')

        # Mover rna.ps a destino
        ps_file = os.path.join(folder, f"{linker}_{tag}_plot.ps")
        png_file = os.path.join(folder, f"{linker}_{tag}_plot.png")

        if os.path.exists("rna.ps"):
            shutil.move("rna.ps", ps_file)
            os.system(f"gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 -dEPSCrop -sOutputFile={png_file} {ps_file}")
            if os.path.exists(png_file):
                generated_images.append((tag, png_file))
        else:
            print(f"[WARN] No se generó 'rna.ps' para {tag}")

    return generated_images
