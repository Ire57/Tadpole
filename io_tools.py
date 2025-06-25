import os
from PIL import Image
from IPython.display import display

import os
import io
from PIL import Image
import streamlit as st

def save_and_plot_structures(seq, structure_unconstr, structure_constr,
                             rna1, linker, rna3, mut1_info, mfe_1, mfe_2,
                             folder_prefix="/tmp/propuestas"):  # Carpeta temporal

    linker_len = len(linker)
    full_len = len(seq)
    rna1_len = len(rna1)

    folder = f"{folder_prefix}/linker_{linker_len}"
    os.makedirs(folder, exist_ok=True)

    for tag, structure, mfe in [
        ("unconstrained", structure_unconstr, mfe_1),
        ("constrained", structure_constr, mfe_2)
    ]:
        fold_file = f"{folder}/{linker}_{tag}.fold"
        with open(fold_file, "w") as f:
            f.write(seq + "\n")
            f.write(structure + "\n")

        pre_command = f"1 {rna1_len} 10 GREEN omark " \
                      f"{rna1_len+1} {rna1_len+linker_len} 10 RED omark " \
                      f"{rna1_len+linker_len+1} {full_len} 10 BLUE omark"

        os.system(f'RNAplot --pre "{pre_command}" < {fold_file}')

        ps_file = f"{folder}/{linker}_{tag}_plot.ps"
        png_file = f"{folder}/{linker}_{tag}_plot.png"

        if os.path.exists("rna.ps"):
            os.rename("rna.ps", ps_file)
            os.system(f"gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r150 "
                      f"-dEPSCrop -sOutputFile={png_file} {ps_file}")

            if os.path.exists(png_file):
                # Abrir la imagen en memoria y mostrarla con Streamlit
                with open(png_file, "rb") as img_file:
                    img_bytes = img_file.read()
                st.image(img_bytes, caption=f"{tag} structure")

                # Opcional: borrar archivos temporales
                os.remove(ps_file)
                os.remove(png_file)
        else:
            st.warning(f"No output from RNAplot for {tag} structure.")
