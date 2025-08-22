FROM continuumio/miniconda3

WORKDIR /app

# Actualiza e instala dependencias del SO necesarias para Viennarna y Ghostscript
RUN apt-get update && apt-get install -y \
    ghostscript \
    build-essential \
    libglib2.0-0 \
    libpango1.0-0 \
    libpangoft2-1.0-0 \
    libjpeg-dev \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Añade canales conda
RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda update -n base -c defaults conda -y

# Instala viennarna, streamlit y demás paquetes Python
RUN conda install viennarna streamlit -y

# Copia requirements y instala con pip
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copia el código de la app
COPY . .

# Expone el puerto de Streamlit
EXPOSE 8501

# Comando para ejecutar la app
CMD ["streamlit", "run", "tadpole_package/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
