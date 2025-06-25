FROM continuumio/miniconda3

WORKDIR /app

# Instala dependencias del sistema necesarias para:
# - Ghostscript
# - WeasyPrint (libpango, libcairo, etc.)
# - Imágenes (libjpeg, libpng)
RUN apt-get update && apt-get install -y \
    ghostscript \
    build-essential \
    libglib2.0-0 \
    libpango1.0-0 \
    libcairo2 \
    libgdk-pixbuf2.0-0 \
    libffi-dev \
    libjpeg-dev \
    libpng-dev \
    libgobject-2.0-0 \
    shared-mime-info \
    fonts-liberation \
    && rm -rf /var/lib/apt/lists/*

# Añade canales conda
RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda update -n base -c defaults conda -y

# Instala ViennaRNA, Streamlit y dependencias base
RUN conda install -y viennarna streamlit

# Instala dependencias Python adicionales
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copia tu código
COPY . .

# Expone el puerto que usa Streamlit
EXPOSE 8501

# Comando para ejecutar la aplicación
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
