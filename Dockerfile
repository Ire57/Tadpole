FROM continuumio/miniconda3

WORKDIR /app

# Add conda-forge and bioconda channels for RNA tools
RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda update -n base -c defaults conda -y

# Install viennarna (latest available), streamlit, and python packages
RUN conda install viennarna streamlit -y

# Copy requirements and install the rest with pip
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app code
COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
