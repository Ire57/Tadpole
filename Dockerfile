# Use miniconda base image
FROM continuumio/miniconda3

# Create a working directory
WORKDIR /app

# Copy your Python requirements file (make sure it’s in your repo)
COPY requirements.txt .

# Install ViennaRNA and Streamlit using conda, plus other pip packages
RUN conda install -c conda-forge viennarna=2.6.4 streamlit -y && \
    pip install --no-cache-dir -r requirements.txt

# Copy all your app files into the container
COPY . .

# Expose the port Streamlit uses
EXPOSE 8501

# Run the Streamlit app
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
