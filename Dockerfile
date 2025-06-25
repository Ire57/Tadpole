FROM python:3.10-slim

# Install system packages including ViennaRNA
RUN apt-get update && apt-get install -y \
    viennarna \
    build-essential \
    libglib2.0-0 \
    libpangocairo-1.0-0 \
    libpangoft2-1.0-0 \
    libjpeg-dev \
    libpng-dev \
    && apt-get clean

# Set working directory
WORKDIR /app

# Copy everything
COPY . /app

# Install Python packages
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Expose Streamlit port
EXPOSE 8501

# Command to run your Streamlit app
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.enableCORS=false"]
