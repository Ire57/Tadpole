# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libglib2.0-0 \
    libpangocairo-1.0-0 \
    libpangoft2-1.0-0 \
    libjpeg-dev \
    libpng-dev \
    curl \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install ViennaRNA from source
RUN mkdir /vienna_tmp \
 && curl -L https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.4.tar.gz | tar xz -C /vienna_tmp \
 && cd /vienna_tmp/ViennaRNA-2.6.4 \
 && mkdir build && cd build \
 && cmake .. \
 && make \
 && make install \
 && cd / \
 && rm -rf /vienna_tmp

# Set working directory in the container
WORKDIR /app

# Copy requirements.txt and install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the app source code
COPY . .

# Expose the port Streamlit runs on
EXPOSE 8501

# Run the Streamlit app
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
