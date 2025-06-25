FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    cmake \
    libglib2.0-0 \
    libpangocairo-1.0-0 \
    libpangoft2-1.0-0 \
    libjpeg-dev \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Install ViennaRNA from source
RUN curl -L https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.4.tar.gz | tar xz \
 && cd ViennaRNA-2.6.4 \
 && mkdir build && cd build \
 && cmake .. \
 && make \
 && make install \
 && cd ../.. \
 && rm -rf ViennaRNA-2.6.4

# Set environment path
ENV PATH="/usr/local/bin:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"

# Install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Add app code
COPY . /app
WORKDIR /app

# Set Streamlit start command
CMD ["streamlit", "run", "app.py"]
