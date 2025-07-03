# Use RunPod-compatible image with CUDA 12.4 and PyTorch
FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

# Install base tools and Miniconda (lightweight)
RUN apt-get update && apt-get install -y \
    wget \
    git \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda (lighter than Miniforge3)
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p $CONDA_DIR && \
    rm miniconda.sh

# Set up Conda and Mamba
RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

# Clone BindCraft repository
RUN git clone https://github.com/A-Yarrow/bindcraft-runpod.git /app/bindcraft

COPY install_bindcraft.sh /app/bindcraft/install_bindcraft.sh
WORKDIR /app/bindcraft
# Install BindCraft (no PyRosetta or weights)
RUN chmod +x install_bindcraft.sh && \
    bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'mamba'

# Set permissions on startup script and notebook
RUN chmod 755 /app/bindcraft/start.sh
RUN chmod 644 /app/bindcraft/bindcraft-runpod-start.ipynb

# Expose Jupyter port
EXPOSE 8888

# Default command
CMD ["/app/bindcraft/start.sh"] 
# For development, launch Jupyter notebook
#EXPOSE 8888
#CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--allow-root", "--NotebookApp.token=", "--NotebookApp.password="]