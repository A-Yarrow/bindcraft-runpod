# Base image for RTX6000 Ada: RunPod PyTorch with CUDA 12.2 for Ada GPUs
FROM runpod/pytorch:2.4.0-py3.11-cuda12.2.0-devel-ubuntu22.04

# Install base tools and libraries
RUN apt-get update && apt-get install -y \
    wget \
    git \
    libgfortran5 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda (lighter than Miniforge3)
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $CONDA_DIR && \
    rm miniforge.sh

# Set up Conda and Mamba
RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

# Clone BindCraft repository
RUN git clone --branch dev --single-branch https://github.com/A-Yarrow/bindcraft-runpod.git /app/bindcraft

WORKDIR /app/bindcraft
# Install BindCraft (no PyRosetta or weights)
RUN chmod +x install_bindcraft.sh && \
    bash install_bindcraft.sh --cuda '12.2' --pkg_manager 'mamba' 

# Set permissions on startup script and notebook
RUN chmod 755 /app/bindcraft/start.sh
RUN chmod 644 /app/bindcraft/bindcraft-runpod-start.ipynb

# Environment variables for JAX / CUDA
ENV XLA_PYTHON_CLIENT_MEM_FRACTION=0.8
ENV XLA_FLAGS="--xla_gpu_enable_command_buffer=false"

# Expose Jupyter port
EXPOSE 8888

# Default command
CMD ["/app/bindcraft/start.sh"]
