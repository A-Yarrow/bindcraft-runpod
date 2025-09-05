# 1 Base image for This image is optimized for JAX development and includes:
#Ubuntu 22.04.2 LTS with Python 3.10.12, NVIDIA CUDA® 12.1.1, NVIDIA cuDNN 8.9.4
#NVIDIA cuBLAS 12.2.5, NVIDIA NCCL 2.18.3, NVIDIA DALI® 1.28.0, TransformerEngine 0.13.0.dev0+03202c3
FROM nvcr.io/nvidia/jax:23.08-py3

# 2 Install base tools and libraries
RUN apt-get update && apt-get install -y \
    wget \
    git \
    libgfortran5 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 3-5 Install Miniconda (lighter than Miniforge3)
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $CONDA_DIR && \
    rm miniforge.sh

# 6 Set up Conda and Mamba
RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

# 7-8 Clone BindCraft repository
RUN git clone --branch dev --single-branch https://github.com/A-Yarrow/bindcraft-runpod.git /app/bindcraft

WORKDIR /app/bindcraft
# 9 Install BindCraft (no PyRosetta or weights)
RUN chmod +x install_bindcraft.sh && \
    bash install_bindcraft.sh 

# 10 Set permissions on startup script and notebook
RUN chmod 755 /app/bindcraft/start.sh
RUN chmod 644 /app/bindcraft/bindcraft-runpod-start.ipynb

# 11-12 Environment variables for JAX / CUDA
ENV XLA_PYTHON_CLIENT_MEM_FRACTION=0.8
ENV XLA_FLAGS="--xla_gpu_enable_command_buffer=false"

# 13 Expose Jupyter port
EXPOSE 8888

# 14 Default command
CMD ["/app/bindcraft/start.sh"]
