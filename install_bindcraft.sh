#!/bin/bash
set -e  # exit on any error

# Build with: docker build --progress=plain -t bindcraft:test .
################## BindCraft installation script optimized for RunPod 12.4 base image
################## Base: FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04
################## Tested on 2025-09-04

# Hardcoded config
pkg_manager="mamba"   # or "conda" if you prefer
cuda="12.1"

############################################################################################################
################## Initialization
SECONDS=0
install_dir=$(pwd)

# Check conda installation
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "Error: conda not installed or not initialised"; exit 1; }
echo "Conda is installed at: $CONDA_BASE"

############################################################################################################
################## Create BindCraft environment
echo "Creating BindCraft environment"
$pkg_manager create --name BindCraft python=3.10 -y || { echo "Error: Failed to create BindCraft env"; exit 1; }

# Activate environment
echo "Activating BindCraft environment"
source ${CONDA_BASE}/bin/activate ${CONDA_BASE}/envs/BindCraft
[ "$CONDA_DEFAULT_ENV" = "BindCraft" ] || { echo "Error: BindCraft environment not active"; exit 1; }
echo "BindCraft environment activated at ${CONDA_BASE}/envs/BindCraft"

############################################################################################################
################## Install conda requirements
echo "Installing conda packages"
$pkg_manager install \
  pip pandas matplotlib 'numpy<2.0.0' biopython scipy pdbfixer seaborn libgfortran5 tqdm \
  jupyter jupyterlab ffmpeg pyrosetta fsspec py3dmol \
  chex dm-haiku 'flax<0.10.0' dm-tree joblib ml-collections immutabledict optax \
  psutil copyparty \
  -c conda-forge --channel https://conda.graylab.jhu.edu -y \
|| { echo "Error: Failed to install conda packages"; exit 1; }

############################################################################################################
################## Install pip requirements
echo "Installing pip packages"
python -m pip install --upgrade pip wheel jupyter-server-proxy

# JAX for CUDA 12.1
python -m pip install --no-cache-dir \
  jax==0.4.28 \
  jaxlib==0.4.28+cuda12.cudnn89 \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

############################################################################################################
################## Install ColabDesign
echo "Installing ColabDesign"
pip install git+https://github.com/sokrypton/ColabDesign.git --no-deps || { echo "Error: ColabDesign install failed"; exit 1; }
python -c "import colabdesign" >/dev/null 2>&1 || { echo "Error: colabdesign not found"; exit 1; }

############################################################################################################
################## Fix permissions for executables
chmod +x "${install_dir}/functions/dssp" || { echo "chmod dssp failed"; exit 1; }
chmod +x "${install_dir}/functions/DAlphaBall.gcc" || { echo "chmod DAlphaBall.gcc failed"; exit 1; }

############################################################################################################
################## Cleanup
conda deactivate
echo "Cleaning up ${pkg_manager} cache"
$pkg_manager clean -a -y

############################################################################################################
################## Finish
t=$SECONDS 
echo -e "\nSuccessfully finished BindCraft installation!"
echo "Activate environment: $pkg_manager activate BindCraft"
echo "Installation took $(($t / 3600))h $((($t / 60) % 60))m $(($t % 60))s"
