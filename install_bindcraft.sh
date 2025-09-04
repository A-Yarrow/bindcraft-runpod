#!/bin/bash
set -e #exit on any error

# Build with: docker build --progress=plain -t bindcraft:test .
################## BindCraft installation script optimized for RunPod 12.4 base image: FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04
################## tested on 2025-09-04
################## specify conda/mamba folder, and installation folder for git repositories, and whether to use mamba or $pkg_manager

# Define the short and long options
OPTIONS=p:c:
LONGOPTIONS=pkg_manager:,cuda:

# wrap getopt in || true so set -e won't exit on non-zero
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@") || true
eval set -- "$PARSED"

# Default value for pkg_manager
pkg_manager='conda'
cuda=''

# Process the command-line options
while true; do
  case "$1" in
    -p|--pkg_manager)
      pkg_manager="$2"
      shift 2
      ;;
    -c|--cuda)
      cuda="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      echo -e "Invalid option $1" >&2
      exit 1
      ;;
  esac
done

# Example usage of the parsed variables
echo -e "Package manager: $pkg_manager"
echo -e "CUDA: $cuda"

############################################################################################################
############################################################################################################
################## initialisation
SECONDS=0

# set paths needed for installation and check for conda installation
install_dir=$(pwd)
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo -e "Error: conda is not installed or cannot be initialised."; exit 1; }
echo -e "Conda is installed at: $CONDA_BASE"

### BindCraft install begin, create base environment
echo -e "Installing BindCraft environment\n"
$pkg_manager create --name BindCraft python=3.10 -y || { echo -e "Error: Failed to create BindCraft conda environment"; exit 1; }
conda env list | grep -w 'BindCraft' >/dev/null 2>&1 || { echo -e "Error: Conda environment 'BindCraft' does not exist after creation."; exit 1; }

# Load newly created BindCraft environment
echo -e "Loading BindCraft environment\n"
source ${CONDA_BASE}/bin/activate ${CONDA_BASE}/envs/BindCraft || { echo -e "Error: Failed to activate the BindCraft environment."; exit 1; }
[ "$CONDA_DEFAULT_ENV" = "BindCraft" ] || { echo -e "Error: The BindCraft environment is not active."; exit 1; }
echo -e "BindCraft environment activated at ${CONDA_BASE}/envs/BindCraft"

# install required conda packages
echo -e "Installing conda requirements\n"

$pkg_manager install \
  pip pandas matplotlib 'numpy<2.0.0' biopython scipy pdbfixer seaborn libgfortran5 tqdm \
  jupyter jupyterlab ffmpeg pyrosetta fsspec py3dmol \
  chex dm-haiku 'flax<0.10.0' dm-tree joblib ml-collections immutabledict optax \
  psutil copyparty \
  -c conda-forge --channel https://conda.graylab.jhu.edu -y \
|| { echo -e "Error: Failed to install conda packages."; exit 1; }


# make sure all required packages were installed
required_packages=(pip pandas libgfortran5 matplotlib numpy biopython scipy pdbfixer seaborn tqdm jupyter jupyterlab ffmpeg \
fsspec py3dmol chex dm-haiku dm-tree joblib ml-collections immutabledict optax psutil copyparty \
)
missing_packages=()

# Check each package
for pkg in "${required_packages[@]}"; do
    conda list "$pkg" | grep -w "$pkg" >/dev/null 2>&1 || missing_packages+=("$pkg")
done

# If any packages are missing, output error and exit
if [ ${#missing_packages[@]} -ne 0 ]; then
    echo -e "Error: The following packages are missing from the environment:"
    for pkg in "${missing_packages[@]}"; do
        echo -e " - $pkg"
    done
    exit 1
fi

# install pip packages
python -m pip install jupyter-server-proxy

# JAX GPU install (pip)
JAX GPU install
python -m pip install --upgrade pip wheel
python -m pip install --no-cache-dir \
  jax==0.4.28 \
  jaxlib==0.4.28+cuda12.cudnn89 \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html


echo -e "GPU-enabled JAX verification complete\n"

# install ColabDesign
echo -e "Installing ColabDesign\n"
pip install git+https://github.com/sokrypton/ColabDesign.git --no-deps || { echo -e "Error: Failed to install ColabDesign"; exit 1; }
python -c "import colabdesign" >/dev/null 2>&1 || { echo -e "Error: colabdesign module not found after installation"; exit 1; }

# chmod executables
echo -e "Changing permissions for executables\n"
chmod +x "${install_dir}/functions/dssp" || { echo -e "Error: Failed to chmod dssp"; exit 1; }
chmod +x "${install_dir}/functions/DAlphaBall.gcc" || { echo -e "Error: Failed to chmod DAlphaBall.gcc"; exit 1; }

# finish
conda deactivate
echo -e "BindCraft environment set up\n"

############################################################################################################
############################################################################################################
################## cleanup
echo -e "Cleaning up ${pkg_manager} temporary files to save space\n"
$pkg_manager clean -a -y
echo -e "$pkg_manager cleaned up\n"

################## finish script
t=$SECONDS 
echo -e "Successfully finished BindCraft installation!\n"
echo -e "Activate environment using command: \"$pkg_manager activate BindCraft\""
echo -e "\n"
echo -e "Installation took $(($t / 3600)) hours, $((($t / 60) % 60)) minutes and $(($t % 60)) seconds."
