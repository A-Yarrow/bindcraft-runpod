#!/bin/bash

LOG_FILE="/workspace/startup.log"
exec > >(tee -a "$LOG_FILE") 2>&1
set -e

echo "=== [STARTUP] $(date) ==="
echo "Log output redirected to $LOG_FILE"

# Activate Conda environment
echo "[STEP] Activating Conda environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft

# Remove any conflicting files where directories are expected
echo "[STEP] Ensuring /workspace subpaths are directories..."
for d in settings_target settings_filters settings_advanced outputs/PDL1 inputs params; do
  path="/workspace/$d"
  if [ -f "$path" ]; then
    echo "[WARN] '$path' is a file — removing it to create directory"
    rm "$path"
  fi
done

# Create required workspace directories
echo "[STEP] Creating /workspace directories if they don't exist..."
mkdir -p /workspace/{settings_target,settings_filters,settings_advanced,outputs/PDL1,inputs,params}

# Copy default configs if missing
for d in settings_target settings_filters settings_advanced inputs; do
  if [ -z "$(ls -A /workspace/$d 2>/dev/null)" ]; then
    echo "[INFO] Copying default $d to /workspace/$d"
    cp -r /app/bindcraft/$d/* /workspace/$d/
  else
    echo "[INFO] $d already exists in /workspace — skipping copy"
  fi
done

# Install PyRosetta from Graylab (offline install)
ENV_PATH="/opt/conda/envs/BindCraft"
PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_NAME="pyrosetta.tar.bz2"
PACKAGE_DIR="/tmp/pyrosetta"

echo "[STEP] Downloading PyRosetta package..."
mkdir -p "$PACKAGE_DIR"
cd "$PACKAGE_DIR"
wget "$PACKAGE_URL" -O "$PACKAGE_NAME"

echo "[STEP] Installing PyRosetta into BindCraft environment (offline)..."
mamba install --offline "$PACKAGE_NAME" --use-local
if [ $? -ne 0 ]; then
  echo "[ERROR] Failed to install PyRosetta"
  exit 1
else
  echo "[INFO] PyRosetta installed successfully"
fi

# Cleanup PyRosetta tarball
rm -rf "$PACKAGE_DIR"

# Download AlphaFold2 weights if missing
WEIGHTS_DIR="/workspace/params"
WEIGHTS_FILE="${WEIGHTS_DIR}/params_model_5_ptm.npz"
if [ ! -f "$WEIGHTS_FILE" ]; then
  echo "[STEP] Downloading AlphaFold2 weights to $WEIGHTS_DIR..."
  cd "$WEIGHTS_DIR"
  wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
  echo "[INFO] Extracting weights..."
  tar -xf alphafold_params_2022-12-06.tar
  rm alphafold_params_2022-12-06.tar
else
  echo "[INFO] AlphaFold2 weights already present"
fi

# Clean up old Jupyter PID files
echo "[STEP] Cleaning up stale Jupyter runtime files..."
rm -f /root/.local/share/jupyter/runtime/*.pid || true

# Launch JupyterLab
echo "[STEP] Launching JupyterLab on port 8888..."
jupyter lab /app/bindcraft-runpod-start.ipynb \
  --ip=0.0.0.0 \
  --port=8888 \
  --allow-root \
  --NotebookApp.token='' \
  --NotebookApp.password=''

echo "=== [FINISHED] Startup complete: $(date) ==="

