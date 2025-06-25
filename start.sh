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

# Create required workspace directories
echo "[STEP] Creating persistent /workspace directories if they don't exist..."
mkdir -p /workspace/{settings_target,settings_filters,settings_advanced,outputs,params}

# Copy default configs if missing
for d in settings_target settings_filters settings_advanced; do
  if [ -z "$(ls -A /workspace/$d 2>/dev/null)" ]; then
    echo "[INFO] Copying default $d to /workspace/$d"
    cp -r /app/bindcraft/$d/* /workspace/$d/
  else
    echo "[INFO] $d already exists in /workspace — skipping copy"
  fi
done

# Install PyRosetta if missing
echo "[STEP] Checking PyRosetta installation..."
if ! python -c "import pyrosetta" &>/dev/null; then
  echo "[INFO] Installing PyRosetta into BindCraft environment..."
  mamba install -y pyrosetta -c https://conda.graylab.jhu.edu
else
  echo "[INFO] PyRosetta already installed"
fi

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
  echo "[INFO] AlphaFold2 weights already present at $WEIGHTS_FILE"
fi

# Clean up old Jupyter PID files
echo "[STEP] Cleaning up old Jupyter runtime files..."
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
