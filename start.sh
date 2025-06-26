#!/bin/bash

LOG_FILE="/workspace/startup.log"
STATUS_FILE="/workspace/status.log"
exec > >(tee -a "$LOG_FILE") 2>&1
# set -e  # Temporarily disabled during testing

echo "=== [STARTUP] $(date) ==="
echo "Log output redirected to $LOG_FILE"

# Activate Conda environment
echo "[STEP] Activating Conda environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft || {
  echo "[FAIL] Failed to activate Conda environment" | tee -a "$STATUS_FILE"
}

# Create required workspace directories
echo "[STEP] Creating persistent /workspace directories if they don't exist..."
mkdir -p /workspace/{settings_target,settings_filters,settings_advanced, outputs/PDL1, inputs, params}

# Copy default configs if missing
for d in settings_target settings_filters settings_advanced inputs; do
  if [ -z "$(ls -A /workspace/$d 2>/dev/null)" ]; then
    echo "[INFO] Copying default $d to /workspace/$d"
    cp -r /app/bindcraft/$d/* /workspace/$d/ || {
      echo "[FAIL] Failed to copy $d configs" | tee -a "$STATUS_FILE"
    }
  else
    echo "[INFO] $d already exists in /workspace — skipping copy"
  fi
done

# Download and install PyRosetta (offline)
PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_NAME="pyrosetta.tar.bz2"
PACKAGE_DIR="/tmp/pyrosetta"
ENV_PATH="/opt/conda/envs/BindCraft"

echo "[STEP] Downloading PyRosetta package..."
mkdir -p "$PACKAGE_DIR"
cd "$PACKAGE_DIR"
wget "$PACKAGE_URL" -O "$PACKAGE_NAME" || {
  echo "[FAIL] Failed to download PyRosetta" | tee -a "$STATUS_FILE"
}

echo "[STEP] Installing PyRosetta into BindCraft environment (offline)..."
mamba install "$PACKAGE_NAME" --offline || {
  echo "[FAIL] Failed to install PyRosetta" | tee -a "$STATUS_FILE"
}
rm -rf "$PACKAGE_DIR"

# Download AlphaFold2 weights if missing
WEIGHTS_DIR="/workspace/params"
WEIGHTS_FILE="${WEIGHTS_DIR}/params_model_5_ptm.npz"
if [ ! -f "$WEIGHTS_FILE" ]; then
  echo "[STEP] Downloading AlphaFold2 weights to $WEIGHTS_DIR..."
  cd "$WEIGHTS_DIR"
  wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar || {
    echo "[FAIL] Failed to download AlphaFold2 weights" | tee -a "$STATUS_FILE"
  }
  echo "[INFO] Extracting weights..."
  tar -xf alphafold_params_2022-12-06.tar || {
    echo "[FAIL] Failed to extract AlphaFold2 weights" | tee -a "$STATUS_FILE"
  }
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
  --NotebookApp.password='' || {
    echo "[FAIL] Failed to launch JupyterLab" | tee -a "$STATUS_FILE"
}

echo "=== [FINISHED] Startup complete: $(date) ==="


