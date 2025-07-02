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
mkdir -p \
  /workspace/settings_target \
  /workspace/settings_filters \
  /workspace/settings_advanced \
  /workspace/outputs/PDL1 \
  /workspace/inputs \
  /workspace/params

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
PACKAGE_NAME=$(basename "$PACKAGE_URL")
PACKAGE_DIR="/tmp/pyrosetta"
STATUS_FILE="/workspace/status.log"

if [ -f "$PACKAGE_DIR/$PACKAGE_NAME" ]; then
  echo "[INFO] PyRosetta package already exists at $PACKAGE_DIR/$PACKAGE_NAME"
else

  echo "[STEP] Downloading PyRosetta package..."
  mkdir -p "$PACKAGE_DIR"
  cd "$PACKAGE_DIR"
  wget "$PACKAGE_URL" -O "$PACKAGE_NAME" || {
  echo "[FAIL] Failed to download PyRosetta" | tee -a "$STATUS_FILE"
  exit 1
}
fi
# Download and install PyRosetta (offline)
PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_NAME=$(basename "$PACKAGE_URL")
PACKAGE_DIR="/tmp/pyrosetta"
STATUS_FILE="/workspace/status.log"

if [ -f "$PACKAGE_DIR/$PACKAGE_NAME" ]; then
  echo "[INFO] PyRosetta package already exists at $PACKAGE_DIR/$PACKAGE_NAME"
else
  echo "[STEP] Downloading PyRosetta package..."
  mkdir -p "$PACKAGE_DIR"
  cd "$PACKAGE_DIR"
  wget "$PACKAGE_URL" -O "$PACKAGE_NAME" || {
    echo "[FAIL] Failed to download PyRosetta" | tee -a "$STATUS_FILE"
    exit 1
  }
fi

# Activate environment again in case subshell lost it
echo "[STEP] Activating BindCraft environment..."

CONDA_BASE=$(dirname $(dirname $(which conda)))
echo "[INFO] Using conda base at $CONDA_BASE"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate BindCraft || {
  echo "[FAIL] Failed to activate Conda environment" | tee -a "$STATUS_FILE"
  exit 1
}

# Install PyRosetta
echo "[STEP] Installing PyRosetta into BindCraft environment (offline)..."
mamba install -y "$PACKAGE_DIR/$PACKAGE_NAME" --offline || {
  echo "[FAIL] Failed to install PyRosetta" | tee -a "$STATUS_FILE"
  exit 1
  }

# Verify PyRosetta import
echo "[STEP] Verifying PyRosetta installation..."
if python -c "import pyrosetta; pyrosetta.init()" >/dev/null 2>&1; then
  echo "[SUCCESS] PyRosetta import and init successful!"
  rm -rf "$PACKAGE_DIR"
else
  echo "[FAIL] PyRosetta import failed" | tee -a "$STATUS_FILE"
  exit 1
fi

# Download AlphaFold2 weights if missing
WEIGHTS_DIR="/workspace/params"
WEIGHTS_FILE="${WEIGHTS_DIR}/params_model_5_ptm.npz"
if [ ! -f "${WEIGHTS_FILE}" ]; then
  echo "[STEP] Downloading AlphaFold2 weights to ${WEIGHTS_DIR}..."
  cd "${WEIGHTS_DIR}"
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

# Install pykernel
mamba install ipykernel -y || {
  echo "[FAIL] Failed to install ipykernel" | tee -a "$STATUS_FILE"
}
# Register the kernel
echo "[STEP] Registering Jupyter kernel..."
python -m ipykernel install --user --name=BindCraft --display-name="Python (BindCraft)" || {
  echo "[FAIL] Failed to register Jupyter kernel" | tee -a "$STATUS_FILE"
}

# Launch JupyterLab
echo "[STEP] Launching JupyterLab on port 8888..."
jupyter lab bindcraft-runpod-start.ipynb \
  --ip=0.0.0.0 \
  --port=8888 \
  --allow-root \
  --NotebookApp.token='' \
  --NotebookApp.password='' || {
    echo "[FAIL] Failed to launch JupyterLab" | tee -a "$STATUS_FILE"
}

echo "=== [FINISHED] Startup complete: $(date) ==="
