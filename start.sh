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
mkdir -p /workspace/{settings_target,settings_filters,settings_advanced, outputs/PDL1, inputs, params}

# Copy default configs if missing
for d in settings_target settings_filters settings_advanced inputs; do
  if [ -z "$(ls -A /workspace/$d 2>/dev/null)" ]; then
    echo "[INFO] Copying default $d to /workspace/$d"
    cp -r /app/bindcraft/$d/* /workspace/$d/
  else
    echo "[INFO] $d already exists in /workspace — skipping copy"
  fi
done

# Install PyRosetta from Graylab (offline method)
ENV_PATH="/opt/conda/envs/BindCraft"
PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_NAME="pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_DIR="/tmp/pyrosetta"

mkdir -p "$PACKAGE_DIR"
cd "$PACKAGE_DIR"
wget "$PACKAGE_URL" -O "$PACKAGE_NAME"

echo "[STEP] Installing PyRosetta locally..."
mamba install --offline "$PACKAGE_NAME" --use-local
if [ $? -ne 0 ]; then
  echo "[ERROR] Failed to install PyRosetta from Graylab"
  exit 1
else
  echo "[INFO] Successfully installed PyRosetta"
fi
rm -rf "$PACKAGE_DIR"  # Clean up temporary directory

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
mkdir -p /workspace/{settings_target,settings_filters,settings_advanced, outputs/PDL1, inputs, params}

# Copy default configs if missing
for d in settings_target settings_filters settings_advanced inputs; do
  if [ -z "$(ls -A /workspace/$d 2>/dev/null)" ]; then
    echo "[INFO] Copying default $d to /workspace/$d"
    cp -r /app/bindcraft/$d/* /workspace/$d/
  else
    echo "[INFO] $d already exists in /workspace — skipping copy"
  fi
done

# Install PyRosetta from graylab channel
ENV_PATH="/opt/conda/envs/BindCraft"
PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_NAME="pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_DIR="/tmp/pyrosetta"

# Create temp dir and download the package
mkdir -p "$PACKAGE_DIR"
cd "$PACKAGE_DIR"
wget "$PACKAGE_URL" -O "$PACKAGE_NAME"

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft

# Use mamba to install from the local file
mamba install --offline "$PACKAGE_NAME" --use-local
if [ $? -ne 0 ]; then
  echo "[ERROR] Failed to install PyRosetta from graylab channel"
  exit 1
else
  echo "[INFO] Successfully installed PyRosetta"
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
