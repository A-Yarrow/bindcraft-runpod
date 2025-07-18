#!/bin/bash

# Logging setup
LOG_FILE="/workspace/startup.log"
chmod 644 "$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== [STARTUP] $(date) ==="
echo "Log output redirected to $LOG_FILE"

# Config variables
BINDCRAFT_DIR="/app/bindcraft"
WORKSPACE_DIR="/workspace"

NOTEBOOK_SRC="$BINDCRAFT_DIR/bindcraft-runpod-start.ipynb"
NOTEBOOK_DEST="$WORKSPACE_DIR/bindcraft-runpod-start.ipynb"

PYROSETTA_PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PYROSETTA_PACKAGE_NAME=$(basename "$PYROSETTA_PACKAGE_URL")
PYROSETTA_PACKAGE_DIR="/tmp/pyrosetta"

ALPHAFOLD_WEIGHTS_ARCHIVE="alphafold_params_2022-12-06.tar"
ALPHAFOLD_WEIGHTS_FILE="$BINDCRAFT_DIR/params/params_model_5_ptm.npz"

JUPYTER_IP="0.0.0.0"
JUPYTER_PORT=8888

# Environment Setup

cd "$BINDCRAFT_DIR" || {
  echo "[FAIL] Cannot change to $BINDCRAFT_DIR" 
  exit 1
}

echo "[STEP] Activating Conda environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft || {
  echo "[FAIL] Failed to activate Conda environment" | tee -a "$STATUS_FILE"
  exit 1
}

# Directory Setup (persistent workspace)

echo "[STEP] Creating persistent /workspace directories..."
mkdir -p \
  "$WORKSPACE_DIR/settings_target" \
  "$WORKSPACE_DIR/settings_filters" \
  "$WORKSPACE_DIR/settings_advanced" \
  "$WORKSPACE_DIR/outputs/PDL1" \
  "$WORKSPACE_DIR/inputs"

if [ ! -d "$BINDCRAFT_DIR/params" ]; then
  echo "[INFO] Creating params directory in $BINDCRAFT_DIR"
  mkdir -p "$BINDCRAFT_DIR/params"
else
  echo "[INFO] Params directory already exists in $BINDCRAFT_DIR"
fi

# Copy Starter Notebook if missing

if [ ! -f "$NOTEBOOK_DEST" ]; then
  echo "[INFO] Copying starter notebook to $WORKSPACE_DIR..."
  cp "$NOTEBOOK_SRC" "$NOTEBOOK_DEST" && chmod 666 "$NOTEBOOK_DEST" || {
    echo "[FAIL] Failed to copy starter notebook" | tee -a "$STATUS_FILE"
    exit 1
  }
else
  echo "[INFO] Notebook already exists in $WORKSPACE_DIR — skipping copy"
fi


# Copy Default Configs if missing

for d in settings_target settings_filters settings_advanced inputs; do
  if [ -z "$(ls -A "$WORKSPACE_DIR/$d" 2>/dev/null)" ]; then
    echo "[INFO] Copying default $d to $WORKSPACE_DIR/$d"
    cp -r "$BINDCRAFT_DIR/$d/"* "$WORKSPACE_DIR/$d/" || {
      echo "[FAIL] Failed to copy $d configs" | tee -a "$STATUS_FILE"
    }
  else
    echo "[INFO] $d already exists in $WORKSPACE_DIR — skipping copy"
  fi
done


# PyRosetta Offline Install

if [ ! -f "$PYROSETTA_PACKAGE_DIR/$PYROSETTA_PACKAGE_NAME" ]; then
  echo "[STEP] Downloading PyRosetta package..."
  mkdir -p "$PYROSETTA_PACKAGE_DIR"
  cd "$PYROSETTA_PACKAGE_DIR"
  wget "$PYROSETTA_PACKAGE_URL" -O "$PYROSETTA_PACKAGE_NAME" || {
    echo "[FAIL] Failed to download PyRosetta" | tee -a "$STATUS_FILE"
    exit 1
  }
else
  echo "[INFO] PyRosetta package already exists at $PYROSETTA_PACKAGE_DIR/$PYROSETTA_PACKAGE_NAME"
fi

echo "[STEP] Installing PyRosetta (offline)..."
mamba install -y "$PYROSETTA_PACKAGE_DIR/$PYROSETTA_PACKAGE_NAME" --offline || {
  echo "[FAIL] Failed to install PyRosetta" | tee -a "$STATUS_FILE"
  exit 1
}

echo "[STEP] Verifying PyRosetta..."
if python -c "import pyrosetta; pyrosetta.init()" >/dev/null 2>&1; then
  echo "[SUCCESS] PyRosetta import and init successful!"
  rm -rf "$PYROSETTA_PACKAGE_DIR"
else
  echo "[FAIL] PyRosetta import failed" | tee -a "$STATUS_FILE"
  exit 1
fi

# AlphaFold Weights Download and Extraction
if [ ! -f "$ALPHAFOLD_WEIGHTS_FILE" ]; then
  echo "[STEP] Downloading AlphaFold2 weights..."
  cd "$BINDCRAFT_DIR/params"
  wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar || {
    echo "[FAIL] Failed to download AlphaFold2 weights" | tee -a "$STATUS_FILE"
  }
  echo "[INFO] Extracting weights..."
  tar -xf "$ALPHAFOLD_WEIGHTS_ARCHIVE" && rm "$ALPHAFOLD_WEIGHTS_ARCHIVE" || {
    echo "[FAIL] Failed to extract weights" | tee -a "$STATUS_FILE"
  }
else
  echo "[INFO] AlphaFold2 weights already present"
fi


# Kernel Registration

echo "[STEP] Installing ipykernel and registering BindCraft kernel..."
mamba install -y ipykernel || {
  echo "[FAIL] Failed to install ipykernel" | tee -a "$STATUS_FILE"
}

python -m ipykernel install --name=BindCraft --display-name="Python (BindCraft)" --sys-prefix || {
  echo "[FAIL] Failed to register Jupyter kernel" | tee -a "$STATUS_FILE"
}

# JupyterLab Launch

echo "[STEP] Cleaning up old Jupyter runtime files..."
rm -f /root/.local/share/jupyter/runtime/*.pid || true

echo "[STEP] Disabling unused Jupyter extensions..."
jupyter server extension disable jupyter_archive
jupyter server extension disable nbclassic

cd "$WORKSPACE_DIR"
echo "[STEP] Launching JupyterLab on port $JUPYTER_PORT..."

if [ -f "$NOTEBOOK_DEST" ]; then
  echo "[INFO] Trusting notebook: $(basename "$NOTEBOOK_DEST")"
  jupyter trust "$(basename "$NOTEBOOK_DEST")"
  START_NOTEBOOK="$(basename "$NOTEBOOK_DEST")"
else
  echo "[WARN] Notebook not found. Launching without a specific notebook."
  START_NOTEBOOK=""
fi

jupyter lab $START_NOTEBOOK \
  --notebook-dir="$WORKSPACE_DIR" \
  --NotebookApp.allow_origin='*' \
  --ip="$JUPYTER_IP" \
  --port="$JUPYTER_PORT" \
  --allow-root \
  --NotebookApp.token='' \
  --NotebookApp.password='' \
  --no-browser || {
    echo "[FAIL] Failed to launch JupyterLab" | tee -a "$STATUS_FILE"
}

echo "=== [FINISHED] Startup complete: $(date) ==="