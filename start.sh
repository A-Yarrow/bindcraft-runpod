#!/bin/bash

LOG_FILE="/workspace/startup.log"
STATUS_FILE="/workspace/status.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== [STARTUP] $(date) ==="
echo "Log output redirected to $LOG_FILE"

# --- Environment Setup ---
cd /app/bindcraft || {
  echo "[FAIL] Cannot change to /app/bindcraft" | tee -a "$STATUS_FILE"
  exit 1
}

echo "[STEP] Activating Conda environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft || {
  echo "[FAIL] Failed to activate Conda environment" | tee -a "$STATUS_FILE"
  exit 1
}

# --- Directory Setup ---
echo "[STEP] Creating persistent /workspace directories..."
mkdir -p \
  /workspace/settings_target \
  /workspace/settings_filters \
  /workspace/settings_advanced \
  /workspace/outputs/PDL1 \
  /workspace/inputs \
  /workspace/params

# --- Copy Starter Notebook ---
NOTEBOOK_SRC="/app/bindcraft/bindcraft-runpod-start.ipynb"
NOTEBOOK_DEST="/workspace/bindcraft-runpod-start.ipynb"

if [ ! -f "$NOTEBOOK_DEST" ]; then
  echo "[INFO] Copying starter notebook to /workspace..."
  cp "$NOTEBOOK_SRC" "$NOTEBOOK_DEST" && chmod 666 "$NOTEBOOK_DEST" || {
    echo "[FAIL] Failed to copy starter notebook" | tee -a "$STATUS_FILE"
    exit 1
  }
else
  echo "[INFO] Notebook already exists in /workspace — skipping copy"
fi

# --- Copy Default Configs ---
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

# --- PyRosetta Offline Install ---
PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PACKAGE_NAME=$(basename "$PACKAGE_URL")
PACKAGE_DIR="/tmp/pyrosetta"

if [ ! -f "$PACKAGE_DIR/$PACKAGE_NAME" ]; then
  echo "[STEP] Downloading PyRosetta package..."
  mkdir -p "$PACKAGE_DIR"
  cd "$PACKAGE_DIR"
  wget "$PACKAGE_URL" -O "$PACKAGE_NAME" || {
    echo "[FAIL] Failed to download PyRosetta" | tee -a "$STATUS_FILE"
    exit 1
  }
else
  echo "[INFO] PyRosetta package already exists at $PACKAGE_DIR/$PACKAGE_NAME"
fi

echo "[STEP] Installing PyRosetta (offline)..."
mamba install -y "$PACKAGE_DIR/$PACKAGE_NAME" --offline || {
  echo "[FAIL] Failed to install PyRosetta" | tee -a "$STATUS_FILE"
  exit 1
}

echo "[STEP] Verifying PyRosetta..."
if python -c "import pyrosetta; pyrosetta.init()" >/dev/null 2>&1; then
  echo "[SUCCESS] PyRosetta import and init successful!"
  rm -rf "$PACKAGE_DIR"
else
  echo "[FAIL] PyRosetta import failed" | tee -a "$STATUS_FILE"
  exit 1
fi

# --- AlphaFold Weights ---
WEIGHTS_FILE="/workspace/params/params_model_5_ptm.npz"
if [ ! -f "$WEIGHTS_FILE" ]; then
  echo "[STEP] Downloading AlphaFold2 weights..."
  cd /workspace/params
  wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar || {
    echo "[FAIL] Failed to download AlphaFold2 weights" | tee -a "$STATUS_FILE"
  }
  echo "[INFO] Extracting weights..."
  tar -xf alphafold_params_2022-12-06.tar && rm alphafold_params_2022-12-06.tar || {
    echo "[FAIL] Failed to extract weights" | tee -a "$STATUS_FILE"
  }
else
  echo "[INFO] AlphaFold2 weights already present"
fi

# --- Kernel Registration ---
echo "[STEP] Installing ipykernel and registering BindCraft kernel..."
mamba install -y ipykernel || {
  echo "[FAIL] Failed to install ipykernel" | tee -a "$STATUS_FILE"
}

python -m ipykernel install --name=BindCraft --display-name="Python (BindCraft)" --sys-prefix || {
  echo "[FAIL] Failed to register Jupyter kernel" | tee -a "$STATUS_FILE"
}

# --- JupyterLab Launch ---
echo "[STEP] Cleaning up old Jupyter runtime files..."
rm -f /root/.local/share/jupyter/runtime/*.pid || true

echo "[STEP] Disabling unused Jupyter extensions..."
jupyter server extension disable jupyter_archive
jupyter server extension disable nbclassic

cd /workspace
echo "[STEP] Launching JupyterLab on port 8888..."

if [ -f "bindcraft-runpod-start.ipynb" ]; then
  jupyter lab bindcraft-runpod-start.ipynb \
    --NotebookApp.allow_origin='*' \
    --ip=0.0.0.0 \
    --port=8888 \
    --allow-root \
    --NotebookApp.token='' \
    --NotebookApp.password='' || {
      echo "[FAIL] Failed to launch JupyterLab" | tee -a "$STATUS_FILE"
    }
else
  jupyter lab \
    --NotebookApp.allow_origin='*' \
    --ip=0.0.0.0 \
    --port=8888 \
    --allow-root \
    --NotebookApp.token='' \
    --NotebookApp.password='' || {
      echo "[FAIL] Failed to launch JupyterLab" | tee -a "$STATUS_FILE"
    }
fi

echo "=== [FINISHED] Startup complete: $(date) ==="

