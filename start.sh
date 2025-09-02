#!/bin/bash
# Config variables
BINDCRAFT_DIR="/app/bindcraft"
WORKSPACE_DIR="/workspace"

NOTEBOOK_SRC="$BINDCRAFT_DIR/bindcraft-runpod-start.ipynb"
NOTEBOOK_DEST="$WORKSPACE_DIR/bindcraft-runpod-start.ipynb"

PYROSETTA_PACKAGE_URL="https://conda.graylab.jhu.edu/linux-64/pyrosetta-2025.03+release.1f5080a079-py310_0.tar.bz2"
PYROSETTA_PACKAGE_NAME=$(basename "$PYROSETTA_PACKAGE_URL")
PYROSETTA_PACKAGE_DIR="$BINDCRAFT_DIR/packages/pyrosetta"
PYROSETTA_PACKAGE_PATH="$PYROSETTA_PACKAGE_DIR/$PYROSETTA_PACKAGE_NAME"

ALPHAFOLD_WEIGHTS_ARCHIVE="alphafold_params_2022-12-06.tar"
ALPHAFOLD_WEIGHTS_FILE="$BINDCRAFT_DIR/params/params_model_5_ptm.npz"


#Jupyter Configeration
JUPYTER_IP="0.0.0.0"
JUPYTER_PORT=8888
JUPYTER_PASS_FILE="$WORKSPACE_DIR/jupyter_password.txt"
JUPYTER_CONFIG_FILE="$HOME/.jupyter/jupyter_server_config.py"
JUPYTER_LAB_CONFIG_FILE="$HOME/.jupyter/labconfig/coparty.json"

# Logging setup
LOG_FILE="$WORKSPACE_DIR/startup.log"
touch "$LOG_FILE"
chmod 644 "$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== [STARTUP] $(date) ==="
echo "Log output redirected to $LOG_FILE"

# Environment Setup
mkdir -p "$PYROSETTA_PACKAGE_DIR"

cd "$BINDCRAFT_DIR" || {
  echo "[FAIL] Cannot change to $BINDCRAFT_DIR" 
  exit 1
}

echo "[STEP] Activating Conda environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft || {
  echo "[FAIL] Failed to activate Conda environment" | tee -a "$LOG_FILE"
  exit 1
}

# Directory Setup (persistent workspace)

echo "[STEP] Creating persistent /workspace directories..."
mkdir -p \
  "$WORKSPACE_DIR/settings_target" \
  "$WORKSPACE_DIR/settings_filters" \
  "$WORKSPACE_DIR/settings_advanced" \
  "$WORKSPACE_DIR/outputs" \
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
    echo "[FAIL] Failed to copy starter notebook" | tee -a "$LOG_FILE"
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
      echo "[FAIL] Failed to copy $d configs" | tee -a "$LOG_FILE"
    }
  else
    echo "[INFO] $d already exists in $WORKSPACE_DIR — skipping copy"
  fi
done



# Checking GPU info
echo "[INFO] NVIDIA-SMI output:"
nvidia-smi

echo "[INFO] CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

# PyTorch check
python - <<END
import torch
print(f"[INFO] PyTorch CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"[INFO] GPU name: {torch.cuda.get_device_name(0)}")
    print(f"[INFO] Total VRAM (GB): {torch.cuda.get_device_properties(0).total_memory / 1024**3:.2f}")
END

# PyRosetta Install

if ! python -c "import pyrosetta" >/dev/null 2>&1; then
    echo "[STEP] PyRosetta not installed. Installing..."

  if [ ! -f "$PYROSETTA_PACKAGE_PATH" ]; then
    echo "[STEP] Downloading PyRosetta package..."
    cd "$PYROSETTA_PACKAGE_DIR"
    wget -nc -O "$PYROSETTA_PACKAGE_NAME" "$PYROSETTA_PACKAGE_URL" || {
      echo "[FAIL] Failed to download PyRosetta" | tee -a "$LOG_FILE"
      exit 1
    }
  else
    echo "[INFO] PyRosetta package already exists at $PYROSETTA_PACKAGE_DIR/$PYROSETTA_PACKAGE_NAME"
  fi

  #Install Package
  echo "[STEP] Installing PyRosetta (offline)..."
  mamba install -y "$PYROSETTA_PACKAGE_DIR/$PYROSETTA_PACKAGE_NAME" --offline || {
  echo "[FAIL] Failed to install PyRosetta" | tee -a "$LOG_FILE"
  exit 1
  }

  echo "[STEP] Verifying PyRosetta..."
  if python -c "import pyrosetta; pyrosetta.init()" >/dev/null 2>&1; then
    echo "[SUCCESS] PyRosetta import and init successful!"
  else
    echo "[FAIL] PyRosetta import failed" | tee -a "$LOG_FILE"
    exit 1
  fi
else 
  echo "[INFO] PyRosetta is allready installed."
fi

# AlphaFold Weights Download and Extraction
if [ ! -f "$ALPHAFOLD_WEIGHTS_FILE" ]; then
  echo "[STEP] Downloading AlphaFold2 weights..."
  cd "$BINDCRAFT_DIR/params"
  wget -nc https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar || {
    echo "[FAIL] Failed to download AlphaFold2 weights" | tee -a "$LOG_FILE"
    exit 1
  }
  echo "[INFO] Extracting weights..."
  tar --no-same-owner -xf "$ALPHAFOLD_WEIGHTS_ARCHIVE" && rm "$ALPHAFOLD_WEIGHTS_ARCHIVE" || {
    echo "[FAIL] Failed to extract weights" | tee -a "$LOG_FILE"
    exit 1
  }
  rm -f "$ALPHAFOLD_WEIGHTS_ARCHIVE"
else
  echo "[INFO] AlphaFold2 weights already present"
fi

# Kernel Registration

echo "[STEP] Installing ipykernel and registering BindCraft kernel..."
mamba install -y ipykernel || {
  echo "[FAIL] Failed to install ipykernel" | tee -a "$LOG_FILE"
}

python -m ipykernel install --name=BindCraft --display-name="Python (BindCraft)" --sys-prefix || {
  echo "[FAIL] Failed to register Jupyter kernel" | tee -a "$LOG_FILE"
}

# Cleanup old Jupyter runtime files
rm -f /root/.local/share/jupyter/runtime/*.pid || true

# JupyterLab Launch
cd "$WORKSPACE_DIR"
echo "[STEP] Launching JupyterLab on port $JUPYTER_PORT..."

if [ -f "$NOTEBOOK_DEST" ]; then
  echo "[INFO] Trusting notebook: $(basename "$NOTEBOOK_DEST")"
  jupyter trust "$NOTEBOOK_DEST"
  START_NOTEBOOK="$(basename "$NOTEBOOK_DEST")"
else
  echo "[WARN] Notebook not found. Launching without a specific notebook."
  START_NOTEBOOK=""
fi

# Write config to disable missing extensions (avoid warnings) + CopyParty setup
mkdir -p ~/.jupyter
cat > ~/.jupyter/jupyter_server_config.py <<'PYCONF'
c = get_config()

# Disable noisy extensions
c.ServerApp.jpserver_extensions = {
    "nbclassic": False,
    "jupyter_nbextensions_configurator": False,
    "jupyter_archive": False,
}

# Copyparty integration
c.ServerProxy.servers = {
    "copyparty": {
        "command": ["coparty", "-d", "/workspace", "--no-auth", "--http-only", "--port={port}"],
        "timeout": 60,
        "launcher_entry": {
            "title": "Copyparty",
            "icon_path": ""
        }
    }
}
PYCONF

PYCONF

# Generate a random password
JUPYTER_PASS=$(openssl rand -hex 16)

# Hash the password for Jupyter
HASHED_PASS=$(python -c "from jupyter_server.auth import passwd; print(passwd('$JUPYTER_PASS'))")

# Write config with hashed password
cat >> "$JUPYTER_CONFIG_FILE" <<EOF
c.ServerApp.identity_provider_class = "jupyter_server.auth.identity.PasswordIdentityProvider"
c.PasswordIdentityProvider.hashed_password = "${HASHED_PASS}"
EOF

#Save password file
echo "$JUPYTER_PASS" | tee "$JUPYTER_PASS_FILE"
chmod 600 "$JUPYTER_PASS_FILE"

echo "=====================================" 
echo "Your Jupyter password is: $JUPYTER_PASS" 
echo "your password is saved in: $JUPYTER_PASS_FILE" 
echo "=====================================" 
echo "=== [FINISHED] Startup complete: $(date) ===" 
echo "Access JupyterLab at http://$JUPYTER_IP:$JUPYTER_PORT" 

jupyter lab $START_NOTEBOOK \
  --ServerApp.root_dir="$WORKSPACE_DIR" \
  --ServerApp.allow_origin='*' \
  --ServerApp.allow_remote_access=True \
  --ip="$JUPYTER_IP" \
  --port="$JUPYTER_PORT" \
  --allow-root \
  --ServerApp.token="" \
  --no-browser || {
    echo "[FAIL] Failed to launch JupyterLab" | tee -a "$LOG_FILE"
}

