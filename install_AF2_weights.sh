#!/bin/bash
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
