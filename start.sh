#!/bin/bash

# Activate Conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate BindCraft

# Ensure persistent workspace directories exist
mkdir -p /workspace/{settings_target,settings_filters,settings_advanced,outputs}

# Copy defaults into /workspace only if missing
for d in settings_target settings_filters settings_advanced; do
  if [ -z "$(ls -A /workspace/$d 2>/dev/null)" ]; then
    echo "Copying default $d to /workspace/$d"
    cp -r /app/bindcraft/$d/* /workspace/$d/
  else
    echo "$d already exists in /workspace — skipping copy"
  fi
done

# Clean up old Jupyter runtime files (optional)
rm -f /root/.local/share/jupyter/runtime/*.pid || true

# Launch JupyterLab on port 8888 with no token/password
echo "Starting JupyterLab..."
jupyter lab /app/bindcraft-runpod-start.ipynb \
  --ip=0.0.0.0 \
  --port=8888 \
  --allow-root \
  --NotebookApp.token='' \
  --NotebookApp.password=''

