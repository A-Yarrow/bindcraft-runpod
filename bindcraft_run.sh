#!/bin/bash
# BindCraft Run Script Template
log_file="/home/yarrow/projects/bindcraft-runpod/outputs/PDL1/PDL1-bindcraft_log.txt"

mkdir -p "/home/yarrow/projects/bindcraft-runpod/outputs/PDL1"

echo "[INFO] Checking target settings file..."
if [ ! -f "/home/yarrow/projects/bindcraft-runpod/settings_target/PDL1.json" ]; then
  echo "[ERROR] Target file 'PDL1.json' not found"
  exit 1
fi

echo "[INFO] Starting BindCraft with target 'PDL1'"
echo "[INFO] Logging output to '/home/yarrow/projects/bindcraft-runpod/outputs/PDL1/PDL1-bindcraft_log.txt'"

nohup python /app/bindcraft/bindcraft.py \
  --settings "/home/yarrow/projects/bindcraft-runpod/settings_target/PDL1.json" \
  --filters "/home/yarrow/projects/bindcraft-runpod/settings_filters/no_filters.json" \
  --advanced "/home/yarrow/projects/bindcraft-runpod/settings_advanced/betasheet_4stage_multimer.json" \
  > "/home/yarrow/projects/bindcraft-runpod/outputs/PDL1/PDL1-bindcraft_log.txt" 2>&1 &

echo "[INFO] BindCraft is now running in the background"