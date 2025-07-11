#!/bin/bash
# BindCraft Run Script Template
log_file="$LOG_DIR/$TARGET_NAME-bindcraft_log.txt"

mkdir -p "$LOG_DIR"

echo "[INFO] Checking target settings file..."
if [ ! -f "$TARGET_FILE_PATH" ]; then
  echo "[ERROR] Target file '$TARGET_FILE_NAME' not found"
  exit 1
fi

echo "[INFO] Starting BindCraft with target '$TARGET_NAME'"
echo "[INFO] Logging output to '$LOG_DIR/$TARGET_NAME-bindcraft_log.txt'"

nohup python /app/bindcraft/bindcraft.py \
  --settings "$TARGET_FILE_PATH" \
  --filters "$FILTERS_FILE_PATH" \
  --advanced "$ADVANCED_FILE_PATH" \
  > "$LOG_DIR/$TARGET_NAME-bindcraft_log.txt" 2>&1 &

echo "[INFO] BindCraft is now running in the background"