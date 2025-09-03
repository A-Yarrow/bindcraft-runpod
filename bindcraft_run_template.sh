#!/bin/bash
# BindCraft Run Script Template

# These variables are substituted by Python when generating the script
# ENV, LOG_DIR, TARGET_NAME, TARGET_FILE_PATH, TARGET_FILE_NAME,
# FILTERS_FILE_PATH, ADVANCED_FILE_PATH

log_file="${LOG_DIR}/${JOB_NAME}-bindcraft_log.txt"
pid_file="${PID_DIR}/${JOB_NAME}-bindcraft_pid.txt"
dummy_run_python_file="/home/yarrow/projects/bindcraft-runpod/functions/dummy_run.py"
#export PID_FILE="$pid_file"

mkdir -p "$LOG_DIR"
mkdir -p "$PID_DIR"

echo "[INFO] Checking target settings file..."
if [ ! -f "$TARGET_FILE_PATH" ]; then
  echo "[ERROR] Target file '$TARGET_FILE_NAME' not found"
  exit 1
fi

echo "[INFO] Starting BindCraft with target '$TARGET_NAME'"
echo "[INFO] Logging output to $log_file" 

# Function to cleanup PID file if script exits
cleanup() {
    if [ -f "$pid_file" ]; then
        rm -f "$pid_file"
        echo "[INFO] PID file $pid_file removed"
    fi
}
trap cleanup EXIT INT TERM

# Determine Python command based on ENV
if [ "$ENV" = "DEV" ]; then
    python_cmd="python -u $dummy_run_python_file"
elif [ "$ENV" = "PROD" ]; then
    python_cmd="python -u /app/bindcraft/bindcraft.py \
    --settings "$TARGET_FILE_PATH" \
    --filters "$FILTERS_FILE_PATH" \
    --advanced "$ADVANCED_FILE_PATH" \
    "
else
    echo "[ERROR] Unknown ENV: $ENV"
    exit 1
fi

start_time=$(date +%s) #Start timestamp

# Start the Python job in the background
echo "[INFO] Starting BindCraft job in background..."
nohup $python_cmd >> "$log_file" 2>&1 &
bindcraft_pid=$!
echo "$bindcraft_pid" > "$pid_file"
echo "[INFO] BindCraft PID: $bindcraft_pid"

# Wait for the background job to finish
wait $bindcraft_pid

# Calculate elapsed time
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))
# Log completion
echo "[INFO] === BindCraft Run Finished ===" >> "$log_file"
echo "[INFO] Total runtime: ${hours}h ${minutes}m ${seconds}s" >> "$log_file"