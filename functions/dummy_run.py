from datetime import datetime
import sys
import time
import itertools
import signal

print("=== Dummy BindCraft Run Started ===", flush=True)
steps = ["Loading data", "Running computations", "Saving results"]

def handle_sigterm(signum, frame):
    print("=== Dummy BindCraft Run Received Termination Signal ===", flush=True)
    sys.exit(0)

signal.signal(signal.SIGTERM, handle_sigterm)
signal.signal(signal.SIGINT, handle_sigterm)

# Loop until killed

starting_time = datetime.now()
for i in itertools.count(1):
    current_time = datetime.now()
    print(f"[Step {i}] {steps[i % len(steps)]}...", flush=True)
    dt = current_time - starting_time
    time.sleep(2)
    print(dt, flush=True)
