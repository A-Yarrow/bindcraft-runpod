import sys
import os
DEV_DIR = '/home/yarrow/projects/bindcraft-runpod/functions/'
PROD_DIR = '/app/bindcraft/functions'

path_used = DEV_DIR if os.path.exists(DEV_DIR) else PROD_DIR
print(f"Using function path: {path_used}")
sys.path.append(path_used)

from main_UI import launch_all_ui

launch_all_ui()