import os
import sys
from pathlib import Path
import ipywidgets as widgets
from IPython.display import display

# Configuration for paths (adjust as needed for cloud/local)
SETTINGS_DIRS = [
    '/workspace/settings_filters',
    '/workspace/settings_advanced'
]
BASE_PATH = '/workspace'  # or os.environ.get('BINDCRAFT_BASE', '/workspace')
BINDCRAFT_TEMPLATE_PATH = '/app/bindcraft/bindcraft_run_template.sh'
BINDCRAFT_OUTPUT_PATH = f'{BASE_PATH}/bindcraft_run.sh'
FUNCTIONS_PATH = '/app/bindcraft/functions'
DEFAULT_TARGET_JSON = f'{BASE_PATH}/settings_target/PDL1.json'

# Add path to import custom modules
sys.path.append(FUNCTIONS_PATH)
sys.path.append(str(Path(FUNCTIONS_PATH).parent))

# Import both UIs
from ui_target_editor import main_launch_target_editor
from ui_bindcraft_launch import main_launch_bindcraft_UI

# Shared widget that will be updated across both steps
selected_json_target_path_widget = widgets.Text(
    value=DEFAULT_TARGET_JSON,
    description='Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='80%')
)

def launch_all_ui():
    """
    Launch the target editor, then the BindCraft UI.
    The JSON path widget is updated after save and passed into BindCraft launcher.
    """
    print("Step 1: Edit your target JSON file.")
    main_launch_target_editor(
        json_target_path=selected_json_target_path_widget.value.strip(),
        path_output_widget=selected_json_target_path_widget  # <-- Live update on save
    )

    print("\nStep 2: Select settings and run BindCraft.")
    main_launch_bindcraft_UI(
        json_target_path=selected_json_target_path_widget.value.strip(),
        settings_dirs=SETTINGS_DIRS,
        base_path=BASE_PATH,
        bindcraft_template_run_file=BINDCRAFT_TEMPLATE_PATH,
        output_dir=BASE_PATH  # same as /workspace
    )
