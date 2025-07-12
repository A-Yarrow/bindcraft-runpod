import os
import sys
from pathlib import Path
import ipywidgets as widgets
from IPython.display import display

SETTINGS_DIRS =
['/workspace/settings_filters',
'/workspace/settings_advanced'
]
BASE_PATH ='/home/yarrow/projects/bindcraft-runpod/workspace'
BINDCRAFT_TEMPLATE_PATH = '/home/yarrow/projects/bindcraft-runpod/workspace/bindcraft_run_template.sh'
BINDCRAFT_OUTPUT_PATH = '/home/yarrow/projects/bindcraft-runpod/workspace/bindcraft_run.sh'
FUNCTIONS_PATH = '/home/yarrow/projects/bindcraft-runpod/functions'
DEFAULT_TARGET_JSON = '/home/yarrow/projects/bindcraft-runpod/workspace/settings_target/PDL1.json'
# Add path to import custom modules
sys.path.append(FUNCTIONS_PATH)

# Import both UIs
from ui_target_editor import main_launch_target_editor
from ui_bindcraft_launch import main_launch_bindcraft_UI

# Shared variable to hold the path to the edited target JSON
selected_json_target_path_widget = widgets.Text(
    value=DEFAULT_TARGET_JSON,
    description='Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='80%')
)
json_target_path=selected_json_target_path_widget.value.strip()

def launch_all_ui():
    """
    Launch the target editor, then settings + bindcraft launch UI.
    The path from the editor is reused as the target input for the BindCraft launcher.
    """
    print("Step 1: Edit your target JSON file.")
    main_launch_target_editor(json_target_path)
    print("\nStep 2: Select settings and run BindCraft.")
    
    # This launches the second UI and passes in the path to the edited JSON
    main_launch_bindcraft_UI(
        json_target_path=json_target_path,
        settings_dirs=SETTINGS_DIRS,
        base_path=BASE_PATH,
        bindcraft_template_path=BINDCRAFT_TEMPLATE_PATH,
        bindcraft_output_path=BINDCRAFT_OUTPUT_PATH,
    )

