import os
import sys
from pathlib import Path
import ipywidgets as widgets
from IPython.display import display

# Add path to import custom modules
sys.path.append('/workspace/../functions')

# Import both UIs
from ui_target_editor import main_launch_editor
from ui_bindcraft_runner import launch_bindcraft_UI

# Shared variable to hold the path to the edited target JSON
selected_target_json_path = widgets.Text(
    value='/workspace/settings_target/PDL1.json',
    description='Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='80%')
)

def launch_all_ui():
    """
    Launch the target editor, then settings + runner UI.
    The path from the editor is reused as the target input for the BindCraft runner.
    """
    print("🔧 Step 1: Edit your target JSON file.")
    main_launch_editor(json_target_path=selected_target_json_path.value)

    print("\n⚙️ Step 2: Select settings and run BindCraft.")
    settings_dirs = [
        '/workspace/settings_filters',
        '/workspace/settings_advanced'
    ]

    # This launches the second UI and passes in the path to the edited JSON
    launch_bindcraft_UI(
        settings_dirs=settings_dirs,
        SETTINGS_FILE_PATH='/workspace',
        BINDCRAFT_TEMPLATE_PATH='/workspace/bindcraft_run_template.sh',
        BINDCRAFT_OUTPUT_PATH='/workspace/bindcraft_run.sh',
        target_path_widget=selected_target_json_path  # shared input
    )

