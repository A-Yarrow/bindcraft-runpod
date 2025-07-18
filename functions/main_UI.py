import os
import sys
import logging
from pathlib import Path
import ipywidgets as widgets
from IPython.display import display
from settings import ENV, SETTINGS  # 👈 import the centralized config

# Load relevant settings from the dictionary
BASE_PATH = SETTINGS[f'{ENV}_RUN_DIR']
FUNCTIONS_PATH = SETTINGS[f'{ENV}_FUNCTIONS_PATH']
SETTINGS_DIRS = SETTINGS[f'{ENV}_SETTINGS_DIRS']
BINDCRAFT_TEMPLATE_PATH = SETTINGS[f'{ENV}_BINDCRAFT_TEMPLATE_PATH']
BINDCRAFT_OUTPUT_PATH = SETTINGS[f'{ENV}_BINDCRAFT_OUTPUT_PATH']
DEFAULT_TARGET_JSON = SETTINGS[f'{ENV}_DEFAULT_TARGET_JSON']

# Logging configuration

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

file_handler = logging.FileHandler('Bindcraft_launch_UI.log')
file_handler.setFormatter(formatter)

logger.addHandler(stream_handler)
logger.addHandler(file_handler)

# Add function path for custom imports
sys.path.append(FUNCTIONS_PATH)
sys.path.append(str(Path(FUNCTIONS_PATH).parent))  # Optional: for modules in parent

# Import the UIs
from ui_target_editor import main_launch_target_editor
from ui_bindcraft_launch import main_launch_bindcraft_UI

# Shared widget that gets updated between steps
selected_json_target_path_widget = widgets.Text(
    value=DEFAULT_TARGET_JSON,
    description='Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='80%')
)


def upload_and_save_file(save_directory: str=BASE_PATH, description: str = "Upload File", filetypes: str = ''):
    """
    Creates a file upload widget and saves the uploaded file to the specified directory.

    Parameters:
        save_directory (str): Directory to save the uploaded file.
        description (str): Description shown next to the upload button.
        filetypes (str): Accepted file extensions (e.g. '.pdb,.json').

    Returns:
        widgets.FileUpload: The file upload widget instance.
    """
    uploader = widgets.FileUpload(
        accept=filetypes,
        multiple=False,
        description=description
    )
    display(uploader)

    def handle_upload(change):
        if not uploader.value:
            print("No file uploaded.")
            return
            
        for file_info in uploader.value.values():
            filename = file_info['metadata']['name']
            content = file_info['content']
            save_path = os.path.join(save_directory, filename)
            with open(save_path, 'wb') as f:
                f.write(content)   
            print(f"File saved to: {save_path}")

    uploader.observe(handle_upload, names='value')
    return uploader


def launch_all_ui():
    """
    Launch the target editor, then the BindCraft UI.
    The JSON path widget is updated after save and passed into BindCraft launcher.
    """
    print("Step 1: Edit your target JSON file.")
    main_launch_target_editor(
        json_target_path=selected_json_target_path_widget.value.strip(),
        path_output_widget=selected_json_target_path_widget
    )

    print("\nStep 2: Select settings and run BindCraft.")
    main_launch_bindcraft_UI(
        json_target_path=selected_json_target_path_widget.value.strip(),
        settings_dirs=SETTINGS_DIRS,
        base_path=BASE_PATH,
        bindcraft_template_run_file=BINDCRAFT_TEMPLATE_PATH,
        output_dir=BASE_PATH
    )
    print("Optional: Upload a JSON file to your settings_target directory.")
    upload_and_save_file(
    save_directory=f'{BASE_PATH}/settings_target',
    description="Upload JSON",
    filetypes='.json'
)
    


