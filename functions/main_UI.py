import os
import sys
import logging
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException
from pathlib import Path
import ipywidgets as widgets
from IPython.display import display
from settings import ENV, SETTINGS  # 👈 import the centralized config
from typing import Callable, Optional

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
target_editor_container = widgets.VBox()
selected_json_target_path_widget = widgets.Text(
    value=DEFAULT_TARGET_JSON,
    description='Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='80%')
)

PDB_UI_OUTPUT_WIDGET = widgets.Output()
JSON_UI_OUTPUT_WIDGET = widgets.Output()

import ipywidgets as widgets
from IPython.display import display

def step_box(text, font_size='16px'):
    return widgets.HTML(
        value=f"<div style='padding: 8px; background-color: #eef; font-size: {font_size}; border: 1px solid #ccd; border-radius: 6px;'><strong>{text}</strong></div>"
    )


uploaded_pdb_path = None
def validate_pdb_file(pdb_path: str) -> bool:
    """
    Validate if the provided PDB file is valid.
    Returns True if valid, False otherwise.
    """
    parser = PDBParser(QUIET=False)
    try:
        structure = parser.get_structure('temp', pdb_path)
        with PDB_UI_OUTPUT_WIDGET:
            print(f"PDB validated and saved to: {pdb_path}")
        logger.info(f"PDB validated and file saved to: {pdb_path}")
        return True
    
    except PDBConstructionException as e:
        with PDB_UI_OUTPUT_WIDGET:
            print(f"Invalid PDB file: {e}")
        logger.error(f"Invalid PDB file: {e}")
        return False    
    
    except Exception  as e:
        with PDB_UI_OUTPUT_WIDGET:
            print(f"Invalid PDB file: {e}")
        logger.error(f"Invalid PDB file: {e}")
        return False    

def upload_and_save_file(
    save_directory: str=BASE_PATH, 
    description: str = "Upload File", 
    filetypes: str = '',
    on_success: Optional[Callable[[str], None]] = None

    ):
    """
    Creates a file upload widget and saves the uploaded file to the specified directory.
    Returns:
        widgets.FileUpload: The file upload widget instance.
    """
    uploader = widgets.FileUpload(
        accept=filetypes,
        multiple=False,
        description=description
    )
    
    #Display the uploader widget
    def on_uploader_change(change):
        global uploaded_pdb_path 
        if uploader.value:
            uploaded_file = uploader.value[0]
            saved_path = os.path.join(save_directory, uploaded_file['name'])

            if uploaded_file['name'].endswith('.json'):
                with open(saved_path, 'wb') as f:
                    f.write(bytes(uploaded_file['content']))
                    with JSON_UI_OUTPUT_WIDGET:
                        print(f"File uploaded and saved as: {saved_path}")
                selected_json_target_path_widget.value = saved_path
                #logger.info(f"File uploaded and saved as: {saved_path}")
                if on_success:
                    with JSON_UI_OUTPUT_WIDGET:
                        print("....Updating Json Fields")
                    on_success(saved_path)
                return

            elif uploaded_file['name'].endswith('.pdb'):
                name = uploaded_file['name']
                pdb_path = os.path.join(BASE_PATH, 'inputs', name)
                os.makedirs(os.path.dirname(pdb_path), exist_ok=True)
                
                with open(pdb_path, 'wb') as f:
                    f.write(uploaded_file['content'])
                validate_pdb_file(pdb_path)
                return
        
    uploader.observe(on_uploader_change, names='value')
    display(uploader)
    return uploader

def refresh_target_editor(new_json_path):
    selected_json_target_path_widget.value = new_json_path

    # Launch a fresh editor UI into the container
    editor_ui = main_launch_target_editor(
        json_target_path = new_json_path,
        path_output_widget = selected_json_target_path_widget
    )
    target_editor_container.children = [editor_ui]


def launch_all_ui():
    """
    Launch the target editor, then the BindCraft UI.
    The JSON path widget is updated after save and passed into BindCraft launcher.
    """
    global target_editor_container

    display(step_box("""
        <strong>Welcome to the BindCraft User Interface</strong><br><br>
        There are three basic steps to running Bindcraft:<br>
        1. Select a target JSON file.<br>
        2. Select a <em>settings filters</em> file (<code>settings_filters.json</code>) — this determines the filters used when selecting binding candidates.<br>
        3. Select an <em>advanced settings</em> file (<code>settings_advanced.json</code>).<br><br>

        The default JSON files are a good place to start.<br>
        To learn more, I encourage you to check out the BindCraft repo README:
        <a href="https://github.com/martinpacesa/BindCraft" target="_blank">GitHub Repo</a><br>
        as well as the wiki:
        <a href="https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft" target="_blank">BindCraft Wiki</a><br><br>

        <strong>The Target JSON file</strong> is the high-level director for the Bindcraft run.<br>

        This is where you can enter:<br>
        1. <strong>File Name</strong> — Name of your JSON file.<br>
        2. <strong>Design Path</strong> — Directory to place designs.<br>
        3. <strong>Binder Name</strong> — Prefix for your output files.<br>
        4. <strong>Starting PDB</strong> — Path to your target PDB file (for Runpod, it will always be in <code>/workspace/inputs</code>).<br>
        5. <strong>Chain</strong> — The chain from the input PDB that you wish to design binders to.<br>
        6. <strong>Target Hotspot Residues</strong> — Area of focus in your PDB (not strictly enforced — Bindcraft may design outside these if metrics are much more favorable).<br>
        7. <strong>Lengths</strong> — Range of amino acid lengths for the designed binders.<br>
        8. <strong>Number of Final Designs</strong> — Number of output binders that pass filters.<br>
        """))
        
    display(step_box("<strong>Step 1:</strong> Optional: Upload a JSON file to your settings_target directory"))
    
    upload_and_save_file(
        save_directory=f'{BASE_PATH}/settings_target',
        description="Upload JSON",
        filetypes='.json',
        on_success = refresh_target_editor
)   
    display(JSON_UI_OUTPUT_WIDGET)

    #Load the default target json at startup
    display(step_box("Optional: If needed, edit or create a new target JSON file."))
    refresh_target_editor(DEFAULT_TARGET_JSON)
    display(target_editor_container)

    display(step_box(f"<strong> Step 2:</strong> If not aleady present, upload a PDB file to your {BASE_PATH}/inputs directory."))
    
    upload_and_save_file(
        save_directory=f'{BASE_PATH}/inputs',
        description="Upload PDB",
        filetypes='.pdb'
)   
    display(PDB_UI_OUTPUT_WIDGET)

    display(step_box("Step 3: Select your settings files and run BindCraft"))
    main_launch_bindcraft_UI(
        json_target_path_widget=selected_json_target_path_widget,
        settings_dirs=SETTINGS_DIRS,
        base_path=BASE_PATH,
        bindcraft_template_run_file=BINDCRAFT_TEMPLATE_PATH,
        output_dir=BASE_PATH
    )


