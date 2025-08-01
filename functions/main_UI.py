import os
import sys
import logging
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException
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

UI_OUTPUT_WIDGET = widgets.Output()
uploaded_pdb_path = None

def validate_pdb_file(pdb_path: str) -> bool:
    """
    Validate if the provided PDB file is valid.
    Returns True if valid, False otherwise.
    """
    parser = PDBParser(QUIET=False)
    try:
        structure = parser.get_structure('temp', pdb_path)
        with UI_OUTPUT_WIDGET:
            print(f"Valid PDB file: {pdb_path}")
        logger.info(f"Valid PDB file: {pdb_path}")
        return True
    
    except PDBConstructionException as e:
        with UI_OUTPUT_WIDGET:
            print(f"Invalid PDB file: {e}")
        logger.error(f"Invalid PDB file: {e}")
        return False    
    
    except Exception  as e:
        with UI_OUTPUT_WIDGET:
            print(f"Invalid PDB file: {e}")
        logger.error(f"Invalid PDB file: {e}")
        return False    

def upload_and_save_file(save_directory: str=BASE_PATH, description: str = "Upload File", filetypes: str = ''):
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
                    with UI_OUTPUT_WIDGET:
                        print(f"File uploaded and saved as: {saved_path}")
                selected_json_target_path_widget.value = saved_path
                #logger.info(f"File uploaded and saved as: {saved_path}")
            
            elif uploaded_file['name'].endswith('.pdb'):
                with open(saved_path, 'wb') as f:
                    f.write(uploaded_file['content'])
                uploaded_pdb_path = saved_path 
                validate_pdb_file(uploaded_pdb_path)

            logger.info(f"File uploaded and saved as: {saved_path}")
            with UI_OUTPUT_WIDGET:
                print(f"File uploaded and saved as: {saved_path}")
        
        
    uploader.observe(on_uploader_change, names='value')
    display(uploader)
    return uploader

def launch_all_ui():
    """
    Launch the target editor, then the BindCraft UI.
    The JSON path widget is updated after save and passed into BindCraft launcher.
    """
    print("Optional: Upload a JSON file to your settings_target directory.")
    upload_and_save_file(
    save_directory=f'{BASE_PATH}/settings_target',
    description="Upload JSON",
    filetypes='.json'
)
    print("Optional: Upload a PDB file to your settings_target directory.")
    upload_and_save_file(
    save_directory=f'{BASE_PATH}/inputs',
    description="Upload PDB",
    filetypes='.pdb'
)   
    display(UI_OUTPUT_WIDGET)

    print("Step 1: Edit or Update your target JSON file.")
    main_launch_target_editor(
        json_target_path_widget=selected_json_target_path_widget,
        path_output_widget=selected_json_target_path_widget
    )

    print("\nStep 2: Select settings and run BindCraft.")
    main_launch_bindcraft_UI(
        json_target_path_widget=selected_json_target_path_widget,
        settings_dirs=SETTINGS_DIRS,
        base_path=BASE_PATH,
        bindcraft_template_run_file=BINDCRAFT_TEMPLATE_PATH,
        output_dir=BASE_PATH
    )
    

