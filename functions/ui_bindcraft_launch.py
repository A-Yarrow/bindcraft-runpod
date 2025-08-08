import os
import logging
import subprocess
from glob import glob
from string import Template
import ipywidgets as widgets
from IPython.display import display
from functools import partial
from settings import ENV, SETTINGS
from main_UI import logger
import threading
import time


# SETTINGS
DEFAULT_SETTINGS_FILTER = os.path.join(SETTINGS[f'{ENV}_SETTINGS_DIRS'][0], 
SETTINGS[f'{ENV}_DEFAULT_SETTINGS_FILTER'])
DEFAULT_SETTINGS_ADVANCED = os.path.join(SETTINGS[f'{ENV}_SETTINGS_DIRS'][1], 
SETTINGS[f'{ENV}_DEFAULT_SETTINGS_ADVANCED'])

# GLOBAL HOLDERS
selected_paths_holder = {}
refresh_dropdown_output_box = widgets.Output()
bindcraft_launch_box = widgets.Output()
LOG_FILE = None


def settings_widget(dirs: list) -> dict:
    """Build dropdown widgets for JSON settings files in specified directories."""
    settings_widget_dict = {}
    for d in dirs:
        file_list = sorted(glob(f'{d}/*.json'))
        filenames = [os.path.basename(file) for file in file_list]
        file_dropdown = widgets.Dropdown(
            options=filenames,
            description=f'{os.path.basename(d)} files:',
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='50%')
        )
        settings_widget_dict[d] = file_dropdown
    
    #Set defaults
    settings_widget_dict[SETTINGS[f'{ENV}_SETTINGS_DIRS'][0]].value = os.path.basename(DEFAULT_SETTINGS_FILTER)
    settings_widget_dict[SETTINGS[f'{ENV}_SETTINGS_DIRS'][1]].value = os.path.basename(DEFAULT_SETTINGS_ADVANCED)
    
    return settings_widget_dict

def on_submit_settings_clicked(button, 
                               base_path:str, 
                               json_target_path_widget:widgets.Text, 
                               bindcraft_template_run_file, 
                               output_dir:str, 
                               settings_widget_dict: dict,
                               json_target_dropdown: widgets.Dropdown = None) -> str:
    
    # Use dropdown if selected (widget is not none and contains a value), else use text widget
    global LOG_FILE
    json_target_path = json_target_dropdown.value.strip() if json_target_dropdown and json_target_dropdown.value else json_target_path_widget.value.strip()

    with bindcraft_launch_box:
        logger.info(f'Json path set to:{json_target_path}')

    selected_paths = {os.path.basename(k): dd.value for k, dd in settings_widget_dict.items()}
    selected_paths_holder.clear()
    selected_paths_holder.update(selected_paths)
    print("Selected settings:")
    for k, v in selected_paths_holder.items():
        with bindcraft_launch_box:
            print(f"  {k}: {v}")

    """Fill template with settings paths and run bindcraft_run.sh."""
    for folder, file in selected_paths_holder.items():
        full_path = os.path.join(base_path, folder, file)
        selected_paths_holder[folder] = full_path

    try:
        filters_path = selected_paths_holder['settings_filters']
        advanced_path = selected_paths_holder['settings_advanced']
        selected_paths_holder['settings_target'] = json_target_path
    except KeyError:
        with bindcraft_launch_box:
            logger.exception("KeyError: Missing one or more settings paths.")
        return

    if not filters_path or not advanced_path or not json_target_path:
        with bindcraft_launch_box:
            logger.debug("Missing required paths in selected_paths_holder.")
        return

    print(f"Using:")
    print(f"  Filters:  {filters_path}")
    print(f"  Advanced: {advanced_path}")
    print(f"  Target:   {json_target_path}")

    with open(bindcraft_template_run_file, 'r') as f:
        template = Template(f.read())

    TARGET_NAME = os.path.splitext(os.path.basename(json_target_path))[0]
    LOG_DIR = f"{output_dir}/outputs/{TARGET_NAME}"
    new_template = template.substitute(
        FILTERS_FILE_PATH=filters_path,
        ADVANCED_FILE_PATH=advanced_path,
        TARGET_FILE_PATH=json_target_path,
        TARGET_FILE_NAME=os.path.basename(json_target_path),
        TARGET_NAME=TARGET_NAME,
        LOG_DIR=LOG_DIR
    )
    LOG_FILE = f"{LOG_DIR}/{TARGET_NAME}-bindcraft_log.txt"

    bindcraft_run_file = bindcraft_template_run_file.replace('_template', '')
    with open(bindcraft_run_file, 'w') as out:
        out.write(new_template)
        logger.debug(f"BindCraft run script written to: {bindcraft_run_file}")
    os.chmod(bindcraft_run_file, 0o755)
    with bindcraft_launch_box:
        logger.info(f"Script written to {bindcraft_run_file}.")
        logger.info(f"To tail log in a terminal, run: tail -f {LOG_FILE}")
    
    return LOG_FILE

def tail_log_widget(log_file: str, N: int = 30, refresh: float = 1.0):
    """Display and live-update the last N lines of a log file in Jupyter."""
    if not os.path.exists(log_file):
        raise FileNotFoundError(f"Log file not found: {log_file}")
    
    # Text area widget
    output_box = widgets.Textarea(
        value="",
        placeholder="Waiting for log output...",
        description="Log:",
        layout=widgets.Layout(width="100%", height="300px"),
        disabled=True
    )
    display(output_box)

    # Function to update log output
    def update_loop():
        while True:
            try:
                with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                    lines = f.readlines()[-N:]
                output_box.value = "".join(lines)
            except Exception as e:
                output_box.value = f"Error reading log: {e}"
            time.sleep(refresh)

    # Start background thread
    t = threading.Thread(target=update_loop, daemon=True)
    t.start()


def run_bindcraft(button=None, bindcraft_run_file: str = None) -> None:
    """Run the updated BindCraft run file, unless in DEV mode."""
    with bindcraft_launch_box:
        logger.info(f"run_bindcraft() triggered with: {bindcraft_run_file}")
    from settings import ENV  # Import here if ENV isn't global

    if ENV == 'DEV':
        with bindcraft_launch_box:
            logger.debug("Skipping BindCraft execution in DEV mode.")
        return

    if not bindcraft_run_file or not os.path.exists(bindcraft_run_file):
        with bindcraft_launch_box:
            logger.error(f"BindCraft run file does not exist or is invalid: {bindcraft_run_file}")
        return

    with bindcraft_launch_box:
        logger.info(f"Running BindCraft with script: {bindcraft_run_file}")
    
    try:
        subprocess.run(["bash", bindcraft_run_file], check=True)
    except subprocess.CalledProcessError as e:
        with bindcraft_launch_box:
            logger.exception(f"BindCraft run launch failed: {e}")
    
    tail_log_widget(LOG_FILE)

def main_launch_bindcraft_UI(
    json_target_path_widget: widgets.Text,
    settings_dirs: list,
    base_path: str,
    bindcraft_template_run_file: str,
    output_dir: str):

    def refresh_dropdowns(b):
        with refresh_dropdown_output_box:
            print("Refresh clicked")
        target_files = sorted(glob(f'{base_path}/settings/*.json'))
        json_target_dropdown.options = [os.path.basename(f) for f in target_files]
    
    """Main UI to select settings and run BindCraft."""
    settings_widget_dict = settings_widget(settings_dirs)
    
    """Json target dropdwon widget"""
    json_target_dropdown = widgets.Dropdown(
    options=sorted(glob(f'{base_path}/settings_target/*.json')),
    description='Select Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='70%')
    )
    refresh_button = widgets.Button(
    description="🔄 Refresh Json Dropdown Menu",
    button_style='warning',
    layout=widgets.Layout(width='30%')
    )
    
    refresh_button.on_click(refresh_dropdowns)
    

    submit_settings_button = widgets.Button(
        description='Generate BindCraft Run Script with Settings',
        button_style='success',
        layout=widgets.Layout(width='50%')
    )

    submit_settings_button.on_click(partial(on_submit_settings_clicked,
        base_path=base_path,
        json_target_path_widget=json_target_path_widget,
        bindcraft_template_run_file=bindcraft_template_run_file,
        output_dir=output_dir,
        settings_widget_dict=settings_widget_dict,
        json_target_dropdown=json_target_dropdown
    ))

    run_bindcraft_button = widgets.Button(
        description='Run BindCraft',
        button_style='primary',
        layout=widgets.Layout(width='30%')
    )

    run_bindcraft_button.on_click(
        partial(run_bindcraft,
        bindcraft_run_file = bindcraft_template_run_file.replace('_template', '')))
    
    display(json_target_dropdown)
    
    display(refresh_button)
    display(refresh_dropdown_output_box)
    display(widgets.Label("Remember to click 'Refresh Json Dropdown' after uploading, editing or saving a new Target Json file"))
    for directory, widget in settings_widget_dict.items():
        display(widget)
    
    display(submit_settings_button)
    display(run_bindcraft_button)
    display(bindcraft_launch_box)

