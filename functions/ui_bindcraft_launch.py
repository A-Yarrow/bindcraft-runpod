import os
import signal
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
PID_DIR = SETTINGS[f'{ENV}_PID_DIR']
os.makedirs(PID_DIR, exist_ok=True)

# GLOBAL HOLDERS
selected_paths_holder = {}
refresh_dropdown_output_box = widgets.Output()
bindcraft_launch_box = widgets.Output()

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
                               json_target_dropdown: widgets.Dropdown = None,
                               env: str = "DEV") -> str:
    
    # Use dropdown if selected (widget is not none and contains a value), else use text widget
    json_target_path = json_target_dropdown.value.strip() if json_target_dropdown and json_target_dropdown.value else json_target_path_widget.value.strip()

    with bindcraft_launch_box:
        logger.info(f'Json path set to:{json_target_path}')

    selected_paths = {os.path.basename(k): dd.value for k, dd in settings_widget_dict.items()}
    selected_paths_holder.clear()
    selected_paths_holder.update(selected_paths)

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
        LOG_DIR=LOG_DIR,
        ENV=env,
        
    )
    LOG_FILE = f"{LOG_DIR}/{TARGET_NAME}-bindcraft_log.txt"
    
    # Ensure log directory exists and create an empty log file
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, 'a') as f:
        pass  # create empty file if it doesn't exist

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

def get_pids(pid_dir:str) -> dict:
    # Load existing PID files
    # Make Dictionary to hold running jobs {job_name: PID}
    running_bindcraft_jobs = {}
    
    job_list = glob(f"{pid_dir}/*-bindcraft_pid.txt")
    for pid_file in job_list:
        job_name = os.path.basename(pid_file).replace("-bindcraft_pid.txt", "")
        with open(pid_file, "r") as f:
            running_bindcraft_jobs[job_name] = int(f.read().strip())
    
    return running_bindcraft_jobs

def job_name_widget(default_job_name):
    return widgets.Text(
        value=f"{default_job_name}",
        description='Enter a job id for the Bindcraft Run (Required)',
        disabled=False,
    )

def job_name_on_clicked(button, job_name_widget):
    job_name = job_name_widget.value.strip()
    with bindcraft_launch_box:
        logger.info(f"Job Name: {job_name} has been registered")

    
def run_bindcraft(button=None, bindcraft_run_file: str = None,
                job_name: str = None,
                log_file: str = None,
                env: str = "DEV") -> None:

    """Run the updated BindCraft run file, unless in DEV mode."""
    running_bindcraft_jobs = get_pids(PID_DIR)
    
    if not bindcraft_run_file or not os.path.exists(bindcraft_run_file):
        with bindcraft_launch_box:
            logger.error(f"BindCraft run file does not exist or is invalid: {bindcraft_run_file}")
        return
    
    if not job_name:
        with bindcraft_launch_box:
            print("Please enter a job name")
        return

    if env == 'DEV':
        with bindcraft_launch_box:
            logger.debug("Skipping BindCraft execution in DEV mode.")
            logger.info(f"DEV mode: Dummy run for job name: {job_name}")

    try:
        proc = subprocess.Popen(["bash", bindcraft_run_file])
        pid = proc.pid
        
        #Store PID in memory and in file
        running_bindcraft_jobs[job_name] = pid
        pid_file_path = os.path.join(PID_DIR, f"{job_name}-bindcraft_pid.txt" )
        
        #Store PID in file
        with open(pid_file_path, "w") as f:
            f.write(str(pid))
            
        with bindcraft_launch_box:
                logger.info(f"Registered BindCraft job Name '{job_name}' with PID {pid}")
        
        with bindcraft_launch_box:
            logger.info(f"run_bindcraft() triggered with: '{bindcraft_run_file}'")

    except subprocess.CalledProcessError as e:
        with bindcraft_launch_box:
            logger.exception(f"BindCraft run launch failed: {e}")
    
    if log_file:
        tail_log_widget(log_file)

def main_launch_bindcraft_UI(
    json_target_path_widget: widgets.Text,
    settings_dirs: list,
    base_path: str,
    bindcraft_template_run_file: str,
    output_dir: str):

    def refresh_dropdowns(b):
        with refresh_dropdown_output_box:
            print("Refresh clicked")
        target_files = sorted(glob(f'{base_path}/settings_target/*.json'))
        json_target_dropdown.options = [os.path.basename(f) for f in target_files]
    
    """Helper function to read in job_name widget value"""
    def on_run_bindcraft_clicked(button):
        job_name = job_name_text_widget.value.strip() # get current widget value
        run_bindcraft(
            bindcraft_run_file = bindcraft_template_run_file.replace('_template', ''),
            log_file=log_file,
            env=env,
            job_name=job_name
        )
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

    # Job Name widget
    default_job_name = os.path.basename(json_target_path_widget.value.strip())[:-5]+'-1'
    job_name_text_widget = job_name_widget(default_job_name)
    
    job_name_button = widgets.Button(
        description='Submit Job Name',
        button_style='primary',
        layout=widgets.Layout(width='30%')
    )
    
    job_name_button.on_click(
        partial(job_name_on_clicked, job_name_widget=job_name_text_widget)
        )


    run_bindcraft_button = widgets.Button(
        description='Run BindCraft',
        button_style='primary',
        layout=widgets.Layout(width='30%')
        )

    env = ENV #From settings
    # Pre-generate run script for logging
    log_file = on_submit_settings_clicked(
    button=None,
    base_path=base_path,
    json_target_path_widget=json_target_path_widget,
    bindcraft_template_run_file=bindcraft_template_run_file,
    output_dir=output_dir,
    settings_widget_dict=settings_widget_dict,
    json_target_dropdown=json_target_dropdown,
    env=env,
        )
    
    run_bindcraft_button.on_click(on_run_bindcraft_clicked)
 
    
    display(json_target_dropdown)
    
    display(refresh_button)
    display(refresh_dropdown_output_box)
    display(widgets.Label("Remember to click 'Refresh Json Dropdown' after uploading, editing or saving a new Target Json file"))
    for directory, widget in settings_widget_dict.items():
        display(widget)
    
    display(submit_settings_button)
    display(job_name_text_widget)
    display(job_name_button)
    display(run_bindcraft_button)
    display(bindcraft_launch_box)

