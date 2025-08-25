import os
import signal
import logging
import psutil
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
JSON_TARGET_DIR = SETTINGS[f'{ENV}_JSON_TARGET_DIR']
PID_DIR = SETTINGS[f'{ENV}_PID_DIR']
os.makedirs(PID_DIR, exist_ok=True)

# GLOBAL HOLDERS
selected_paths_holder = {}
refresh_dropdown_output_box = widgets.Output()
job_registration_box = widgets.Output()
bindcraft_launch_box = widgets.Output()
json_path_output_box = widgets.Output()
generate_bindcraft_run_script_box = widgets.Output()
refresh_jobs_output_box = widgets.Output()
pid_output_box = widgets.Output()

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
                               env: str = "DEV",
                               job_name_widget:widgets.Text = None,
                               pid_dir:str = None) -> str:
    
    job_name = job_name_widget.value.strip() if job_name_widget else "Default Job"
    # Use dropdown if selected (widget is not none and contains a value), else use text widget
    # Get value from dropdown if available, else from text widget
    # Here the user enters the whole path

    if json_target_dropdown and json_target_dropdown.value:
        filename = json_target_dropdown.value.strip()
    else:
        filename = json_target_path_widget.value.strip()
    
    json_target_path = os.path.join(JSON_TARGET_DIR, filename)
    
    with json_path_output_box:
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
    new_template = template.safe_substitute(
        FILTERS_FILE_PATH=filters_path,
        ADVANCED_FILE_PATH=advanced_path,
        TARGET_FILE_PATH=json_target_path,
        TARGET_FILE_NAME=os.path.basename(json_target_path),
        TARGET_NAME=TARGET_NAME,
        LOG_DIR=LOG_DIR,
        ENV=env,
        JOB_NAME=job_name,
        PID_DIR=pid_dir
        
    )
    LOG_FILE = f"{LOG_DIR}/{TARGET_NAME}-bindcraft_log.txt"
    
    # Ensure log directory exists and create an empty log file
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, 'a') as f:
        pass  # create empty file if it doesn't exist

    bindcraft_run_file = bindcraft_template_run_file.replace('_template', '')
    with open(bindcraft_run_file, 'w') as out:
        out.write(new_template)
    os.chmod(bindcraft_run_file, 0o755)
    with generate_bindcraft_run_script_box:
        logger.info(f"Script written to {bindcraft_run_file}.")
        logger.info(f"To tail log in a terminal, run: tail -f {LOG_FILE}")
    
    return LOG_FILE

def tail_log_widget(log_file: str, pid_file: str, N: int = 30, refresh: float = 1.0):
    """Display and live-update the last N lines of a log file in Jupyter.
       If pid_file is provided, stop tailing when the PID disappears.
       Handles cases where the PID file may appear a little late.
    """
    
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

    def update_loop():
        pid = None
        pid_checked = False  # flag to avoid spamming errors

        while True:
            # Read last N lines of log
            try:
                with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                    lines = f.readlines()[-N:]
                output_box.value = "".join(lines)
            except Exception as e:
                output_box.value = f"Error reading log: {e}"

            job_finished = False

            if pid_file:
                # Try to read PID if it exists
                if os.path.exists(pid_file):
                    if pid is None:
                        for _ in range(5):  # retry loop
                            try:
                                with open(pid_file, "r") as f:
                                    pid = int(f.read().strip())
                                break
                            except Exception:
                                time.sleep(0.2)
                        if pid is None and not pid_checked:
                            output_box.value += f"\n[WARN] Unable to read PID file yet: {pid_file}"
                            pid_checked = True  # only warn once
                    # Check if PID exists. bindcraft_run.sh will clean it up when it exits
                    if pid is not None and not psutil.pid_exists(pid):
                        output_box.value += "\n=== Job Finished ==="
                        job_finished = True
                else:
                    # PID file doesn't exist yet
                    if pid_checked:
                        # already warned, just wait
                        pass
                    else:
                        output_box.value += f"\n[INFO] Waiting for PID file to appear: {pid_file}"
                        pid_checked = True
            # Stop loop if job finished
            if job_finished:
                break

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
    
    #rebuild dictionary using only valid running pids
    running_bindcraft_jobs = {
        job_name:pid
        for job_name, pid in running_bindcraft_jobs.items()
        if psutil.pid_exists(pid)
    }
    
    return running_bindcraft_jobs

def clean_stale_pids(pid_dir):
    for pid_file in glob(f"{pid_dir}/*-bindcraft_pid.txt"):
        try:
            with open(pid_file, "r") as f:
                pid = int(f.read().strip())
            if not psutil.pid_exists(pid):
                os.remove(pid_file)
                with pid_output_box:
                    print(f"Removed stale PID file: {pid_file}")
        except Exception:
            # bad file, just remove it
            os.remove(pid_file)

def abort_run_widget(pid_dir:str) -> widgets.Dropdown:
    running_bindcraft_jobs = get_pids(pid_dir)
    options = sorted(running_bindcraft_jobs.keys())
    value = options[0] if options else None #Empty list evaluates to False
    disabled = len(options) == 0 #True if options is empty list 
    
    return widgets.Dropdown(
        options=options,
        value=value,
        description="Select a job to abort",
        style={'description_width': 'initial'},
        disabled=disabled
    )

def on_abort_run_widget_click(button, 
                              pid_dir: str, 
                              abort_dropdown:widgets.Dropdown, 
                              ):
    running_bindcraft_jobs = get_pids(pid_dir) #refresh each click
    abort_job = abort_dropdown.value
    
    if not abort_job:
        with pid_output_box:
            print ('There are no running jobs to abort')
        return
    
    try:
        abort_pid = running_bindcraft_jobs[abort_job]
        process = psutil.Process(abort_pid)
        process.terminate()
        with pid_output_box:
            print(f'job ID:{abort_job}, with PID:{abort_pid} has been terminated')
        #refresh the abort dropdown
        
        def wait_for_termination():
            while psutil.pid_exists(abort_pid):
                time.sleep(1)
            refresh_abort_dropdown(abort_dropdown)
            with pid_output_box:
                logger.info(f"Job '{abort_job}' has been terminated and dropdown refreshed.")
    
        threading.Thread(target=wait_for_termination, daemon=True).start()
    
    except psutil.NoSuchProcess:
        with pid_output_box:
            print(f"Job '{abort_job}' was not running. Refreshing dropdown.")
            refresh_abort_dropdown(abort_dropdown)
    except Exception as e:
        with pid_output_box:
            logger.warning(f"Error terminating {abort_job}: {e}")

def refresh_jobs(button, abort_dropdown:widgets.Dropdown):
    abort_dropdown.options = list(get_pids(pid_dir=PID_DIR).keys())
    with refresh_jobs_output_box:
        print("Running Jobs List Refreshed")

def refresh_abort_dropdown(abort_dropdown: widgets.Dropdown):
    clean_stale_pids(PID_DIR)
    jobs = list(get_pids(PID_DIR).keys())
    abort_dropdown.options = jobs
    if jobs:
        abort_dropdown.values = jobs[0]
        abort_dropdown.disabled = False
    else:
        abort_dropdown.value = None
        abort_dropdown.disabled = True

def job_name_widget():
    return widgets.Text(
        value=None,
        description='Enter a job name for the Bindcraft Run (Required)',
        disabled=False,
        layout=widgets.Layout(width='50%'),
        style={'description_width':'initial'}
    )

def job_name_on_clicked(button, job_name_widget):
    job_name = job_name_widget.value.strip()
    with job_registration_box:
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
            print("Please enter a job name (required)")
        return

    if env == 'DEV':
        with bindcraft_launch_box:
            logger.debug("Skipping BindCraft execution in DEV mode.")
            logger.info(f"DEV mode: Dummy run for job name: {job_name}")

    try:
        proc = subprocess.Popen(["bash", bindcraft_run_file])
        shell_pid = proc.pid
        # Note: The actual bindcraft PID is written by the .sh script itself
        # into $PID_DIR/$JOB_NAME-bindcraft_pid.txt

        #path to pid file
        pid_file_path = os.path.join(PID_DIR, f"{job_name}-bindcraft_pid.txt")
        with bindcraft_launch_box:
            print(f'pid_file_path: {pid_file_path}')
        #Wait until PID file exists
        for i in range(120):
            if os.path.exists(pid_file_path):
                with bindcraft_launch_box:
                    logger.debug(f"Found pid file: {os.path.basename(pid_file_path)}")
                break
            time.sleep(1)
        else:
            with bindcraft_launch_box:
                logger.error(f"PID file was not created: {pid_file_path}")
                return

        #Read in PID bindcraft_run.sh wrote
        with open(pid_file_path) as f:
            pid = int(f.read().strip())

        #Store PID in memory and in file
        running_bindcraft_jobs[job_name] = pid
            
        with bindcraft_launch_box:
                logger.info(f"Registered BindCraft job Name '{job_name}' with PID {pid}")
                logger.info(f"run_bindcraft() triggered with: '{bindcraft_run_file}'")

    except Exception as e:
        with bindcraft_launch_box:
            logger.exception(f"BindCraft run launch failed: {e}")
    
    if log_file:
        tail_log_widget(log_file, pid_file_path)
    
def main_launch_bindcraft_UI(
    json_target_path_widget: widgets.Text,
    settings_dirs: list,
    base_path: str,
    bindcraft_template_run_file: str,
    output_dir: str):

    
    # Job Name widget
    #default_job_name = os.path.basename(json_target_path_widget.value.strip())[:-5]+'-1'
    job_name_text_widget = job_name_widget()
    
    job_name_button = widgets.Button(
        description='Submit Job Name',
        button_style='primary',
        style={'description_width':'initial'}
    )
    
    job_name_button.on_click(
        partial(job_name_on_clicked, job_name_widget=job_name_text_widget)
        )
    
    def refresh_dropdowns(b):
        with refresh_dropdown_output_box:
            print("Refresh clicked")
        target_files = sorted(glob(f'{base_path}/settings_target/*.json'))
        json_target_dropdown.options = [os.path.basename(f) for f in target_files]
    
    """Helper function to read in job_name widget value"""
    def on_run_bindcraft_clicked(button, abort_dropdown):
        job_name = job_name_text_widget.value.strip() # get current widget value
        #call on_submit_settings_clicked again so you get the new log file name.
        
        log_file = on_submit_settings_clicked(
            button=None,
            base_path=base_path,
            json_target_path_widget=json_target_path_widget,
            bindcraft_template_run_file=bindcraft_template_run_file,
            output_dir=output_dir,
            settings_widget_dict=settings_widget_dict,
            json_target_dropdown=json_target_dropdown,
            env=env,
            job_name_widget=job_name_text_widget,
            pid_dir=PID_DIR
        )
        run_bindcraft(
            bindcraft_run_file = bindcraft_template_run_file.replace('_template', ''),
            log_file=log_file,
            env=env,
            job_name=job_name
        )
        # refresh abort dropdown after job launch
        refresh_abort_dropdown(abort_dropdown)
    
    # JSON TARGET DROPDOWN WIDGET
    
    """Main UI to select settings and run BindCraft."""
    settings_widget_dict = settings_widget(settings_dirs)
    
    """Json target dropdwon widget"""
    json_target_dropdown = widgets.Dropdown(
    options=sorted(glob(f'{base_path}/settings_target/*.json')),
    value=None,
    description='Select Target JSON:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='70%')
    )
    refresh_button = widgets.Button(
    description="🔄 Refresh Json Dropdown Menu",
    button_style='warning',
    layout=widgets.Layout(width='30%'),
    style={'description_width': 'initial'}
    )
    # Status Message
    json_target_status_label = widgets.HTML(
        value = f"<b>Selected Target JSON: </b>{json_target_dropdown.value}")
    
    # Udpate JSON TARGET DROPDOWN CALLBACK
    def update_json_target_dropdown_status(change):
        if change['type'] == 'change' and change['name'] == 'value':
            json_target_status_label.value = f"<b>Selected Target JSON: </b>{change['new']}"
    
    json_target_dropdown.observe(update_json_target_dropdown_status, names='value')
    
    refresh_button.on_click(refresh_dropdowns)
    
    submit_settings_button = widgets.Button(
        description='Generate BindCraft Run Script with Settings',
        button_style='success',
        layout=widgets.Layout(width='50%'),
        style={'description_width': 'initial'}
    )

    submit_settings_button.on_click(partial(on_submit_settings_clicked,
        base_path=base_path,
        json_target_path_widget=json_target_path_widget,
        bindcraft_template_run_file=bindcraft_template_run_file,
        output_dir=output_dir,
        settings_widget_dict=settings_widget_dict,
        json_target_dropdown=json_target_dropdown,
        job_name_widget=job_name_text_widget,
        pid_dir=PID_DIR 
        ))

    run_bindcraft_button = widgets.Button(
        description='Run BindCraft',
        button_style='primary',
        layout=widgets.Layout(width='30%'),
        style={'description_width':'initial'}
        )
    
    clean_stale_pids(pid_dir=PID_DIR)
    abort_dropdown = abort_run_widget(pid_dir=PID_DIR)
    abort_jobs_button = widgets.Button(
        description='Abort a running job',
        button_style='warning',
        style={'description_width':'initial'}
        )
    abort_jobs_button.on_click(
        partial(on_abort_run_widget_click,
                pid_dir=PID_DIR,
                abort_dropdown=abort_dropdown
        )
    )

    refresh_jobs_button = widgets.Button(
    description="🔄 Refresh Jobs Dropdown Menu",
    button_style='warning',
    style={'description_width': 'initial'},
    Layout=widgets.Layout(width = '35%')
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
    
    run_bindcraft_button.on_click(
        partial(on_run_bindcraft_clicked, abort_dropdown=abort_dropdown))
    
    refresh_jobs_button.on_click(lambda b:refresh_abort_dropdown(abort_dropdown))
                                
    #-- Display Widgets --
    display(job_name_text_widget, job_name_button, job_registration_box)
    
    
    display(json_target_dropdown, json_target_status_label, refresh_button, refresh_dropdown_output_box)
  
    display(widgets.Label("Remember to click 'Refresh Json Dropdown' after uploading, editing or saving a new Target Json file"))
    for widget in settings_widget_dict.values():
        display(widget)
    display(submit_settings_button, generate_bindcraft_run_script_box, run_bindcraft_button, bindcraft_launch_box)
    
    #Abort Jobs
    display(abort_dropdown, abort_jobs_button, refresh_jobs_button, pid_output_box)