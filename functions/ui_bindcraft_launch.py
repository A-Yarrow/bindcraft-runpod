import os
import subprocess
from glob import glob
from string import Template
import ipywidgets as widgets
from IPython.display import display
from functools import partial

# GLOBAL HOLDERS
selected_paths_holder = {}

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
    return settings_widget_dict


def on_submit_settings_clicked(button, 
                               base_path:str, 
                               target_json_path:str, 
                               bindcraft_template_run_file, 
                               output_dir:str, 
                               settings_widget_dict: dict) -> str:
    """Save selected paths from dropdowns to global holder."""
    selected_paths = {os.path.basename(k): dd.value for k, dd in settings_widget_dict.items()}
    selected_paths_holder.clear()
    selected_paths_holder.update(selected_paths)
    print("Selected settings:")
    for k, v in selected_paths_holder.items():
        print(f"  {k}: {v}")

    """Fill template with settings paths and run bindcraft_run.sh."""
    for folder, file in selected_paths_holder.items():
        full_path = os.path.join(base_path, folder, file)
        selected_paths_holder[folder] = full_path

    try:
        filters_path = selected_paths_holder['settings_filters']
        advanced_path = selected_paths_holder['settings_advanced']
        selected_paths_holder['settings_target'] = target_json_path
    except KeyError:
        print("Please click the Submit button to select settings files.")
        return

    if not filters_path or not advanced_path or not target_json_path:
        print("Missing required paths. Make sure all fields are selected.")
        return

    print(f"Using:")
    print(f"  Filters:  {filters_path}")
    print(f"  Advanced: {advanced_path}")
    print(f"  Target:   {target_json_path}")

    with open(bindcraft_template_run_file, 'r') as f:
        template = Template(f.read())

    new_template = template.substitute(
        FILTERS_FILE_PATH=filters_path,
        ADVANCED_FILE_PATH=advanced_path,
        TARGET_FILE_PATH=target_json_path,
        TARGET_FILE_NAME=os.path.basename(target_json_path),
        TARGET_NAME=os.path.splitext(os.path.basename(target_json_path))[0],
        LOG_DIR=f"{output_path}/{os.path.splitext(os.path.basename(target_json_path))[0]}"
    )

    bindcraft_run_file = bindcraft_template_run_file.replace('template', '')
    with open(bindcraft_run_file, 'w') as out:
        out.write(new_template)
    
    os.chmod(bindcraft_run_file, 0o755)
    print(f"🚀 Script written to {bindcraft_run_file}. Now launching...")
    run_bindcraft(bindcraft_run_file)

def run_bindcraft(bindcraft_run_file:str=None) -> None:
    """Run the updated bincraft run file."""
    if not os.path.exists(bindcraft_run_file):
        print(f"Template file {bindcraft_run_file} does not exist.")
        return
    try:
        subprocess.run(["bash", bindcraft_run_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"BindCraft failed: {e}")


def main_launch_bindcraft_UI(
    settings_dirs: list,
    base_path: str,
    bindcraft_template_run_file: str,
    output_dir: str,
    target_json_path_widget_input: widgets.Text = None):
    
    """Main UI to select settings and run BindCraft."""
    settings_widget_dict = settings_widget(settings_dirs)

    submit_button = widgets.Button(
        description='Run BindCraft With Selected Settings',
        button_style='success',
        layout=widgets.Layout(width='30%')
    )

    submit_button.on_click(partial(
        on_submit_settings_clicked,
        base_path=base_path,
        target_json_path=target_json_path_widget_input.value if target_json_path_widget_input else '',
        bindcraft_template_run_file=bindcraft_template_run_file,
        output_dir=output_dir,
        settings_widget_dict=settings_widget_dict
    ))

    run_button = widgets.Button(
        description='Run BindCraft',
        button_style='primary',
        layout=widgets.Layout(width='30%')
    )

    if target_json_path_widget_input is None:
        target_json_path_widget_input = widgets.Text(
            value='',
            description='Target JSON path:',
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='80%')
        )

    run_button.on_click(lambda b: run_bindcraft(
        base_path,
        target_json_path_widget_input.value,
        bindcraft_template_run_file,
        output_dir
    ))

    print("�� Select your settings files:")
    for directory, widget in settings_widget_dict.items():
        display(widget)

    display(target_json_path_widget_input)
    display(submit_button)
    display(run_button)

