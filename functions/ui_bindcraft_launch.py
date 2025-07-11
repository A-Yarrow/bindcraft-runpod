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


def on_submit_clicked(button, settings_widget_dict: dict) -> None:
    """Save selected paths from dropdowns to global holder."""
    selected_paths = {os.path.basename(k): dd.value for k, dd in settings_widget_dict.items()}
    selected_paths_holder.clear()
    selected_paths_holder.update(selected_paths)
    print("✅ Selected settings:")
    for k, v in selected_paths_holder.items():
        print(f"  {k}: {v}")


def run_bindcraft(
    SETTINGS_FILE_PATH: str,
    target_path: str,
    BINDCRAFT_TEMPLATE_PATH: str,
    BINDCRAFT_OUTPUT_PATH: str
) -> None:
    """Fill template with settings paths and run bindcraft_run.sh."""
    for folder, file in selected_paths_holder.items():
        full_path = os.path.join(SETTINGS_FILE_PATH, folder, file)
        selected_paths_holder[folder] = full_path

    try:
        filters_path = selected_paths_holder['settings_filters']
        advanced_path = selected_paths_holder['settings_advanced']
        selected_paths_holder['settings_target'] = target_path
    except KeyError:
        print("❌ Please click the Submit button to select settings files.")
        return

    if not filters_path or not advanced_path or not target_path:
        print("❌ Missing required paths. Make sure all fields are selected.")
        return

    print(f"✅ Using:")
    print(f"  Filters:  {filters_path}")
    print(f"  Advanced: {advanced_path}")
    print(f"  Target:   {target_path}")

    with open(BINDCRAFT_TEMPLATE_PATH, 'r') as f:
        template = Template(f.read())

    filled = template.substitute(
        FILTERS_FILE_PATH=filters_path,
        ADVANCED_FILE_PATH=advanced_path,
        TARGET_FILE_PATH=target_path,
        TARGET_FILE_NAME=os.path.basename(target_path),
        TARGET_NAME=os.path.splitext(os.path.basename(target_path))[0],
        LOG_DIR=f"/workspace/outputs/{os.path.splitext(os.path.basename(target_path))[0]}"
    )

    with open(BINDCRAFT_OUTPUT_PATH, 'w') as out:
        out.write(filled)
    os.chmod(BINDCRAFT_OUTPUT_PATH, 0o755)
    print(f"🚀 Script written to {BINDCRAFT_OUTPUT_PATH}. Now launching...")

    try:
        subprocess.run(["bash", BINDCRAFT_OUTPUT_PATH], check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ BindCraft failed: {e}")


def launch_bindcraft_UI(
    settings_dirs: list = ['/workspace/settings_filters', '/workspace/settings_advanced'],
    SETTINGS_FILE_PATH: str = '/workspace',
    BINDCRAFT_TEMPLATE_PATH: str = '/workspace/bindcraft_run_template.sh',
    BINDCRAFT_OUTPUT_PATH: str = '/workspace/bindcraft_run.sh',
    target_path_widget: widgets.Text = None
):
    """Main UI to select settings and run BindCraft."""

    settings_widget_dict = settings_widget(settings_dirs)

    submit_button = widgets.Button(
        description='Submit Selected Settings',
        button_style='success',
        layout=widgets.Layout(width='30%')
    )

    submit_button.on_click(partial(
        on_submit_clicked,
        settings_widget_dict=settings_widget_dict
    ))

    run_button = widgets.Button(
        description='Run BindCraft',
        button_style='primary',
        layout=widgets.Layout(width='30%')
    )

    if target_path_widget is None:
        target_path_widget = widgets.Text(
            value='',
            description='Target JSON path:',
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='80%')
        )

    run_button.on_click(lambda b: run_bindcraft(
        SETTINGS_FILE_PATH,
        target_path_widget.value,
        BINDCRAFT_TEMPLATE_PATH,
        BINDCRAFT_OUTPUT_PATH
    ))

    print("�� Select your settings files:")
    for directory, widget in settings_widget_dict.items():
        display(widget)

    display(target_path_widget)
    display(submit_button)
    display(run_button)

