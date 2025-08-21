# file: functions/ui_target_editor.py

import json
import logging
import os
from pathlib import Path
from functools import partial
import ipywidgets as widgets
from IPython.display import display
from main_UI import logger

editor_output_box = widgets.Output()

def load_target_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def save_target_json(filename, json_target_path, data):
    base_dir = os.path.dirname(json_target_path)
    target_filepath = os.path.join(base_dir)
    os.makedirs(target_filepath, exist_ok=True)
    full_path = os.path.join(target_filepath, filename)

    with open(full_path, 'w') as f:
        json.dump(data, f, indent=4)
    with editor_output_box:
        logger.info(f"Data saved to {full_path}")

def update_data(data, widget_dict):
    updated_data = {}
    for key, widget in widget_dict.items():
        val = widget.value
        if isinstance(data[key], list):
            try:
                val = json.loads(val)
            except json.JSONDecodeError:
                with editor_output_box:
                    logger.error(f"Error decoding JSON for key '{key}', keeping original value.")
                val = data[key]
        updated_data[key] = val
    return updated_data

def on_save_clicked(button, 
                    data, widget_dict, 
                    filename_widget, 
                    json_target_path,
                    path_output_widget):
    
    filename = filename_widget.value.strip()
    updated_data = update_data(data, widget_dict)
    save_target_json(filename, json_target_path, updated_data)
    
    new_path = os.path.join(os.path.dirname(json_target_path),filename)
    path_output_widget.value = new_path 

def target_editor_widget(data):
    widget_dict = {}
    for k, v in data.items():
        if isinstance(v, list):
            widget_dict[k] = widgets.Text(value=str(v), description=f'{k}:',
                                          style={'description_width': 'initial'})
        elif isinstance(v, int):
            widget_dict[k] = widgets.IntText(value=v, description=f'{k}:',
                                             style={'description_width': 'initial'})
        else:
            widget_dict[k] = widgets.Text(value=v, description=f'{k}:',
                                          style={'description_width': 'initial'})
    return widget_dict

def main_launch_target_editor(json_target_path: str, path_output_widget: widgets.Text) -> None:
    """
    Launch the target editor UI to edit a JSON file.
    """
    data = load_target_json(json_target_path)

    filename_widget = widgets.Text(
        value=os.path.basename(json_target_path),
        description='File Name:',
        layout=widgets.Layout(width='50%'),
        style={'description_width':'initial'}
    )

    widget_dict = target_editor_widget(data)

    save_button = widgets.Button(description="Save Changes", button_style='success')
    
    
    save_button.on_click(partial(
        on_save_clicked,
        data=data,
        widget_dict=widget_dict,
        filename_widget=filename_widget,
        json_target_path=json_target_path,
        path_output_widget=path_output_widget
    ))

    
    all_widgets = [filename_widget] + list(widget_dict.values()) + [save_button, editor_output_box]
    return widgets.VBox(all_widgets)
     
