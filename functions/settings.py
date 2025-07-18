#ENV = 'DEV'  # Change to 'PROD' for development settings
ENV = 'PROD'  # Change to 'DEV' for development settings
SETTINGS = {
    # Production
    'PROD_RUN_DIR': '/workspace',
    'PROD_BINDCRAFT_DIR': '/app/bindcraft',

    'PROD_SETTINGS_DIRS': [
        '/workspace/settings_filters',
        '/workspace/settings_advanced'
    ],
    'PROD_BINDCRAFT_TEMPLATE_PATH': '/app/bindcraft/bindcraft_run_template.sh',
    'PROD_BINDCRAFT_OUTPUT_PATH': '/workspace/bindcraft_run.sh',
    'PROD_FUNCTIONS_PATH': '/app/bindcraft/functions',
    'PROD_DEFAULT_TARGET_JSON': '/workspace/settings_target/PDL1.json',

    # Development
    'DEV_RUN_DIR': '/home/yarrow/projects/bindcraft-runpod',
    'DEV_BINDCRAFT_DIR': '/home/yarrow/projects/bindcraft-runpod',

    'DEV_SETTINGS_DIRS': [
        '/home/yarrow/projects/bindcraft-runpod/settings_filters',
        '/home/yarrow/projects/bindcraft-runpod/settings_advanced'
    ],
    'DEV_BINDCRAFT_TEMPLATE_PATH': '/home/yarrow/projects/bindcraft-runpod/bindcraft_run_template.sh',
    'DEV_BINDCRAFT_OUTPUT_PATH': '/home/yarrow/projects/bindcraft-runpod/bindcraft_run.sh',
    'DEV_FUNCTIONS_PATH': '/home/yarrow/projects/bindcraft-runpod/functions',
    'DEV_DEFAULT_TARGET_JSON': '/home/yarrow/projects/bindcraft-runpod/settings_target/PDL1.json'
}
