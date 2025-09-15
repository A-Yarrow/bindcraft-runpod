BindCraft-RunPod

Streamline Your Protein Engineering with automated cloud deployment and a user-friendly Interface. For detailed instructions on set-up please see this blog post.

Quick Start for Protein Engineers: For those eager to begin, I have developed a specialized tool to simplify your workflow. 

You can immediately access and utilize BindCraft with a pre-configured environment on RunPod (plus an intuitive GUI). This allows you to bypass complex setup and focus directly on your design tasks. Click here to try the BindCraft RunPod template: https://console.runpod.io/deploy?template=sfpbu5gqqw&ref=tt816r6b


This repository is not a fork of the original BindCraft repo, but a repackaged, cloud-optimized version that:

	Runs in a RunPod container with a lightweight Docker image.

	PyRosetta and AF2 weights and PyRosetta are auto-installed after container mounting.
	
 	Offers a Jupyter-based GUI for editing targets and launching jobs

	Supports user-provided design inputs via mounted volumes


⚠️ This repo reuses key components from the original BindCraft repository (MIT License), including settings templates and backend logic. 

Full credit and detailed documentation for BindCraft is available at: https://github.com/martinpacesa/BindCraft Preprint


Quick Start on RunPod
Clone the latest development repository:

git clone -b dev https://github.com/A-Yarrow/bindcraft-runpod.git

cd bindcraft-runpod

Build and Push the Docker image:

docker build -t your-dockerhub-login-name/bindcraft-cloud:yourtag.

docker push your-dockerhub-login-name/bindcraft-cloud:yourtag

Choose an appropriate Cloud GPU. Alternatively, use my RunPod Template: https://console.runpod.io/deploy?template=sfpbu5gqqw&ref=tt816r6b

	This Template has been tested and works with RunPods: A100 SXM, H100 SXM, L40, L40S

Open Jupyter: JupyterLab will be available on port 8888. Use these detailed  instructions to launch your first job.


Jupyter Notebook and GUI
A starter notebook bindcraft-runpod-start.ipynb is provided to:

	Help users mount settings and parameter folders

	Launch designs using bindcraft.py

	View outputs and logs

	Upload custom target files

	GUI Workflow

The user-friendly GUI treamlines the design process:

	Uploads:
	
 		Upload PDB or Create Target JSON:

		Upload a PDB file directly. The GUI will check for chain breaks and other inconsistencies.

		Alternatively, create a new target JSON file by filling in the required fields (design path, starting PDB, hotspot residues, binder lengths, and number of designs).

		Submit Job:

			Enter a unique job name.
			
   			Select the target JSON file from the dropdown menu.

			Click "Submit Job" to start the design process.

	Monitor Progress:

		Logs are stored in a subdirectory within the outputs folder, labeled with the job name.

What’s Inside This Repo

	bindcraft-runpod/

		├── Dockerfile                   # Main build script

		├── install_bindcraft.sh         # Custom install script for BindCraft with mamba/conda

		├── start.sh                     # Launch script that copies settings and starts Jupyter

		├── bindcraft-runpod-start.ipynb # Starter notebook for launching designs

	├── bindcraft.py                 # Main binder design script (from upstream)

	├── settings_target/             # Example target config

	├── settings_filters/            # Default filter settings

	├── settings_advanced/           # Default design settings

	└── README.md

Why I Built This

	As a protein design engineer, I needed a reproducible and cloud-friendly setup to:

	Run high-throughput binder design without local GPU constraints

	Use RunPod or other cloud GPU services with persistent state
	
 	Share a ready-to-use Docker image for others to reproduce my workflows

	This repo lets users skip complex local setup and run everything via a containerized workflow, while still preserving full flexibility of BindCraft’s powerful design stack.

Acknowledgments

This project reuses and repackages tools developed by:

Martin Pacesa et al. — BindCraft
Sergey Ovchinnikov — ColabDesign
Justas Dauparas — ProteinMPNN
RosettaCommons — PyRosetta

BindCraft is licensed under the MIT License. ColabDesign and ProteinMPNN are licensed under Apache 2.0.

Feel free to reach out via GitHub Issues if you encounter bugs or would like to contribute.

