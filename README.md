# BindCraft-Cloud

**Cloud-native Docker deployment of the [BindCraft](https://github.com/martinpacesa/BindCraft) binder design pipeline**, with support for RunPod, GPU acceleration, and persistent volume-based workflows.

This repository is **not a fork** of the original BindCraft repo, but a repackaged, cloud-optimized version that:
- Installs BindCraft and its dependencies (including ColabDesign, MPNN, JAX, and optionally PyRosetta)
- Runs in a RunPod container with a lightweight Docker image
- Supports user-provided design inputs via mounted volumes
- Offers Jupyter-based interaction for editing targets and launching jobs
- Allows persistent storage of AlphaFold2 weights and PyRosetta installs

> ⚠️ This repo reuses key components from the original BindCraft repository (MIT License), including settings templates and backend logic.
> Full credit and detailed documentation for BindCraft is available at:  
>  https://github.com/martinpacesa/BindCraft  
>  [Preprint](https://www.biorxiv.org/content/10.1101/2024.09.30.615802)

---

## What This Project Does

This repo provides:
- A **Dockerfile** to install BindCraft with all core dependencies, excluding weights
- A **Jupyter-based entrypoint** to design binders on RunPod (or other GPU cloud platforms)
- Scripts to install AlphaFold2 weights and PyRosetta **after launch**, allowing reuse across sessions via mounted volumes
- Volume-mount support for persistent configuration and output directories

---

## Quick Start on RunPod

### 1. Clone this repository:
```bash
git clone https://github.com/YOUR_USERNAME/BindCraft-Cloud
cd BindCraft-Cloud

2. Build the Docker image:
bash
Copy
Edit
docker build -t bindcraft-cloud:latest .
3. Launch container (with GPU and mounted workspace):
bash
Copy
Edit
docker run --gpus all -dit --name bindcraft_session \
  -v /your/local/workspace:/workspace \
  bindcraft-cloud:latest
4. Open Jupyter:
JupyterLab will be available on port 8888 with no password or token.

Installing AlphaFold2 Weights and PyRosetta
These are not included in the image (due to size and licensing), but you can install them inside the container using provided instructions in the notebook or script.

AlphaFold2 weights (~5.3 GB)

PyRosetta (requires a license for commercial use)

Persistent Volume Layout
Volume Path	Description
/workspace/settings_target/	User-defined binder target configs (JSON)
/workspace/settings_filters/	Design filter configs
/workspace/settings_advanced/	Advanced model config
/workspace/outputs/	Output designs and logs
/workspace/params/	AlphaFold2 weights & PyRosetta install (optional)

Jupyter Notebook
A starter notebook bindcraft-runpod-start.ipynb is provided to:

Help users mount settings and parameter folders

Launch designs using bindcraft.py

View outputs and logs

Upload custom target files

What’s Inside This Repo
bash
Copy
Edit
bindcraft-cloud/
├── Dockerfile                   # Main build script
├── install_bindcraft.sh         # Custom install script for BindCraft with mamba/conda
├── start.sh                     # Launch script that copies settings and starts Jupyter
├── bindcraft-runpod-start.ipynb # Starter notebook for launching designs
├── bindcraft.py                 # Main binder design script (from upstream)
├── functions/                   # ColabDesign utilities and backend logic
│   └── colabdesign_utils.py
├── settings_target/             # Example target config
├── settings_filters/            # Default filter settings
├── settings_advanced/           # Default design settings
└── README.md
🧠 Why I Built This
As a protein design engineer, I needed a reproducible and cloud-friendly setup to:

Run high-throughput binder design without local GPU constraints

Use RunPod or other cloud GPU services with persistent state

Share a ready-to-use Docker image for others to reproduce my workflows

This repo lets users skip complex local setup and run everything via a containerized workflow, while still preserving full flexibility of BindCraft’s powerful design stack.

📜 Acknowledgments
This project reuses and repackages tools developed by:
Martin Pacesa et al. — BindCraft
Sergey Ovchinnikov — ColabDesign
Justas Dauparas — ProteinMPNN
RosettaCommons — PyRosetta
BindCraft is licensed under the MIT License.
ColabDesign and ProteinMPNN are licensed under Apache 2.0.
Feel free to reach out via GitHub Issues if you encounter bugs or would like to contribute.