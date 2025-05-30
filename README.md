# CLDeVa
Pipeline for cluster detection and validation in LSST-DESC.

# Conda Environments Setup

This repository includes pre-configured conda environment files located in the `conda_envs/` directory.
Follow the instructions below to set up an environment.

## Prerequisites

- If at CC-IN2P3 or NERSC, load conda via `module load conda`.
- Else, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/).

## Setting Up the Conda Environments

1. **Navigate to the repository directory:**

   ```bash
   cd path/to/this/repo
   ```
   
2. **Create the environment from the `.yaml` file:**
   ## Option 1: Create by custom path (using --prefix)
   This method creates the environment in a specific directory without adding it to the global environment list.
   This is the advised method as it allows you to build the environment in a place of your choosing (i.e. a place with more memory allocation).
   ```bash
   conda env create -f conda_envs/<environment_name>.yaml --prefix <custom_env_path>
   conda activate <custom_env_path>
   ```
   Replace <environment_name> with the YAML filename, and <custom_env_path> with the name or path of the folder where you'd like the environment to be created.

   ## Option 2: Create by environment name
   This method registers the environment globally under the name defined in the YAML file.
   This location may not always be optimal for memory usage.
   ```bash
   conda env create -f conda_envs/<environment_name>.yaml
   conda activate <environment_name>
   ```
   Replace `<environment_name>` with the name of the YAML file (without the .yaml extension).

