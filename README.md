# D<sup>3</sup>TaLES High-throughput Workflow Management
This repository contains [Fireworks](https://materialsproject.github.io/fireworks/)-based code for implementing 
high-throughput calculations for the [D<sup>3</sup>TaLES database](https://d3tales.as.uky.edu/database/). The
contents of this repository include

- `ex_config`: Directory containing example [Fireworks](https://materialsproject.github.io/fireworks/) 
configuration files. 
- `management`: Directory containing scripts for workflow management including the `populate_lpad` file 
for launching high-throughput quantum chemical workflows
- `misc_scripts`: Directory containing miscellaneous scripts for database management (not as useful outside 
of the D<sup>3</sup>TaLES project)
- `parameters`: Directory containing parameter files that set the parameters (e.g., functional, basis set, etc.)
for high-throughput calculations 
- `workflows`: Directory containing Fireworks workflows for high-throughput calculations including Firetask,
Firework, and workflow classes. 

Most of the most useful code in this repository will be found in the `workflows` directory. 


## Environment
### Install
The primary package this repo requires is the [D<sup>3</sup>TaLES API](https://github.com/D3TaLES/d3tales_api). 
It is recommended that you create an environment to host the required packages. 
```bash
conda create --name d3tales_fw --file d3tales_fw.yml
conda activate d3tales_fw
pip install git+https://github.com/d3tales/d3tales_api.git
```

### Activate
Note that you must set the `DB_INFO_FILE` environment variable as stipulated in the
[D<sup>3</sup>TaLES API Docs](https://github.com/D3TaLES/d3tales_api).You will also
need to set the fireworks variable `FW_CONFIG_FILE` (see the [Fireworks website](https://materialsproject.github.io/fireworks/) 
for more details) and a `PYTHONPATH`. For example, a user might set the environment variables as follows: 
```bash
conda activate d3tales_fw
export PYTHONPATH='C:/Users/Lab/d3tales_fw'
export FW_CONFIG_FILE='C:/Users/Lab/d3tales_fw/ex_config/FW_config.yaml'
export DB_INFO_FILE='C:/Users/Lab/db_infos.json'
```


## Running and Viewing Jobs

To view jobs run: 
```bash
lpad webgui
```

To launch a D3TaLES workflow or a series of D3TaLES workflows for given molecule(s), use the 
`management populate_lpad.py` script. 