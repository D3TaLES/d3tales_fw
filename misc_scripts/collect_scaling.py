import os
import re
import json


def get_runtime(log_path):
        """
        Collects runtime in core hours from a logfile
        """
        time_patt = re.compile(r"\d+\.d+|\d+")
        time_data = []
        with open(log_path, "r") as f:
            line = f.readline()
            while line != "":
                if re.match("Job cpu time", line.strip()):
                    time_data.extend(time_patt.findall(line))
                line = f.readline()
        if time_data:
            time_data = [float(time) for time in time_data]
            runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
        else:
            runtime = 0
        return round(runtime, 3)

home = os.getcwd()
data_dict = {}
mols = [m for m in os.listdir(home) if m.startswith('mol')]
for mol in mols:
    mol_path = os.path.join(home, mol)
    for cpus in os.listdir(mol_path):
        cpus_path = os.path.join(mol_path, cpus)
        log_files = [l for l in os.listdir(cpus_path) if l.endswith('.log')]
        for log in log_files:
            calc_type = log.split('.')[0]
            log_path = os.path.join(cpus_path, log)
            runtime = get_runtime(log_path)
            if not data_dict.get(calc_type, ):
                data_dict[calc_type] = {}
            if not data_dict[calc_type].get(mol):
                data_dict[calc_type][mol] = {}
            data_dict[calc_type][mol][cpus] = runtime

with open(os.path.join(home, 'scaling_data.json'), 'w') as fn:
    json.dump(data_dict, fn, indent=2)
