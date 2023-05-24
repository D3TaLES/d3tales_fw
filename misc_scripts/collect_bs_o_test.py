import os
import re
import json
from pymatgen.io.gaussian import GaussianOutput


home = os.getcwd()
data_dict = {}
tests = [t for t in os.listdir(home) if t.endswith('sTZ') or t.endswith('sDZ')]

for test in tests:
    structure = test.split('_')[-1]
    omega = test.split('_')[-2]
    badO = True if test.split('_')[0] == 'bado' else False
    test_dir = os.path.join(home, test)
    try:
        log_file = [f for f in os.listdir(test_dir) if f.endswith('.log')][0]
    except IndexError:
        continue
    log_path = os.path.join(test_dir, log_file)
    mol = GaussianOutput(log_path)
    if mol.properly_terminated:
        energy = mol.final_energy
        basis_set = mol.basis_set + "_badO" if badO else mol.basis_set
        route_params = {k.lower(): v for k, v in mol.route_parameters.items()}
        omega_string = re.sub("[^0-9]", "", route_params['3108'])
        if not data_dict.get(basis_set, ):
            data_dict[basis_set] = {}
        if not data_dict[basis_set].get(omega):
            data_dict[basis_set][omega] = {}
        data_dict[basis_set][omega][structure] = energy

with open(os.path.join(home, 'bs_omega_test_data.json'), 'w') as fn:
    json.dump(data_dict, fn, indent=2)

