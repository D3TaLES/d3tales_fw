import os
from d3tales_api.D3database.restapi import RESTAPI

data_home = os.path.join(os.getcwd(), 'd3tales')
mols = [m for m in os.listdir(data_home) if m.startswith('05')]
for mol_id in mols:
    mol_path = os.path.join(data_home, mol_id, 'computation', 'gaussian')
    calcs = []
    for zip_file in os.listdir(mol_path):
        zip_path = os.path.join(mol_path, zip_file)
        calc_type = "_".join(zip_file.split("_")[1:-1])
        if calc_type in calcs:
            continue
        RESTAPI(method='post', endpoint='tools/upload/computation-gaussian',
                url="https://d3tales.as.uky.edu", expected_endpoint="tools/user_uploads",
                upload_file=zip_path, params=dict(molecule_id=mol_id, calculation_type=calc_type))
        calcs.append(calc_type)
    break

# import os
# from d3tales_api.D3database.restapi import RESTAPI
#
# data_home = os.path.join(os.getcwd())
# for zip_file in os.listdir(data_home):
#     zip_path = os.path.join(data_home, zip_file)
#     mol_id = zip_file.split("_")[0]
#     calc_type = "_".join(zip_file.split("_")[1:]).split('.')[0]
#     print(mol_id, calc_type)
#     RESTAPI(method='post', endpoint='tools/upload/computation-gaussian',
#             url="https://d3tales.as.uky.edu", expected_endpoint="tools/user_uploads",
#             upload_file=zip_path, params=dict(molecule_id=mol_id, calculation_type=calc_type))
