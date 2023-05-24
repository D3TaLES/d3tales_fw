import os
import argparse
from d3tales_api.D3database.restapi import RESTAPI
from d3tales_api.D3database.d3database import FrontDB

calc_loc_dict = {
    "tddft": "tddft/energy",
    "opt": "gas_phase/opt",
    "freq": "gas_phase/opt",
    "solv": "solv_acetonitrile/energy",
    "wtuning": "wtuning",
}

chg_type_dict = {
    "c1c1": "cation1",
    "gsgs": "groundState",
    "c2c2": "cation2",
}


def resubmit(zip_path, runfile_path="/scratch/rdu230/d3tales/high_throughput/runfiles"):
    mol_id = zip_path.split('/')[0]
    calc_type = zip_path.split('/')[-1].strip(mol_id + '_').strip('.zip').replace("hf_", "").replace("freq_opt", "freq")
    upload_path = os.path.join(runfile_path, zip_path)
    RESTAPI(method='post', endpoint='tools/upload/computation-gaussian',
            url="https://d3tales.as.uky.edu", expected_endpoint="tools/user_uploads",
            upload_file=upload_path, params=dict(molecule_id=mol_id, calculation_type=calc_type))
    print("File resubmitted: {} calculation for molecule {}.".format(calc_type, mol_id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Resubmit zip files that were not successfully processed and approved.')
    parser.add_argument('-f', '--filename', type=str, help='Filepath to a .txt file containing filepaths to resubmit')
    parser.add_argument('-z', '--zip_file', type=str, help='Filepath to zipfile to resubmit')
    parser.add_argument('-w', '--write', action='store_true', help='write a textfile for all failed resubmissions')
    parser.add_argument('-hf', '--hf', action='store_true', help='resubmit all HF files')
    parser.add_argument('-solv', '--solv', action='store_true', help='resubmit all "solv" files')
    args = parser.parse_args()

    if args.write:
        frontend_db = FrontDB(collection_name="submission")
        if args.hf:
            unprocessed_cursor = frontend_db.coll.find({"data_type": "gaussian", "metadata.calculation_type": {"$regex": "solv_.*"}})
        elif args.solv:
            unprocessed_cursor = frontend_db.coll.find({"data_type": "gaussian", "file": {"$regex": ".*_hf_.*"}})
        else:
            unprocessed_cursor = frontend_db.coll.find({"data_type": "gaussian", "approved": False, "processed_data": {}})
        path_list = []
        for mol in unprocessed_cursor:
            try:
                file_name = mol.get('file')
                mol_id = mol.get('molecule_id')
                data_type = mol.get('data_type')
                calc_type = file_name.replace("_hf_", "_").split("_")[1].strip(".zip")
                calc_location = calc_loc_dict[calc_type]
                if calc_type == "solv":
                    chg_type = file_name.replace("_hf_", "_").split("_")[3].split('.')[0]
                    charge_state = chg_type_dict[chg_type]
                    zip_path = "{}/{}/{}/{}/{}".format(mol_id, data_type, calc_location, charge_state, file_name)
                else:
                    zip_path = "{}/{}/{}/{}".format(mol_id, data_type, calc_location, file_name)
                path_list.append(zip_path)
            except:
                pass
        with open("paths_to_resubmit.txt", 'w') as fn:
            fn.writelines("%s\n" % path for path in path_list)
    elif args.filename:
        with open(args.filename) as fn:
            path_list = [path.rstrip() for path in fn.readlines()]
            for path in path_list:
                try:
                    resubmit(path, runfile_path=os.getcwd())
                except FileNotFoundError:
                    print(path, " not found")
    elif args.zip_file:
        resubmit(args.zip_file)
