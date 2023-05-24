import os
import json
import glob
import argparse
from d3tales_api.D3database.restapi import RESTAPI


def resubmit(zip_path):
    # 06OPWX_hf_opt_cation2_R17.zip
    file_name = zip_path.split('/')[-1]
    mol_id = file_name.split('_')[0]
    calc_type = "_".join(file_name.split('_')[:-1]).strip(mol_id + '_').replace("hf_", "").replace("freq_opt", "freq")
    RESTAPI(method='post', endpoint='tools/upload/computation-gaussian',
            url="https://d3tales.as.uky.edu", expected_endpoint="tools/user_uploads",
            upload_file=zip_path, params=dict(molecule_id=mol_id, calculation_type=calc_type))
    print("File resubmitted: {} calculation for molecule {}.".format(calc_type, mol_id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Resubmit files from storage from processing.')
    parser.add_argument('-f', '--filename', type=str, help='Filepath to a JSON file containing IDs to resubmit')
    parser.add_argument('-s', '--single_file', type=str, help='Specific filepath to resubmit')
    args = parser.parse_args()

    if args.single_file:
        resubmit(args.single_file)

    home = os.getcwd()
    if args.filename:
        with open(args.filename) as fn:
            data = json.load(fn)
        ids = [i for i in data.keys()]
        for i in ids:
            try:
                files_path = os.path.join(home, 'd3tales', i, "computation", "gaussian")
                for z in glob.glob(files_path + "/*.zip"):
                    zip_path = os.path.join(files_path, z)
                    resubmit(zip_path)
                print("Files for {} resubmitted.".format(i))
            except:
                pass

