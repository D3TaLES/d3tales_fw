import os
import shutil
import argparse
from d3tales_api.D3database.d3database import FrontDB


calc_loc_dict = {
    "tddft": "tddft/energy",
    "opt": "gas_phase/opt",
    "freq": "gas_phase/opt",
    "solv": "solv_acetonitrile/energy",
}

chg_type_dict = {
    "c1c1": "cation1",
    "gsgs": "groundState",
    "c2c2": "cation2",
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate a CSV file containing given properties (default: smiles) for all D3TaLES molecules.')
    parser.add_argument('-f', '--filename', type=str, help='filepath for a .txt file with paths to delete', default=False)
    parser.add_argument('-l', '--logfile', type=str, help='filepath for the logfile with paths to delete', default=False)
    parser.add_argument('-d', '--delete', action='store_true', help='delete file paths')
    parser.add_argument('-fd', '--filename_dirs', type=str, help='filepath for a .txt file with a list of directories to delete', default=False)
    parser.add_argument('-dd', '--delete_dirs', action='store_true', help='delete directory paths')
    parser.add_argument('-gd', '--get_dirs', action='store_true', help='get directories to delete')
    args = parser.parse_args()

    home = os.getcwd()
    if args.filename:
        with open(args.filename) as fn:
            delete_paths = fn.readlines()
        paths_to_delete = [path.replace('\n', '').strip(" ") for path in delete_paths if os.path.isfile(path.replace('\n', ''))]
    elif args.logfile:
        with open(args.logfile) as fn:
            delete_paths = fn.readlines()
        paths_to_delete = [path.split(",")[-1].replace('\n', '').strip(" ") for path in delete_paths]
    elif args.filename_dirs:
        with open(args.filename_dirs) as fn:
            delete_paths = fn.readlines()
        paths_to_delete = [path.replace('\n', '') for path in delete_paths]
    elif args.get_dirs or args.delete_dirs:
        frontend_db = FrontDB(collection_name="base")
        completed_cursor = frontend_db.coll.find({'$and': [
            {"mol_characterization.oxidation_potential": {'$exists': True}},
            {"species_characterization.cation1.spectra": {'$exists': True}},
            {"species_characterization.cation2.spectra": {'$exists': True}},
        ]
        })
        with open("paths_to_delete.txt", 'w') as fn:
            fn.writelines("%s\n" % mol.get("_id") for mol in completed_cursor)
        paths_to_delete = [mol.get("_id") for mol in completed_cursor]
        print(len(paths_to_delete), " mol directories gathered: ")

    else:
        submission_db = FrontDB(collection_name="submission")
        unprocessed_cursor = submission_db.coll.find({'$and': [
            {"approved": True},
            {"processed_data": {'$ne': {}}}
        ]
        })
        path_list = []
        for mol in unprocessed_cursor:
            try:
                file = mol.get('file')
                mol_id = mol.get('molecule_id')
                data_type = mol.get('data_type')
                calc_type = file.split("_")[1]
                calc_location = calc_loc_dict[calc_type]
                files_to_delete = '_'.join(file.split("_")[1:]).split('.')[0] + '*'
                if calc_type =="solv":
                    chg_type = file.split("_")[3].split('.')[0]
                    charge_state = chg_type_dict[chg_type]
                    zip_path1 = "{}/{}/{}/{}/{}".format(mol_id, data_type, calc_location, charge_state, file)
                    zip_path2 = "{}/{}/{}/{}/{}".format(mol_id, data_type, calc_location, charge_state, files_to_delete)
                else:
                    zip_path1 = "{}/{}/{}/{}".format(mol_id, data_type, calc_location, file)
                    zip_path2 = "{}/{}/{}/{}".format(mol_id, data_type, calc_location, files_to_delete)
            except:
                continue
            try:
                if args.delete:
                    os.remove(zip_path1)
                    print("Deleted: ", zip_path1)
                    os.remove(zip_path2)
                    print("Deleted: ", zip_path2)
                path_list.extend([zip_path1, zip_path2])
            except FileNotFoundError:
                print("File Not Found.")
        with open("paths_to_delete.txt", 'w') as fn:
            fn.writelines("%s\n" % path for path in path_list)
        paths_to_delete = []

    count = 0
    if args.delete_dirs:
        existing_dirs = os.listdir(home)
        mols_to_delete = [mol for mol in paths_to_delete if mol in existing_dirs]
        print("Deleting {} directories...".format(len(mols_to_delete)))
        for mol_id in mols_to_delete:
            try:
                shutil.rmtree(os.path.join(home, mol_id))
                count += 1
            except OSError as e:
                print("Error: %s : %s" % (os.path.join(home, mol_id), e.strerror))
        print("{} directories deleted.".format(count))
    if args.delete:
        paths_f = [p for p in paths_to_delete if os.path.isfile(p)]
        print("Deleting {} files...".format(len(paths_f)))
        for path in paths_f:
            try:
                os.remove(path)
                count += 1
            except OSError as e:
                print("Error: %s : %s" % (os.path.join(home, path), e.strerror))
        paths_d = [p for p in paths_to_delete if os.path.isdir(p)]
        print("Deleting {} directories...".format(len(paths_d)))
        for path in paths_d:
            try:
                shutil.rmtree(path)
                count += 1
            except OSError as e:
                print("Error: %s : %s" % (os.path.join(home, path), e.strerror))
        print("{} files deleted.".format(count))
