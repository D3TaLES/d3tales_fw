import os
import argparse

PREFIX = ''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clean storge for D3TaLES')
    parser.add_argument('-cc', '--check_calc_type', type=str, help='check if a given calculation type is in storage', default=None)
    parser.add_argument('-d', '--delete', action='store_true', help='delete files with too long names')
    parser.add_argument('-dd', '--delete_duplicates', action='store_true', help='delete oldest files for duplicate calculations')
    parser.add_argument('-sd', '--storage_directory', type=str, help='storage directory to clean', default='d3tales')
    args = parser.parse_args()

    data_home = os.path.join(os.getcwd(), args.storage_directory)
    mols = [m for m in os.listdir(data_home) if m.startswith(PREFIX)]
    for mol_id in mols:
        calcs = []
        mol_path = os.path.join(data_home, mol_id, 'computation', 'gaussian')
        if not os.path.isdir(mol_path):
            continue
        zip_files = sorted(os.listdir(mol_path), key=lambda x: os.path.getmtime(os.path.join(mol_path, x)), reverse=True)
        for zip_file in zip_files:
            calc_type = "_".join(zip_file.split("_")[1:-1])
            zip_path = os.path.join(mol_path, zip_file)
            # Remove names tags
            name_tags = ['hf_', 'nmr_', 'hfgeom_', 'dftgeom_']
            zip_name = zip_file
            for t in name_tags:
                zip_name = zip_name.replace(t, '')
            name_parts = zip_name.split("_")

            try:
                if name_parts[1] == 'solv' or name_parts[1] == 'freq':
                    if len(name_parts) > 5:
                        print(zip_path)
                        if args.delete:
                            os.remove(zip_path)
                    else:
                        pass
                elif name_parts[1] == 'wtuning' and len(name_parts) == 4:
                    print(zip_path)
                    if args.delete:
                        os.remove(zip_path)
                elif len(name_parts) > 4:
                    print(zip_path)
                    if args.delete:
                        os.remove(zip_path)
                if calc_type in calcs:
                    print(zip_path)
                    if args.delete_duplicates:
                        os.remove(zip_path)
                calcs.append(calc_type)

            except:
                continue
        if args.check_calc_type:
            if args.check_calc_type not in calcs:
                print("{}/gaussian/gas_phase/opt/{}_{}.zip".format(mol_id, mol_id, args.check_calc_type))
    print("Cleaned {} molecules.".format(len(mols)))
