import os
import json
import argparse
from d3tales_api.D3database.d3database import FrontDB
from rdkit.Chem import MolFromSmiles, MolToSmiles

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get D3TaLES IDs from SMILES JSON file.')
    parser.add_argument('-f', '--filename', type=str, help='filepath for a JSON file from which to get SMILES', default='')
    parser.add_argument('-i', '--ids_list', type=str, help='list of ids from which to get smiles', default='')
    args = parser.parse_args()

    if os.path.isfile(args.filename):
        with open(args.filename, 'r') as f:
            smiles_dict = json.load(f)
        origin = smiles_dict.pop("origin") if smiles_dict.get("origin") else 'Risko'

        id_data = {}
        for name, smiles in smiles_dict.items():
            try:
                rdkmol = MolFromSmiles(smiles)
                clean_smiles = MolToSmiles(rdkmol)
            except:
                continue
            _id = FrontDB(smiles=clean_smiles, group=origin).generate_id()
            id_data[_id] = name
        print(id_data)
        print(','.join(id_data.keys()))

    if args.ids_list:
        return_dict = {}
        for i in args.ids_list.split(','):
            query = FrontDB().make_query({"_id": i}, {"mol_info.smiles": 1})
            return_dict[i] = query.mol_info[0].get("smiles")
        with open("structure_data/ids_smiles.json", 'w') as f:
            json.dump(return_dict, f, indent=2)





