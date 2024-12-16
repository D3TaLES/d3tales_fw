import json
import argparse
import functools
import numpy as np
import pandas as pd
from d3tales_api.D3database.d3database import FrontDB


LIMIT = 0
SOLVENT = {"name": "Acetonitrile", "model": "implicit_solvent", "dielectric_constant": 35.688}
DISPLAY_CONDITIONS = {"data_source": ['dft'], 'functional': ['LC-wHPBE', 'ulc-whpbe'], 'basis_set': ['Def2SVP']}

mol_i_props = ["smiles", "source_group", "groundState_charge", "number_of_atoms", "molecular_weight", "sa_score"]
mol_props = ["hole_reorganization_energy", "electron_reorganization_energy", "relaxation_groundState_cation1",
             "relaxation_cation1_groundState", "relaxation_groundState_anion1", "relaxation_anion1_groundState",
             "vertical_ionization_energy", "vertical_electron_affinity", "adiabatic_ionization_energy",
             "adiabatic_electron_affinity", "oxidation_potential", "reduction_potential", "rmsd_groundState_cation1",
             "rmsd_cation1_cation2", "rmsd_groundState_anion1", "rmsd_anion1_anion2", "omega"]
mol_solv_props = ["oxidation_potential", "reduction_potential"]
species_props = ["charge", "globular_volume", "radical_stability_score", "homo_lumo_gap", "dipole_moment",
                 "solvation_energy", "homo", "lumo", "singlet_states", "triplet_states"]
species_solv_props = ["globular_volume", "radical_stability_score", "homo_lumo_gap", "dipole_moment",
                      "solvation_energy", "homo", "lumo"]
species = ["groundState", "cation1", "cation2", "anion1", "anion2"]

MOL_INFO = ",".join([f"mol_info.{p}" for p in mol_i_props])
MOL_CHAR = ",".join([f"mol_characterization.{p}" for p in mol_props])
SPECIES_CHAR = ",".join([f"species_characterization.{s}.{p}" for p in species_props for s in species])
SPECIES_SOLV_CHAR = ",".join([f"{s}.{p}" for p in species_solv_props for s in species] + mol_solv_props)
ALL_PROPERTIES = ",".join(["public", MOL_INFO, MOL_CHAR, SPECIES_CHAR])


def rgetkeys(_dict, keys, **kwargs):
    def _getkey(_dict, key):
        _dict = _dict or {}
        if isinstance(_dict, dict):
            return _dict.get(key, **kwargs)
        if isinstance(_dict, list) and key.isdigit():
            return _dict[int(key)]

    return functools.reduce(_getkey, [_dict] + keys.split('.'))


def get_value(x, display_conditions: dict = DISPLAY_CONDITIONS, solv=None):
    display_conditions['solvent'] = [solv]
    if isinstance(x, list):
        for prop_item in x:
            conditions = prop_item.get('conditions')
            display = True
            for display_cond, display_values in display_conditions.items():
                if conditions.get(display_cond) not in display_values:
                    display = False
                    break
            if display:
                return prop_item.get("value", prop_item.get("excitations", [[None]])[0][0])
        return None
    if isinstance(x, dict):
        return x.get("value", x.get("excitations", [[None]])[0][0])
    return x

def generate_csv(collect_properties, out_file="database.csv", ids_filename=None, public=True, standard_deviations=True,
                 std_props=mol_props):
    # Database search
    frontend_db = FrontDB(collection_name="base")
    projection = {p: 1 for p in collect_properties}
    if ids_filename:
        with open(ids_filename) as f:
            ids_data = json.load(f)
        ids = list(ids_data.keys())
        cursor = frontend_db.coll.find({'_id': {"$in": ids}}, {projection: 1}).limit(LIMIT)
    else:
        cursor = frontend_db.coll.find({}, projection).limit(LIMIT)

    # Data cleaning with pandas
    master_data = pd.DataFrame.from_records(cursor)
    master_data.set_index('_id', inplace=True)
    print("PULLED DATA SHAPE: ", master_data.shape)
    prop_paths = {
        p.replace('mol_info.', "").replace('species_characterization.', "").replace('mol_characterization.', ""): p for
        p in collect_properties}
    columns = []
    for prop_name, prop_path in prop_paths.items():
        prop_df = pd.DataFrame()
        prop_df[prop_name] = master_data.apply(lambda x: rgetkeys(x.to_dict(), prop_path), axis=1)
        prop_df[prop_name] = prop_df[prop_name].apply(lambda x: get_value(x, solv=None))
        if prop_name in SPECIES_SOLV_CHAR:
            prop_df["solv_" + prop_name] = master_data.apply(lambda x: rgetkeys(x.to_dict(), prop_path), axis=1)
            prop_df["solv_" + prop_name] = prop_df["solv_" + prop_name].apply(lambda x: get_value(x, solv=SOLVENT))
        columns.append(prop_df)
    final_data = pd.concat(columns, axis=1)

    final_data.dropna(axis=1, how="all", inplace=True)
    if public:
        final_data.public.fillna(False, inplace=True)
        final_data = final_data[final_data.public]
    final_data.to_csv(out_file.replace(".csv", "_raw.csv"))
    if standard_deviations:
        numeric_data = final_data[[c for c in final_data.select_dtypes(include=[np.number]).columns if c in std_props]]
        normal_rows = ((numeric_data - numeric_data.mean()) / numeric_data.std(ddof=0)).fillna(0).abs().le(standard_deviations).all(axis=1)
        final_data = final_data[normal_rows]
    final_data.to_csv(out_file)
    print("FINAL DATA SHAPE: ", master_data.shape)
    return final_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate a CSV file containing given properties (default: smiles) for all D3TaLES molecules.')
    parser.add_argument('-f', '--filename', type=str, help='filepath for a CSV file to which to save the data',
                        default='structure_data/d3tales.csv')
    
    parser.add_argument('-p', '--collect_properties', type=str, help='properties to include in CSV',
                        default=ALL_PROPERTIES)
    parser.add_argument('-ids', '--ids_filename', type=str,
                        help='filepath for JSON file containing ids (as keys) to search')
    parser.add_argument('-std', '--standard_deviations', type=int, default=3,
                        help='the number of standard deviations to keep. If 0, not std. dev. filter will be applied')
    parser.add_argument('-pub', '--public', action='store_true', help='save only public data to csv')
    args = parser.parse_args()

    # Generate query CSV
    coll_properties = args.collect_properties.split(',')
    final_data_df = generate_csv(coll_properties, out_file=args.filename, standard_deviations=args.standard_deviations,
                                 ids_filename=args.ids_filename, std_props=mol_props[:-1], public=args.public)
    print("Success! D3TaLES molecules and {} property were saved to {}".format(",".join(list(final_data_df.columns)), args.filename))
