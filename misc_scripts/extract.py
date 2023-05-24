import os
from zipfile import ZipFile

names_dict = {'05QDYA': 'mol_NN_cp_frag', '05DIRJ': 'mol_ta_cp_whol', '05RPNH': 'mol_bp_ma_whol', '05MXLE': 'mol_ta_cp_half', '05AWCY': 'mol_ta_pa_whol', '05OKDT': 'mol_ta_pa_half', '05JZJR': 'mol_bp_cp_whol', '05POEN': 'mol_bp_ma_half', '06ESQE': 'mol_bp_NN_frag', '05CKZO': 'mol_ta_ma_half', '05XJXO': 'mol_bp_pa_half', '05WOBI': 'mol_bp_cp_half', '05SMPN': 'mol_ta_ma_whol', '05JNUA': 'mol_bp_pa_whol', '05MNXL': 'mol_ta_NN_frag', '05DUQU': 'mol_NN_pa_frag', '05VEUX': 'mol_NN_ma_frag', '90PSTI': 'mol_ph_ma_whol', '90TRQT': 'mol_ph_ma_half', '90LPMI': 'mol_ph_fa_whol', '90QDBN': 'mol_bp_fa_whol', '90WROC': 'mol_bp_fa_half', '90MOFW': 'mol_bp_fm_whol', '90JCLK': 'mol_bp_fm_half', '90PNAB': 'mol_ph_fa_half', '90NCLQ': 'mol_ph_fm_whol', '90MOMH': 'mol_ph_fm_half', '05SADF': 'mol_ph_ef_whol', '05XQSP': 'mol_bpo_fm_whol', '05TWGB': 'mol_bpo_fm_half', '05WKJG': 'mol_NN_fm_frag', '05VXVT': 'mol_ph_pa_whol'}

home = os.getcwd()
id_path = os.path.join(home, 'ids')
name_path = os.path.join(home, 'names')
for mol in [m for m in os.listdir(id_path) if os.path.isdir(os.path.join(id_path, m))]:
    mol_path = os.path.join(id_path, mol)
    name = names_dict.get(mol)
    if not name:
        continue
    print("Extracting data for: ", mol, name)
    mol_name_path = os.path.join(name_path, name)
    if not os.path.isdir(mol_name_path): os.mkdir(mol_name_path)
    gas_opt_dir = os.path.join(mol_path, 'computation', 'gaussian')
    if os.path.isdir(gas_opt_dir):
        for zip_file in [f for f in os.listdir(gas_opt_dir) if f.endswith('.zip')]:
            zip_path = os.path.join(gas_opt_dir, zip_file)
            with ZipFile(zip_path, "r") as target_zip:
                target_zip.extractall(mol_name_path)
