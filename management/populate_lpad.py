import argparse
from pathlib import Path
from fireworks import LaunchPad
from d3tales_fw.workflows.wf_writer import *
from d3tales_fw.workflows.ParamSet import GausParamSet

BASE_DIR = Path(__file__).resolve().parent.parent


def populate_d3tales_lpad(
    filename=None,
    workflow_function=d3tales_wf,
    public=False,
    priority=4,
    check_if_already_run=True,  # norm: True
    solvent='acetonitrile',  # norm: acetonitrile  ( DiMethylSulfoxide )
    param_tag='gaus_',  # norm: gaus_  ( huckaba_ dihedral_ gaus_singlets_)
    hf_mol_opt=False,  # norm: False
    restricted=True,  # norm: True
    use_iop=True,  # norm: True
    wtune=True,  # norm: True
    run_nto=False,  # norm: False

    default_origin='Zinc',
    set_smiles_list=[],  # norm: []  (must be empty for json to be read)
    set_smiles_json={},  # norm: {}
    json_of_ids=False,  # norm: False
    get_name=True,  # norm: True
    name_string='',  # norm: ''
    submit=True,  # norm: True
    lpad_file=None,  # norm: None

    # gaussian gaussian_special
    email=None,  # "rdu230@uky.edu"
    username=None,  # "Rebekah Duke"
    wf_tag=""  # ""
):
    """
    Function to initiate a series of D3TaLES workflows for a given set of molecules.
    """
    gaussian_file_name = 'gaussian_ur' if not restricted else 'gaussian_noTune' if not use_iop else 'gaussian_hf' if hf_mol_opt else 'gaussian'
    lpad_file = os.path.join(BASE_DIR.parent, 'config', 'my_launchpad.yaml')
    param_file = os.path.join(BASE_DIR, 'parameters', param_tag + 'parameter_file.json')

    if filename and os.path.isfile(filename):
        with open(args.filename, 'r') as f:
            smiles_dict = json.load(f)
    else:
        smiles_dict = set_smiles_json

    origin = smiles_dict.pop("origin") if smiles_dict.get("origin") else default_origin
    paramset = GausParamSet().from_json(param_file)
    smiles_list = [s for s in smiles_dict.values()]
    names_dict = {smiles: name for name, smiles in smiles_dict.items()}

    if set_smiles_list:
        smiles_list = set_smiles_list

    fw_id_list = []
    kwarg_dict = dict(paramset=paramset, public=public, hf_mol_opt=hf_mol_opt, submit=submit, use_iop=use_iop,
                      check_if_already_run=check_if_already_run, solvent=solvent, priority=priority,
                      gaussian_file_name=gaussian_file_name, email=email, username=username, restricted=restricted,
                      wtune=wtune, run_nto=run_nto, origin_group=origin, wf_tag=wf_tag)
    for smiles in smiles_list:
        if set_smiles_list or json_of_ids:
            kwarg_dict.update(identifier=smiles)
        else:
            name = names_dict.get(smiles) if get_name else None
            if name_string not in name:
                continue
            kwarg_dict.update(smiles=smiles, mol_name=name)
        wf = workflow_function(**kwarg_dict)
        lpad = LaunchPad().from_file(lpad_file)
        info = lpad.add_wf(wf)
        fw_id = list(info.values())[0]
        fw_id_list.append(fw_id)

    print("Added {} workflows to the launchpad".format(len(smiles_list)))
    print("FW ID List: ", ''.join(str(fw_id_list).split(' ')))
    return fw_id_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Launch D3TaLES FireWorks workflow.')
    parser.add_argument('filename', metavar='filename', type=str, help='filepath for a JSON molecule file',
                        default=False)
    args = parser.parse_args()

    populate_d3tales_lpad(args.filename)

