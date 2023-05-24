### generates the cube files for every calculation for every mol in the given directory
import os
import glob
import subprocess
from sys import argv

print_cmds = True
gaussian_cmd = "module load ccs/gaussian/g16-A.03/g16-haswell"


def homo_lumo_cubes(log_file, num_orbitals=2, print_cmds=False):
    from pymatgen.io.gaussian import GaussianOutput
    name = log_file[:-4]
    gout = GaussianOutput(log_file)
    n_electrons = gout.electrons[0]
    mos = list(gout.eigenvalues.values())[0]
    data_dict = {"name": name}
    for i in range(num_orbitals):
        orb = n_electrons - 1 - i
        lumo_orb = n_electrons - 1 - orb + n_electrons
        state = "homo_{}".format(orb - n_electrons + 1)
        lumo_state = "lumo_{}".format(lumo_orb - n_electrons)
        data_dict.update({state: mos[orb] * 27.2114, lumo_state: mos[lumo_orb] * 27.2114})
        if not os.path.isfile(name + "_" + state):
            command = "cubegen 32 MO={} {}.fchk {}.cube 80".format(orb + 1, name, name + "_" + state)
            if print_cmds:
                print(command)
            else:
                os.system(command)
        if not os.path.isfile(name + "_" + lumo_state):
            command = "cubegen 32 MO={} {}.fchk {}.cube 80".format(lumo_orb + 1, name, name + "_" + lumo_state)
            if print_cmds:
                print(command)
            else:
                os.system(command)
        orb -= 1
    return data_dict


def electro_potent_cubes(log_file, print_cmds=False):
    name = log_file[:-4]
    if not os.path.isfile(name + "{}_potential.cube 80"):
        command = "cubegen 0 potential {}.fchk {}_potential.cube 80".format(name, name)
        if print_cmds:
            print(command)
        else:
            os.system(command)
    if not os.path.isfile(name + "{}_density.cube 80"):
        command = "cubegen 0 density {}.fchk {}_density.cube 80".format(name, name)
        if print_cmds:
                print(command)
        else:
            os.system(command)


if __name__ == "__main__":
    try:
        arguments = argv[1:]
        kwargs = [k for k in arguments if not k.startswith('-')]
        tags = [k for k in arguments if k.startswith('-')]
    except IndexError:
        kwargs = []
        tags = []
    data_dir = kwargs[0] if kwargs else ''
    tags = tags or ['-hl', '-ep']

    subprocess.call(gaussian_cmd, shell=True)
    data = []
    data_home = os.path.join(os.getcwd(), data_dir)
    mols = [m for m in os.listdir(data_home) if os.path.isdir(os.path.join(data_home, m))]
    for mol_id in mols:
        mol_path = os.path.join(data_home, mol_id)
        if not os.path.isdir(mol_path):
            continue
        os.chdir(mol_path)
        for log_file in glob.glob("*opt.log"):
            if '-hl' in tags:
                d = homo_lumo_cubes(log_file, print_cmds=print_cmds)
                data.append(d)
            if '-ep' in tags:
                electro_potent_cubes(log_file, print_cmds=print_cmds)

    os.chdir(data_home)
    if '-hl' in tags:
        import pandas as pd
        pd.DataFrame(data).to_csv("mo.csv", index=False)
