import os
import subprocess
import time
from rdkit import Chem
from rdkit.Chem import AllChem
class lig:
    def __init__(self, smiles, regular_name, molecule, charge, dir):
        self.PDBMAKER(regular_name,smiles,molecule)
        self.dir = dir
        self.smiles = smiles
        self.mol = molecule
        self.charge=int(charge)
        subprocess.run([f'mkdir {self.dir}/{self.mol}script'], shell=True)
        conda_activate = f"source /project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/miniconda3/bin/activate && conda activate ligpg"
        export_bossdir = f" export BOSSdir=/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/boss"
        #igpargen_cmd = f"ligpargen -s '{self.smiles}' -n {self.mol} -p /mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs/{self.mol} -r {self.mol} -c {self.charge} -o 0 -cgen CM1A"
        ligpargen_cmd = f"ligpargen -i /project/cmri235_uksr/shasanka_conda_boss/launch/{regular_name}/gaussian/gas_phase/opt/{self.mol}.pdb -n {self.mol} -p /mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs/{self.mol} -r {regular_name[:3]} -c 0 -o 0 -cgen CM1A"

        singularity_container = f"/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/Fast/f.sif"

        cmd = ["singularity", "exec", singularity_container, "bash", "-c",
               f'{conda_activate} && {export_bossdir} && {ligpargen_cmd}']
        subprocess.Popen(cmd).wait()

        while os.path.isfile(f'{self.dir}/{self.mol}/{self.mol}-debug.pdb') == False:
            print("still making pdb")
            time.sleep(1)
        pdb_location = Chem.MolFromPDBFile(f'{self.dir}/{self.mol}/{self.mol}-debug.pdb')
        pdb_location = Chem.AddHs(pdb_location)
        AllChem.EmbedMolecule(pdb_location)
        Chem.MolToMolFile(pdb_location, f'{self.dir}/{self.mol}/{self.mol}.mol')
        comand2 = f'obabel -imol {self.dir}/{self.mol}/{self.mol}.mol -ogjf -O {self.dir}/{self.mol}script/{self.mol}.gjf --gen3D'
        comand3 = f'obabel -igro {self.dir}/{self.mol}/{self.mol}.gmx.gro -opdb -O {self.dir}/{self.mol}/{self.mol}.pdb '
        print(comand3)
        subprocess.run([comand2], shell=True)
        subprocess.run([comand3], shell=True)

    def PDBMAKER(self, name, smiles, type):
        subprocess.run([
                           f"obabel /project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/freq_opt_groundState.log -opdb -O /project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/out.pdb"],
                       shell=True)
        subprocess.run(
            [f"touch /project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/out.pdb"],
            shell=True)
        with open(
                f"/project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/out.pdb") as file, open(
                f"/project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/{type}.pdb",
                'a') as new:
            f = file.readlines()
            line_to_print = []
            for iteams in f:
                line_to_print.append(iteams.strip("\n").split(' '))
            for lines in line_to_print:
                if 'UNL' in lines:
                    lines[lines.index('UNL')] = "MET"
            for lines in line_to_print:
                for line in lines:
                    if line == "":
                        new.write(f'')
                    new.write(f'{line} ')
                new.write("\n")