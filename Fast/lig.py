import os
import subprocess
import time
from rdkit import Chem
from rdkit.Chem import AllChem
class lig:
    def __init__(self, sm, m,charge,dir):
        self.dir = dir
        self.smiles = sm
        self.mol = m
        self.charge=int(charge)
        subprocess.run([f'mkdir {self.dir}/{self.mol}script'], shell=True)
        conda_activate = f"source /project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/miniconda3/bin/activate && conda activate ligpg"
        export_bossdir = f" export BOSSdir=/mnt/gpfs2_4m/scratch/sla296/test_run/py/d3tales_fw/boss"
        ligpargen_cmd = f"ligpargen -s '{self.smiles}' -n {self.mol} -p /mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs/{self.mol} -r {self.mol} -c {self.charge} -o 0 -cgen CM1A"
        singularity_container = f"/project/cmri235_uksr/shasanka_conda_boss/sla296/Desktop/scratch/test_run/py/d3tales_fw/Fast/f.sif"

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
