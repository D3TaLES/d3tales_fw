import os
import subprocess
import time
import rdkit as rd

class lig:
    def __init__(self, smiles, regular_name, molecule, charge, dir, own, own_path):
        self.dir = dir
        self.smiles = smiles
        self.mol = molecule
        self.charge = int(charge) or 0
        path_to_pdb= os.path.join(self.dir, f"{self.mol}.pdb")
        print(path_to_pdb)
        path_to_output=os.path.join(self.dir,self.mol)
        if own == True:
            self.own(smiles, regular_name, molecule, charge, dir, own, own_path)
        else:
            tr_mol = rd.Chem.MolFromSmiles(smiles)
            b = rd.Chem.AddHs(tr_mol)
            rd.Chem.AllChem.EmbedMolecule(b)
            rd.Chem.MolToPDBFile(b,path_to_pdb)
            conda_activate = f"source /project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/miniconda3/bin/activate && conda activate ligpg"
            export_bossdir = f" export BOSSdir=/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/boss"
            ligpargen_cmd = f"ligpargen -i {path_to_pdb} -n {self.mol} -p {path_to_output} -r {regular_name[:3]} -c {self.charge} -o 0 -cgen CM1A"
            singularity_container = f"/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/Fast/f.sif"
            moving_pdb_command= f"mv {path_to_pdb} {os.path.join(self.dir, self.mol )}"

            cmd = ["singularity", "exec", singularity_container, "bash", "-c",
                   f'{conda_activate} && {export_bossdir} && {ligpargen_cmd} && {moving_pdb_command}']
            try:
                print("in try")
                print(cmd)
                subprocess.Popen(cmd).wait()
                while os.path.isfile(f'{self.dir}/{self.mol}/{self.mol}-debug.pdb') == False:
                    print("still making pdb")
                    time.sleep(1)
            except:
                print(f'Ligpargen was not able to find a parameter, user input files is being used. Please rerun with your own itp and pdb files. This is passed for smiles, regular_name, molecule, charge, dir {(smiles, regular_name, molecule, charge, dir)}')


    def PDBMAKER(self, name, smiles, type):
        if name == "MET":
           name = "MeT"
        

        with open(
                f"/project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/out.pdb") as fie, open(
                f"/project/cmri235_uksr/shasanka_conda_boss/launch/{name}/gaussian/gas_phase/opt/{type}.pdb",
                'a') as new:
            f = fie.readlines()
            line_to_print = []
            for iteams in f:
                line_to_print.append(iteams.strip("\n").split(' '))
            for lines in line_to_print:
                if 'UNL' in lines:
                    lines[lines.index('UNL')] = f"{name[:3]}"
            for lines in line_to_print:
                for line in lines:
                    if line == "":
                        new.write(f'')
                    new.write(f'{line} ')
                new.write("\n")
    def own(self,smiles, regular_name, molecule, charge, dir, own, own_path):
        subprocess.run([f'mkdir {self.dir}/{self.mol}script'], shell=True)
        subprocess.run([f'mkdir {self.dir}/{self.mol}'], shell=True)
        path1=os.path.join(own_path, f"{regular_name}.pdb")
        path2=os.path.join(own_path, f"{regular_name}.itp")
        cmd1= f"cp {path1} {self.dir}/{self.mol}/{molecule}.pdb"
        cmd2 = f"cp {path2} {self.dir}/{self.mol}/{molecule}.gmx.itp"

        subprocess.run([cmd1], shell=True)
        subprocess.run([cmd2], shell=True)



