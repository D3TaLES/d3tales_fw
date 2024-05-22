import subprocess
import d3tales_fw.Fast.reorg as g
import d3tales_fw.Fast.topology as top


class gro:
    def __init__(self, solvent, solute, solvent2, di, x, y, z, key, initla_system,path_to_file):
        self.dir = di
        self.x = x
        self.y = y
        self.z = z
        self.fullname_solvent=solvent
        self.solvent = solvent[:3]
        self.key=key
        self.command12 = f'cp -r /project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/Fast/MDP {self.dir}/InputGrofiles{key}'

        self.solute = solute
        self.solvent2 = solvent2
        self.initla_system = initla_system
        self.path_to_file = path_to_file
        if self.initla_system:
            self.inital()

        else:
            print(f"in the gro solutes: {self.solute}")
            command = f'gmx_mpi editconf -f {self.dir}/Packmol{key}/Mixture.pdb -box {self.x} {self.y} {self.z} -o {self.dir}/Packmol{key}/solvated.gro'
            subprocess.run(command, shell=True, check=True)
            print("gro file made")

            command4 = f'mv {self.dir}/{solvent}_Solvent/{solvent}_Solvent.pdb {self.dir}/InputGrofiles{key} && mv {self.dir}/{solvent}_Solventscript {self.dir}/InputGrofiles{key}'

            command5 = f'mv {self.dir}/Packmol{key}/solvated.gro {self.dir}/InputGrofiles{key}'
            command6 = f'rm -r {self.dir}/{solvent}_Solvent && rm -r {self.dir}/Packmol{key}'

            command8 = f'mv {self.dir}/{solvent2}_Solvent2/{solvent2}_Solvent2.pdb {self.dir}/InputGrofiles && mv {self.dir}/{solvent2}_Solvent2/{solvent2}_Solvent2.gmx.itp {self.dir}/InputGrofiles{key} && mv {self.dir}/{solvent2}_Solvent2script {self.dir}/InputGrofiles{key}'
            command10 = f'rm -r {self.dir}/{solvent2}_Solvent2'
            print("cleaning up")

            subprocess.run(command4, shell=True, check=True)
            subprocess.run(command5, shell=True, check=True)
            subprocess.run(command6, shell=True, check=True)
            if len(self.solvent2) != 0:
                subprocess.run(command8, shell=True, check=True)
                subprocess.run(command10, shell=True, check=True)

            for i in range(len(self.solute)):
                com1 = (
                    f"mv {self.dir}/{self.solute[i].strip()[:3]}_Solute{1}/"
                    f"{self.solute[i].strip()[:3]}_Solute{1}.pdb {self.dir}/InputGrofiles{key} && "
                    f"mv {self.dir}/{self.solute[i].strip()[:3]}_Solute{1}script {self.dir}/InputGrofiles{key}")
                com2 = f'rm -r {self.dir}/{self.solute[i].strip()[:3]}_Solute{1}'
                subprocess.run(com1, shell=True, check=True)
                subprocess.run(com2, shell=True, check=True)

            g.reorg(self.solvent + "_Solvent", self.dir + f"/InputGrofiles{key}", key)

            if len(self.solvent2) != 0:
                g.reorg(self.solvent2 + "_Solvent2", self.dir + f"/InputGrofiles{key}")

            for j in range(len(self.solute)):
                print(f'ran soulte{j+1}')
                g.reorg(self.solute[j].strip()[:3] + f"_Solute{j + 1}", self.dir + f"/InputGrofiles{key}", key)

            subprocess.run(self.command12, shell=True, check=True)

            print("Starting the simulation")
    def initla(self):
        top.toopol(currentdir=self.dir,inital=True)
        subprocess.run(self.command12, shell=True, check=True)
        move_the_inital_system=f"cp {self.path_to_file}/solvated.gro {self.dir}/InputGrofiles{self.key}"
        move_topol=f"cp {self.path_to_file}/topol.top {self.dir}/InputGrofiles{self.key}"
        move_nmol=f"cp {self.path_to_file}/nmol.itp {self.dir}/InputGrofiles{self.key}"
        move_solvent_itp=  f"{self.path_to_file}/{self.fullname_solvent}.itp {self.dir}/InputGrofiles{self.key}"
        move_solvent_atomtype= f"{self.path_to_file}/{self.fullname_solvent}_atomtypes.itp {self.dir}/InputGrofiles{self.key}"

        subprocess.run(move_the_inital_system, shell=True, check=True)
        subprocess.run(move_topol, shell=True, check=True)
        subprocess.run(move_nmol, shell=True, check=True)
        subprocess.run(move_solvent_itp, shell=True, check=True)
        try:
            subprocess.run(move_solvent_atomtype, shell=True, check=True)
        except FileNotFoundError:
            print("Looks like you did not provide the atomtypes for your solvent. The system continued with out it.")
        for i in self.solute:
            move_solute_itp = f"{self.path_to_file}/{i}.itp {self.dir}/InputGrofiles{self.key}"
            move_solute_atomtype = f"{self.path_to_file}/{i}_atomtypes.itp {self.dir}/InputGrofiles{self.key}"
            subprocess.run(move_solute_itp, shell=True, check=True)
            try:
                subprocess.run(move_solute_atomtype, shell=True, check=True)
            except FileNotFoundError:
                print(
                    f"Looks like you did not provide the atomtypes for {i}. The system continued with out it.")












