
import subprocess
import os
import multiprocessing
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import d3tales_fw.workflows.ASMD as run

from atomate.utils.utils import get_logger, env_chk
from fireworks import FiretaskBase, explicit_serialize, FWAction


logger = get_logger(__name__)
cpus = [multiprocessing.cpu_count() if multiprocessing.cpu_count() < 16 else 16]
nprocs = str(cpus[0])


# Copyright 2021, University of Kentucky


@explicit_serialize
class MDInit(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameter
        return FWAction()


@explicit_serialize
class Ligpargen(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameter
        self.ligpargen_cmd = env_chk(self.get("ligpargen_cmd"), fw_spec)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.full_name = name_tag + self['name']
        self.mol = self['name']
        self.BOSS= self.get("boss") or fw_spec.get("boss")
        self.paramset = self["paramset"]
        self.charge= self.get("charge") or fw_spec.get("charge") or "0"
        self.identifier = fw_spec.get("identifier", ) or self.get("identifier")
        self.ligpargen_fn = fw_spec.get("ligpargen_fn") or self.get("ligpargen_fn") or "ligpargen"
        self.dir = self.get("dir") or fw_spec.get("dir")
        self.smiles = fw_spec.get("smiles") or self.get("smiles")
        subprocess.run([f'mkdir {self.dir}/{self.mol}script'], shell=True)
        conda_activate = f"source {self.dir}/miniconda3/bin/activate && conda activate ligpg"
        export_bossdir = f"export BOSSdir={self.BOSS}"
        ligpargen_cmd = f"ligpargen -s '{self.smiles}' -n {self.mol} -p {self.mol} -r {self.mol} -c 0 -o 0 -cgen CM1A"
        self.singularity_container = self.get("singu") or fw_spec.get("singu")

        cmd = ["singularity", "exec", self.singularity_container, "bash", "-c",
               f'{conda_activate} && {export_bossdir} && {ligpargen_cmd}']
        subprocess.Popen(cmd).wait()

        while os.path.isfile(f'{self.dir}/{self.mol}/{self.mol}-debug.pdb') == False:
            print("still making pdb")
            time.sleep(1)
        pdb_location = Chem.MolFromPDBFile(f'{self.dir}/{self.mol}/{self.mol}-debug.pdb')
        pdb_location = Chem.AddHs(pdb_location)
        AllChem.EmbedMolecule(pdb_location)
        Chem.MolToMolFile(pdb_location, f'{self.mol}/{self.mol}.mol')
        comand2 = f'obabel -imol {self.dir}/{self.mol}/{self.mol}.mol -ogjf -O {self.mol}script/{self.mol}.gjf --gen3D'
        comand3 = f'obabel -igro {self.dir}/ {self.mol}/{self.mol}.gmx.gro -opdb -O {self.mol}/{self.mol}.pdb '
        subprocess.run([comand2], shell=True)
        subprocess.run([comand3], shell=True)

        return FWAction(update_spec={})


@explicit_serialize
class LigParGenInit(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameters
        self.ligpargen_cmd = env_chk(self.get("ligpargen_cmd"), fw_spec)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.full_name = name_tag + self['name']
        self.calc_name = self['name']
        self.paramset = self["paramset"]
        self.identifier = fw_spec.get("identifier", ) or self.get("identifier")
        self.gaussian_file_name = fw_spec.get("ligpargen_fn") or self.get("ligpargen_fn") or "ligpargen"
    ##Not sure what this is
        return FWAction(update_spec={})

##we need to run the dft here for the charge deravation

@explicit_serialize
class LigParGenUpdate(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameters
        self.ligpargen_cmd = env_chk(self.get("ligpargen_cmd"), fw_spec)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.full_name = name_tag + self['name']
        self.name = self['name']
        self.paramset = self["paramset"]
        self.identifier = fw_spec.get("identifier", ) or self.get("identifier")
        self.gaussian_file_name = fw_spec.get("ligpargen_fn") or self.get("ligpargen_fn") or "ligpargen"
        self.dir = self.get("dir") or fw_spec.get("dir")
        chagreMatrix = []
        with open(f'{self.dir}/{self.name}script/{self.name}-final.log', 'r') as log, open(
                f'{self.dir}/{self.name}/{self.name}.gmx.itp', 'r') as gmx, open(
                f'{self.dir}/{self.name}script/{self.name}.itp', 'a') as itp:
            lines = log.readlines()
            gmx_lines = gmx.readlines()
            if len(gmx_lines) == 0:
                print("gmx is empty, try again")
                subprocess.run(["ls"], shell=True)
            index_ESP = 0
            index_lastEsp = 0
            index_GMXESP = 0
            index_lastGMXESP = 0
            for line in gmx_lines:
                if line.strip() == "[ atoms ]":
                    index_GMXESP = gmx_lines.index(line) + 2
                    break
            for line in gmx_lines:
                if line.strip() == "[ bonds ]":
                    index_lastGMXESP = gmx_lines.index(line) - 2
            for line in lines:
                if line.strip() == "ESP charges:":
                    index_ESP = lines.index(line) + 2
                    break
            for line in lines:
                if line.strip()[0:18] == "Sum of ESP charges":
                    index_lastEsp = lines.index(line) - 1
                    break

            for a in range(index_lastEsp - index_ESP + 1):
                charge_int = float(lines[index_ESP + a].split()[2])
                charge = f'{charge_int:.4f}'.rjust(11)

                chagreMatrix.append(charge)
            for char in range(index_lastGMXESP - index_GMXESP + 1):
                chargeIndexGmx = gmx_lines[index_GMXESP - 1].index("charge") - 6

                gmx_lines[index_GMXESP + char] = gmx_lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                    chagreMatrix[char]) + gmx_lines[index_GMXESP + char][chargeIndexGmx + 11:]
            for line in gmx_lines:
                itp.writelines(line)
        print("Itp file is made")
        return FWAction(update_spec={})


@explicit_serialize
class MDPrep(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameters
        self.ligpargen_cmd = env_chk(self.get("ligpargen_cmd"), fw_spec)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.full_name = name_tag + self['name']
        self.calc_name = self['name']
        self.paramset = self["paramset"]
        self.identifier = fw_spec.get("identifier", ) or self.get("identifier")
        self.gaussian_file_name = fw_spec.get("ligpargen_fn") or self.get("ligpargen_fn") or "ligpargen"

        return FWAction(update_spec={})

@explicit_serialize
class Packer(FiretaskBase):

    def run_task(self, fw_spec):
        self.dir = self.get("dir") or fw_spec.get("dir")
        self.solvent = self.get("solvent1") or fw_spec.get("solvetn1")
        self.solvent2 =  self.get("solvent2") or fw_spec.get("solvetn2") or None
        self.solutes = self.get("solute_matrix") or fw_spec.get("solute_matrix")
        self.con = self.get("con_matrix") or fw_spec.get("con_matrix")
        self.den = float(self.get("density1") or fw_spec.get("density1"))
        try:
            self.den2 =  float(self.get("density2") or fw_spec.get("density2"))
        except:
            self.den2 = ''
        try:
            self.ratio =  float(self.get("ratio") or fw_spec.get("ratio"))
        except:
            self.ratio = 1
        self.x = 10 * float(self.get("x-dim") or fw_spec.get("x-dim"))
        self.y = 10 * float(self.get("y-dim") or fw_spec.get("y-dim"))
        self.z = 10 * float(self.get("y-dim") or fw_spec.get("y-dim"))

        self.numberSolvent = int(
            self.den * (6.4e-20) * 6.02214e23 * (((self.x / 10) * (self.y / 10) * (self.z / 10)) / (40 * 40 * 40)))
        try:
            self.numberSolvent2 = int(self.den2 * (6.4e-20) * 6.02214e23 * (
                        ((self.x / 10) * (self.y / 10) * (self.z / 10)) / (40 * 40 * 40)) * self.ratio)
        except:
            self.numberSolvent2 = 0

        subprocess.run(f"mkdir {self.dir}/Packmol", shell=True)

        with open(f'{self.dir}/Packmol/mixture.inp', 'a') as m:
            info_all = ["#", f"# A mixture of {self.solvent},{self.solvent2} and solutes", "#", "", "tolerance 2.0",
                        "filetype pdb", f"output Packmol/Mixture.pdb",
                        "", f"structure {self.solvent}_Solvent/{self.solvent}_Solvent.pdb",
                        f"  number {self.numberSolvent}", f"  inside box 0. 0. 0. {self.x} {self.y} {self.z}",
                        "end structure", "",
                        f"structure {self.solvent2}_Solvent2/{self.solvent2}_Solvent2.pdb",
                        f"  number {self.numberSolvent2}",
                        f"  inside box 0. 0. 0. {self.x} {self.y} {self.z}", "end structure",
                        ]
            for i in range(len(self.con)):
                self.mol = int(
                    float(self.con[i].strip()) * (6.4e-20) * 6.02214e23 * (
                                ((self.x / 10) * (self.y / 10) * (self.z / 10)) / (40 * 40 * 40)))
                info_all.append(" ")
                info_all.append(
                    f"structure {self.solute[i].strip()[:3]}_Solute{i + 1}/{self.solute[i].strip()[:3]}_Solute{i + 1}.pdb")
                info_all.append(f"  number {self.mol}")
                info_all.append(f"  inside box 0. 0. 0. {self.x} {self.y} {self.z}")
                info_all.append("end structure")

            info_none = ["#", f"# A mixture of {self.solvent},and solutes", "#", "", "tolerance 2.0", "filetype pdb",
                         f"output Packmol/Mixture.pdb",
                         "", f"structure {self.solvent}_Solvent/{self.solvent}_Solvent.pdb",
                         f"  number {self.numberSolvent}", f"  inside box 0. 0. 0. {self.x} {self.y} {self.z}",
                         "end structure"]

            for i in range(len(self.con)):
                self.molnum = int(
                    float(self.con[i].strip()) * (6.4e-20) * 6.02214e23 * (
                                ((self.x / 10) * (self.y / 10) * (self.z / 10)) / (40 * 40 * 40)))
                info_none.append(" ")
                info_none.append(
                    f"structure {self.solute[i].strip()[:3]}_Solute{i + 1}/{self.solute[i].strip()[:3]}_Solute{i + 1}.pdb")
                info_none.append(f"  number {self.molnum}")
                info_none.append(f"  inside box 0. 0. 0. {self.x} {self.y} {self.z}")
                info_none.append("end structure")

            if len(self.solvent2) != 0:
                for lines in info_all:
                    m.writelines(lines + "\n")
            if len(self.solvent2) == 0:
                for lines in info_none:
                    m.writelines(lines + "\n")

        while not os.path.isfile(f"{self.dir}/Packmol/mixture.inp"):
            time.sleep(1)
        conda_activate = f"source {self.dir}/miniconda3/bin/activate"
        packmol_cmd = f"packmol < {self.dir}/Packmol/mixture.inp"

        cmd = ["bash", "-c",
               f'{conda_activate} && {packmol_cmd}']
        subprocess.Popen(cmd).wait()
        self.to()
        return FWAction(update_spec={})

    def to(self, ):
        self.toopol(self.solvent, self.numberSolvent, self.solvent2, self.numberSolvent2, self.solutes, self.con,
                   self.dir, self.x / 10, self.y / 10, self.z / 10)
    def toopol(self, solvent1, solvent1N, solvent2, solvent2N, solute, con, currentdir, x,y,z):
        self.solvent = solvent1
        self.solvent2 = solvent2
        self.solute1 = solute
        self.con = con
        self.dir = currentdir
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        command2 = f'mkdir {self.dir}/InputGrofiles'
        if os.path.isdir(f'{self.dir}/InputGrofiles') == True:
            folder_number=0
            if os.path.isdir(f'{self.dir}/Old_InputGrofiles_Output')==False:
                make_old= f'mkdir {self.dir}/Old_InputGrofiles_Output'
                subprocess.run(make_old, shell=True, check=True)
            while True:
                if os.path.isdir(f'{self.dir}/Old_InputGrofiles_Output/InputGrofiles{folder_number+1}'):
                    folder_number +=1
                else:
                    move = f'mv {self.dir}/InputGrofiles {self.dir}/Old_InputGrofiles_Output/InputGrofiles{folder_number+1}'
                    subprocess.run(move, shell=True, check=True)
                    break
        if os.path.isdir(f'{self.dir}/Output') == True:
            folder_number = 0
            while True:
                if os.path.isdir(f'{self.dir}/Old_InputGrofiles_Output/Output{folder_number + 1}'):
                    folder_number += 1
                else:
                    move = f'mv {self.dir}/Output {self.dir}/Old_InputGrofiles_Output/Output{folder_number + 1}'
                    subprocess.run(move, shell=True, check=True)
                    break
        subprocess.run(command2, shell=True, check=True)
        cmd = f"touch {self.dir}/InputGrofiles/topol.top && touch InputGrofiles/nmol.itp"
        subprocess.run(cmd, shell=True, check=True)
        with open(f"{self.dir}/InputGrofiles/topol.top", 'a') as top, open(
                f"{self.dir}/InputGrofiles/nmol.itp", 'a') as nmol:

            lines_atomtypes = ['#include "oplsaa.ff/forcefield.itp"', f'#include "{self.solvent}_Solvent_atomtype.itp"'
                          ]
            lines_itp=[f'#include "{self.solvent}_Solvent.itp"']

            last_lines=["", '[system]',
                          f'{solvent1}', "", '#include "nmol.itp"']

            if len(self.solvent2) !=0:
                lines_atomtypes.append(f'#include "{self.solvent2}_Solvent2_atomtype.itp"')
                lines_itp.append(f'#include "{self.solvent2}_Solvent2.itp"')
            for i in range(len(self.solute1)):
                lines_atomtypes.append(f'#include "{self.solute1[i].strip()[:3]}_Solute{i+1}_atomtype.itp"')
                lines_itp.append(f'#include "{self.solute1[i].strip()[:3]}_Solute{i+1}.itp"')




            for l in lines_atomtypes:
                top.write(l + "\n")
            for sha in lines_itp:
                top.write(sha + "\n")
            for lami in last_lines:
                top.write(lami + "\n")


            linesNmol = ['[ molecules ]', ';molecules #molecules', "", f'{self.solvent}   {solvent1N}'
                            ]
            if len(self.solvent2) != 0:
                linesNmol.append(f'{self.solvent2}   {solvent2N}')

            for i in range(len(self.solute1)):
                self.mol = int(
                    float(self.con[i].strip()) * (6.4e-20) * 6.02214e23 * ((self.x * self.y * self.z) / (40 * 40 * 40)))
                linesNmol.append(f'{self.solute1[i].strip()[:3]}   {self.mol}')


            for li in linesNmol:
                nmol.write(li + "\n")







@explicit_serialize
class EnergyMinimization(FiretaskBase):

    def run_task(self, fw_spec):
        runer = run.ASMD()
        a = runer.EnergyMin()
        return FWAction(update_spec={"outputEM":a, "runner":runer})

@explicit_serialize
class NVT(FiretaskBase):

    def run_task(self, fw_spec):
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output_NVT= b.NVT()

        return FWAction(update_spec={"outputNVT":output_NVT})


@explicit_serialize
class NTP(FiretaskBase):

    def run_task(self, fw_spec):
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output_NPT = b.NPT()

        return FWAction(update_spec={"outputNPT": output_NPT})


@explicit_serialize
class Density(FiretaskBase):

    def run_task(self, fw_spec):
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output_density=b.calculate_density()

        return FWAction(update_spec={"outputden": output_density})


@explicit_serialize
class Den_checker(FiretaskBase):

    def run_task(self, fw_spec):
        den= float(self.get("density1") or fw_spec.get("density1"))
        molarMass= self.get("MM") or fw_spec.get("MM")
        Density= self.get("outputden") or fw_spec.get("outputden")
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output = b.check_density_accuracy( float(den)* float(molarMass),Density)


        return FWAction(update_spec={"outputdenCheck": output})

@explicit_serialize
class trj_corrector(FiretaskBase):

    def run_task(self, fw_spec):
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output = b.correct()

        return FWAction(update_spec={"outputTrjCor": output})
@explicit_serialize
class Index(FiretaskBase):

    def run_task(self, fw_spec):
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output = b.index_file()

        return FWAction(update_spec={"outputIndex": output})
@explicit_serialize
class residue(FiretaskBase):

    def run_task(self, fw_spec):
        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        output = b.extract_residues_from_itp()

        return FWAction(update_spec={"outputextract": output})

@explicit_serialize
class rdf(FiretaskBase):

    def run_task(self, fw_spec):
        residue =fw_spec.get("outputextract") or self.get("outputextract")

        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        out=b.rdf(residue,5)

        return FWAction(update_spec={"rdf":out})


@explicit_serialize
class cord(FiretaskBase):

    def run_task(self, fw_spec):
        cord =fw_spec.get("rdf") or self.get("rdf")

        b = fw_spec.get("runner") or self.get("runner") or run.ASMD()
        b.cordination_number(cord)

        return FWAction(update_spec={})


























