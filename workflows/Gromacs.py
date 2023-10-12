import subprocess
import os
import multiprocessing
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import d3tales_fw.Fast.ASMD_1 as run
import d3tales_fw.Fast.scriptMaker as m
import d3tales_fw.Fast.lig as l
import d3tales_fw.Fast.dft as dft
import d3tales_fw.Fast.waiting as wait
import d3tales_fw.Fast.chargeTrasnfer as transfer
import d3tales_fw.Fast.packmol as pack
import d3tales_fw.Fast.make_gro as gro

from atomate.utils.utils import get_logger, env_chk
from fireworks import FiretaskBase, explicit_serialize, FWAction
logger = get_logger(__name__)
cpus = [multiprocessing.cpu_count() if multiprocessing.cpu_count() < 16 else 16]
nprocs = str(cpus[0])


# Copyright 2021, University of Kentucky


@explicit_serialize
class Ligpargen(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameter

        return FWAction()


@explicit_serialize
class Ligpargen(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameter
        self.dir = fw_spec.get("dir")
        self.smiles = self.get("smile") or fw_spec.get("smile")
        self.charge = self.get("charge") or fw_spec.get("charge")
        self.names = self.get("name") or fw_spec.get("name")
        self.type = self.get("Type") or fw_spec.get("TYPE") or ""

        l.lig(self.smiles, self.names + f"_{self.type}", self.charge, self.dir)

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
        key = self.get("key")
        chagreMatrix = []
        with open(f'{self.dir}/{self.name}script{key}/{self.name}-final.log', 'r') as log, open(
                f'{self.dir}/{self.name}{key}/{self.name}.gmx.itp', 'r') as gmx, open(
            f'{self.dir}/{self.name}script{key}/{self.name}.itp', 'a') as itp:
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
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.dir = self.get("dir") or fw_spec.get("dir")
        self.charge = self.get("charge") or fw_spec.get("charge")
        self.solvent_name = self.get("solvent_name") or fw_spec.get("solvent_name")
        self.Solname = self.solvent_name[0]
        self.solute_smiles= self.get("solute_smiles")
        self.solvent_smiles= self.get("solvent_smiles")
        self.solute_name = self.get("solute_name") or fw_spec.get("solute_name")
        self.xdim = float(self.get("x") or fw_spec.get("x"))
        self.ydim = float(self.get("y") or fw_spec.get("y"))
        self.zdim = float(self.get("z") or fw_spec.get("z"))
        self.conmatrix = self.get("conmatrix") or fw_spec.get("conmatrix")
        key = self.get("key")
        if self.conmatrix == None:
            prtin("shit failed")
            exit()

        self.Density = self.get("den") or fw_spec.get("den")
        pack.Solvate(self.Solname[:3], self.solute_name, self.conmatrix, self.Density, '', None, self.xdim, self.ydim,
                     self.zdim, self.dir, None, key)
        names=[]+ self.solvent_name+self.solute_name 
        smiles=[] + self.solvent_smiles + self.solute_smiles
        key=self.get("key")
        dft_folder="/project/cmri235_uksr/shasanka_conda_boss/dft_folder"
        i = 0
        for iteams, name in zip(smiles,names):
            i+=1
            print(f'{dft_folder}/{iteams}')
            if os.path.isfile(f'{dft_folder}/{iteams}'):
                if i >1:
                    transfer.trans(f"{name[:3]}_Solute1",iteams,key,1,self.dir,dft_folder)
                else:
                    transfer.trans(f"{name[:3]}_Solvent",iteams,key,1,self.dir,dft_folder)
            else:
                if i >1:
                    transfer.trans(f"{name[:3]}_Solute1",iteams,key,0,self.dir,dft_folder)
                else:
                    transfer.trans(f"{name[:3]}_Solvent",iteams,key,0,self.dir,dft_folder)
        gro.gro(self.Solname[:3], self.solute_name, '', self.dir, self.xdim, self.ydim, self.zdim, key)
        return FWAction(update_spec={})


@explicit_serialize
class EnergyMinimization(FiretaskBase):

    def run_task(self, fw_spec):
        runer = run.ASMD(self.get("key"))
        a = runer.EnergyMin()
        return FWAction(update_spec={"outputEM": a, "runner": runer})


@explicit_serialize
class NVT(FiretaskBase):

    def run_task(self, fw_spec):
        b = run.ASMD(self.get("key"))
        output_NVT = b.NVT()

        return FWAction(update_spec={"outputNVT": output_NVT})


@explicit_serialize
class NPT(FiretaskBase):

    def run_task(self, fw_spec):
        b = run.ASMD(self.get("key"))
        output_NPT = b.NPT()

        return FWAction(update_spec={"outputNPT": output_NPT})


@explicit_serialize
class Density(FiretaskBase):

    def run_task(self, fw_spec):
        b = run.ASMD(self.get("key"))
        output_density = b.calculate_density()

        return FWAction(update_spec={"outputden": output_density})


@explicit_serialize
class Den_checker(FiretaskBase):

    def run_task(self, fw_spec):
        den = float(self.get("den") or fw_spec.get("den"))
        molarMass = float(self.get("MM") or fw_spec.get("MM"))
        Density = self.get("outputden") or fw_spec.get("outputden")
        b = run.ASMD(self.get("key"))
        output = b.check_density_accuracy(float(den * molarMass), Density)

        return FWAction(update_spec={"outputdenCheck": output})


@explicit_serialize
class trj_corrector(FiretaskBase):

    def run_task(self, fw_spec):
        b = run.ASMD(self.get("key"))
        output = b.correct()

        return FWAction(update_spec={"outputTrjCor": output})


@explicit_serialize
class Index(FiretaskBase):

    def run_task(self, fw_spec):
        b = run.ASMD(self.get("key"))
        output = b.index_file()

        return FWAction(update_spec={"outputIndex": output})


@explicit_serialize
class residue(FiretaskBase):

    def run_task(self, fw_spec):
        b = run.ASMD(self.get("key"))
        output = b.extract_residues_from_itp()

        return FWAction(update_spec={"outputextract": output})


@explicit_serialize
class rdf(FiretaskBase):

    def run_task(self, fw_spec):
        residue = fw_spec.get("outputextract") or self.get("outputextract")

        b = run.ASMD(self.get("key"))
        out = b.rdf(residue, 5)

        return FWAction(update_spec={"rdf": out})


@explicit_serialize
class cord(FiretaskBase):

    def run_task(self, fw_spec):
        cord = fw_spec.get("rdf") or self.get("rdf")

        b = run.ASMD(self.get("key"))
        b.cordination_number(cord)

        return FWAction(update_spec={})


@explicit_serialize


class key_gen(FiretaskBase):

    def run_task(self, fw_spec):
        key = fw_spec.get("key_dic")
        self.dir = self.get("dir") or fw_spec.get("dir")
        with open(f"{self.dir}/key", 'a') as k:
            for i, j in key.items():
                k.write(f'{i}: {j} \n')

        return FWAction(update_spec={})

@explicit_serialize

class dft_checker(FiretaskBase):
    def run_task(self, fw_spec):
        smiles= fw_spec.get("smiles_list")
        names=fw_spec.get("name_list")
        key=self.get("key")
        dft_folder="/project/cmri235_uksr/shasanka_conda_boss/dft_folder"
        for iteams, name in zip(smiles,names):
            if os.path.isfile(dft_folder+iteams):
                transfer.trans(name[:3],iteams,key,1)
            else:
                transfer.trans(name[:3], iteams, key, 0)

        return FWAction(update_spec={})
























