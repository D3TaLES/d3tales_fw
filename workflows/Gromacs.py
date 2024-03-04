import subprocess
import os
import multiprocessing
import datetime
import sys
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
import d3tales_fw.Fast.titration as titrate
import d3tales_fw.Fast.titrationPlotter as plotter


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
        self.dir = fw_spec.get("dir")
        self.smiles = self.get("smile") or fw_spec.get("smile")
        self.charge = self.get("charge") or fw_spec.get("charge")
        self.names = self.get("name") or fw_spec.get("name")
        self.type = self.get("Type") or fw_spec.get("TYPE","")

        l.lig(self.smiles, self.names[:3] ,self.names[:3] + f"_{self.type}", self.charge, self.dir)

        return FWAction(update_spec={})





@explicit_serialize
class MDPrep(FiretaskBase):

    def run_task(self, fw_spec):
        self.dir = self.get("dir") or fw_spec.get("dir")
        self.titration = self.get("titration_constant")
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
            print("did not work")
            exit()

        self.Density = self.get("den") or fw_spec.get("den")
        print(self.Solname or "did not work")

        pack.Solvate(self.Solname[:3], self.solute_name, self.conmatrix, self.Density, '', None, self.xdim, self.ydim,
                     self.zdim, self.dir, None, key)
        names=[]+ self.solvent_name+self.solute_name 
        smiles=[] + self.solvent_smiles + self.solute_smiles
        key=self.get("key")
        dft_folder=f"/project/cmri235_uksr/shasanka_conda_boss/launch"
        print(smiles)
        print(names)
        i = 0
        for iteams, name in zip(smiles,names):
            print(i)
            print(f'{dft_folder}/{iteams}/gaussian/gas_phase/opt/opt_groundState.log')
            if os.path.isfile(f'{dft_folder}/{iteams}/gaussian/gas_phase/opt/opt_groundState.log'):
                print("found dft")
                if i >=1:
                    transfer.trans(f"{name[:3]}_Solute1",iteams,key,1,self.dir,dft_folder,self.titration)
                else:
                    transfer.trans(f"{name[:3]}_Solvent",iteams,key,1,self.dir,dft_folder, self.titration)
            else:
                print("no dft")
                if i >=1:
                    transfer.trans(f"{name[:3]}_Solute1",iteams,key,0,self.dir,dft_folder, self.titration)
                else:
                    transfer.trans(f"{name[:3]}_Solvent",iteams,key,0,self.dir,dft_folder, self.titration)
            i+=1
        gro.gro(self.Solname[:3], self.solute_name, '', self.dir, self.xdim, self.ydim, self.zdim, key)
        return FWAction(update_spec={})

@explicit_serialize
class TitrationChargeScaler(FiretaskBase):

    def run_task(self, fw_spec):
        self.titration_matrix=fw_spec.get("titration_list") or self.get("titration_list")
        self.dir = self.get("dir") or fw_spec.get("dir")
        self.key = self.get("key")
        self.solvent_name = self.get("solvent_name") or fw_spec.get("solvent_name")
        self.solvent = self.solvent_name[0][:3]
        self.solute_name = self.get("solute_name") or fw_spec.get("solute_name")
        titrate.titration(self.titration_matrix, self.key, self.dir, self.solvent, self.solute_name)


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
        key=self.get("key")
        output_density = b.calculate_density()
        average=0.0
        for i in output_density:
            average = average+float(i)
        average= average/(float(len(output_density)))
        density_mat=fw_spec.get(f"Average_den{key}")
        density_mat.append(average)

        return FWAction(update_spec={f"Average_den{key}": density_mat})


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
                k.writelines(f'{i}: {j} \n   ')
        print(f'{key},{fw_spec.get("date_sumbit")},{self.dir}')
        self.orgainze(key,fw_spec.get("date_sumbit"),self.dir)

        return FWAction(update_spec={})

    def orgainze(self,key, date_sumbited, dirs):
        dir = dirs
        name = f'{str(date_sumbited)}_'
        for i in (str(datetime.datetime.now()).split()[1]).split(":"):
            name += f"_{str(i)[:2]}"
        subprocess.run([f"mkdir {dir}/run_{name}"], shell=True, check=True)
        subprocess.run([f"mv {dir}/key {dir}/run_{name}"], shell=True, check=True)
        print(key)
        for iteams in key:
            subprocess.run([f"mv {dir}/InputGrofiles{key[iteams]} {dir}/run_{name}"], shell=True)
            subprocess.run([f"mv {dir}/Output{key[iteams]} {dir}/run_{name}"], shell=True)


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



@explicit_serialize
class Graph_plotter(FiretaskBase):

    def run_task(self, fw_spec):
        key=self.get("This_key")
        titration_list=fw_spec.get("titration_list")
        Average_simulation_density=fw_spec.get(f"Average_den{key}")
        Average_simulation_density_dic={i:j for i,j in zip(titration_list,Average_simulation_density) }
        plotter.TitrationPlotter( Average_simulation_density_dic, key)

        return FWAction(update_spec={})