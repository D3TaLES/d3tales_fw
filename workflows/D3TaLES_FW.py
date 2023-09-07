from fireworks import Firework
from d3tales_fw.workflows.Gromacs import *
# from d3tales_fw.workflows.Gaussian import *
# from d3tales_fw.workflows.Initialize import *
from d3tales_fw.workflows.envwf import G16_CMD, PATH, PROC_VM_KEY, RUNFILE_LOG

# Copyright 2021-2022, University of Kentucky

CALC_CATEGORY = 'gaussian'  # gaussian gaussian_mcc


class InitializeMolecule(Firework):
    def __init__(self, name="initial_molecule_data", parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [MoleculeInit(**kwargs)]
        super(InitializeMolecule, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class MolOpt(Firework):
    def __init__(self, name="opt_mol", parents=None, priority=None, name_tag='', **kwargs):
        kwargs.pop("use_iop") if "use_iop" in kwargs.keys() else None
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        t = [RunGaussianOpt(name=name, name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY,
                            path=PATH, runfile_log=RUNFILE_LOG, skip_freq=True, use_iop=False, **kwargs)]
        super(MolOpt, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class WTuning(Firework):
    def __init__(self, name="wtuning", parents=None,  priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        t = [RunWtuning(path=PATH, runfile_log=RUNFILE_LOG, name=name, name_tag=name_tag, geometry="opt_mol", **kwargs)]
        super(WTuning, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class Optimization(Firework):
    def __init__(self, species="groundState", geometry="opt_mol", solvent=None, parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        _type = "solv_" + solvent if solvent else "gas_phase"
        name = "solv_opt_" + species if solvent else "opt_" + species
        t = [RunGaussianOpt(name=name, name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH, runfile_log=RUNFILE_LOG,
                            geometry=geometry, type=_type, subtype=species, solvent=solvent, **kwargs)]
        super(Optimization, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class Energy(Firework):
    def __init__(self, species="groundState", geometry="groundState", parents=None, priority=None, name_tag='', solvent=None, **kwargs):
        abbrev_dict = {"groundState": "gs", "cation1": "c1", "anion1": "a1", "cation2": "c2", "anion2": "a2", }
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        _type = "solv_" + solvent if solvent else "gas_phase"
        prefix = "solv_energy" if solvent else "energy"
        name = "{}_{}{}".format(prefix, abbrev_dict[geometry.split("_")[-1]], abbrev_dict[species])
        t = [RunGaussianEnergy(name=name, name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH, runfile_log=RUNFILE_LOG,
                               geometry=geometry, type=_type, subtype=species, solvent=solvent, **kwargs)]
        super(Energy, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class TDDFT(Firework):
    def __init__(self, species='groundState', solvent=None, parents=None, priority=None, name_tag='', solv_geom=False, **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        name = "tddft_{}".format(species)
        geom_prefix = "solv_opt_" if solv_geom else "opt_"
        route_keys = list(map(str.lower, kwargs.get("paramset").route_parameters.get("td", [])))
        tddft_prefix = "singlet_" if "singlets" in route_keys else "triplet_" if "triplets" in route_keys else "singlet_triplet_" if "50-50" in route_keys else "tddft_"
        t = [RunGaussianTDDFT(name=name, name_tag=name_tag, path=PATH, runfile_log=RUNFILE_LOG, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, type='tddft',
                              geometry="{}{}".format(geom_prefix, species), prefix=tddft_prefix+species, solvent=solvent, **kwargs)]
        super(TDDFT, self).__init__(t, parents=parents, spec=spec, name="{}tddft_{}".format(name_tag, species))


class LowestEConformer(Firework):
    def __init__(self, parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        t = [GetLowestEConformer(name="conformers", name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH,
                                 runfile_log=RUNFILE_LOG, type="gas_phase", **kwargs)]
        super(LowestEConformer, self).__init__(t, parents=parents, spec=spec, name="{}conformers".format(name_tag))


class DihedRot(Firework):
    def __init__(self, geometry="groundState", dihed_degree=0, parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        name = "dihedrot_" + str(dihed_degree).zfill(3)
        t = [RunGaussianDihedRot(name=name, name_tag=name_tag, dihed_degree=dihed_degree, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH, runfile_log=RUNFILE_LOG,
                                 geometry=geometry, type="gas_phase", subtype=str(dihed_degree).zfill(3), **kwargs)]
        super(DihedRot, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class EmailStart(Firework):
    def __init__(self, parents=None, priority=None, name_tag='', identifier="", email="", username="", **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [EmailStarting(identifier=identifier, email=email, username=username)]
        super(EmailStart, self).__init__(t, parents=parents, spec=spec, name="{}email_starting".format(name_tag))


class EmailEnd(Firework):
    def __init__(self, parents=None, priority=None, name_tag='', identifier="", email="", username="", **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [EmailFinished(identifier=identifier, email=email, username=username)]
        super(EmailEnd, self).__init__(t, parents=parents, spec=spec, name="{}email_finished".format(name_tag))


# ---------------- Molecular Dynamics FWs ----------------
class InitializeMD(Firework):
    def __init__(self, name="initial_MD_data", parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [MDInit(**kwargs)]
        super(InitializeMD, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class Ligpargen_FW(Firework):
    def __init__(self, name=None, parents=None, priority=None, name_tag='', smiles=None, con=None,Type=None, di=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority, 'smiles':smiles, 'name':name, 'charge':con, "Type":Type} if priority else {'_category': 'gromacs','smiles':smiles, 'name':name, 'charge':con, "Type":Type, "dir":di, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [Ligpargen(name=name[:3],smile=smiles, charge=con,Type=Type,**kwargs)] ## this is passed for self, self.get() gets this
        super(Ligpargen_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class Pack_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [MDPrep(solute_name=solute_name, solvent_name=solvent_name, x=x,y=y,z=z,di=di,conmatrix=conmatrix,den=den,key=key )] ## this is passed for self, self.get() gets this
        super(Pack_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class EM_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [EnergyMinimization(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(EM_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class NVT_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [NVT(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(NVT_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class NPT_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [NPT(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(NPT_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class Density_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [Density(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(Density_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class TR_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [trj_corrector(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(TR_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class Index_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [Index(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(Index_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class RES_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [residue(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(RES_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class RDF_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [rdf(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(RDF_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class CORD_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [cord(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(CORD_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class Check_FW(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,mm=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den,"MM":mm, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [Den_checker(**kwargs,key=key)] ## this is passed for self, self.get() gets this
        super(Check_FW, self).__init__(t, parents=parents, spec=spec, name=name)

class key_GEN(Firework):
    def __init__(self,name=None, parents=None, priority=None, name_tag='',solute_name=[], solvent_name=[], x=None,y=None,z=None, di=None,conmatrix=None,den=None,key=None,**kwargs):
        spec = {'_category': 'processing', '_priority': priority,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"dir":di, **kwargs} if priority else {'_category': 'gromacs', "dir":di,"solute_name":solute_name, "solvent_name":solvent_name, "x":x,"y":y,"z":z,"conmatrix":conmatrix,"den":den, **kwargs} ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [key_gen(**kwargs)] ## this is passed for self, self.get() gets this
        super(key_GEN, self).__init__(t, parents=parents, spec=spec, name="key_gen")

