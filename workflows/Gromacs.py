import subprocess
import multiprocessing

from atomate.utils.utils import get_logger, env_chk
from fireworks import FiretaskBase, explicit_serialize, FWAction


logger = get_logger(__name__)
cpus = [multiprocessing.cpu_count() if multiprocessing.cpu_count() < 16 else 16]
nprocs = str(cpus[0])


# Copyright 2021, University of Kentucky


@explicit_serialize
class MDInit(FiretaskBase):

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


@explicit_serialize
class LigParGenUpdate(FiretaskBase):

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
class EnergyMinimization(FiretaskBase):

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
