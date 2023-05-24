import os
import re
import subprocess
from pathlib import Path
from pymatgen.io.gaussian import GaussianOutput


class HFCalcs:
    """
    For a mixed valence molecule, this sets up Hartree Fock calculations
    Copyright 2021, University of Kentucky
    """

    def __init__(self, mol_id, run_dir=None, collect_dir=None, still_running=True):
        self.mol_id = mol_id
        self.still_running = still_running
        self.run_dir = run_dir or os.getcwd()
        self.collect_dir = collect_dir or os.path.join(self.run_dir, 'hf_calculations/results')
        calc_dir = os.path.join(self.run_dir, mol_id, 'gaussian')
        self.gas_phase_dir = os.path.join(calc_dir, 'gas_phase')
        self.charges_dict = {'groundState': 0, 'cation1': 1, 'cation2': 2}
        self.multiplicity_dict = {'groundState': 1, 'cation1': 2, 'cation2': 1}

        if not os.path.isdir(self.collect_dir): os.mkdir(self.collect_dir)

    def generate_gjf(self, orig_gaus, out_file, charge='groundState', calculation='opt', functional='LC-wHPBE'):
        # Convert an individually inputted gaussian file to a gjf Gaussian input file

        calculation_name = 'hf_{}_{}'.format(calculation, charge) if calculation == 'hf' else '{}_{}'.format(
            calculation, charge)
        charge_num = self.charges_dict.get(charge)

        original_mol = GaussianOutput(orig_gaus)
        link0_parameters = {'%mem': '5GB', '%chk': '{}.chk'.format(calculation_name)}
        route_parameters = {calculation: '', 'SCF': '(MaxCycle=512)', 'Int': '(Grid=SuperFine)'}
        if calculation == 'energy':
            del route_parameters[calculation]
            del route_parameters['SCF']
            del route_parameters['Int']
        # get omega value
        orig_route_params = original_mol.route_parameters
        orig_route_params = {k.lower(): v for k, v in orig_route_params.items()}
        try:
            omega_string = re.sub("[^0-9]", "", orig_route_params['3108'])
        except KeyError:
            try:
                omega_string = re.sub("[^0-9]", "", orig_route_params['iop(3107'].split(',')[0])
            except KeyError:
                omega_string = None
        if omega_string:
            route_parameters.update({"iop(3/107={}, 3/108={})".format(omega_string, omega_string): ""})
        # generate gjf
        gau = original_mol.to_input(charge=charge_num, functional=functional, spin_multiplicity=self.multiplicity_dict[charge],
                                    route_parameters=route_parameters, link0_parameters=link0_parameters)

        gjf_file = gau.write_file(out_file)
        return gjf_file

    def check_if_run_finished(self, log_path, runfile_addon=None, opt_calc=True):
        txt_file = os.path.join(self.run_dir, 'hf_calculations', 'folders_to_run{}.txt'.format(runfile_addon or ''))
        dir_path = '/'.join(log_path.split('/')[:-1])
        file_name = log_path.split('/')[-1]

        if not os.path.isfile(log_path):
            print("No log file for {}".format(log_path))
        elif opt_calc:
            with open(log_path) as fn:
                log_fn = fn.readlines()
            last_line = log_fn[-1]
            if re.search('Normal termination', last_line):
                for line in log_fn:
                    if re.search(' Optimized Parameters', line):
                        return True
                print(
                    "Error. {} for {} optimization terminated normally, but not optimized parameters were found.".format(
                        file_name, self.mol_id))
            error = False
            for line in log_fn:
                if re.search('Error termination', line):
                    print("Error. {} for {} optimization did NOT terminated normally".format(file_name, self.mol_id))
                    error = True
                    break
            if not error and self.still_running:
                print("{} for {} optimization may still be running.".format(file_name, self.mol_id))
                return False
        else:
            mol = GaussianOutput(log_path)
            if mol.properly_terminated:
                return True
            print("Error. {} for {} did NOT terminated normally".format(file_name, self.mol_id))

        with open(txt_file, "a+") as fn:
            fn.write(log_path.replace(".log", ".gjf") + "\n")
        return False

    def dft_opt(self, charge):
        out_file = os.path.join(self.gas_phase_dir, 'opt', 'opt_{}.log'.format(charge))
        if os.path.isfile(out_file):
            self.collect(out_file)
            return out_file
        else:
            print("Error. No {} DFT optimization log file for {}".format(charge, self.mol_id))

    def hf_opt(self, charge):
        orig_gaus = self.dft_opt(charge)
        if orig_gaus:
            in_file = os.path.join(self.gas_phase_dir, 'opt', 'hf_opt_{}.gjf'.format(charge))
            out_file = os.path.join(self.gas_phase_dir, 'opt', 'hf_opt_{}.log'.format(charge))
            if os.path.isfile(in_file):
                if self.check_if_run_finished(out_file, runfile_addon='_'+charge):
                    self.collect(out_file)
                    return out_file
            else:
                self.generate_gjf(orig_gaus, in_file, charge=charge, functional='hf')
                self.check_if_run_finished(out_file, runfile_addon='_'+charge)

    def hfgeom_eng(self, charge):
        orig_gaus = self.hf_opt(charge)
        if orig_gaus:
            in_file = os.path.join(self.gas_phase_dir, 'energy', 'hfgeom_energy_{}.gjf'.format(charge))
            out_file = os.path.join(self.gas_phase_dir, 'energy', 'hfgeom_energy_{}.log'.format(charge))
            if os.path.isfile(in_file):
                if self.check_if_run_finished(out_file, runfile_addon='_'+charge, opt_calc=False):
                    self.collect(out_file)
                    return out_file
            else:
                self.generate_gjf(orig_gaus, in_file, charge=charge, calculation='energy', functional='LC-wHPBE')
                self.check_if_run_finished(out_file, runfile_addon='_'+charge, opt_calc=False)

    def collect(self, log_file, load_gaussian="module load gaussian/16rA.03"):
        # convert chk file
        calc_run_dir = '/'.join(log_file.split('/')[:-1])
        os.chdir(calc_run_dir)
        chk_file = log_file.replace(".log", ".chk")
        fchk_file = log_file.replace(".log", ".fchk")
        if not os.path.isfile(fchk_file):
            subprocess.call(load_gaussian, shell=True)
            subprocess.call("formchk {} {}".format(chk_file, fchk_file), shell=True)
        os.chdir(self.run_dir)

        # transfer files
        mol_dir = os.path.join(self.collect_dir, self.mol_id)
        if not os.path.isdir(mol_dir): os.mkdir(mol_dir)
        subprocess.call("cv {} {} {}".format(log_file, fchk_file, mol_dir), shell=True)

    def run_all_charges(self):
        for charge in self.charges_dict.keys():
            # self.hf_opt(charge)
            self.hfgeom_eng(charge)


mol_ids = ['90LPMI', '90NCLQ', '90QDBN', '90MOFW', '90PSTI', '05UXMH', '05JNUA', ]
# scratch = "05UXMH 90LPMI 90NCLQ 90QDBN 90MOFW 90PSTI 05JNUA"

run_dir = Path(os.getcwd()).parent.absolute()
print("Parent Directory: ", run_dir)
for mol_id in mol_ids:
    setup = HFCalcs(mol_id, run_dir=run_dir, still_running=False)
    setup.hfgeom_eng('groundState')
    # setup.hfgeom_eng('cation1')
    # setup.hfgeom_eng('cation2')
    print("Setup complete for ", mol_id)
    # break
