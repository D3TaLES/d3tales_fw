import rdkit.Chem as Chem
from rdkit.Chem import AllChem, Descriptors
import rdkit as rd
import veloxchem as vlx
import subprocess



def get_charge(name, q, smiles, dir):
    tr_mol = rd.Chem.MolFromSmiles(smiles) # I found that if I just use the xyz file RDKIT
    # cannot count num of unpaired e, so I am defining the molecule using smiles first
    b = rd.Chem.AddHs(tr_mol) # adding H so i can later turn XYZ to PDB
    rd.Chem.AllChem.EmbedMolecule(b) # Generating 3D rep of molecule
    rd.Chem.MolToXYZFile(b, f"{name}.xyz") # Generating XYZ, so I can turn this into PDB for Ligpargen
    molecule= vlx.Molecule.read_xyz_file(f"{name}.xyz") # Loading Molecule into veloxchem
    basis = vlx.MolecularBasis.read(molecule, "cc-pVDZ") # setting up for calc
    unpaired=Descriptors.NumRadicalElectrons(tr_mol) # finding num of unpaired e for multiplicity calc
    multiplicity = unpaired+1 + int(q) # 2S+1+charge
    print(f"the spin multi: {multiplicity}")
    vlx.Molecule.set_multiplicity(molecule,multiplicity) # setting multi
    vlx.Molecule.set_charge(molecule, q) # Setting Charge
    scf_drv= vlx.ScfUnrestrictedDriver() # Unrestricted for sate>singlet
    if multiplicity==1:
        scf_drv= vlx.ScfRestrictedDriver() # for singlet
    scf_drv.ostream.mute() # shutting off the output stream for the calc

    scf_drv.xcfun = "B3LYP"
    # avaible functioals: ['SLATER', 'SLDA', 'B88X', 'BLYP', 'B3LYP', 'BHANDH', 'BHANDHLYP', 'PBE', 'PBE0', 'REVPBE',
    # 'BP86', 'PW91', 'MPW1K', 'OLYP', 'O3LYP', 'X3LYP', 'B97', 'B97-1', 'B97-2', 'B97-3', 'TPSS', 'TPSSH',
    # 'REVTPSS', 'PKZB', 'SCAN', 'RSCAN', 'R2SCAN', 'M05', 'M05-2X', 'M06', 'M06-2X', 'M06-HF', 'M06-L', 'M11-L',
    # 'MPW1B95', 'MPWB1K', 'PW6B95', 'PWB6K']

    scf_drv.grid_level = 4  # grid spacing for numerical calc, i dont know how to pick this but this seems to be the
    # defult
    scf_drv.conv_thresh = 1.0e-6  # setting criteria of deltaE for convergence


    try:
        scf_results = scf_drv.compute(molecule, basis)
        esp_drv = vlx.RespChargesDriver() #restraied ESP calc driver
        esp_charges = esp_drv.compute(molecule, basis, scf_results, "esp") # calc
    except AssertionError:
        print(f" Error in DFT clac ")
        return None
    subprocess.run(f"obabel -ixyz {name}.xyz -opdb -O {name}.pdb")
    return esp_charges  #returns an np array

if __name__=="__main__":
    print( get_charge("try", 1, "CCCC"))