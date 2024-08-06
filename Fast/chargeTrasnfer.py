import subprocess
from new_dft import get_charge
import os
class trans:
    def __init__(self,name,smiles,key,re,di,charg,Titration, mult=None):
        charges = get_charge(name, charg, smiles,direc=di,mul=mult)
        self.name = name
        self.titrat=Titration

        self.dir = di
        self.smiles = smiles
        self.file_to_make = os.path.join(self.dir, f"InputGrofiles{key}", f"{name}f.itp")
        self.orignal_gro_outputfile = os.path.join(self.dir, name, f"{name}.gmx.itp")
        chagreMatrix = []
        if re == 1:
            print("it is ran")

            subprocess.run(f"touch {self.file_to_make}", shell=True)
            with open(f'{self.orignal_gro_outputfile}','r') as gmx,   open(f'{self.file_to_make}','a') as itp:
                gmx_lines= gmx.readlines()
                if len(gmx_lines) ==0:
                    print("gmx is empty, try again")
                    subprocess.run(["ls"], shell= True)
                index_GMXESP=0
                index_lastGMXESP=0
                for line in gmx_lines:
                    if line.strip() =="[ atoms ]":
                        index_GMXESP = gmx_lines.index(line) + 2
                        break
                for line in gmx_lines:
                    if line.strip()=="[ bonds ]":
                        index_lastGMXESP = gmx_lines.index(line) -2

                for a in charges:
                    charge_int= float(a)*self.titrat
                    charge = f'{charge_int:.4f}'.rjust(11)

                    chagreMatrix.append(charge)
                for char in range(index_lastGMXESP-index_GMXESP+1):
                    chargeIndexGmx= gmx_lines[index_GMXESP-1].index("charge")-6

                    gmx_lines[index_GMXESP+char]= gmx_lines[index_GMXESP+char][0:chargeIndexGmx] + str(chagreMatrix[char]) + gmx_lines[index_GMXESP+char][chargeIndexGmx+11:]
                for line in gmx_lines:
                    itp.writelines(line)
            print("Itp file is made")
        else:
            print("else is ran")
            subprocess.run(f"mv {self.orignal_gro_outputfile} {self.file_to_make}",shell=True)
