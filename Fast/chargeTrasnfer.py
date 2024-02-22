import subprocess
class trans:
    def __init__(self,name,smiles,key,re,di,dft):

        self.name = name

        self.dir = di
        self.dft=dft
        self.smiles = smiles
        chagreMatrix = []
        if re == 1:
            print("it is ran")
            subprocess.run(f"touch {self.dir}/InputGrofiles{key}/{name}f.itp", shell=True)
            with open (f'{self.dft}/{self.smiles}/gaussian/gas_phase/opt/','r') as log, open(f'{self.dir}/{name}/{name}.gmx.itp','r') as gmx,   open(f'{self.dir}/InputGrofiles{key}/{name}f.itp','a') as itp:
                lines= log.readlines()
                gmx_lines= gmx.readlines()
                if len(gmx_lines) ==0:
                    print("gmx is empty, try again")
                    subprocess.run(["ls"], shell= True)
                index_ESP=0
                index_lastEsp=0
                index_GMXESP=0
                index_lastGMXESP=0
                for line in gmx_lines:
                    if line.strip() =="[ atoms ]":
                        index_GMXESP = gmx_lines.index(line) + 2
                        break
                for line in gmx_lines:
                    if line.strip()=="[ bonds ]":
                        index_lastGMXESP = gmx_lines.index(line) -2
                for line in lines:
                    if line.strip() == "ESP charges:":
                        index_ESP= lines.index(line) +2
                        break
                for line in lines:
                    if line.strip()[0:18] =="Sum of ESP charges":
                        index_lastEsp= lines.index(line) -1
                        break

                for a in range(index_lastEsp - index_ESP + 1):
                    charge_int= float(lines[index_ESP+a].split()[2])
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
            subprocess.run(f"mv {self.dir}/{name}/{name}.gmx.itp {self.dir}/InputGrofiles{key}/{name}f.itp",shell=True)
