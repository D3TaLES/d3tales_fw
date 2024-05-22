import subprocess


class titration:
    def __init__(self, titration_matrix, key, outputdir, solvent, solute, intial, path_inital):
        titration_list = []
        titration_list[:] = titration_matrix
        titration_list.pop(titration_matrix.index(1.0))
        print(intial)
        intial = True

        if not intial:
            for i in range(len(titration_list)):
                chargeMatrix = []
                subprocess.run(
                    [f"cp -r {outputdir}/InputGrofiles{key} {outputdir}/InputGrofiles{key}_{titration_list[i]}"],
                    shell=True)
                with open(
                        f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{solvent[:3]}_Solvent.itp', 'r') as gmx:
                    lines = gmx.readlines()

                    index_GMXESP = 0
                    index_lastGMXESP = 0
                    for line in lines:
                        if line.strip() == "[ atoms ]":
                            index_GMXESP = lines.index(line) + 2
                            break
                    for line in lines:
                        if line.strip() == "[ bonds ]":
                            index_lastGMXESP = lines.index(line) - 2

                    for a in range(index_lastGMXESP - index_GMXESP + 1):
                        charge_int = float(lines[index_GMXESP + a].split()[6]) * float(titration_list[i])
                        charge = f'{charge_int:.4f}'.rjust(11)

                        chargeMatrix.append(charge)
                    for char in range(index_lastGMXESP - index_GMXESP + 1):
                        chargeIndexGmx = lines[index_GMXESP - 1].index("charge") - 6

                        lines[index_GMXESP + char] = lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                            chargeMatrix[char]) + lines[index_GMXESP + char][chargeIndexGmx + 11:]
                    subprocess.run([f"rm {outputdir}/InputGrofiles{key}_{titration_list[i]}/{solvent[:3]}_Solvent.itp"],
                                   shell=True)

                with open(f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{solvent[:3]}_Solvent.itp', 'a') as itp:
                    for line in lines:
                        itp.writelines(line)

                for iteams in solute:
                    chargeMatrix = []
                    with open(
                            f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{iteams[:3]}_Solute1.itp', 'r') as gmx:
                        lines = gmx.readlines()

                        index_GMXESP = 0
                        index_lastGMXESP = 0
                        for line in lines:
                            if line.strip() == "[ atoms ]":
                                index_GMXESP = lines.index(line) + 2
                                break
                        for line in lines:
                            if line.strip() == "[ bonds ]":
                                index_lastGMXESP = lines.index(line) - 2

                        for a in range(index_lastGMXESP - index_GMXESP + 1):
                            charge_int = float(lines[index_GMXESP + a].split()[6]) * float(titration_list[i])
                            charge = f'{charge_int:.4f}'.rjust(11)

                            chargeMatrix.append(charge)
                        for char in range(index_lastGMXESP - index_GMXESP + 1):
                            chargeIndexGmx = lines[index_GMXESP - 1].index("charge") - 6

                            lines[index_GMXESP + char] = lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                                chargeMatrix[char]) + lines[index_GMXESP + char][chargeIndexGmx + 11:]
                        subprocess.run(
                            [f"rm {outputdir}/InputGrofiles{key}_{titration_list[i]}/{iteams[:3]}_Solute1.itp"],
                            shell=True)

                    with open(f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{iteams[:3]}_Solute1.itp',
                              'a') as itp:
                        for line in lines:
                            itp.writelines(line)

        else:
            for i in range(len(titration_list)):
                chargeMatrix = []
                subprocess.run(
                    [f"cp -r {outputdir}/InputGrofiles{key} {outputdir}/InputGrofiles{key}_{titration_list[i]}"],
                    shell=True)
                with open(
                        f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{solvent}.itp', 'r') as gmx:
                    lines = gmx.readlines()

                    index_GMXESP = 0
                    index_lastGMXESP = 0
                    for line in lines:
                        if line.strip() == "[ atoms ]" or line.strip() == "[  atoms  ]":
                            index_GMXESP = lines.index(line) + 2
                            break
                    for line in lines:
                        if line.strip() == "[ bonds ]" or line.strip() == "[  bonds  ]":
                            index_lastGMXESP = lines.index(line) - 2

                    for a in range(index_lastGMXESP - index_GMXESP + 1):
                        charge_int = float(lines[index_GMXESP + a].split()[6]) * float(titration_list[i])
                        charge = f'{charge_int:.4f}'.rjust(11)

                        chargeMatrix.append(charge)
                    for char in range(index_lastGMXESP - index_GMXESP + 1):
                        chargeIndexGmx = lines[index_GMXESP - 1].index("charge") - 6

                        lines[index_GMXESP + char] = lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                            chargeMatrix[char]) + lines[index_GMXESP + char][chargeIndexGmx + 11:]
                    subprocess.run([f"rm {outputdir}/InputGrofiles{key}_{titration_list[i]}/{solvent}.itp"],
                                   shell=True)

                with open(f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{solvent}.itp', 'a') as itp:
                    for line in lines:
                        itp.writelines(line)

                for iteams in solute:
                    chargeMatrix = []
                    with open(
                            f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{iteams}.itp', 'r') as gmx:
                        lines = gmx.readlines()

                        index_GMXESP = 0
                        index_lastGMXESP = 0
                        for line in lines:
                            if line.strip() == "[ atoms ]" or line.strip() == "[  atoms  ]":
                                index_GMXESP = lines.index(line) + 2
                                break
                        for line in lines:
                            if line.strip() == "[ bonds ]" or line.strip() == "[  bonds  ]":
                                index_lastGMXESP = lines.index(line) - 2

                        for a in range(index_lastGMXESP - index_GMXESP + 1):
                            print(lines[index_GMXESP + a].split()[6])
                            print(f"the titration list {titration_list[i]}")
                            charge_int = float(lines[index_GMXESP + a].split()[6]) * float(titration_list[i])
                            charge = f'{charge_int:.4f}'.rjust(11)

                            chargeMatrix.append(charge)
                        for char in range(index_lastGMXESP - index_GMXESP + 1):
                            chargeIndexGmx = lines[index_GMXESP - 1].index("charge") - 6

                            lines[index_GMXESP + char] = lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                                chargeMatrix[char]) + lines[index_GMXESP + char][chargeIndexGmx + 11:]
                        subprocess.run(
                            [f"rm {outputdir}/InputGrofiles{key}_{titration_list[i]}/{iteams}.itp"],
                            shell=True)

                    with open(f'{outputdir}/InputGrofiles{key}_{titration_list[i]}/{iteams}.itp',
                              'a') as itp:
                        for line in lines:
                            itp.writelines(line)