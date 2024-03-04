systems = []
subnamemat = []
solute_titrants = []
solvent_titratnts = []
check_titration=1
number = 2
entries={f"solventname1":"O","solutename1":"C", f"solventname2":"CC","solutename2":"CCCCCC" }
typematrix=[]
nameMatrix=[]
titration_list=[1.0,0.9,0.8]
systemNamemat=[]
smilesMatrix=[]
systemsmilesmat=[]

if check_titration != 0:
    print("in the titration")

    for i in range(number):
        print(" in loop 1")
        solute_matrix = entries[f"solutename{i + 1}"].split(",")
        print(solute_matrix)
        solvent = entries[f"solventname{i + 1}"].strip()
        typematrix.append("Solvent")
        nameMatrix.append(solvent)
        subnamemat.append(solvent)
        print(titration_list)
        for _ in titration_list:
            print("in loop 2")
            titration_name = solvent
            solvent_titratnts.append(titration_name)

        for iteams in solute_matrix:
            print(solute_matrix)
            print("in loop 3")
            typematrix.append("Solute1")
            nameMatrix.append(iteams.strip())
            subnamemat.append(iteams.strip())
        systemNamemat.append(subnamemat)
        print(f"this is the subname mat at {i + 1} iteration : {subnamemat}")

        for _ in titration_list:
            print("in loop 4")
            solute_titrants.append(iteams.strip() + str(_))
    for i, iteams in enumerate(solvent_titratnts):
        print("in loop 5")
        systems.append(f"{iteams}_{solute_titrants[i]}")
    for _ in range(number):
        submatsmiles = []
        a = entries[f"solvetsmiles{_ + 1}"]().strip()
        submatsmiles.append(a)
        smilesMatrix.append(a)
        b = entries[f"solutesmiles{_ + 1}"]().split(",")
        for iteams in b:
            submatsmiles.append(iteams.strip())
            smilesMatrix.append(iteams.strip())
        systemsmilesmat.append(submatsmiles)

    print(systems)
    print(nameMatrix)