import sys
import tkinter as tk
from tkinter import ttk
import argparse
import random
from pathlib import Path
from fireworks import LaunchPad
import datetime as dat
from monty.serialization import dumpfn, loadfn
from d3tales_fw.workflows.wf_writer import *


class GUI:

    """
    A graphical user interface for setting up and launching Molecular Dynamics (MD) calculations.

    Attributes:
    - window (tk.Tk): The main application window.
    - nameMatrix (list): List to store system names.
    - systemNamemat (list): List to store system names in a matrix.
    - systemsmilesmat (list): List to store SMILES codes for systems in a matrix.
    - entries (dict): Dictionary to store tkinter Entry widgets for user input.
    - smilesMatrix (list): List to store SMILES codes for all components.
    - typematrix (list): List to store the types of components (Solvent/Solute).
    - canvas (tk.Canvas): Canvas for scrollable frame.
    - scrollbar (ttk.Scrollbar): Scrollbar for canvas.
    - frame (ttk.Frame): Frame containing the widgets.
    - frame_id (int): ID of the frame in the canvas.
    - numsys (tk.Entry): Entry widget for the number of systems.
    - numsysL (tk.Label): Label for the number of systems.
    - numsysB (tk.Button): Button to create systems based on user input.
    - email (tk.Entry): Entry widget for user email.
    - emailLable (tk.Label): Label for user email.

    Methods:
    - on_frame_configure(event): Adjusts the canvas scroll region based on the frame size.
    - on_canvas_configure(event): Adjusts the canvas width based on the window size.
    - add_widgets(): Adds widgets (Entry, Label, Button) to the frame.
    - make(): Dynamically creates Entry widgets based on the number of systems entered by the user.
    - run(): Initiates MD calculations using user-inputted data.

    Note: The GUI is currently in development and new features may be added while the old ones are dropped
    """
    def __init__(self):
        self.window = tk.Tk()
        self.window.title('ASMD')
        self.window.geometry('850x650')
        self.nameMatrix = []
        self.systemNamemat = []
        self.systemsmilesmat = []
        self.entries = {}
        self.smilesMatrix = []
        self.typematrix = []
        self.canvas = tk.Canvas(self.window)
        self.canvas.grid(row=0, column=0, sticky=tk.W + tk.E + tk.N + tk.S)

        self.window.grid_rowconfigure(0, weight=1)
        self.window.grid_columnconfigure(0, weight=1)

        self.scrollbar = ttk.Scrollbar(self.window, command=self.canvas.yview)
        self.scrollbar.grid(row=0, column=1, sticky=tk.N + tk.S)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.frame = ttk.Frame(self.canvas)

        self.add_widgets()

        self.frame_id = self.canvas.create_window((0, 0), window=self.frame, anchor=tk.NW)

        self.frame.bind("<Configure>", self.on_frame_configure)
        self.canvas.bind("<Configure>", self.on_canvas_configure)

        self.window.mainloop()

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def on_canvas_configure(self, event):
        canvas_width = event.width
        self.canvas.itemconfig(self.frame_id, width=canvas_width)

    def add_widgets(self):
        self.numsys = tk.Entry(self.frame, fg="black", bg="white", width=10)
        self.numsysL = tk.Label(self.frame, text="number of systems")
        self.sumbit_button = tk.Button(self.frame, text="create sys", command= self.charge_titration_or_not)
        self.email = tk.Entry(self.frame, fg="black", bg="white", width=50)
        self.check_titration = tk.IntVar()
        self.charge_tittration_promt=tk.Checkbutton(self.frame, text=f"Is this a charge titration?",variable=self.check_titration)
        self.emailLable = tk.Label(self.frame, text="Email: ".ljust(20))
        self.emailLable.grid(row=2, column=0)
        self.email.grid(row=2, column=1)
        self.numsysL.grid(row=0, column=0)
        self.numsys.grid(row=0, column=1)
        self.sumbit_button.grid(row=1, column=0)
        self.charge_tittration_promt.grid(row=3, column=0)

    def charge_titration_or_not(self):
        if self.check_titration.get() !=0:
            self.charge_titartion_setup()
        else:
            self.make()
    def charge_titartion_setup(self):

        self.Titration_start= tk.Entry(self.frame, fg="black", bg="white", width=50)
        self.Titration_start_promt=tk.Label(self.frame, text="Titration start")

        self.Titration_finish= tk.Entry(self.frame, fg="black", bg="white", width=50)
        self.Titration_finish_promt=tk.Label(self.frame, text="Titration finish")

        self.Titration_steps= tk.Entry(self.frame, fg="black", bg="white", width=50)
        self.Titration_steps_promt=tk.Label(self.frame, text="delta_steps")

        self.charge_sumbit_button=tk.Button(self.frame, text="Next", command= self.charge_titration_maker)

        self.Titration_start.grid(row=4, column=1)
        self.Titration_start_promt.grid(row=4, column=0)

        self.Titration_finish.grid(row=5, column=1)
        self.Titration_finish_promt.grid(row=5, column=0)

        self.Titration_steps.grid(row=6, column=1)
        self.Titration_steps_promt.grid(row=6, column=0)

        self.charge_sumbit_button.grid(row=7, column=0)
    def charge_titration_maker(self):
        self.number_of_titrationsneeded= 1 + ((1*float(self.Titration_finish.get()) - float(self.Titration_start.get()))/float(self.Titration_steps.get()))
        self.iterator = int((round(self.number_of_titrationsneeded, 4)))
        self.titration_list = [float(self.Titration_finish.get()) - (i * float(self.Titration_steps.get())) for i in
                          range(self.iterator)]

        print(round(self.number_of_titrationsneeded,4))
        if round(self.number_of_titrationsneeded,4)%1 != 0.0:
            raise(ValueError("titration inputs are not valid"))

        for i in range(int(self.numsys.get())):
            global number
            global Tnumber
            number = int(self.numsys.get())
            Tnumber = int(self.numsys.get()) * self.iterator

            self.button = tk.Button(self.frame,
                                    text="Submit",
                                    width=5,
                                    height=1,
                                    bg="gray",
                                    fg="black", command=self.run)
            self.button.grid(row=17 + (17 * int(self.numsys.get()) - 1), column=0)
            self.SoluteSmiles = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.Density = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.xdim = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.ydim = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.zdim = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.xdimLabel = tk.Label(self.frame, text=f"x{i + 1} dimensions: ".ljust(20))
            self.ydimLabel = tk.Label(self.frame, text=f"y{i + 1} dimensions: ".ljust(20))
            self.zdimLabel = tk.Label(self.frame, text=f"z{i + 1} dimenisons: ".ljust(20))

            self.labelDensity = tk.Label(self.frame, text="Desity of solvent(M)".ljust(20))

            self.SoluteName = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.SolventSmiles = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.Solname = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.Concentrations = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.labelSolute = tk.Label(self.frame,
                                        text=f"Solute{i + 1} names (3 letters, if using multiple, seprate with commas): ".ljust(
                                            20))

            self.labelSolvent = tk.Label(self.frame, text=f"Solvent{i + 1} name (3 letters): ".ljust(20))

            self.labelConcentration = tk.Label(self.frame,
                                               text=f"Concentration(M) of the Solutes{i + 1}, seprate with commas: ".ljust(
                                                   20))

            self.lableSolventSmiles = tk.Label(self.frame, text=f"Solvent{i + 1} SMILES code:".ljust(20))

            self.labelSoluteSmiles = tk.Label(self.frame,
                                              text=f"Solute SMILES{i + 1} code(separate with commas): ".ljust(20))

            self.system = tk.Label(self.frame, text=f"System_{i + 1}-----------------------------------")

            self.chargeV = tk.IntVar()

            self.molarMass = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.molarlabel = tk.Label(self.frame, text=f"Molarmass{i + 1} of Solvents")

            self.chargeCheck = tk.Checkbutton(self.frame, text=f"Charge on the Solutes{i + 1}", variable=self.chargeV)
            self.chargeMatrix = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.chargeEntry = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.system.grid(row=10 + (13 * i), column=0, padx=50)
            self.labelSolvent.grid(row=11 + (13 * i), column=0)
            self.Solname.grid(row=11 + (13 * i), column=1)
            self.entries[f"solventname{i + 1}"] = self.Solname

            self.lableSolventSmiles.grid(row=12 + (13 * i), column=0)
            self.SolventSmiles.grid(row=12 + (13 * i), column=1)
            self.entries[f"solvetsmiles{i + 1}"] = self.SolventSmiles

            self.labelSolute.grid(row=13+ (13 * i), column=0)
            self.SoluteName.grid(row=13 + (13 * i), column=1)
            self.entries[f"solutename{i + 1}"] = self.SoluteName

            self.labelSoluteSmiles.grid(row=14 + (13 * i), column=0)
            self.SoluteSmiles.grid(row=14 + (13 * i), column=1)
            self.entries[f"solutesmiles{i + 1}"] = self.SoluteSmiles

            self.labelConcentration.grid(row=15 + (13 * i), column=0)
            self.Concentrations.grid(row=15 + (13 * i), column=1)
            self.entries[f"concentration{i + 1}"] = self.Concentrations

            self.labelDensity.grid(row=16 + (13 * i), column=0)
            self.Density.grid(row=16 + (13 * i), column=1)
            self.entries[f"Density{i + 1}"] = self.Density

            self.xdimLabel.grid(row=17 + (13 * i), column=0)
            self.xdim.grid(row=17 + (13 * i), column=1)
            self.entries[f"xdim{i + 1}"] = self.xdim

            self.ydimLabel.grid(row=18 + (13 * i), column=0)
            self.ydim.grid(row=18 + (13 * i), column=1)
            self.entries[f"ydim{i + 1}"] = self.ydim

            self.zdimLabel.grid(row=19 + (13 * i), column=0)
            self.zdim.grid(row=19 + (13 * i), column=1)
            self.entries[f"zdim{i + 1}"] = self.zdim

            self.chargeCheck.grid(row=20 + (13 * i), column=0)
            self.chargeMatrix.grid(row=20+ (13 * i), column=1)
            self.entries[f"charges{i + 1}"] = self.chargeMatrix

            self.molarMass.grid(row=21 + (13 * i), column=1)
            self.molarlabel.grid(row=21 + (13 * i), column=0)
            self.entries[f"molarmass{i + 1}"] = self.molarMass

        return




    def make(self):
        global number
        number = int(self.numsys.get())
        self.button = tk.Button(self.frame,
                                text="Submit",
                                width=5,
                                height=1,
                                bg="gray",
                                fg="black", command=self.run)
        self.button.grid(row=17 + (17 * number - 1), column=0)
        for i in range(number):
            self.SoluteSmiles = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.Density = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.xdim = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.ydim = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.zdim = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.xdimLabel = tk.Label(self.frame, text=f"x{i + 1} dimensions: ".ljust(20))
            self.ydimLabel = tk.Label(self.frame, text=f"y{i + 1} dimensions: ".ljust(20))
            self.zdimLabel = tk.Label(self.frame, text=f"z{i + 1} dimenisons: ".ljust(20))

            self.labelDensity = tk.Label(self.frame, text="Desity of solvent(M)".ljust(20))

            self.SoluteName = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.SolventSmiles = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.Solname = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.Concentrations = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.labelSolute = tk.Label(self.frame,
                                        text=f"Solute{i + 1} names (3 letters, if using multiple, seprate with commas): ".ljust(
                                            20))

            self.labelSolvent = tk.Label(self.frame, text=f"Solvent{i + 1} name (3 letters): ".ljust(20))

            self.labelConcentration = tk.Label(self.frame,
                                               text=f"Concentration(M) of the Solutes{i + 1}, seprate with commas: ".ljust(
                                                   20))

            self.lableSolventSmiles = tk.Label(self.frame, text=f"Solvent{i + 1} SMILES code:".ljust(20))

            self.labelSoluteSmiles = tk.Label(self.frame,
                                              text=f"Solute SMILES{i + 1} code(separate with commas): ".ljust(20))

            self.system = tk.Label(self.frame, text=f"System_{i + 1}-----------------------------------------------")

            self.chargeV = tk.IntVar()

            self.molarMass = tk.Entry(self.frame, fg="black", bg="white", width=50)
            self.molarlabel = tk.Label(self.frame, text=f"Molarmass{i + 1} of Solvents")

            self.chargeCheck = tk.Checkbutton(self.frame, text=f"Charge on the Solutes{i + 1}", variable=self.chargeV)
            self.chargeMatrix = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.chargeEntry = tk.Entry(self.frame, fg="black", bg="white", width=50)

            self.system.grid(row=5 + (13 * i), column=0, padx=50)
            self.labelSolvent.grid(row=6 + (13 * i), column=0)
            self.Solname.grid(row=6 + (13 * i), column=1)
            self.entries[f"solventname{i + 1}"] = self.Solname

            self.lableSolventSmiles.grid(row=7 + (13 * i), column=0)
            self.SolventSmiles.grid(row=7 + (13 * i), column=1)
            self.entries[f"solvetsmiles{i + 1}"] = self.SolventSmiles

            self.labelSolute.grid(row=8 + (13 * i), column=0)
            self.SoluteName.grid(row=8 + (13 * i), column=1)
            self.entries[f"solutename{i + 1}"] = self.SoluteName

            self.labelSoluteSmiles.grid(row=9 + (13 * i), column=0)
            self.SoluteSmiles.grid(row=9 + (13 * i), column=1)
            self.entries[f"solutesmiles{i + 1}"] = self.SoluteSmiles

            self.labelConcentration.grid(row=10 + (13 * i), column=0)
            self.Concentrations.grid(row=10 + (13 * i), column=1)
            self.entries[f"concentration{i + 1}"] = self.Concentrations

            self.labelDensity.grid(row=12 + (13 * i), column=0)
            self.Density.grid(row=12 + (13 * i), column=1)
            self.entries[f"Density{i + 1}"] = self.Density

            self.xdimLabel.grid(row=13 + (13 * i), column=0)
            self.xdim.grid(row=13 + (13 * i), column=1)
            self.entries[f"xdim{i + 1}"] = self.xdim

            self.ydimLabel.grid(row=14 + (13 * i), column=0)
            self.ydim.grid(row=14 + (13 * i), column=1)
            self.entries[f"ydim{i + 1}"] = self.ydim

            self.zdimLabel.grid(row=15 + (13 * i), column=0)
            self.zdim.grid(row=15 + (13 * i), column=1)
            self.entries[f"zdim{i + 1}"] = self.zdim

            self.chargeCheck.grid(row=16 + (13 * i), column=0)
            self.chargeMatrix.grid(row=16 + (13 * i), column=1)
            self.entries[f"charges{i + 1}"] = self.chargeMatrix

            self.molarMass.grid(row=17 + (13 * i), column=1)
            self.molarlabel.grid(row=17 + (13 * i), column=0)
            self.entries[f"molarmass{i + 1}"] = self.molarMass


    def run(self):
        BASE_DIR = Path(__file__).resolve().parent.parent

        parser = argparse.ArgumentParser(description='Launch MD calculations')
        parser.add_argument('-f', '--filename', type=str, help='filepath for a JSON nlp data file', default="")
        parser.add_argument('-p', '--priority', type=int, help='jobs priority', default=5)
        args = parser.parse_args()

        def populate_md_wf(**kwargs):

            lpad_file = os.path.join(BASE_DIR.parent, 'launch', 'md_launchpad.yaml')
            wf = d3tales_md_wf(**kwargs)
            info = LaunchPad().from_file(lpad_file).add_wf(wf)
            fw_id = list(info.values())[0]
            return fw_id

        if not os.path.isfile(args.filename):
            self.systems = []
            self.solute_titrants=[]
            self.solvent_titratnts=[]

            if self.check_titration.get() !=0:
                print("in the titration")

                for i in range(number):
                    self.subnamemat = []

                    print(" in loop 1")
                    self.solute_matrix = self.entries[f"solutename{i + 1}"].get().split(",")
                    print(self.solute_matrix)
                    solvent = self.entries[f"solventname{i + 1}"].get().strip()
                    self.typematrix.append("Solvent")
                    self.nameMatrix.append(solvent)
                    self.subnamemat.append(solvent)
                    print(self.titration_list)
                    for _ in self.titration_list:
                        print("in loop 2")
                        titration_name= solvent
                        self.solvent_titratnts.append(titration_name)

                    for iteams in self.solute_matrix:
                        print(self.solute_matrix)
                        print("in loop 3")
                        self.typematrix.append("Solute1")
                        self.nameMatrix.append(iteams.strip())
                        self.subnamemat.append(iteams.strip())

                        for _ in self.titration_list:
                            print("in loop 4")
                            self.solute_titrants.append(iteams.strip() + str(_))
                    self.systemNamemat.append(self.subnamemat)
                    print(f"this is the subname mat at {i+1} iteration : {self.subnamemat}")

                for i, iteams in enumerate(self.solvent_titratnts):
                    print("in loop 5")
                    self.systems.append(f"{iteams}_{self.solute_titrants[i]}")
                for _ in range(number):
                    self.submatsmiles = []
                    a = self.entries[f"solvetsmiles{_ + 1}"].get().strip()
                    self.submatsmiles.append(a)
                    self.smilesMatrix.append(a)
                    b = self.entries[f"solutesmiles{_ + 1}"].get().split(",")
                    for iteams in b:
                        self.submatsmiles.append(iteams.strip())
                        self.smilesMatrix.append(iteams.strip())
                    self.systemsmilesmat.append(self.submatsmiles)


                print(self.systems)
                print(self.nameMatrix)



            else:
                for _ in range(number):
                    self.subnamemat = []
                    string_to_append = ""
                    a = self.entries[f"solventname{_ + 1}"].get().strip()
                    self.typematrix.append("Solvent")
                    self.subnamemat.append(a)
                    self.nameMatrix.append(a)
                    string_to_append += f'{a}'
                    b = self.entries[f"solutename{_ + 1}"].get().split(",")
                    for iteams in b:
                        self.subnamemat.append(iteams.strip())
                        string_to_append += f'_{iteams}'
                        self.typematrix.append("Solute1")
                        self.nameMatrix.append(iteams.strip())
                    self.systems.append(string_to_append)
                    self.systemNamemat.append(self.subnamemat)
                for _ in range(number):
                    self.submatsmiles = []
                    a = self.entries[f"solvetsmiles{_ + 1}"].get().strip()
                    self.submatsmiles.append(a)
                    self.smilesMatrix.append(a)
                    b = self.entries[f"solutesmiles{_ + 1}"].get().split(",")
                    for iteams in b:
                        self.submatsmiles.append(iteams.strip())
                        self.smilesMatrix.append(iteams.strip())
                    self.systemsmilesmat.append(self.submatsmiles)
                    print(f"submat: {self.submatsmiles}")
                    print(f'systemmat: {self.systemsmilesmat}')



            if self.check_titration.get() !=0:

                number_sys = len(self.systems)
                key_dic = {}
                darte = (str(dat.datetime.now()).split()[0]).split("-")
                date = ""
                for iteams in darte:
                    date += f"_{str(iteams)}"

                titration_list = []
                titration_list[:] = self.titration_list
                titration_list.pop(self.titration_list.index(1.0))

                global number_of_titration
                number_of_titration = len(self.titration_list)
                outer_system=number_sys/(number_of_titration)
                print(int(outer_system))
                for j in range(int(outer_system)):
                    current_index=j*number_of_titration
                    print(current_index)
                    key_number=random.randint(1, 30000000)
                    key_dic[self.systems[current_index]] = key_number
                    for i in range(number_of_titration-1):
                        key_dic[self.systems[current_index+(i+1)]] = f'{key_number}_{titration_list[i]}'

                print(f"the key looks like this:{key_dic} ")

                md_kwargs = {"date_sumbit":date,
                    "smiles_list": self.smilesMatrix,
                    "name_list": self.nameMatrix,  # "cons": 0,
                    "type_list": self.typematrix,
                    "dir": "/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs",
                    "num_systems": f"{number_sys}", "titartion_list":self.titration_list, "populate_name": "MD_FIREWORK", "key_dic": key_dic,"is_titration":True}
                print(f"The titrationlist:{self.titration_list}")
            else:
                number_sys = number

                key_dic = {}
                darte = (str(dat.datetime.now()).split()[0]).split("-")
                date = ""
                for iteams in darte:
                    date += f"_{str(iteams)}"
                for _ in range(number_sys):
                     key_dic[self.systems[_]] = random.randint(1, 3000000000)




                md_kwargs = {"date_sumbit": date,
                             "smiles_list": self.smilesMatrix,
                             "name_list": self.nameMatrix,  # "cons": 0,
                             "type_list": self.typematrix,
                             "dir": "/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs",
                             "num_systems": f"{number_sys}", "populate_name": "MD_FIREWORK", "key_dic": key_dic,"is_titration":False}
            if self.check_titration.get() !=0:

                self.iterator_for_wf = Tnumber

            else:


                self.iterator_for_wf=  number

            for _ in range(self.iterator_for_wf):
                md_kwargs[f"WF_name{_ + 1}"] = self.systems[_]
            if self.check_titration.get() !=0:
                for a in range(number):
                    md_kwargs[f"Average_den{list(key_dic.values())[a*number_of_titration]}"] = []
            if self.check_titration.get() ==0:
                for a in range(number):
                    md_kwargs[f"Average_den{list(key_dic.values())[a]}"] = []
            for _ in range(number):
                md_kwargs[f"den{_ + 1}"] = float(self.entries[f'Density{_ + 1}'].get().strip())
            for _ in range(number):
                try:
                    md_kwargs[f'MM{_ + 1}'] = float(self.entries[f'molarmass{_ + 1}'].get().strip())
                except ValueError:
                    md_kwargs[f'MM{_ + 1}'] = int(self.entries[f'molarmass{_ + 1}'].get().strip())
            for _ in range(number):
                md_kwargs[f'x{_ + 1}'] = float(self.entries[f'xdim{_ + 1}'].get().strip())
                md_kwargs[f'y{_ + 1}'] = float(self.entries[f'ydim{_ + 1}'].get().strip())
                md_kwargs[f'z{_ + 1}'] = float(self.entries[f'zdim{_ + 1}'].get().strip())
            for _ in range(number):
                conamt = []
                preprocessd = self.entries[f'concentration{_ + 1}'].get().split(",")
                for iteams in preprocessd:
                    conamt.append(iteams.strip())
                md_kwargs[f'conmatrix{_ + 1}'] = conamt
            index = 0
            for j in self.systemNamemat:

                solventmat = []
                solutemat = j[1:]
                solventmat.append(j[0])
                md_kwargs[f'solvent_name{index + 1}'] = solventmat
                md_kwargs[f'solute_name{index + 1}'] = solutemat
                index += 1
            index=0
            for j in self.systemsmilesmat:
                solventmat = []
                solutemat = j[1:]
                solventmat.append(j[0])
                md_kwargs[f'solvent_smiles{index + 1}'] = solventmat
                md_kwargs[f'solute_smiles{index + 1}'] = solutemat
                print(index,solutemat,solventmat)
                index +=1
            self.charge = []

            for i in range(number):
                self.charge.append(
                    self.entries[f"charges{i + 1}"].get() if len(self.entries[f"charges{i + 1}"].get()) != 0 else 0)
            md_kwargs["charge_list"] = self.charge
            print(md_kwargs)

            # md_kwargs = {
            #     "smiles_list": ["CO", "CO", "CCO", "CCO"],
            #     "con_list": [0, 0, 0, 0], "WF_name1": "Test_run1", "WF_name2": "Test_run2",
            #     "name_list": ["wa1ter", "wa2ter", "methane1", "ethane2"], "cons": 0,
            #     "type_list": ["Solvent", "Solvent", "Solute1", "Solute1"],
            #     "dir": "/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs", "solvent_name1": ["wa1ter"],
            #     "solute_name1": ["methane"], "solvent_name2": ["wa2ter"], "solute_name2": ["ethane"],
            #     "solvent_smiles1": ["CO"], "solute_smiles1": ["CCO"], "solvent_smiles2": ["CO"],
            #     "solute_smiles2": ["CCO"], "x1": 10, "y1": 10, "z1": 10, "conmatrix1": ["0.1"], "den1": "10",
            #     "MM1": "5", "x2": 10, "y2": 10, "z2": 10, "conmatrix2": ["0.1"], "den2": "10", "MM2": "5",
            #     "num_systems": "2", "populate_name": "test", "key_dic": key_dic}
            # md_kwargs = {
            # "smiles_list": ["O","C"],
            # "con_list": [0,0], "WF_name1":"Test_run1","WF_name2":"Test_run2",  "name_list":["water","methane1"], "cons":0, "type_list":["Solvent","Solute1"], "dir":"/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs", "solvent_name1":["water"], "solute_name1":["methane"],"solvent_name2":["wa2ter"], "solute_name2":["ethane"], "x1":10, "y1":10, "z1":10, "conmatrix1":["0.1"], "den1":"10", "MM1":"10", "x2":10, "y2":10, "z2":10, "conmatrix2":["0.1"], "den2":"10", "MM2":"18.012","num_systems":"1", "populate_name":"test","key_dic":key_dic}
            all_ids = {"test_md_fw": populate_md_wf(**md_kwargs)}
            sys.exit()

        else:
            all_md_data = loadfn(args.filename)
            all_ids = {}
            for mol_name, md_kwargs in all_md_data.items():
                all_ids[mol_name] = populate_md_wf(path="/project/cmri235_uksr/shasanka_conda_boss/d3tales_fw/parameters/md_gaus_parameter_file.json",**md_kwargs)


GUI()
