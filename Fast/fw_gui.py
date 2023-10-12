import tkinter as tk
from tkinter import ttk
import argparse
import random
from pathlib import Path
from fireworks import LaunchPad
from monty.serialization import dumpfn, loadfn
from d3tales_fw.workflows.wf_writer import *

class GUI:
    def __init__(self):
        self.window = tk.Tk()
        self.window.title('ASMD')
        self.window.geometry('850x650')
        self.nameMatrix=[]
        self.systemNamemat=[]
        self.systemsmilesmat=[]
        self.entries={}
        self.smilesMatrix=[]
        self.canvas = tk.Canvas(self.window)
        self.canvas.grid(row=0, column=0, sticky=tk.W + tk.E + tk.N + tk.S)

        self.window.grid_rowconfigure(0, weight=1)
        self.window.grid_columnconfigure(0, weight=1)

        self.scrollbar = ttk.Scrollbar(self.window, command=self.canvas.yview)
        self.scrollbar.grid(row=0, column=1, sticky=tk.N + tk.S)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.frame = ttk.Frame(self.canvas)

        self.add_widgets()

        self.frame_id= self.canvas.create_window((0, 0), window=self.frame, anchor=tk.NW)

        self.frame.bind("<Configure>", self.on_frame_configure)
        self.canvas.bind("<Configure>", self.on_canvas_configure)

        self.window.mainloop()

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def on_canvas_configure(self, event):
        canvas_width = event.width
        self.canvas.itemconfig( self.frame_id,width=canvas_width)

    def add_widgets(self):
        self.numsys = tk.Entry(self.frame, fg="black", bg="white", width=10)
        self.numsysL = tk.Label(self.frame, text="number of systems")
        self.numsysB = tk.Button(self.frame, text="create sys", command=self.make)
        self.email = tk.Entry(self.frame, fg="black", bg="white", width=50)

        self.emailLable = tk.Label(self.frame, text="Email: ".ljust(20))
        self.emailLable.grid(row=2, column=0)
        self.email.grid(row=2, column=1)
        self.numsysL.grid(row=0, column=0)
        self.numsys.grid(row=0, column=1)
        self.numsysB.grid(row=1, column=0)

    def make(self):
        global number
        number = int(self.numsys.get())
        self.button = tk.Button(self.frame,
                                text="Submit",
                                width=5,
                                height=1,
                                bg="gray",
                                fg="black", command=self.run)
        self.button.grid(row=17+(17*number-1), column=0)
        for i in range(number):
                    self.SoluteSmiles = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.Density = tk.Entry(self.frame, fg="black", bg="white", width=50)
                    self.xdim = tk.Entry(self.frame, fg="black", bg="white", width=50)
                    self.ydim = tk.Entry(self.frame, fg="black", bg="white", width=50)
                    self.zdim = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.xdimLabel = tk.Label(self.frame, text=f"x{i+1} dimensions: ".ljust(20))
                    self.ydimLabel = tk.Label(self.frame, text=f"y{i+1} dimensions: ".ljust(20))
                    self.zdimLabel = tk.Label(self.frame, text=f"z{i+1} dimenisons: ".ljust(20))


                    self.labelDensity = tk.Label(self.frame, text="Desity of solvent(M)".ljust(20))

                    self.SoluteName = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.SolventSmiles = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.Solname = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.Concentrations = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.labelSolute = tk.Label(self.frame,
                                                text=f"Solute{i+1} names (3 letters, if using multiple, seprate with commas): ".ljust(
                                                    20))

                    self.labelSolvent = tk.Label(self.frame, text=f"Solvent{i+1} name (3 letters): ".ljust(20))






                    self.labelConcentration = tk.Label(self.frame,
                                                       text=f"Concentration(M) of the Solutes{i+1}, seprate with commas: ".ljust(20))

                    self.lableSolventSmiles = tk.Label(self.frame, text=f"Solvent{i+1} SMILES code:".ljust(20))

                    self.labelSoluteSmiles = tk.Label(self.frame, text=f"Solute SMILES{i+1} code(separate with commas): ".ljust(20))


                    self.system = tk.Label(self.frame, text=f"System_{i+1}-----------------------------------------------")

                    self.chargeV = tk.IntVar()


                    self.molarMass = tk.Entry(self.frame, fg="black", bg="white", width=50)
                    self.molarlabel = tk.Label(self.frame, text=f"Molarmass{i+1} of Solvents")

                    self.chargeCheck = tk.Checkbutton(self.frame, text=f"Charge on the Solutes{i+1}", variable=self.chargeV)
                    self.chargeMatrix = tk.Entry(self.frame, fg="black", bg="white", width=50)

                    self.chargeEntry = tk.Entry(self.frame, fg="black", bg="white", width=50)




                    self.system.grid(row=5+(13*i), column=0, padx=50)
                    self.labelSolvent.grid(row=6+(13*i), column=0)
                    self.Solname.grid(row=6+(13*i), column=1)
                    self.entries[f"solventname{i+1}"]=self.Solname

                    self.lableSolventSmiles.grid(row=7+(13*i), column=0)
                    self.SolventSmiles.grid(row=7+(13*i), column=1)
                    self.entries[f"solvetsmiles{i+1}"]=self.SolventSmiles

                    self.labelSolute.grid(row=8+(13*i), column=0)
                    self.SoluteName.grid(row=8+(13*i), column=1)
                    self.entries[f"solutename{i + 1}"]=self.SoluteName


                    self.labelSoluteSmiles.grid(row=9+(13*i), column=0)
                    self.SoluteSmiles.grid(row=9+(13*i), column=1)
                    self.entries[f"solutesmiles{i + 1}"]= self.SoluteSmiles

                    self.labelConcentration.grid(row=10+(13*i), column=0)
                    self.Concentrations.grid(row=10+(13*i), column=1)
                    self.entries[f"concentration{i + 1}"] = self.Concentrations


                    self.labelDensity.grid(row=12+(13*i), column=0)
                    self.Density.grid(row=12+(13*i), column=1)
                    self.entries[f"Density{i + 1}"] = self.Density


                    self.xdimLabel.grid(row=13+(13*i), column=0)
                    self.xdim.grid(row=13+(13*i), column=1)
                    self.entries[f"xdim{i + 1}"] = self.xdim

                    self.ydimLabel.grid(row=14+(13*i), column=0)
                    self.ydim.grid(row=14+(13*i), column=1)
                    self.entries[f"ydim{i + 1}"] = self.ydim


                    self.zdimLabel.grid(row=15+(13*i), column=0)
                    self.zdim.grid(row=15+(13*i), column=1)
                    self.entries[f"zdim{i + 1}"] = self.zdim


                    self.chargeCheck.grid(row=16+(13*i), column=0)
                    self.chargeMatrix.grid(row=16+(13*i), column=1)
                    self.entries[f"charges{i + 1}"] = self.chargeMatrix


                    self.molarMass.grid(row=17+(13*i), column=1)
                    self.molarlabel.grid(row=17+(13*i), column=0)
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

            for _ in range (number):
                self.subnamemat = []

                string_to_append=""
                a=self.entries[f"solventname{_+1}"].get().strip()
                self.subnamemat.append(a)
                self.nameMatrix.append(a)
                string_to_append+=f'{a}'
                b=self.entries[f"solutename{_ + 1}"].get().split(",")
                for iteams in b:
                    self.subnamemat.append(iteams.strip())
                    string_to_append += f'_{iteams}'
                    self.nameMatrix.append(iteams.strip())
                self.systems.append(string_to_append)
                self.systemNamemat.append(self.subnamemat)
            for _ in range(number):
                self.submatsmiles = []
                a = self.entries[f"solvetsmiles{_ + 1}"].get().strip()
                self.submatsmiles.append(a)
                self.smilesMatrix.append(a)
                b=self.entries[f"solutesmiles{_ + 1}"].get().split(",")
                for iteams in b:
                    self.submatsmiles.append(iteams.strip())
                    self.smilesMatrix.append(iteams.strip())
                self.systemsmilesmat.append(self.submatsmiles)

            number_sys=number

            key_dic = {}
            for _ in range(number_sys):  ###chnage this later to a var
                key_dic[self.systems[_]] = random.randint(1, 30000000)

            md_kwargs = {
                "smiles_list": self.smilesMatrix,
                "con_list": [0, 0, 0, 0],
                "name_list": self.nameMatrix, "cons": 0,
                "type_list": ["Solvent", "Solute1", "Solvent", "Solute1"],
                "dir": "/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs", "solvent_name1": ["wa1ter"],
                "solute_name1": ["methane"], "solvent_name2": ["wa2ter"], "solute_name2": ["ethane"],
                "solvent_smiles1": ["CO"], "solute_smiles1": ["CCO"], "solvent_smiles2": ["CO"],
                "solute_smiles2": ["CCO"], "conmatrix1": ["0.1"],
               "conmatrix2": ["0.1"],
                "num_systems": f"{number_sys}", "populate_name": "test", "key_dic": key_dic}
            for _ in range(number_sys):
                md_kwargs[f"WF_name{_+1}"]= self.systems[_]
            for _ in range(number_sys):
                md_kwargs[f"den{_+1}"] = float(self.entries[f'Density{_+1}'].get().strip())
            for _ in range(number_sys):
                md_kwargs[f'MM{_+1}']=float(self.entries[f'molarmass{_ + 1}'].get().strip())
            for _ in range(number_sys):
                md_kwargs[f'x{_+1}']=float(self.entries[f'x{_ + 1}'].get().strip())
                md_kwargs[f'y{_ + 1}'] = float(self.entries[f'y{_ + 1}'].get().strip())
                md_kwargs[f'z{_ + 1}'] = float(self.entries[f'z{_ + 1}'].get().strip())
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

        else:
            all_md_data = loadfn(args.filename)
            all_ids = {}
            for mol_name, md_kwargs in all_md_data.items():
                all_ids[mol_name] = populate_md_wf(**md_kwargs)

GUI()
