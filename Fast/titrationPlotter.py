import matplotlib.pyplot as plt

def TitrationPlotter( Average_simulation_density_dic, key):
    difrence=[float(i) for i in Average_simulation_density_dic.values() ]
    x_axis=[i for i in Average_simulation_density_dic]

    plt.plot(x_axis,difrence,'o')
    plt.xlabel("System")
    plt.ylabel("Simulation density")
    plt.title(f"Titration Result{key}")
    plt.savefig(f"/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs/Output{key}/TitrationGraph_{key}.pdf")
