import subprocess

# Command to run the shell script and gmx_mpi together
command = 'source /project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/Fast/ASMD1/Input/gromacs_gpu.sh && gmx_mpi'

# Run the command in a subprocess with shell=True
subprocess.run(command, shell=True)

