import os
import subprocess
from sys import argv
from pymatgen.io.gaussian import GaussianOutput


if __name__ == "__main__":
    log_path = argv[1]

    home = '/'.join(log_path.split('/')[:-1])
    os.chdir(home)
    log_file = log_path.split('/')[-1]
    name = log_file.split('.')[0]
    chk_file = name+'.chk'
    gjf_file = name+'.gjf'
    if not os.path.isfile(chk_file):
        chk_file = [f for f in os.listdir(home) if '.chk' in f][0]
    subprocess.call('mv {} old_{}'.format(chk_file, chk_file), shell=True)

    gout = GaussianOutput(log_file)
    gin = gout.to_input()
    gin.link0_parameters.update({'%oldchk': 'old_{}.chk'.format(name), 'Geom': 'AllCheck'})
    gin.route_parameters.update({'Geom': 'AllCheck'})

    gin.write_file(gjf_file)
    print("cd ", home)
    print("g16", gjf_file)
