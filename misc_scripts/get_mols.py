import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Gather data file directories for certain molecules ')
    parser.add_argument('-i', '--id_list', type=str, help='comma seperated list of ids')
    parser.add_argument('-r', '--rclone', action='store_true', help='use rclone to search cloud storage')
    parser.add_argument('--hussein_mols', action='store_true', help='all Hussein mols')
    parser.add_argument('--hussein_whol', action='store_true', help='all Hussein whole mols')
    parser.add_argument('--dihedral', action='store_true', help='all dihedral ring mols')
    parser.add_argument('--huckaba_mols', action='store_true', help='all Huckaba mols')
    parser.add_argument('--database_demo', action='store_true', help='all database demo mols')
    parser.add_argument('--sutton_repliation', action='store_true', help='all Sutton replication mols')
    args = parser.parse_args()

    id_list = []
    if args.id_list:
        id_list.extend(args.ids_list.split(','))
    if args.hussein_mols:
        id_list.extend(['06ESQE', '06ASFG', '05MOUY', '05JCLA', '05YZEU', '05NBNL', '05APQE', '05RPNH', '05POEN', '05XJXO', '05JNUA', '05GDCT', '05ZVEH', '05GMFA', '05NVRP', '05TSYK', '05RCMG', '90PSTI', '90TRQT', '05VXVT', '05DUQU']) # all Hussein mols
    if args.hussein_whol:
        id_list.extend(['05JNUA', '05RPNH', '05JCLA', '05NBNL', '05ACSI', '05KMLD', '05GDCT', '05ZURH', '05VXVT', '90PSTI', '05GMFA', '05TSYK']) #Hussein mols whole
    if args.dihedral:
        id_list.extend(['06YIEK', '06XZRY', '80IKJA', '06ASFG', '80GSXO', '06XWRM', '90GOIB', '90XPEY', '06XYPW', '06SGLH', '06ULAN', '06HEDT', '06INFH', '06FWHE']) # dihedral rings
    if args.huckaba_mols:
        id_list.extend(['10FXVB', '10ZMEC', '10ADWR', '10CGFY', '10IYDE', '10LYCW', '10NUOE', '10CLVI', '10VRYW', '10SBZT']) # Huckaba mols
    if args.database_demo:
        id_list.extend(['90ZYYG', '06ASFG', '06RACX', '06CCVA']) # database_demo
    if args.sutton_repliation:
        id_list.extend(['06LDFR', '06IMCA', '06KKFK']) # sutton replication

    home = os.getcwd()
    copy_dir = os.path.join(home, 'stored', 'gather_dir')
    for _id in id_list:
        print(_id)
        data_dir = os.path.join(home, 'd3tales', _id)
        if args.rclone:
            os.system("rclone copy onedrive:D3talesStorage/d3tales/{} {}/{} --bwlimit 8.6M -vv --s3-upload-cutoff=0".format(_id, copy_dir, _id))
            os.system("rclone copy onedrive:D3talesStorage/new/{} {}/{} --bwlimit 8.6M -vv --s3-upload-cutoff=0".format(_id, copy_dir, _id))
        os.system("cv -r {} {}".format(data_dir, copy_dir))
