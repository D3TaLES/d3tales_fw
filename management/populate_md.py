"""
Parse and push data from NLP parsing and submit calculations for extracted molecules
"""
import argparse
from pathlib import Path
from fireworks import LaunchPad
from monty.serialization import dumpfn
from d3tales_fw.workflows.wf_writer import *

BASE_DIR = Path(__file__).resolve().parent.parent

parser = argparse.ArgumentParser(description='Launch MD calculations')
parser.add_argument('filename', metavar='filename', type=str, help='filepath for a JSON nlp data file', default="")
parser.add_argument('-p', '--priority', type=int, help='jobs priority', default=5)
args = parser.parse_args()


def populate_md_wf(**kwargs):
    lpad_file = os.path.join(BASE_DIR.parent, 'config', 'md_launchpad.yaml')
    wf = d3tales_md_wf(**kwargs)
    info = LaunchPad().from_file(lpad_file).add_wf(wf)
    fw_id = list(info.values())[0]
    return fw_id


if not os.path.isfile(args.filename):
    md_kwargs = {

    }
    all_ids = {"test_md_fw": populate_md_wf(**md_kwargs)}

else:
    all_md_data = loadfn(args.filename)
    all_ids = {}
    for mol_name, md_kwargs in all_md_data.items():
        all_ids[mol_name] = populate_md_wf(**md_kwargs)

dumpfn(all_ids, "md_data/md_{}_molIDs.json".format(args.filename.split("/")[-1]), indent=2)
