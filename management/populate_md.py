"""
Parse and push data from NLP parsing and submit calculations for extracted molecules
"""
import argparse
import os.path

from monty.serialization import loadfn, dumpfn
from d3tales_api.D3database.d3database import *
from d3tales_api.Processors.back2front import *
from d3tales_api.Processors.d3tales_parser import *
from d3tales_fw.management.populate_lpad import populate_d3tales_lpad

parser = argparse.ArgumentParser(description='Launch MD calculations')
parser.add_argument('filename', metavar='filename', type=str, help='filepath for a JSON nlp data file', default="")
parser.add_argument('-p', '--priority', type=int, help='jobs priority', default=5)
args = parser.parse_args()

if os.path.isfile(args.filename):
    md_kwargs = {}


else:
    all_md_data = loadfn(args.filename)
    all_ids = {}
    for mol_name, md_kwargs in all_md_data.items():
        pass

    dumpfn(all_ids, "md_data/md_{}_molIDs.json".format(args.filename.split("/")[-1]), indent=2)
