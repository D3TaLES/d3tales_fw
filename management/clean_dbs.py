import pymongo
import argparse
from d3tales_api.D3database.d3database import D3Database, FrontDB


def get_id_list(query, database=FrontDB()):
    cursor = database.coll.find(query)
    return set([data['_id'] for data in cursor])


def wipe_property(prop_query, id_list):
    for i in id_list:
        FrontDB().coll.update_one({"_id": i}, {"$unset": {prop_query: ""}})
        print("Property {} wiped from {}".format(prop_query, i))


def remove_duplicate_prop(_id, prop, precision=1):
    """
db.getCollection('base').aggregate([
  {
    $addFields: {
      "mol_characterization.omega": {
        $reduce: {
          input: "$mol_characterization.omega",
          initialValue: [],
          in: {
            $cond: [
              {
                $in: [
                  "$$this.value",
                  "$$value.value"
                ]
              },
              "$$value",
              {
                $concatArrays: [
                  "$$value",
                  [
                    "$$this"
                  ]
                ]
              }
            ]
          }
        }
      }
    }
  }
])
    """
    molecule = FrontDB().coll.find({"_id": _id})[0]
    query_terms = prop.split('.')
    queried_prop = molecule
    for term in query_terms:
        queried_prop = queried_prop.get(term, {})
    new_prop_items, values = [], []
    for prop_item in queried_prop:
        value = round(prop_item.get('value', 0), precision)
        if value in values:
            continue
        values.append(value)
        new_prop_items.append(prop_item)
    FrontDB().coll.update_one({"_id": _id}, {"$set": {prop: new_prop_items}}, upsert=True)
    return _id


def check_conditions_omega(_id, prop, correct_omega):
    molecule = FrontDB().coll.find({"_id": _id})[0]
    query_terms = prop.split('.')
    queried_prop = molecule
    for term in query_terms:
        queried_prop = queried_prop.get(term, {})
    if not queried_prop:
        return None
    new_prop_items = []
    for prop_item in queried_prop:
        item_omega = round(prop_item.get('conditions', {}).get('tuning_parameter', 0), 2)
        if item_omega == round(correct_omega, 2):
            new_prop_items.append(prop_item)
    if len(new_prop_items) == 0:
        print("Error! No {} items for {} have the correct omega".format(prop, _id))
        return _id
    FrontDB().coll.update_one({"_id": _id}, {"$set": {prop: new_prop_items}}, upsert=True)


def get_correct_omega(_id):
    molecule = FrontDB().coll.find({"_id": _id})[0]
    omega_list = molecule.get("mol_characterization", {}).get("omega")
    if not omega_list:
        print("Error. There are NO omega values for {}.".format(_id))
        return None
    if len(omega_list) == 1:
        return omega_list[0].get('value')
    print("Error. There are {} omega values for {}.".format(len(omega_list), _id))


def query_and_update(mol_prop_list, species_prop_list, species_list=None, functional_check=False, omega_check=False):
    """
    Function to identify and remove properties with (1) the wrong functional (not LC-wHPBE) if functional_check is True, (2) duplicate
    properties, and (3) properties with the incorrect tuning parameter (i.e., not the tuning parameter listed in
    mol_characteristics) if omega_check is true.
    """
    species_list = species_list or ["groundState", "cation1", "cation2"]
    two_items_props = ["homo_lumo_gap", "dipole_moment", "geometry", ]

    query_list = ["mol_characterization.{}".format(prop) for prop in mol_prop_list]
    for prop in species_prop_list:
        for species in species_list:
            query_list.append("species_characterization.{}.{}".format(species, prop))

    error_ids = []
    print("Getting no omega ids...")
    no_omega_ids = get_id_list({"mol_characterization.omega.1": {'$exists': True}})
    print("starting clean...")
    for query in query_list:
        query_field = "{}.1".format(query) if query.split('.')[-1] not in two_items_props else "{}.2".format(query)
        unapproved_cursor = FrontDB().coll.find({query_field: {'$exists': True}})
        ids_set = set([data['_id'] for data in unapproved_cursor if data['_id'] not in no_omega_ids])
        print("Cleaning {} ids with too many {}".format(len(ids_set), query))
        for _id in ids_set:
            # remove items with a wrong functional
            if functional_check:
                FrontDB().coll.update_one(
                    {"_id": _id},
                    {"$pull": {query: {"conditions.functional": {"$ne": "LC-wHPBE"}}}},
                )
            # remove duplicates
            error_id = remove_duplicate_prop(_id, query)
            if error_id:
                error_ids.append(error_id)

            # remove items with incorrect tuning parameter
            if omega_check:
                correct_omega = get_correct_omega(_id)
                if not correct_omega:
                    error_ids.append(_id)
                    continue
                error_id = check_conditions_omega(_id, query, correct_omega)
                if error_id:
                    error_ids.append(error_id)

    print("----------  Error IDs  ------------")
    print(len(set(error_ids)))
    print(set(error_ids))


if __name__ == "__main__":
    mol_props = ["hole_reorganization_energy", "electron_reorganization_energy", "relaxation_groundState_cation1",
                 "relaxation_cation1_groundState", "relaxation_groundState_anion1", "relaxation_anion1_groundState",
                 "vertical_ionization_energy", "vertical_electron_affinity", "adiabatic_ionization_energy",
                 "adiabatic_electron_affinity", "oxidation_potential", "reduction_potential",
                 "rmsd_groundState_cation1", "rmsd_cation1_cation2", "rmsd_groundState_anion1", "rmsd_anion1_anion2",
                 "omega"]
    species_props = ["charge", "globular_volume", "radical_stability_score", "homo_lumo_gap", "dipole_moment",
                     "solvation_energy", "radical_spin", "radical_buried_vol"]

    MOL_CHARACTERIZATION = ",".join(mol_props)
    SPECIES_CHARACTERIZATION = ",".join(species_props)
    SPECIES = "groundState,cation1,cation2,anion1,anion2"

    parser = argparse.ArgumentParser(description='Clean frontend database for D3TaLES')
    parser.add_argument('-mc', '--mol_props', type=str, help='mol characterization properties to be cleaned', default=MOL_CHARACTERIZATION)
    parser.add_argument('-sc', '--species_props', type=str, help='species characterization properties to be cleaned', default=SPECIES_CHARACTERIZATION)
    parser.add_argument('-s', '--species', type=str, help='list of species to be cleaned', default=SPECIES)
    parser.add_argument('-fc', '--functional_check', action='store_true', help='check if the wrong functional (not LC-wHPBE)')
    parser.add_argument('-oc', '--omega_check', action='store_true', help='check if properties have the incorrect tuning parameter (i.e., not the tuning parameter listed in mol_characteristics)')
    parser.add_argument('-wp', '--wipe_prop', type=str, help='property to wipe clean from database', default='')
    parser.add_argument('-wi', '--wipe_ids', type=str, help='ids from which to wipe wipe_prop', default='')
    args = parser.parse_args()

    if args.wipe_prop and args.wipe_ids:
        # ids = get_id_list({'mol_characterization.omega.1': {'$exists': True}})
        wipe_property(args.wipe_prop, args.wipe_ids.split(','))
    else:
        query_and_update(mol_prop_list=args.mol_props.split(','), species_prop_list=args.species_props.split(','),
                          species_list=args.species.split(','), functional_check=args.functional_check, omega_check=args.omega_check)
