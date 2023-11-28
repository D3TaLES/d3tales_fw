import json

with open("db_unfinished_ids_20231031-133808.json") as fn:
    unfinished = json.load(fn).keys()
with open("db_archived_ids_20230620-100803.json") as fn:
    archived = json.load(fn).keys()

run_ids = [i for i in unfinished if i not in archived]
print(",".join((run_ids)))
