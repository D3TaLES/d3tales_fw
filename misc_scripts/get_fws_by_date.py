from datetime import datetime, timedelta
from d3tales_api.D3database.d3database import D3Database

STATE = "RUNNING"
DELTA_TIME = timedelta(hours=0, minutes=0)
CUTOFF_DATE = datetime(2022, 7, 12)

if __name__ == "__main__":
    if DELTA_TIME:
        CUTOFF_DATE = datetime.now() - DELTA_TIME
    defuse_ids = []
    fireworks_db = D3Database(database='fireworks', collection_name='fireworks')
    launch_db = D3Database(database='fireworks', collection_name='launches')

    fw_cursor = fireworks_db.coll.find({'$and': [
                {"state": STATE},
                # {"updated_on": {'$lt': CUTOFF_DATE.isoformat()}},
            ]})

    fw_ids = [l.get("fw_id") for l in fw_cursor]
    print('lpad rerun_fws -i ', ''.join(str(fw_ids).split(' ')).replace('[', '').replace(']', ''))

