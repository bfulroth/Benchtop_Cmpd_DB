from sqlalchemy import engine, create_engine, Table, MetaData, insert
from datetime import datetime

engine = create_engine('sqlite:////Users/bfulroth/PycharmProjects/CompoundDB/181018_bench_cmpds_db.sqlite')

conn = engine.connect()

# Reflect the table
metadata = MetaData()
cmpd_tbl = Table('cmpd_tbl', metadata, autoload=True, autoload_with=engine)

# Prepare insert statement
stmt = insert(cmpd_tbl)

# Open a connection to the file you want to import
import csv

csv_file_cpd_tbl = '/Users/bfulroth/PycharmProjects/CompoundDB/190503_db_cmpd_import_file_.csv'

with open(file=csv_file_cpd_tbl, newline='') as file:
    csv_reader = csv.reader(file)

    values_list = []
    row_count = 0

    # read past first line containing the header
    next(csv_reader)

    for idx, row in enumerate(csv_reader):
        records = {'Broad_ID': row[0], 'Barcode': row[1], 'Well': row[2], 'Plate': row[3], 'Conc_mM': row[4],
                   'Ori_Vol_uL': row[5], 'Rem_Vol_uL': row[6], 'MW': row[7],
                   'Date_Aliquoted': datetime.strptime(row[8],'%m/%d/%y')}

        values_list.append(records)
        row_count += 1

    results = conn.execute(stmt, values_list)

    print(str(row_count) + ' Rows were inserted into table *cmpd_tbl*')
