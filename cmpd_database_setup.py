from sqlalchemy import Table, Column, String, Float, Date, MetaData, create_engine, insert
from datetime import datetime

# engine is a common interface to a database.
engine = create_engine('sqlite:////Users/bfulroth/PycharmProjects/CompoundDB/181018_bench_cmpds_db.sqlite')

# Create a connection to the database.
conn = engine.connect()

# Create metadata object
metadata = MetaData()

# Create a Table for Compounds
cmpd_tbl = Table('cmpd_tbl', metadata,
                    Column ('Broad_ID', String(22), nullable=False),
                    Column ('Barcode', String(10)),
                    Column ('Well', String(3)),
                    Column ('Plate', String(7)),
                    Column ('Conc_mM', Float()),
                    Column ('Ori_Vol_uL', Float()),
                    Column ('Rem_Vol_uL', Float()),
                    Column ('MW', Float()),
                    Column ('Date_Aliquoted', Date))

# Create seperate solubility table.
sol_tbl = Table('cmpd_sol_tbl', metadata,
                    Column('Broad_ID', String(22)),
                    Column('Buffer', String()),
                    Column('Sol_uM', Float()),
                    Column('Exp_Date', Date),
                    Column('Source', String())
                   )


# Create the metadata tables.
metadata.create_all(engine)

# Print Results
print('Compound Tbl:')
print('')
print(repr(cmpd_tbl))
print('Solubility Tbl:')
print('')
print(repr(sol_tbl))
print('')
print('printing metadata table keys...')
print('')
print(metadata.tables.keys())

# Create the SQlite Database from the csv files and for compound data and solubility data.
import csv

# Create insert statment for compound table
stmt = insert(cmpd_tbl)

csv_file_cpd_tbl = '/Users/bfulroth/PycharmProjects/CompoundDB/181018_cmpd_tbl_FINAL.csv'

with open(csv_file_cpd_tbl, newline='') as file:
    csv_reader = csv.reader(file)

    values_list = []
    row_count = 0

    # Read past header
    next(csv_reader)

    for idx, row in enumerate(csv_reader):
        records = {'Broad_ID': row[0], 'Barcode': row[1], 'Well': row[2], 'Plate': row[3], 'Conc_mM': row[4],
                   'Ori_Vol_uL': row[5], 'Rem_Vol_uL': row[6], 'MW': row[7],
                   'Date_Aliquoted': datetime.strptime(row[8],'%m/%d/%y')}
        values_list.append(records)
        row_count += 1

    results = conn.execute(stmt, values_list)

print('Row count for compounds table: ' + str(row_count))

# Create insert statement for solubility table
stmt_sol = insert(sol_tbl)

csv_file_sol_tbl = '/Users/bfulroth/PycharmProjects/CompoundDB/181018_solubility_data.csv'

with open(csv_file_sol_tbl, newline='') as sol_file:
    csv_reader_sol = csv.reader(sol_file)

    sol_values_list = []
    row_count = 0

    # Read past header
    next(csv_reader_sol)

    for idx, row in enumerate(csv_reader_sol):
        records = {'Broad_ID': row[0], 'Buffer': row[1], 'Sol_uM': row[2],
                   'Exp_Date': datetime.strptime(row[3],'%m/%d/%y'),'Source': row[4]}
        sol_values_list.append(records)

        row_count += 1

    results_sol = conn.execute(stmt_sol, sol_values_list)

print('Row count for Solubility table: ' + str(row_count))


