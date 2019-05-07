import pandas as pd
from sqlalchemy import Table, engine, create_engine, MetaData, select

engine = create_engine('sqlite:////Users/bfulroth/PycharmProjects/CompoundDB/181018_bench_cmpds_db.sqlite')

conn = engine.connect()

metadata = MetaData()

cmpds = Table('cmpd_tbl', metadata, autoload=True, autoload_with=engine)

# Compounds to look query
df_brds = pd.read_clipboard(header=None)

df_brds.columns = ['BRD']
df_brd_ls = df_brds['BRD'].tolist()
print("List of BRD's: ")
print(df_brd_ls)


stmt = select([cmpds]).where(cmpds.columns.Broad_ID.in_(df_brd_ls))

results = conn.execute(stmt).fetchall()

df_results = pd.DataFrame(results)

df_results.columns = cmpds.columns.keys()

print(df_results)

df_results.to_clipboard()

conn.close()
