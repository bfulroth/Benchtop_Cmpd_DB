from sqlalchemy import engine, create_engine, Table, MetaData, select, and_
import my_utilities
import pandas as pd

engine = create_engine('sqlite:////Users/bfulroth/PycharmProjects/CompoundDB/181018_bench_cmpds_db.sqlite')

conn = engine.connect()

metadata = MetaData()

cmpd_tbl = Table('cmpd_tbl', metadata, autoload=True, autoload_with=engine)

bar_df = pd.read_clipboard(header=None)
bar_df.columns = ['Bar']
bar_ls = bar_df['Bar'].tolist()

stmt = select([cmpd_tbl]).where(cmpd_tbl.c.Barcode.in_(bar_ls))

results = conn.execute(stmt).fetchall()

df_results = pd.DataFrame(results)

df_results.columns = cmpd_tbl.columns.keys()

print(df_results)

df_setup_tbl = my_utilities.spr_from_db_to_setup_tbl(df_results)

print(df_setup_tbl)

df_setup_tbl.to_clipboard(index=None)

conn.close()
