import pandas as pd
"""
This script takes a compound managment file and puts it into the correct format for inserting the compounds
into the compound table in my compound database.
"""

# Global variables
path_cm_file = '/Users/bfulroth/Library/Mobile Documents/com~apple~CloudDocs/Broad Files 2/CM_Orders/190503_made_up_plate_map_from_cm_for_Amgen_covalent_modifer_BRD-K96889205-001-01-7.xlsx'

file_xlsx = True

# Date aliquoted in 'dd/mm/yy' format.
date_aliquoted = '05/03/19'


# Remaining volume
remain_vol_uL = 20

# Read in original CM file into a data frame
if file_xlsx:
    df_cm_ori = pd.read_excel(path_cm_file)

else:
    df_cm_ori = pd.read_csv(path_cm_file, sep='\t')

# Reformat the the compound management file for insertion into my compound database.
# Add new columns
df_cm_ori['Date_Aliquoted'] = date_aliquoted
df_cm_ori['Rem_Vol_uL'] = remain_vol_uL

df_for_db = df_cm_ori.loc[:, ['Broad ID', 'Tube Barcode', 'Position', 'Plate', 'Concentration (mM)', 'Volume (ul)'
    ,'Rem_Vol_uL', 'Molecular Weight', 'Date_Aliquoted']]

# Rename the columns
df_for_db.columns = ['Broad_ID', 'Barcode', 'Well', 'Plate', 'Conc_mM', 'Ori_Vol_uL', 'Rem_Vol_uL', 'MW',
                     'Date_Aliquoted']

df_for_db = df_for_db.dropna(subset=['Broad_ID', 'Barcode', 'Well', 'Conc_mM', 'Ori_Vol_uL', 'Rem_Vol_uL', 'MW',
                     'Date_Aliquoted'])

parsed_date = date_aliquoted.split('/')

df_for_db.to_csv('/Users/bfulroth/PycharmProjects/CompoundDB/' + parsed_date[2] + parsed_date[0] + parsed_date[1] +
                 '_db_cmpd_import_file_.csv', index=None)