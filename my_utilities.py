
import pandas as pd
from glob import glob
import re
import numpy as np
from numpy import NaN
import matplotlib.pyplot as plt
import os
import functools
import sys
from datetime import datetime


def dup_item(df, col_name, times_dup=3, sort=False):
    """
    Takes a DataFrame and a column name with items to be replicated. Sorts the list and replicates the number of
    times specified by the parameter times_dup. Copies the replicated values to the clip board.

    :param df: A DataFrame containing the column of values to be replicated.
    :param col_name: Name of the column containing values to replicate.
    :param times_dup: Number of times to replicate each value in the specified column.
    :param sort: Boolean to sort the replicated values.
    :type sort: bool
    """

    dup_list = []

    try:
        for item in df[col_name]:
            for i in range(times_dup):
                dup_list.append(item)

        a = pd.Series(dup_list)

        if sort == True:
             b = a.sort_values()
             print(b)
             b.to_clipboard(index=False)
        else:
            print(a)
            a.to_clipboard(index=False)
        print("Done!")
    except:
        print("The DataFrame does not have a " + col_name + " column.")


def count_entries(df, col_name = 'BRD'):
    """
    Accepts a DataFrame as a parameter and a column name in the DataFrame and returns a dict with unique items as keys
    and their counts as values.
    :param df: DataFrame containing column of items we want to count.
    :arg col_name: Name of the column to count unique items.
    :return Dict: containing
    occurrences as value for each key."""

    # Initialize an empty dictionary: item_occur
    item_occur = {}
    try:
        # Extract column from DataFrame
        col = df[col_name]

        # Iterate over column in DataFrame
        for entry in col:

            # If the item is in item_occur, add 1
            if entry in item_occur.keys():
                item_occur[entry] += 1
            # Else add the item to item_occur, set the value to 1
            else:
                item_occur[entry] = 1

        # Return the item_occur dictionary
        return item_occur

    except:
        print("The DataFrame does not have a " + col_name + " column.")


def agg_plate_data(dir, file_extension, skip_rows=1, nrows=17, sep=','):
    """
    This function aggregates plate data from Envision Plate Reader

    :param dir: Directory of where the Envison files are stored.
    :param file_extension:
    :param skip_rows: Number of rows to skip before the start of the data.
    :param nrows: Number of rows to read that has data.
    :param sep: deliminator of the data.
    :returns A new data DataFrame composed of all of the Envision files concatenated together.

    """

    # Change the current working directory.
    os.chdir(dir)

    if file_extension == 'txt':
        pattern = '*.txt'
    else:
        pattern = '*.csv'

    frames = []  # initialize an empty list that will contain DataFrames.

    text_files = glob.glob(pattern)

    for txt in text_files:
        df = pd.read_csv(txt, skip_rows, nrows, sep)
        frames.append(df)

    all_frames = pd.concat(frames)
    return all_frames

# Note when calling this method on a DataFrame using the apply method df.apply(clean_two_cols, axis = 1,
# pattern=pattern)
# with axis = 1 the method will work row wise.
pattern = re.compile("^\$\d*\.\d{2}$")
def clean_two_cols(row, pattern, col1, col2):
    """Cleans up two columns of a DataFrame.  Creates a new subtracted column"""

    col1 = row[col1]
    col2 = row[col2]

    if bool(pattern.match(col1)) and bool(pattern.match(col2)):
        col1 = col1.replace('$', '')
        col2 = col2.replace('$', '')

        col1 = float(col1)
        col2 = float(col2)

        return col1 - col2
    else:
        return NaN

def clean_one_col(row, rep_val, col):
    """Cleans up a columns of a data frame"""
    col = row[col]

    col = col.replace(rep_val, '')

    col = float(col)

    return col

# Define recode row
def recode_col(col_value, replace1, replace2):

    # Return 1 if col_value is replace1
    if col_value == replace1:
        return 1

    # Return 0 if col_value is replace2
    elif col_value == replace2:
        return 0

    # Return np.nan
    else:
        return np.nan

# Apply the function to the sex column
# tips['sex_recode'] = tips.sex.apply(recode_sex, replace1='Male', replace2='Female')
# tips['total_dollar_replace'] = tips.total_dollar.apply(lambda x: x.replace('$', ''))

# multiply two columns togeather
#df['Multiply MW and G12D'] = df.apply(lambda x: x['MW'] * x['G12D'], axis=1)



def list_files(dir, fileEx='pdf'):
    """This function takes a directory and a file extension and
    :returns a list of files in that directory with that extension.  Also saves the list to the clipboard."""

    #Change the current working directory.
    os.chdir(dir)

    pattern = '*.' + fileEx

    textFiles = glob(pattern)

    file_list = pd.Series(textFiles)

    file_list.to_clipboard(index=False)
    return file_list


def create_spr_summary(df, ass_funct=False, create_protein_col=True, groupby_brd=False, group_by_protein = False):
    """This method takes an SPR results DataFrame and trims it down to a summary table with a focus on
     KD value and percent binding information. """

    #Trim the DataFrame to the revelent columns
    df_trim = df.loc[:, ['broad_id', 'Percent_Binding_at_top_cmpd', 'KD_steady_state_affinity_uM', 'Flow_cell',
                         'Date_']]

    # Rename columns
    df_trim.columns = ['BRD', '%_Bind_Top', 'KD_uM', 'FC', 'Date']

    # coherce KD column to numeric value.
    #df_trim['KD_uM'] = pd.to_numeric(df_trim['KD_uM'], errors='coerce')

    if create_protein_col:

        if ass_funct == False:
            # Affinity Assay
            # Create a dictionary of proteins on each flow cell
            protein = {'FC2-1corr' : 'G12D', 'FC3-1corr' : 'G12D/T35S', 'FC4-1corr' : 'WT'}
        else:
            # Functional assay 6.2
            # Create a dictionary of proteins on each flow cell
            protein = {'FC2-1corr': 'BRAF-RBD', 'FC3-1corr': 'CRAF-RBD', 'FC4-1corr': 'ARAF-RBD'}

        #  Add a protein column using the protein dictionary above.
        df_trim['Protein'] = df_trim['FC'].map(protein)

        # reorder the colums
        df_trim_order = df_trim.loc[:, ['BRD', 'KD_uM', '%_Bind_Top', 'FC', 'Protein', 'Date']]

        if group_by_protein == True:
            try:
                df_trim_order = df_trim_order.groupby(['Protein', 'BRD'])
                [['KD_uM', '%_Bind_Top']].agg([np.mean, np.std, np.count_nonzero])
            except:
                print("Error in performing grouby operaton.")
    else:
        df_trim_order = df_trim.loc[:, ['BRD', 'KD_uM', '%_Bind_Top', 'FC', 'Date']]

    if groupby_brd == True:

        try:
            df_trim_order = df_trim_order.groupby(['BRD'])[['KD_uM', '%_Bind_Top']].agg([np.mean, np.std,
                                                                                         np.count_nonzero])
            # Add in a % CV calculation for KD.
            df_trim_order['KD', '%_COV'] = df_trim_order.apply(lambda row: (row['KD_uM', 'std']/ row['KD_uM', 'mean'])*100,
                                                             axis=1)
        except:
            print("Error in performing groupby operation.")

    # return a new DataFrame that can be further processed.
    return df_trim_order


def create_spr_summary_v2(df, ass_funct=False, create_protein_col=True, groupby_brd=False, group_by_protein = False):
    """
    This method takes an SPR ADLP upload file as a DataFrame and returns a summary of the KD and % binding for each
    compound tested.

    Note: Version 2 of this method to be used for SPR ADLP Excel upload forms created after 11/13/2018.

    :param df: ADLP upload form in a DataFrame
    :param ass_funct:
    :param create_protein_col: Boolean to create a new column of KRAS protein names.
    :param groupby_brd: Boolean to group the new DataFrame by BRD.
    :param group_by_protein: Boolean to group the new DataFrame by protein.
    :return: New DataFrame that is a summary of the original ADLP SPR upload file.
    """

    #Trim the DataFrame to the revelent columns
    df_trim = df.loc[:, ['BROAD_ID', '%_BINDING_TOP', 'KD_SS_UM', 'FC', 'EXP_DATE']]

    # Rename columns
    df_trim.columns = ['BRD', '%_Bind_Top', 'KD_uM', 'FC', 'Date']

    # coherce KD column to numeric value.
    df_trim['KD_uM'] = pd.to_numeric(df_trim['KD_uM'], errors='coerce')

    if create_protein_col:

        if ass_funct == False:
            # Affinity Assay
            # Create a dictionary of proteins on each flow cell
            protein = {'FC2-1Corr': 'G12D', 'FC3-1Corr': 'G12D/T35S', 'FC4-1Corr': 'WT'}
        else:
            # Functional assay 6.2
            # Create a dictionary of proteins on each flow cell
            protein = {'FC2-1Corr': 'BRAF-RBD', 'FC3-1Corr': 'CRAF-RBD', 'FC4-1Corr': 'ARAF-RBD'}

        #  Add a protein column using the protein dictionary above.
        df_trim['Protein'] = df_trim['FC'].map(protein)

        # reorder the colums
        df_trim_order = df_trim.loc[:, ['BRD', 'KD_uM', '%_Bind_Top', 'FC', 'Protein', 'Date']]

        if group_by_protein == True:
            try:
                df_trim_order = df_trim_order.groupby(['Protein', 'BRD'])
                [['KD_uM', '%_Bind_Top']].agg([np.mean, np.std, np.count_nonzero])
            except:
                print("Error in performing grouby operaton.")
    else:
        df_trim_order = df_trim.loc[:, ['BRD', 'KD_uM', '%_Bind_Top', 'FC', 'Date']]

    if groupby_brd == True:

        try:
            df_trim_order = df_trim_order.groupby(['BRD'])[['KD_uM', '%_Bind_Top']].agg([np.mean, np.std,
                                                                                         np.count_nonzero])
            # Add in a % CV calculation for KD.
            df_trim_order['KD', '%_COV'] = df_trim_order.apply(lambda row: (row['KD_uM', 'std']/ row['KD_uM', 'mean'])*100,
                                                             axis=1)
        except:
            print("Error in performing groupby operation.")

    # return a new DataFrame that can be further processed.
    return df_trim_order


def plot_spr_kd(df, index_col_name='BRD', results='KD_uM', tight_layout=True):
    """This method takes a DataFrame, averages the result values passed in as a column name and plots the mean of the
    values and error bars for each BRD.
    :return plot object."""

    #Create a new data frame with the KD values averaged.
    df.set_index(index_col_name)
    df_pivot = pd.pivot_table(df, values=[results], index=[index_col_name], aggfunc=[np.mean, np.std])
    df_pivot.columns = df_pivot.columns.droplevel(1)

    #Create the graph
    graph = df_pivot.plot(kind='bar', y='mean', yerr='std', capsize=2)
    graph.set_ylabel('KD_uM', fontsize=12)

    if tight_layout:
        plt.tight_layout()

    return graph


def match(row, pattern, colName='BRD'):
    """When placed inside the pandas apply method, this method accepts a row in a DataFrame and a regex expression
    as a compiled pattern and matches the pattern to each element of the specified column."""
    BRD = row[colName]

    if bool(pattern.match(BRD)):
        return "T"
    else:
        return "F"


def rename_cols_spr_agg_curves(df):
        """Method that renames the columns of the agg_spr_func_curves() DataFrames and puts
        the concentration into the correct units."""

        # Determine how many columns in the DataFrame
        x, numCols = df.shape

        # empty list to store replicate headers
        rep_list = []

        # Variable to store the number of replicates
        rep_num = 1

        for col in range(0, numCols - 1):
            rep_list.append('N' + str(rep_num))
            rep_num += 1

        rep_list.insert(0, "Conc. (uM)")

        # Change all the header names
        df.columns = rep_list

        # Change concentration in M to uM.
        df['Conc. (uM)'] = df.apply(lambda x: df['Conc. (uM)'] * 1000000)

        df = df.round(2)

        return df


def merge_dfs(file_dir, df_list, merge_all_dfs_in_df_list=False, mrg_col_name='X'):
    """
    This function takes a list of text files and merges the files based on a single column.

    :param file_dir: path of the directory of the main folder where the files to aggregate are located.
    :param df_list: a list of DataFrames to aggregate.
    :param merge_all_dfs_in_df_list: Booolean to merge all files.
    :param mrg_col_name: Name of the column to merge the data on.
    :returns the merged files as a new DataFrame.
    """

    if merge_all_dfs_in_df_list == False:
        # Change the current working directory.
        os.chdir(file_dir)

        pattern = '*.txt'
        frames = []

        text_files = glob(pattern)
        text_files.sort()

        for text in text_files:
            df = pd.read_csv(text, sep='\t')
            df = df.loc[:, ['X', 'Y']]

            df = df.dropna(how='all', axis='rows')
            frames.append(df)

        df_all = functools.reduce(lambda x, y: pd.merge(x, y, on=mrg_col_name), frames)

        return df_all

    else:
        df_all = functools.reduce(lambda x, y: pd.merge(x, y, on=mrg_col_name), df_list)

        return df_all




def get_list_folder_names(path):
    """This function returns a sorted list of folder names in the file path passed as an argument to this method."""

    dirs = os.listdir(path)

    # Need to remove this hidden file as it causes issues when looping through the folder names.
    dirs.remove('.DS_Store')
    dirs.sort()

    return dirs


def extract_rmax_theory_spr(txtFilePath, sort=False):
    """This function uses a table of calculated RMAX values in SPR Protocol Version 3 (3 Proteins)
    and creates a list of RMAX value's in the order the compounds were run.  This list is copied to the clipboard and can
    be pasted directly into the SPR data analysis sheet."""

    # Read data into a DataFrame
    df_rmax_ori = pd.read_csv(txtFilePath, sep='\t')

    # Delete the molecular weight columns
    del df_rmax_ori['MW']

    if sort:
        df_rmax_ori = df_rmax_ori.sort_values(by=['Broad ID'])

    # create a new transposed DataFrame
    df_rmax_ori_trans = df_rmax_ori.transpose()

    # Slice out the columns you want
    df_rmax_ori_trans = df_rmax_ori_trans.iloc[1:4, :]

    # Loop through each column and add each value in order to a new list.

    list = []

    for col in df_rmax_ori_trans:
        x = df_rmax_ori_trans[col].tolist()
        for num in x:
            list.append(num)

    # convert the list to a DataFrame
        rmax = pd.DataFrame(list)
        rmax.to_clipboard(index=False)

    return rmax


def spr_setup_sheet(df, path='', from_clip=True, to_excel=True):
    """Creates the setup file necessary to run a protocol on a Biacore instrument.  Takes a reference to a text file
    that is a copied from the setup compound plate table used to prepare the compound plate in my notebook. E180601-1
    is and example of the table used.

    :param df: DataFrame containing the data used as a template in my notebook to setup KRAS Biacore binding exps.
    :param path: Directory path of where the DataFrame is located if not using the clipboard.
    :type from_clip: bool
    :param from_clip: Wheather the DataFrame exists in the clipboard.
    :type to_excel: bool
    :param to_excel: Write the file to an Excel Workbook.
    """
    try:

        if from_clip:
            df_setup_ori = df

        else:
            # Read in the DataFrame.
            df_setup_ori = pd.read_csv(path, sep='\t')

        # Trim the sheet down to only the columns we need for the SPR setup sheet.
        df_setup_trim = df_setup_ori.loc[:, ['Broad ID','MW', 'Barcode',
           'Test [Cpd] uM', 'fold_dil', 'num_pts']]

        # Start building the setup sheet.
        # Store the dimensions of the DataFrame in variables that are used later in the method.
        nRows, nCols = df_setup_trim.shape

        # Create empty list used to build up the final DataFrame.
        brd_list = []
        mw_list = []
        bar_list = []
        conc_list = []

        # Inner method uses the original DataFrame to construct each column of the setup sheet.
        def create_lists(header, list):

            if header == 'Broad ID':
                unique_brd = 1
                for cmpd in range(nRows):
                    value = df_setup_trim.iloc[cmpd][header]

                    for i in range(int(df_setup_trim.iloc[cmpd]['num_pts']) + 2):

                        # As the SPR field limit is only 15 characters trim the BRD's
                        if len(value) == 22:
                            v = value[:3] + '-' + value[9:13] + '_' + str(unique_brd)
                            list.append(v)
                        else:
                            v = value + '_' + str(unique_brd)
                            list.append(v)
                    unique_brd += 1
            else:
                for cmpd in range(nRows):
                    value = df_setup_trim.iloc[cmpd][header]
                    for i in range(int(df_setup_trim.iloc[cmpd]['num_pts']) + 2):
                        list.append(value)

        # Inner method needed to create the dose response column.
        def dose_conc_list():

            for cmpd in range(nRows):
                # empty list to store each concentration in the dose response
                dose_list = []

                # Include two blank injections for each compound
                dose_list.append(0)
                dose_list.append(0)

                top = df_setup_trim.iloc[cmpd]['Test [Cpd] uM']  #store top dose in a variable.

                for i in range(int(df_setup_trim.iloc[cmpd]['num_pts'])):
                    dose_list.append(top)
                    top = top / int(df_setup_ori.iloc[cmpd]['fold_dil']) # use value in original DataFrame to determine
                    # the fold of dilution.
                dose_list.sort(reverse=False)

                # Want one final concentration list.  So add each concentration in the dose_list to the final conc_list.
                for c in dose_list:
                    conc_list.append(c)

        # Create the columns in the setup sheet
        create_lists(header='Broad ID', list=brd_list)
        create_lists(header='MW', list=mw_list)
        dose_conc_list()
        create_lists(header='Barcode', list=bar_list)

        # Create the final DataFrame from all of the lists.
        final_df = pd.DataFrame({'BRD': brd_list, 'MW': mw_list, 'CONC': conc_list, 'BAR': bar_list})

        # Reorder the columns
        final_df = final_df.loc[:, ['BRD', 'MW', 'CONC', 'BAR']]

        now = datetime.now()

        # Truncate the year in the file name.
        now = datetime.now()
        now = now.strftime('%y%m%d')

        if to_excel:
            try:
                final_df.to_excel('/Volumes/tdts_users/BFULROTH/' + now + '_spr_setup_affinity.xlsx')
            except:
                print('Issue connecting to Flynn. Mount drive and try again.')
                print('')

                final_df.to_excel('/Users/bfulroth/Desktop/' + now + '_spr_setup_affinity.xlsx')
                print('File created on desktop.')

        return final_df

    except:
        print("Something is wrong. Check the original file.")
        raise


def spr_setup_sheet_funct_aba(df, path='', from_clip=True, to_excel=True):
    """Creates the setup file necessary to run and an aba functional assay protocol on a Biacore instrument.
    Takes a reference to a text file that is a copied from the setup compound plate table used to prepare the compound
    plate in my notebook. E180830-1 is and example of the table used.

     :param df: DataFrame containing the data used as a template in my notebook to setup KRAS Biacore binding
    experiments.
    :param path: Directory path of where the DataFrame is located if not using the clipboard.
    :type from_clip: bool
    :param from_clip: Wheather the DataFrame exists in the clipboard.
    :type to_excel: bool
    :param to_excel: Write the file to an Excel Workbook.
    """

    try:

        if from_clip:
            df_setup_ori = df

        else:
            # Read in the DataFrame.
            df_setup_ori = pd.read_csv(path, sep='\t')

        # Trim the sheet down to only the columns we need for the SPR setup sheet.
        df_setup_trim = df_setup_ori.loc[:, ['Broad ID','MW', 'Barcode',
           'Test [Cpd] uM', 'fold_dil', 'num_pts']]



        # Start building the setup sheet.
        # Store the dimensions of the DataFrame in variables that are used later in the method.
        nRows, nCols = df_setup_trim.shape

        # Create empty list used to build up the final DataFrame.
        sample_sol_list = []
        flank_sol_list = []
        mw_list = []
        bar_list = []
        conc_list = []

        # Inner method uses the original DataFrame to construct the sample solution column.
        def create_sample_sol_list(header, list):

            for cmpd in range(nRows):
                value = df_setup_trim.iloc[cmpd][header]

                for i in range(int(df_setup_trim.iloc[cmpd]['num_pts']) + 2):

                    # As the SPR field limit is only 15 characters trim the BRD's
                    if len(value) == 22:
                        v = value[:3] + '-' + value[9:13] + '_' + 's'
                        list.append(v)
                    else:
                        v = value + '_' + 's'
                        list.append(v)

        # Inner method uses the original DataFrame to construct the flanking solution column.
        def create_flank_sol_list(header, list):

            unique_brd = 1
            for cmpd in range(nRows):
                value = df_setup_trim.iloc[cmpd][header]

                for i in range(int(df_setup_trim.iloc[cmpd]['num_pts']) + 2):

                    # As the SPR field limit is only 15 characters trim the BRD's
                    if len(value) == 22:
                        v = value[:3] + '-' + value[9:13] + '_' + str(unique_brd)
                        list.append(v)
                    else:
                        v = value + '_' + str(unique_brd)
                        list.append(v)
                    unique_brd += 1
                unique_brd = 1

        # Inner method needed to create the dose response column.
        def create_dose_conc_list():

            for cmpd in range(nRows):
                # empty list to store each concentration in the dose response
                dose_list = []

                # Include two blank injections for each compound
                dose_list.append(0)
                dose_list.append(0)

                top = df_setup_trim.iloc[cmpd]['Test [Cpd] uM'] #store top dose in a variable.

                for i in range(int(df_setup_trim.iloc[cmpd]['num_pts'])):
                    dose_list.append(top)
                    top = top / int(df_setup_ori.iloc[cmpd]['fold_dil']) # use value in original DataFrame to determine
                    # the fold of dilution.
                dose_list.sort(reverse=False)

                # Want one final concentration list.  So add each concentration in the dose_list to the final conc_list.
                for c in dose_list:
                    conc_list.append(c)

            # Inner method used to create the mw or bc columns.
        def create_mw_or_bc_list(header, list):

            for cmpd in range(nRows):
                value = df_setup_trim.iloc[cmpd][header]
                for i in range(int(df_setup_trim.iloc[cmpd]['num_pts']) + 2):
                    list.append(value)

        # Create the columns in the setup sheet by calling the inner methods above.
        create_sample_sol_list(header='Broad ID', list=sample_sol_list)
        create_flank_sol_list(header='Broad ID', list=flank_sol_list)
        create_dose_conc_list()
        create_mw_or_bc_list(header='MW', list=mw_list)
        create_mw_or_bc_list(header='Barcode', list=bar_list)

        # Create the final DataFrame from all of the lists.
        final_df = pd.DataFrame({'Sample Solution': sample_sol_list, 'Flanking Solution (A)': flank_sol_list,
                                 'CONC': conc_list, 'MW': mw_list, 'BAR': bar_list})

        # Reorder the columns
        final_df = final_df.loc[:, ['Sample Solution', 'Flanking Solution (A)', 'CONC', 'MW', 'BAR']]

        now = datetime.now()

        # Truncate the year in the file name.
        now = now.strftime('%y%m%d')

        if to_excel:
            try:
                final_df.to_excel('/Volumes/tdts_users/BFULROTH/' + now + '_spr_setup_ABA.xlsx')
            except:
                print('Issue connecting to Flynn. Mount drive and try again.')
                print('')

                final_df.to_excel('/Users/bfulroth/Desktop/' + now + '_spr_setup_ABA.xlsx')
                print('File created on desktop.')

        return final_df

    except:
        print("Something is wrong. Check the original file.")


def norm_spr_funct_curve_zero(df):

    # Extract the first row and convert to list
    df_fstRow = df.iloc[0]
    lst_fstRow = df_fstRow.tolist()

    # Extract current headers and save the Conc. Header
    col = df.columns
    col = col.tolist()
    col_zero = col[0]

    # Rename the columns
    nRows, nCols = df.shape
    col_num = 1
    col_list = []
    col_list.append(col_zero)

    for col in range(1, nCols):
        col_list.append('N' + str(col_num))
        col_num += 1

    # Replace old headers with new headers
    df.columns = col_list

    # Skip the first columns by storing the column headers as in an iterable object and calling next() on this iterable.
    iter_c = iter(df.columns)
    next(iter_c)

    # Main loop that normalizes each value in a column such that the first value is 0.
    i = 1
    for c in iter_c:
        df[c] = df.apply(lambda x: df[c] - lst_fstRow[i])
        i += 1

    return df


def agg_spr_func_curves(raw_folder_path, normalize=True, return_df_list=False):
    """This method aggregates all SPR functional raw curve data in a the specified raw data folder, the path of which
    is passed as a parameter to this method."""

    folder_list = os.listdir(raw_folder_path)

    if '.DS_Store' in folder_list:
        folder_list.remove('.DS_Store')

    # Need to create a new list of file names as integers so they can be sorted correctly.
    # Strings may not be sorted correctly
    folder_list = [int(i) for i in folder_list]
    folder_list.sort()

        # empty list of all the data frames. Each data frame are the aggregated results of 1 compound.
    df_list = []

    print('Order of the folder names:\n')

    for folder in folder_list:
        # need to convert integer names in folder list back to strings as strings are needed as a folder path.
        print(folder)
        active_folder = raw_folder_path + '/' + str(folder)
        active_file_data = merge_dfs(active_folder, df_list, merge_all_dfs_in_df_list=False)

        if normalize & return_df_list:
            active_file_data = df_col_renamed = rename_cols_spr_agg_curves(active_file_data)

        elif return_df_list:
            df_col_renamed = rename_cols_spr_agg_curves(active_file_data)

        else:
            df_col_renamed = active_file_data

        df_list.append(df_col_renamed)

    if return_df_list:
        return df_list

    # Merge all data if not returning a list of DataFrames.
    else:
        all_data = merge_dfs(raw_folder_path, df_list, merge_all_dfs_in_df_list=True)
        df_col_renamed_all_data = rename_cols_spr_agg_curves(all_data)

        if normalize:
            all_data_norm = norm_spr_funct_curve_zero(df_col_renamed_all_data)
            return all_data_norm

        else:
            return df_col_renamed_all_data


def spr_funct_norm_to_blank(report_pt_file, df_norm_curves):
    """This method calculates the percent displacement of our KRAS functional assay.

    :param report_pt_file: Reference to a file that contains the SPR runs report point information.
    :param df_norm_curves: DataFrame containing all the normalized curves df_norm_curves.
    :returns New DataFrame with all values normalized to the blank injection.

    Additional Info:
    Method extracts the blank injections (protein only) for each compound and then calculates the average signal for
    the blank. The method then calculates the percent displacement for each of the values in the df_norm_curves
    DataFrame and returns a new DataFrame.
    """

    try:
        # Read in data
        df_rpt_pts_all = pd.read_excel(report_pt_file, sheetname='Report Point Table', skiprows=3)

    except:
        raise FileNotFoundError('The files could not be imported please check.')

    # Check that the columns in the report point file match the expected values.
    expected_cols = ['Unnamed: 0', 'Cycle','Fc','Report Point','Time [s]','Window [s]','AbsResp [RU]','SD',
                     'Slope [RU/s]','LRSD','RelResp [RU]',	'Baseline',	'AssayStep','Assay Step Purpose',
                    'Buffer','Cycle Type','Temp','Sample_1_Conc [µM]','Sample_1_Ligand','Sample_1_MW [Da]',
                     'Sample_1_Sample']

    if df_rpt_pts_all.columns.tolist() != expected_cols:
        raise ValueError('The columns in the report point file do not match the expected names.')

    # Exclude first column as it is blank
    df_rpt_pts_trim = df_rpt_pts_all.iloc[:, 1:]

    # Remove other not needed columns
    df_rpt_pts_trim = df_rpt_pts_trim.loc[:,
                      ['Cycle', 'Fc', 'Report Point', 'Time [s]', 'RelResp [RU]', 'AssayStep', 'Cycle Type',
                       'Sample_1_Conc [µM]',
                       'Sample_1_Sample']]

    # Remove not needed rows.
    df_rpt_pts_trim = df_rpt_pts_trim[(df_rpt_pts_trim['Report Point'] == 'binding') &
                                      (df_rpt_pts_trim['AssayStep'] != 'Startup') &
                                      (df_rpt_pts_trim['AssayStep'] != 'Solvent correction') &
                                      (df_rpt_pts_trim['Sample_1_Conc [µM]'] == 0)]

    # Filter out non-corrected data.
    df_rpt_pts_trim['FC_Type'] = df_rpt_pts_trim['Fc'].str.split(' ', expand=True)[1]
    df_rpt_pts_trim = df_rpt_pts_trim[df_rpt_pts_trim['FC_Type'] == 'corr']

    # Store the mean of the RelResp [RU] column.
    mean = round(df_rpt_pts_trim['RelResp [RU]'].mean(), 3)
    print('Binding Stats for the blank injection of effector to KRAS')
    print('The Mean: ' + str(mean))

    stdev = np.std(df_rpt_pts_trim['RelResp [RU]'])
    stdev = round(stdev, 3)
    print('The Stdev: ' + str(stdev))

    # Calculate the percent displacement across the DataFrame using the mean of the blank.
    df_norm_percent_dis = df_norm_curves.iloc[:, 1:].apply(lambda x: -1*round((x/mean) * 100, 2))

    return df_norm_percent_dis


def spr_binding_top(report_pt_file, instrument='BiacoreS200'):
    """This method calculates the binding in RU at the top concentration.

        :param report_pt_file: reference to the report point file exported from the Biacore Instrument.
        :param instrument: The instrument as a string. (e.g. 'BiacoreS200', 'BiacoreT200')
        :returns DataFrame with the percent binding INCLUDING concentrations 50uM and 100uM.
        For compounds run at 100uM the 50uM concentration will need to be removed manually.
        """
    if (instrument != 'BiacoreS200') & (instrument != 'BiacoreT200'):
        raise ValueError('Instrument argument must be BiacoreS200 or BiacoreT200')

    try:
        # Read in data
        df_rpt_pts_all = pd.read_excel(report_pt_file, sheetname='Report Point Table', skiprows=3)
    except:
        raise FileNotFoundError('The files could not be imported please check.')

    # Biacore instrument software for the S200 and T200 instruments exports different column names.
    # Check that the columns in the report point file match the expected values.
    if instrument=='BiacoreT200':
        expected_cols = ['Cycle', 'Fc', 'Time', 'Window', 'AbsResp', 'SD', 'Slope', 'LRSD', 'Baseline', 'RelResp',
                         'Report Point', 'AssayStep', 'AssayStepPurpose', 'Buffer', 'CycleType', 'Temp',
                         'Sample_1_Sample', 'Sample_1_Ligand', 'Sample_1_Conc', 'Sample_1_MW', 'General_1_Solution']

    # Check that the columns in the report point file match the expected values.
    if instrument == 'BiacoreS200':
        expected_cols = ['Unnamed: 0', 'Cycle','Fc','Report Point','Time [s]','Window [s]','AbsResp [RU]','SD',
                     'Slope [RU/s]','LRSD','RelResp [RU]',	'Baseline',	'AssayStep','Assay Step Purpose',
                    'Buffer','Cycle Type','Temp','Sample_1_Barcode','Sample_1_Conc [µM]','Sample_1_Ligand',
                         'Sample_1_MW [Da]', 'Sample_1_Sample', 'General_1_Solution']

    if df_rpt_pts_all.columns.tolist() != expected_cols:
        raise ValueError('The columns in the report point file do not match the expected names.')

    # For BiacoreS200
    # Remove first column
    if instrument == 'BiacoreS200':
        df_rpt_pts_trim = df_rpt_pts_all.iloc[:, 1:]

    # Remove other not needed columns
        df_rpt_pts_trim = df_rpt_pts_trim.loc[:,
                      ['Cycle', 'Fc', 'Report Point', 'Time [s]', 'RelResp [RU]', 'AssayStep', 'Cycle Type',
                       'Sample_1_Conc [µM]',
                       'Sample_1_Sample']]

        # Remove not needed rows.
        df_rpt_pts_trim = df_rpt_pts_trim[df_rpt_pts_trim['Report Point'] == 'binding']
        df_rpt_pts_trim = df_rpt_pts_trim[(df_rpt_pts_trim['AssayStep'] != 'Startup') &
                                          (df_rpt_pts_trim['AssayStep'] != 'Solvent correction') &
                                          (df_rpt_pts_trim['Sample_1_Conc [µM]'] == 25) |
                                          (df_rpt_pts_trim['Sample_1_Conc [µM]'] == 50) |
                                          (df_rpt_pts_trim['Sample_1_Conc [µM]'] == 100)]

    # For Biacore T200 and T100
    else:
        # Remove other not needed columns
        df_rpt_pts_trim = df_rpt_pts_all.loc[:,
                          ['Cycle', 'Fc', 'Report Point', 'Time', 'RelResp', 'AssayStep', 'CycleType',
                           'Sample_1_Conc',
                           'Sample_1_Sample']]

        # Remove not needed rows.
        df_rpt_pts_trim = df_rpt_pts_trim[df_rpt_pts_trim['Report Point'] == 'binding']
        df_rpt_pts_trim = df_rpt_pts_trim[(df_rpt_pts_trim['AssayStep'] != 'Startup') &
                                          (df_rpt_pts_trim['AssayStep'] != 'Solvent correction') &
                                          (df_rpt_pts_trim['Sample_1_Conc'] == 25) |
                                          (df_rpt_pts_trim['Sample_1_Conc'] == 50) |
                                          (df_rpt_pts_trim['Sample_1_Conc'] == 100)]

    # Filter out non-corrected data.
    df_rpt_pts_trim['FC_Type'] = df_rpt_pts_trim['Fc'].str.split(' ', expand=True)[1]
    df_rpt_pts_trim = df_rpt_pts_trim[df_rpt_pts_trim['FC_Type'] == 'corr']

    df_final = df_rpt_pts_trim

    return df_final


def spr_cm_sheet_to_setup_tbl(df, df_brds):
    """
    This method takes a compound management compound distribution sheet and puts it in the correct format for my SPR
    setup table.

    :param df: The original tracking sheet as a DataFrame.
    :param df_brds: a list of BRD's that will be run.
    :return: New DataFrame in the same format as the SPR setup template.
    """

    # Convert DataFrame of Broad ID's to a list.
    df_brds.columns = ['Broad ID']
    list_brds = df_brds['Broad ID'].tolist()

    # Filter out compounds in the rack map that will not be run.
    df_cpds_to_run = df[df['Broad ID'].isin(list_brds)]

    # Add columns not present in the rack map but that are needed in the template.
    df_cpds_to_run.loc[:,'Comment'] = 'Test'
    df_cpds_to_run.loc[:,'BC Added'] = ''
    df_cpds_to_run.loc[:,'Sol. (uM)'] = ''

    # Rearrange columns so that they match the SPR setup table.
    df_cpds_to_run = df_cpds_to_run.loc[:, ['Broad ID', 'Comment', 'Molecular Weight', 'Sol. (uM)',
                                            'Tube Rack Barcode', 'Tube Barcode', 'BC Added', 'Position',
                                            'Concentration (mM)']]

    # Rename the columns to match the setup sheet.
    df_cpds_to_run.columns = ['Broad ID', 'Comment', 'MW', 'Sol. (uM)',
                                            'Plate', 'Barcode', 'BC Added',
                                            'Well', 'Conc (uM)']

    return df_cpds_to_run


def spr_from_db_to_setup_tbl(df_db):
    """
    This method takes a Compound Database result as a DataFrame argument and return a new DataFrame in the format
    needed to setup an SPR run (Setup Table in Notebook).
    :param df_db: DataFrame of a database compound query.
    :return: New Dataframe in the format needed for setting up an SPR binding experiment.
    """

    # Trim df_db for setup_tbl
    df_setup_tbl = df_db.loc[:, ['Broad_ID', 'Comment', 'MW', 'Sol. (uM)', 'Plate', 'Barcode', 'BC Added', 'Well',
                                'Conc_mM']]

    df_setup_tbl.columns = ['Broad ID', 'Comment', 'MW', 'Sol. (uM)', 'Plate', 'Barcode', 'BC Added', 'Well',
                                'Conc. (mM)']

    return df_setup_tbl









