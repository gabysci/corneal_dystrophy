import glob
import pandas as pd
import os
import numpy as np
from copy import deepcopy

# only the files that have _S*.xlsx in their name will be included i.e. the 5 sample files
file_list = glob.glob("./patient_samples/*_S*.xlsx")

# creating an empty dataframe and an empty list to place the filenames of the samples
list_dataframes = list()
list_samples = list()


for file_name in file_list:
    data1 = pd.read_excel(file_name)
    list_samples.append(os.path.basename(file_name).split('.')[0])

    # to exclude syn and intron variants that are found in consequence column
    condition1 = data1["Consequence"] != "synonymous_variant"
    condition2 = data1["Consequence"] != "intron_variant"
    filtered = data1[condition1 & condition2]
    columns_dataset = list(filtered.columns)
    list_dataframes.append(filtered)


final_dataframe = pd.DataFrame(columns=list_dataframes[0].columns)

sample_specific_columns_old = ['Type', 'Genotype', 'Filters', 'Quality', 'GQX', 'Alt Variant Freq', 'Read Depth',
                           'Alt Read Depth', 'Allelic Depths', 'Num Transcripts']
final_dataframe = final_dataframe.drop(sample_specific_columns_old, axis=1)
variant_specific_columns = final_dataframe.columns

for x in range(len(list_dataframes)):
    sample_specific_columns = ['Type', 'Genotype', 'Filters', 'Quality', 'GQX', 'Alt Variant Freq', 'Read Depth',
                               'Alt Read Depth', 'Allelic Depths', 'Num Transcripts']
    for y in range(len(sample_specific_columns)):
        name_col = '_' + list_samples[x]
        sample_specific_columns[y] = sample_specific_columns[y] + name_col
    sample_dataframe = pd.DataFrame(columns=sample_specific_columns)

    final_dataframe = pd.concat([final_dataframe, sample_dataframe])

sample_specific_column_list = final_dataframe.columns[45:]

# Commence the cross-comparison between different samples
# independent copy of filtered data is made
dataframe_1 = deepcopy(list_dataframes[0])

# a df is created containing a column named ID_variant where data of the 3 columns that
# make it possible to identify a variant are combined into one
dataframe_1['ID_variant'] = (dataframe_1['Chr'].astype(str) + '_'
                             + dataframe_1['Coordinate'].astype(str) + '_' + dataframe_1['Variant'].astype(str))
list_vars_1 = (list(dataframe_1['ID_variant']))

for i, var_id in enumerate(list_vars_1):
    print('Running variant number ' + str(i) + ' of ' + str(len(list_vars_1)) + ' total variants')
    temp_variant_list = []
    continue_variant = True
    # Add variant specific columns from sample 1, irrespective of if it is present in other samples
    for column_name in variant_specific_columns:
        temp_variant_list.append(list(dataframe_1[column_name])[i])
    # Add sample specific columns from sample 1, irrespective of if it is present in other samples
    for column_name in sample_specific_columns_old:
        temp_variant_list.append(list(dataframe_1[column_name])[i])
    # Commence check against other samples
    for dataframe_index2 in range(1, len(list_dataframes)):
        dataframe_2 = list_dataframes[dataframe_index2]
        dataframe_2['ID_variant'] = (dataframe_2['Chr'].astype(str) + '_'
                                     + dataframe_2['Coordinate'].astype(str) + '_' + dataframe_2['Variant'].astype(str))
        list_vars_2 = np.array(list(dataframe_2['ID_variant']))
        index_var = np.where(list_vars_2 == var_id)
        if len(index_var[0]) > 0:
            for column_name in sample_specific_columns_old:
                temp_variant_list.append(list(dataframe_2[column_name])[index_var[0][0]])
        else:
            break

        if dataframe_index2 == len(list_dataframes)-1:
            final_dataframe.loc[len(final_dataframe.index)] = temp_variant_list

final_dataframe.to_excel('filtered_variants_patients.xlsx')

