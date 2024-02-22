import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.preprocessing import deseq2_norm


data_in = pd.read_csv('./results_new/combined_raw_reads.csv')

data_in = data_in.set_index('Gene_ID')

samples_to_process = ['LB2B51488', 'LB2B51495', 'LB2B51501', 'LB2B51518', './PURO_Scra_Cntrl.bam', 'LB2B51525', './PURO_Neg_Cntrl.bam', './Drug_Negative.bam', 'LB2B51532', 'LB2B51549', 'LB2B51556', 'LB2B51563']

sample_info = ['kd1', 'kd2', 'kd3', 'scr_control_1', 'scr_control_2', 'neg_control_1', 'neg_control_2', 'drug_neg', 'low_li', 'high_li', 'low_mito', 'high_mito']

samples_required = ['low_li', 'high_li', 'neg_control_1', 'neg_control_2', 'drug_neg']
sample_vals = []
for i, x in enumerate(sample_info):
    if x in samples_required:
        sample_vals.append(samples_to_process[i])

data_in = data_in[sample_vals]

data_in = data_in[data_in.sum(axis=1) > 0]
data_in = data_in.T


condition_list = []
for index in data_in.index:
    val = np.where(np.array(samples_to_process) == index)[0][0]
    index = sample_info[val]
    if 'neg' in index:
        condition_list.append('C')
    else:
        condition_list.append('RS')


metadata_ = pd.DataFrame(zip(data_in.index, condition_list), columns=['Sample', 'Condition'])
metadata_ = metadata_.set_index('Sample')


dds = DeseqDataSet(counts=data_in, metadata=metadata_, design_factors='Condition')

dds.deseq2()
stat_res = DeseqStats(dds, n_cpus=8, contrast=['Condition', 'RS', 'C'])
stat_res.summary()
res = stat_res.results_df

# print(res)

res.to_csv('./results_new/deseq_licl_vs_3_control.tsv', sep='\t')





