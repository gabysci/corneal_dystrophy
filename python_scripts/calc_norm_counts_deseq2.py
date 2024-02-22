import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.preprocessing import deseq2_norm


# data_in = pd.read_csv('./results_dante/RAW_COUNTS/RAW_counts_matrix.csv')
data_in = pd.read_csv('./results_new/combined_raw_reads_2.csv')

data_in = data_in.set_index('Gene_ID')

data_in = data_in[data_in.sum(axis=1) > 0]
data_in = data_in.T

normalized_counts = deseq2_norm(data_in)[0]

normalized_counts.to_csv('./results_new/normalized_counts_all_2.tsv', sep='\t')