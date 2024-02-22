import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns


samples_to_process = ['LB2B51488', 'LB2B51495', 'LB2B51501', 'LB2B51518', './PURO_Scra_Cntrl.bam', 'LB2B51525', './PURO_Neg_Cntrl.bam', './Drug_Negative.bam', 'LB2B51532', 'LB2B51549', 'LB2B51556', 'LB2B51563']

sample_info = ['kd1', 'kd2', 'kd3', 'scr_control_1', 'scr_control_2', 'neg_control_1', 'neg_control_2', 'drug_neg', 'low_li', 'high_li', 'low_mito', 'high_mito']


corr_matrix = np.ones((len(sample_info), len(sample_info)))

for x in range(len(sample_info)):
    for y in range(x+1, len(sample_info)):
        data_in = pd.read_csv('./results_new/normalized_counts.tsv', sep='\t', header=0, index_col=0)
        data_in = data_in.reindex(samples_to_process)

        data_in = data_in.set_axis(sample_info, axis='rows')
    
        print(str(sample_info[x]), str(sample_info[y]))
        samples_required = []
        samples_required.append(sample_info[x])
        samples_required.append(sample_info[y])
        sort_by_val = samples_required[0]
        other_val = samples_required[1]

        data_in = data_in[data_in.sum(axis=1) > 0]
        data_in = data_in.T[samples_required]
        
        sorted_df = data_in.sort_values(by=sort_by_val)
        
        std_1 = np.std(sorted_df[sort_by_val])        
        std_2 = np.std(sorted_df[other_val])
        mean_1 = np.mean(sorted_df[sort_by_val])        
        mean_2 = np.mean(sorted_df[other_val])
        upper_1 = mean_1 + (2*std_1)
        lower_1 = mean_1 - (2*std_1)
        upper_2 = mean_2 + (2*std_2)
        lower_2 = mean_2 - (2*std_2)
        
        sorted_df = sorted_df[sorted_df[sort_by_val] < upper_1]
        sorted_df = sorted_df[sorted_df[sort_by_val] > lower_1]
        sorted_df = sorted_df[sorted_df[other_val] < upper_2]
        sorted_df = sorted_df[sorted_df[other_val] > lower_2]
        
        pearson_ = stats.pearsonr(sorted_df[sort_by_val], sorted_df[other_val])
        
        plt.figure(figsize=(15, 7))
        plt.scatter(sorted_df[sort_by_val], sorted_df[other_val], alpha=0.3, s=5)
        a, b = np.polyfit(sorted_df[sort_by_val], sorted_df[other_val], 1)
        plt.plot(sorted_df[sort_by_val], a*sorted_df[sort_by_val]+b)
        plt.title(str(sample_info[x]) + ' vs ' + str(sample_info[y]) + ' - Pearson: ' + str(np.round(pearson_.statistic, 3)) + ', p-val=' + str(pearson_.pvalue))    
        plt.savefig('corr_plots/corr_plot_' + str(sample_info[x]) + ' vs ' + str(sample_info[y]) + '.png', dpi=300)
        
        corr_matrix[x, y] = np.round(pearson_.statistic, 4)
        corr_matrix[y, x] = np.round(pearson_.statistic, 4)

plt.figure(figsize=(10, 10))
sns.heatmap(corr_matrix, annot=True, annot_kws={"fontsize":6}, fmt='.3f', xticklabels=sample_info, yticklabels=sample_info)

plt.savefig('corr_plots/corr_plot_pearson_all.png', dpi=600)
        

