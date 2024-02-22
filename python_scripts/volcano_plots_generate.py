import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


data_lfc = pd.read_csv('./results_new/deseq_kd_1_2_3_vs_5_control_all.tsv', sep='\t')
pval = np.array(data_lfc['pvalue'])
lfc = np.array(data_lfc['log2FoldChange'])

p_vals_no_sig = np.where(pval > 0.05)
p_vals_sig = np.where(pval <= 0.05)
  
log2fc_vals_pos = np.where(lfc >= 1.5)
log2fc_vals_neg = np.where(lfc <= -1.5)
 
pos_sig = []
pos_sig_log2fc = []
pos_sig_pval = []
neg_sig = []
neg_sig_log2fc = []
neg_sig_pval = []
data_lfc['Gene_ID'] = data_lfc['Gene_ID'].fillna(0)
  
for val in np.array(p_vals_sig[0]):
    if val in np.array(log2fc_vals_pos):
        pos_sig_log2fc.append(data_lfc['log2FoldChange'][val])
        pos_sig_pval.append(data_lfc['pvalue'][val])
        pos_sig.append(data_lfc['Gene_ID'][val])

    if val in np.array(log2fc_vals_neg):
        neg_sig_log2fc.append(data_lfc['log2FoldChange'][val])
        neg_sig_pval.append(data_lfc['pvalue'][val])
        neg_sig.append(data_lfc['Gene_ID'][val])

data_out_pos = pd.DataFrame()
data_out_pos['Gene'] = pos_sig
data_out_pos['Log2FC'] = pos_sig_log2fc
data_out_pos['PVal'] = pos_sig_pval
data_out_pos.to_csv('./de_gene_lists/deseq_kd_1_2_3_vs_5_control_sig_upreg_genes.csv', index=False)
data_out_neg = pd.DataFrame()
data_out_neg['Gene'] = neg_sig
data_out_neg['Log2FC'] = neg_sig_log2fc
data_out_neg['PVal'] = neg_sig_pval
data_out_neg.to_csv('./de_gene_lists/deseq_kd_1_2_3_vs_5_control_sig_downreg_genes.csv', index=False)
  
plt.figure(figsize=(10, 10))
plt.title('Volcano Plot for KD (1,2,3) vs Control Exp (log2FC cutoff 1.5, p-val<0.05)')
plt.xlabel('log2FC(KD_vs_Control)')
plt.ylabel('-log10(P-val)')
plt.scatter(lfc[log2fc_vals_pos], -1*(np.log10(pval[log2fc_vals_pos])),
            c='red', s=4, label='upreg')
plt.scatter(lfc[log2fc_vals_neg], -1*(np.log10(pval[log2fc_vals_neg])),
            c='blue', s=4, label='downreg')
plt.scatter(lfc[p_vals_no_sig], -1*(np.log10(pval[p_vals_no_sig])),
            c='grey', s=4, label='not sig')
plt.legend()
plt.savefig('./volcano_plots/volcano_plot_deseq_kd_1_2_3_vs_5_control_all.png', dpi=300)
