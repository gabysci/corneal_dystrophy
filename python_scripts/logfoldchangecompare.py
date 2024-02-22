import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


log2fc_threshold = 1.5

knock_data_deseq = pd.read_csv('./results_new/deseq_kd_1_2_3_vs_5_control_all.tsv', sep='\t')

treat_data_deseq = pd.read_csv('./results_new/deseq_treat_vs_3_control.tsv', sep='\t')

knock_data = knock_data_deseq
treat_data = treat_data_deseq

treat_data_df = pd.DataFrame()
treat_data_df['Gene_ID_TR'] = treat_data['Gene_ID']
treat_data_df['PVAL_TR'] = treat_data['pvalue']
treat_data_df['Log2FC_TR'] = treat_data['log2FoldChange']


combined_data = pd.DataFrame()
combined_data['Gene_ID'] = knock_data['Gene_ID']
combined_data['Log2FC_KD'] = knock_data['log2FoldChange']
combined_data['PVAL_KD'] = knock_data['pvalue']

combined_data = combined_data.merge(treat_data_df, left_on='Gene_ID', right_on='Gene_ID_TR', how='inner')

combined_data = combined_data.drop('Gene_ID_TR', axis=1)

# filter data on pvals
combined_data = combined_data[(combined_data["PVAL_KD"]<0.05) & (combined_data["PVAL_TR"]<0.05)]

plt.figure(figsize=(10, 10))

upreg = []
downreg = []
notsig = []
upreg_kd_downreg_tr = []
upreg_tr_downreg_kd = []

list_genes = list(combined_data["Gene_ID"])
list_pval_kd = list(combined_data["Log2FC_KD"])
list_pval_tr = list(combined_data["Log2FC_TR"])


for x in range(len(list_pval_kd)):
    if list_pval_kd[x] <= -1*log2fc_threshold and list_pval_tr[x] <= -1*log2fc_threshold:
        downreg.append(x)
    elif list_pval_kd[x] >= log2fc_threshold and list_pval_tr[x] >= log2fc_threshold:
        upreg.append(x)
    elif list_pval_kd[x] >= log2fc_threshold and list_pval_tr[x] <= -1*log2fc_threshold:
        upreg_kd_downreg_tr.append(x)
    elif list_pval_kd[x] <= -1*log2fc_threshold and list_pval_tr[x] >= log2fc_threshold:
        upreg_tr_downreg_kd.append(x)
    else:
        notsig.append(x)
        
list_pval_kd = np.array(list_pval_kd)
list_pval_tr = np.array(list_pval_tr)

plt.scatter(list_pval_kd[upreg], list_pval_tr[upreg], label='log2fc > ' + str(log2fc_threshold) + ' (KD and TR)', c='green')
plt.scatter(list_pval_kd[downreg], list_pval_tr[downreg], label='log2fc < -'  + str(log2fc_threshold) +  ' (KD and TR)', c='red')
plt.scatter(list_pval_kd[upreg_kd_downreg_tr], list_pval_tr[upreg_kd_downreg_tr], label='log2fc > ' + str(log2fc_threshold) + ' (KD) and log2fc < -' + str(log2fc_threshold) + ' (TR)', c='blue')
plt.scatter(list_pval_kd[upreg_tr_downreg_kd], list_pval_tr[upreg_tr_downreg_kd], label='log2fc < -' + str(log2fc_threshold) + ' (KD) and log2fc > ' + str(log2fc_threshold) + ' (TR)', c='orange')
plt.scatter(list_pval_kd[notsig], list_pval_tr[notsig], label='-' + str(log2fc_threshold) + ' < log2fc < ' + str(log2fc_threshold) + ' (KD or treatment)', c='gray')

plt.legend()
plt.title('Log2FC in Knockdown (KD) and Treatment (TR) Samples (pval<0.05)')
plt.xlabel('Log2FC_KD_vs_control')
plt.ylabel('Log2FC_TR_vs_control')
plt.savefig('./log2fc_compare/deseq_log2fc_compare_thr' + str(log2fc_threshold).replace('.','_') + '.png', dpi=600)

upreg_gene_list = []
for upr in upreg:
    upreg_gene_list.append(list_genes[upr])

downreg_gene_list = []
for dwr in downreg:
    downreg_gene_list.append(list_genes[dwr])

upr_downreg_gene_list = []
for upr_dwr in upreg_kd_downreg_tr:
    upr_downreg_gene_list.append(list_genes[upr_dwr])

dwr_upreg_gene_list = []
for dwr_upr in upreg_tr_downreg_kd:
    dwr_upreg_gene_list.append(list_genes[dwr_upr])


file = open('./log2fc_compare/deseq_upreg_genes_KD_upreg_genes_TR_thr' + str(log2fc_threshold).replace('.','_') + '_2.txt', 'w')
for gene_ in upreg_gene_list:
    file.write(gene_ + "\n")
file.close()

file = open('./log2fc_compare/deseq_downreg_genes_KD_downreg_genes_TR_thr' + str(log2fc_threshold).replace('.','_') + '_2.txt', 'w')
for gene_ in downreg_gene_list:
    file.write(gene_ + "\n")
file.close()

file = open('./log2fc_compare/deseq_upreg_genes_KD_downreg_genes_TR_thr' + str(log2fc_threshold).replace('.','_') + '_2.txt', 'w')
for gene_ in upr_downreg_gene_list:
    file.write(gene_ + "\n")
file.close()

file = open('./log2fc_compare/deseq_downreg_genes_KD_upreg_genes_TR_thr' + str(log2fc_threshold).replace('.','_') + '_2.txt', 'w')
for gene_ in dwr_upreg_gene_list:
    file.write(gene_ + "\n")
file.close()
