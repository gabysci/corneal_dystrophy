import pandas as pd
import numpy as np

all_variants = pd.read_excel('./filtered_variants_patients.xlsx')
genes_o_i = pd.read_excel('./fullList_Literature_TGFBI_associatedGenes.xlsx')

genes_o_i = list(genes_o_i['Gene (orange-altered in KD)'])

# Filter by genes of interest list
filtered_variants_genes_o_i = all_variants[all_variants['Gene'].isin(genes_o_i)]
filtered_variants_genes_o_i.to_csv('filtered_variants_patients_genes_o_i.csv')

gene_o_i_variants = pd.DataFrame()
gene_o_i_variants['Gene'], gene_o_i_variants['Number of Variants'] = np.unique(np.array(list(filtered_variants_genes_o_i['Gene'])), return_counts=True)

# print(gene_o_i_variants)
gene_o_i_variants.to_csv('genes_o_i_containing_gcd_patient_variants.csv')
