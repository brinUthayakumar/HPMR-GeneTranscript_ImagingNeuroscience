import pandas as pd
from gseapy import prerank
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import numpy as np

# Use perm from 
gene_expression_file = 'LH_Schaefer200GeneExpression.csv'  
brain_map_file = 'BPSurrMap_forOligoSurrMapWorkflowDiag.csv'            
gene_sets_file = 'geneset_LAKE.gmt'                       

gene_expr = pd.read_csv(gene_expression_file)

phenotype = pd.read_csv(brain_map_file)


common_regions = gene_expr.index.intersection(phenotype.index)
gene_expr = gene_expr.loc[common_regions]
phenotype = phenotype.loc[common_regions]

# generate correlations using the surrogate map and genes to run prerank on 
gene_scores = gene_expr.apply(lambda col: spearmanr(col, phenotype).correlation, axis=0)


prerank_df = pd.DataFrame({
    'Gene': gene_scores.index,
    'Score': gene_scores.values
}).sort_values('Score', ascending=False)


gsea_res = prerank(
    rnk=prerank_df,
    gene_sets=gene_sets_file,
    outdir=None,        
    verbose=False
)


term = 'Oli'  
ep = gsea_res.plot(terms=term, figsize=(6,5))

# Save the cleaned plot
ep.savefig('GSEA_Oli.pdf', bbox_inches='tight')
