from brainsmash.mapgen.base import Base 
from brainsmash.workbench.geo import cortex
from brainsmash.workbench.geo import parcellate 
from brainsmash.mapgen.eval import base_fit
from brainsmash.mapgen.stats import spearmanr, pairwise_r, pearsonr
import numpy as np
from scipy import stats
from brainsmash.mapgen.stats import nonparp
import pandas as pd
import gseapy as gp
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import SpatialAutocorrelationGSEA
import glob




# BICARB -------------------------------
ES_nullDistribution_temp = np.array([])

files = sorted(glob.glob("nullDistSchaefer100/NullDistSchaefer100_Bic_nperm_1303*.csv"))
arrays = [pd.read_csv(f) for f in files]
ES_nullDistribution_temp = pd.concat(arrays,ignore_index=True).to_numpy()[:100000]
print(ES_nullDistribution_temp.shape)

celltypes = ('Ast','End','Ex1','Ex2','Ex3a','Ex3b','Ex3c','Ex3d','Ex3e','Ex4','Ex5a','Ex5b','Ex6a','Ex6b','Ex8','In1a','In1b',
             'In1c','In2','In3','In4a','In4b','In6a','In6b','In7','In8','Mic','OPC','Oli','Per')
ES_nullDistribution = pd.DataFrame(ES_nullDistribution_temp,columns = list(celltypes))

lake_gmt = "geneset_LAKE.gmt"
geneExpressionCSV = "Schaefer100_lh_expression.csv"
brain_phenotype_file = "formatted_bstemNormBicSchaefer100.txt"

Bic_empiricalSpearman = SpatialAutocorrelationGSEA.empiricalSpearman(brain_phenotype_file,geneExpressionCSV)

rnk = Bic_empiricalSpearman.sort_values(by = 'EmpiricalSpearman',ascending = False)

empiricalPrerank = gp.prerank(rnk = rnk,
                     gene_sets = "geneset_LAKE.gmt",
                     ascending = False)


correctedPvalues_temp = []
correctedPvalues = []
for geneSetID in range(len(celltypes)):
    n = len(ES_nullDistribution)
    correctedPvalues_temp = np.sum(np.abs(ES_nullDistribution.iloc[:,geneSetID]) > abs(
                               empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[geneSetID,'ES']))/n
    correctedPvalues.append(correctedPvalues_temp)

geneSetTerm = empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'Term']
enrichmentScore = empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'ES']
df_correctedPvalues = pd.DataFrame(correctedPvalues)
# FDR correction 
FDRCorrectedPvals = pd.Series(stats.false_discovery_control(df_correctedPvalues[0],method = 'BH'))
df_correctedLabelledPvalues = pd.concat([geneSetTerm,enrichmentScore,df_correctedPvalues,FDRCorrectedPvals],ignore_index=True,axis=1)
df_correctedLabelledPvalues = df_correctedLabelledPvalues.set_axis(['geneSetTerm','EmpiricalES','SA_CorrectedPval','FDRandSA_Corrected'],
                                                                   axis = 1)
print(df_correctedLabelledPvalues)



# BICARB -------------------------------

significantGenesets = [i for i, j in zip(df_correctedLabelledPvalues['geneSetTerm'],
                                         df_correctedLabelledPvalues['SA_CorrectedPval']) if j <0.05]
print(significantGenesets)
EP = empiricalPrerank.plot(terms = significantGenesets,
                           figsize = (3,4))
EP.savefig('OneMill_Schaefer100-bstemNormalized_Bic_SA_0.05.pdf',format ='pdf',bbox_inches="tight")


# LACTATE -------------------------------
ES_nullDistribution_temp = np.array([])

files = sorted(glob.glob("nullDistSchaefer100/NullDistSchaefer100_Lac_nperm_1303*.csv"))
arrays = [pd.read_csv(f) for f in files]
ES_nullDistribution_temp = pd.concat(arrays,ignore_index=True).to_numpy()[:100000]
print(ES_nullDistribution_temp.shape)

celltypes = ('Ast','End','Ex1','Ex2','Ex3a','Ex3b','Ex3c','Ex3d','Ex3e','Ex4','Ex5a','Ex5b','Ex6a','Ex6b','Ex8','In1a','In1b',
             'In1c','In2','In3','In4a','In4b','In6a','In6b','In7','In8','Mic','OPC','Oli','Per')
ES_nullDistribution = pd.DataFrame(ES_nullDistribution_temp,columns = list(celltypes))

lake_gmt = "geneset_LAKE.gmt"
geneExpressionCSV = "Schaefer100_lh_expression.csv"
brain_phenotype_file = "formatted_bstemNormLacSchaefer100.txt"

Lac_empiricalSpearman = SpatialAutocorrelationGSEA.empiricalSpearman(brain_phenotype_file,geneExpressionCSV)

rnk = Lac_empiricalSpearman.sort_values(by = 'EmpiricalSpearman',ascending = False)

empiricalPrerank = gp.prerank(rnk = rnk,
                     gene_sets = "geneset_LAKE.gmt",
                     ascending = False)


correctedPvalues_temp = []
correctedPvalues = []
for geneSetID in range(len(celltypes)):
    n = len(ES_nullDistribution)
    correctedPvalues_temp = np.sum(np.abs(ES_nullDistribution.iloc[:,geneSetID]) > abs(
                               empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[geneSetID,'ES']))/n
    correctedPvalues.append(correctedPvalues_temp)

geneSetTerm = empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'Term']
enrichmentScore = empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'ES']
df_correctedPvalues = pd.DataFrame(correctedPvalues)
# FDR correction 
FDRCorrectedPvals = pd.Series(stats.false_discovery_control(df_correctedPvalues[0],method = 'BH'))
df_correctedLabelledPvalues = pd.concat([geneSetTerm,enrichmentScore,df_correctedPvalues,FDRCorrectedPvals],ignore_index=True,axis=1)
df_correctedLabelledPvalues = df_correctedLabelledPvalues.set_axis(['geneSetTerm','EmpiricalES','SA_CorrectedPval','FDRandSA_Corrected'],
                                                                   axis = 1)
print(df_correctedLabelledPvalues)



# LACTATE -------------------------------

significantGenesets = [i for i, j in zip(df_correctedLabelledPvalues['geneSetTerm'],
                                         df_correctedLabelledPvalues['SA_CorrectedPval']) if j <0.05]
print(significantGenesets)
EP = empiricalPrerank.plot(terms = significantGenesets,
                           figsize = (3,4))
EP.savefig('OneMill_Schaefer100-bstemNormalized_Lac_SA_0.05.pdf',format ='pdf',bbox_inches="tight")




# PYRUVATE  -------------------------------

ES_nullDistribution_temp = np.array([])

files = sorted(glob.glob("nullDistSchaefer100/NullDistSchaefer100_Pyr_nperm_1303*.csv"))
arrays = [pd.read_csv(f) for f in files]
ES_nullDistribution_temp = pd.concat(arrays,ignore_index=True).to_numpy()[:100000]
print(ES_nullDistribution_temp.shape)

celltypes = ('Ast','End','Ex1','Ex2','Ex3a','Ex3b','Ex3c','Ex3d','Ex3e','Ex4','Ex5a','Ex5b','Ex6a','Ex6b','Ex8','In1a','In1b',
             'In1c','In2','In3','In4a','In4b','In6a','In6b','In7','In8','Mic','OPC','Oli','Per')
ES_nullDistribution = pd.DataFrame(ES_nullDistribution_temp,columns = list(celltypes))

lake_gmt = "geneset_LAKE.gmt"
geneExpressionCSV = "Schaefer100_lh_expression.csv"
brain_phenotype_file = "formatted_bstemNormPyrSchaefer100.txt"

Pyr_empiricalSpearman = SpatialAutocorrelationGSEA.empiricalSpearman(brain_phenotype_file,geneExpressionCSV)

rnk = Pyr_empiricalSpearman.sort_values(by = 'EmpiricalSpearman',ascending = False)

empiricalPrerank = gp.prerank(rnk = rnk,
                     gene_sets = "geneset_LAKE.gmt",
                     ascending = False)


correctedPvalues_temp = []
correctedPvalues = []
for geneSetID in range(len(celltypes)):
    n = len(ES_nullDistribution)
    correctedPvalues_temp = np.sum(np.abs(ES_nullDistribution.iloc[:,geneSetID]) > abs(
                               empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[geneSetID,'ES']))/n
    correctedPvalues.append(correctedPvalues_temp)

geneSetTerm = empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'Term']
enrichmentScore = empiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'ES']
df_correctedPvalues = pd.DataFrame(correctedPvalues)
# FDR correction 
FDRCorrectedPvals = pd.Series(stats.false_discovery_control(df_correctedPvalues[0],method = 'BH'))
df_correctedLabelledPvalues = pd.concat([geneSetTerm,enrichmentScore,df_correctedPvalues,FDRCorrectedPvals],ignore_index=True,axis=1)
df_correctedLabelledPvalues = df_correctedLabelledPvalues.set_axis(['geneSetTerm','EmpiricalES','SA_CorrectedPval','FDRandSA_Corrected'],
                                                                   axis = 1)
print(df_correctedLabelledPvalues)



# PYRUVATE -------------------------------

significantGenesets = [i for i, j in zip(df_correctedLabelledPvalues['geneSetTerm'],
                                         df_correctedLabelledPvalues['SA_CorrectedPval']) if j <0.05]
print(significantGenesets)
EP = empiricalPrerank.plot(terms = significantGenesets,
                           figsize = (3,4))
EP.savefig('OneMill_Schaefer100-bstemNormalized_Pyr_SA_0.05.pdf',format ='pdf',bbox_inches="tight")




