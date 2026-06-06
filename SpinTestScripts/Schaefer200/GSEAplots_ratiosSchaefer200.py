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
import glob #For bash like path filenames




# BICARB TO PYRUVATE RATIO -------------------------------


files = sorted(glob.glob("nullDistSchaefer200ratios/NullDistSchaefer200_BP_nperm_1310*.csv"))
arrays = [pd.read_csv(f) for f in files]
ES_nullDistribution_temp = pd.concat(arrays,ignore_index=True).to_numpy()[:100000]
print(ES_nullDistribution_temp.shape)

celltypes = ('Ast','End','Ex1','Ex2','Ex3a','Ex3b','Ex3c','Ex3d','Ex3e','Ex4','Ex5a','Ex5b','Ex6a','Ex6b','Ex8','In1a','In1b',
             'In1c','In2','In3','In4a','In4b','In6a','In6b','In7','In8','Mic','OPC','Oli','Per')
ES_nullDistribution = pd.DataFrame(ES_nullDistribution_temp,columns = list(celltypes))

lake_gmt = "geneset_LAKE.gmt"
geneExpressionCSV = "LH_Schaefer200GeneExpression.csv"
brain_phenotype_file = "formatted_segstatsFilesSchaefer200/formatted_BP-Schaefer200.txt"

BP_empiricalSpearman = SpatialAutocorrelationGSEA.empiricalSpearman(brain_phenotype_file,geneExpressionCSV)

rnk = BP_empiricalSpearman.sort_values(by = 'EmpiricalSpearman',ascending = False)

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



# BICARB TO PYRUVATE RATIO -------------------------------

significantGenesets = [i for i, j in zip(df_correctedLabelledPvalues['geneSetTerm'],
                                         df_correctedLabelledPvalues['SA_CorrectedPval']) if j <0.05]
print(significantGenesets)
EP = empiricalPrerank.plot(terms = significantGenesets,
                           figsize = (3,4))
EP.savefig('OneMill_Schaefer200-BP_SA_0.05.pdf',format ='pdf',bbox_inches="tight")



# LACTATE TO PYRUVATE RATIO ES plot  -------------------------------

files = sorted(glob.glob("nullDistSchaefer200ratios/NullDistSchaefer200_LP_nperm_1310*.csv"))
arrays = [pd.read_csv(f) for f in files]
ES_nullDistribution_temp = pd.concat(arrays,ignore_index=True).to_numpy()[:100000]
print(ES_nullDistribution_temp.shape)

celltypes = ('Ast','End','Ex1','Ex2','Ex3a','Ex3b','Ex3c','Ex3d','Ex3e','Ex4','Ex5a','Ex5b','Ex6a','Ex6b','Ex8','In1a','In1b',
             'In1c','In2','In3','In4a','In4b','In6a','In6b','In7','In8','Mic','OPC','Oli','Per')
ES_nullDistribution = pd.DataFrame(ES_nullDistribution_temp,columns = list(celltypes))

lake_gmt = "geneset_LAKE.gmt"
geneExpressionCSV = "LH_Schaefer200GeneExpression.csv"
brain_phenotype_file = "formatted_segstatsFilesSchaefer200/formatted_LP-Schaefer200.txt"

LP_empiricalSpearman = SpatialAutocorrelationGSEA.empiricalSpearman(brain_phenotype_file,geneExpressionCSV)

rnk = LP_empiricalSpearman.sort_values(by = 'EmpiricalSpearman',ascending = False)

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



# LACTATE TO PYRUVATE RATIO ES plot  -------------------------------

significantGenesets = [i for i, j in zip(df_correctedLabelledPvalues['geneSetTerm'],
                                         df_correctedLabelledPvalues['FDRandSA_Corrected']) if j <0.25]

terms = empiricalPrerank.res2d.Term

EP = empiricalPrerank.plot(terms = significantGenesets,
                           figsize = (3,4))
EP.savefig('OneMill_Schaefer200-LP_fdr0.25.pdf',format ='pdf',bbox_inches="tight")





# LACTATE TO BICARBONATE RATIO ES plot  -------------------------------

files = sorted(glob.glob("nullDistSchaefer200ratios/NullDistSchaefer200_LB_nperm_1310*.csv"))
arrays = [pd.read_csv(f) for f in files]
ES_nullDistribution_temp = pd.concat(arrays,ignore_index=True).to_numpy()[:100000]
print(ES_nullDistribution_temp.shape)


celltypes = ('Ast','End','Ex1','Ex2','Ex3a','Ex3b','Ex3c','Ex3d','Ex3e','Ex4','Ex5a','Ex5b','Ex6a','Ex6b','Ex8','In1a','In1b',
             'In1c','In2','In3','In4a','In4b','In6a','In6b','In7','In8','Mic','OPC','Oli','Per')
ES_nullDistribution = pd.DataFrame(ES_nullDistribution_temp,columns = list(celltypes))

lake_gmt = "geneset_LAKE.gmt"
geneExpressionCSV = "LH_Schaefer200GeneExpression.csv"
brain_phenotype_file = "formatted_segstatsFilesSchaefer200/formatted_LB-Schaefer200.txt"

LB_empiricalSpearman = SpatialAutocorrelationGSEA.empiricalSpearman(brain_phenotype_file,geneExpressionCSV)

rnk = LB_empiricalSpearman.sort_values(by = 'EmpiricalSpearman',ascending = False)

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



# LACTATE TO BICARBONATE RATIO ES plot  -------------------------------

significantGenesets = [i for i, j in zip(df_correctedLabelledPvalues['geneSetTerm'],
                                         df_correctedLabelledPvalues['FDRandSA_Corrected']) if j <0.25]

terms = empiricalPrerank.res2d.Term

EP = empiricalPrerank.plot(terms = significantGenesets,
                           figsize = (3,4))
EP.savefig('OneMill_Schaefer200-LB_fdr0.25.pdf',format ='pdf',bbox_inches="tight")




