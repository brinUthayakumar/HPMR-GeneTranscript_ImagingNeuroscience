#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy import stats
import gseapy as gp
from gseapy import barplot, dotplot
import time 
import matplotlib.pyplot as plt
from collections import Counter 
# Bring in the spatial autocorrelation jupyter notebook to be able to use its functions 
get_ipython().run_line_magic('run', 'SpatialAutocorrelationCorrelationsFunctions.ipynb')


# In[2]:


################################### FOR BP RATIO ##########################################

geneExpressionCSV = "LH_Schaefer200GeneExpression.csv"
brain_phenotype_file = "C13Data/formatted_BP-Schaefer200.txt"
dist_mat_file = "Schaefer200_L_GeodesicParcelDistMat_midThickness_Left.txt"
nh = 10
pv = 70
npermutations = 1000

SurrMapsBP = generateSurrogateMaps(brain_phenotype_file,dist_mat_file,npermutations,
                                   nh,pv)
geneWiseOutputBP = geneWiseNull(geneExpressionCSV,SurrMapsBP)
SAcorrectedSpearmanBP = spatialAutocorrCorrectedPval(geneWiseOutputBP,geneExpressionCSV,brain_phenotype_file)

# Uncomment below to write to csv
#SAcorrectedSpearmanBP.to_csv('SAcorrectedSpearmanBP-Schaefer200_1000perm.csv', index = False)
################################### FOR LP RATIO ##########################################

geneExpressionCSV = "LH_Schaefer200GeneExpression.csv"
brain_phenotype_file = "C13Data/formatted_LP-Schaefer200.txt"
dist_mat_file = "Schaefer200_L_GeodesicParcelDistMat_midThickness_Left.txt"
nh = 10
pv = 70
npermutations = 1000

SurrMapsLP = generateSurrogateMaps(brain_phenotype_file,dist_mat_file,npermutations,
                                   nh,pv)
geneWiseOutputLP = geneWiseNull(geneExpressionCSV,SurrMapsLP)
SAcorrectedSpearmanLP = spatialAutocorrCorrectedPval(geneWiseOutputLP,geneExpressionCSV,brain_phenotype_file)

# Uncomment below to write to csv
#SAcorrectedSpearmanLP.to_csv('SAcorrectedSpearmanLP-Schaefer200_1000perm.csv', index = False)


################################### FOR LB RATIO ##########################################

geneExpressionCSV = "LH_Schaefer200GeneExpression.csv"
brain_phenotype_file = "C13Data/formatted_LB-Schaefer200.txt"
dist_mat_file = "Schaefer200_L_GeodesicParcelDistMat_midThickness_Left.txt"
nh = 10
pv = 70
npermutations = 1000

SurrMapsLB = generateSurrogateMaps(brain_phenotype_file,dist_mat_file,npermutations,
                                   nh,pv)
geneWiseOutputLB = geneWiseNull(geneExpressionCSV,SurrMapsLB)
SAcorrectedSpearmanLB = spatialAutocorrCorrectedPval(geneWiseOutputLB,geneExpressionCSV,brain_phenotype_file)

# Uncomment below to write to csv
#SAcorrectedSpearmanLB.to_csv('SAcorrectedSpearmanLB-Schaefer200_1000perm.csv', index = False)


# In[11]:


# Results of GOI for BP
df_significantCorrelation_BP = df_correctedLabelledPvalues_BP[df_correctedLabelledPvalues_BP['GeneLabel'].isin(df_correctedLabelledPvalues_BP)]
df_significantCorrelation_BP.sort_values(by = 'EmpiricalSpearman',ascending = False)
genesOfInterest = ['SLC16A1','SLC16A7','LDHA','LDHB','HIF1A']
print(df_correctedLabelledPvalues_BP[df_correctedLabelledPvalues_BP['GeneLabel'].isin(genesOfInterest)])


# In[12]:


# Results of GOI for LP
df_significantCorrelation_LP = df_correctedLabelledPvalues_LP[df_correctedLabelledPvalues_LP['GeneLabel'].isin(df_correctedLabelledPvalues_LP)]
df_significantCorrelation_LP.sort_values(by = 'EmpiricalSpearman',ascending = False)
genesOfInterest = ['SLC16A1','SLC16A7','LDHA','LDHB','HIF1A']
print(df_correctedLabelledPvalues_LP[df_correctedLabelledPvalues_LP['GeneLabel'].isin(genesOfInterest)])


# In[13]:


# Results of GOI for LB
df_significantCorrelation_LB = df_correctedLabelledPvalues_LB[df_correctedLabelledPvalues_LB['GeneLabel'].isin(df_correctedLabelledPvalues_LB)]
df_significantCorrelation_LB.sort_values(by = 'EmpiricalSpearman',ascending = False)
genesOfInterest = ['SLC16A1','SLC16A7','LDHA','LDHB','HIF1A']
print(df_correctedLabelledPvalues_LB[df_correctedLabelledPvalues_LB['GeneLabel'].isin(genesOfInterest)])


# In[14]:


# Top 5 highest and lowest Empirical spearmans that are also significant with FDR and SA correction for BP
significantGenes = [i for i, j in zip(df_correctedLabelledPvalues_BP['GeneLabel'],
                                     df_correctedLabelledPvalues_BP['FDRandSA_Corrected']) if j <0.10]
df_SignificantCorrectedGenes_BP = (df_correctedLabelledPvalues_BP[df_correctedLabelledPvalues_BP['GeneLabel'].isin(significantGenes)])

print(df_SignificantCorrectedGenes_BP.sort_values('EmpiricalSpearman',ascending = False).head(5))
print(df_SignificantCorrectedGenes_BP.sort_values('EmpiricalSpearman',ascending = True).head(5))


# In[15]:


# Top 5 highest and lowest Empirical spearmans that are also significant with FDR and SA correction for LP
significantGenes = [i for i, j in zip(df_correctedLabelledPvalues_LP['GeneLabel'],
                                     df_correctedLabelledPvalues_LP['FDRandSA_Corrected']) if j <0.1]
df_SignificantCorrectedGenes_LP = (df_correctedLabelledPvalues_LP[df_correctedLabelledPvalues_LP['GeneLabel'].isin(significantGenes)])

print(df_SignificantCorrectedGenes_LP.sort_values('EmpiricalSpearman',ascending = False).head(5))
print(df_SignificantCorrectedGenes_LP.sort_values('EmpiricalSpearman',ascending = True).head(5))


# In[12]:


# Top 5 highest and lowest Empirical spearmans that are also significant with FDR and SA correction for LB
significantGenes = [i for i, j in zip(df_correctedLabelledPvalues_LB['GeneLabel'],
                                     df_correctedLabelledPvalues_LB['FDRandSA_Corrected']) if j <0.05]
df_SignificantCorrectedGenes_LB = (df_correctedLabelledPvalues_LB[df_correctedLabelledPvalues_LB['GeneLabel'].isin(significantGenes)])

print(df_SignificantCorrectedGenes_LB.sort_values('EmpiricalSpearman',ascending = False).head(5))
print(df_SignificantCorrectedGenes_LB.sort_values('EmpiricalSpearman',ascending = True).head(5))


# In[7]:


################################### RESULTS SUMMARY FOR BP RATIO ###########################

# CORRECTED PVALS = SA CORRECTED
# FDR CORRECTED PVALS = FDR CORRECTION OF SA CORRECTED PVALS
BPCorrectedPvals_Lake = pd.read_csv('SpatialAutocorrelationPermutation_CSV/SAcorrectedSpearmanBP-Schaefer200_1000perm.csv')
FDRCorrectedPvals = pd.DataFrame(stats.false_discovery_control(BPCorrectedPvals_Lake['CorrectedPval'],method = 'BH'))

df_correctedLabelledPvalues_BP = pd.concat([BPCorrectedPvals_Lake,FDRCorrectedPvals],ignore_index=True,axis=1)
df_correctedLabelledPvalues_BP
df_correctedLabelledPvalues_BP = df_correctedLabelledPvalues_BP.set_axis(['GeneLabel',
                                                            'EmpiricalSpearman',
                                                            'SA_CorrectedPval',
                                                            'FDRandSA_Corrected'],
                                                           axis = 1)

print(df_correctedLabelledPvalues_BP[df_correctedLabelledPvalues_BP['FDRandSA_Corrected']<0.05].sort_values
      ('EmpiricalSpearman',ascending=False).head(20))


# In[8]:


################################### RESULTS SUMMARY FOR LP RATIO ###########################

# CORRECTED PVALS = SA CORRECTED
# FDR CORRECTED PVALS = FDR CORRECTION OF SA CORRECTED PVALS
LPCorrectedPvals_Lake = pd.read_csv('SpatialAutocorrelationPermutation_CSV/SAcorrectedSpearmanLP-Schaefer200_1000perm.csv')
FDRCorrectedPvals = pd.DataFrame(stats.false_discovery_control(LPCorrectedPvals_Lake['CorrectedPval'],method = 'BH'))

df_correctedLabelledPvalues_LP = pd.concat([LPCorrectedPvals_Lake,FDRCorrectedPvals],ignore_index=True,axis=1)
df_correctedLabelledPvalues_LP
df_correctedLabelledPvalues_LP = df_correctedLabelledPvalues_LP.set_axis(['GeneLabel',
                                                            'EmpiricalSpearman',
                                                            'SA_CorrectedPval',
                                                            'FDRandSA_Corrected'],
                                                           axis = 1)

print(df_correctedLabelledPvalues_LP[df_correctedLabelledPvalues_LP['FDRandSA_Corrected']<0.00005].sort_values
      ('EmpiricalSpearman',ascending=False).head(20))
print(df_correctedLabelledPvalues_LP[df_correctedLabelledPvalues_LP['FDRandSA_Corrected']<0.00005].sort_values
      ('EmpiricalSpearman',ascending=True).head(20))


# In[10]:


################################### RESULTS SUMMARY FOR LB RATIO ###########################

# CORRECTED PVALS = SA CORRECTED
# FDR CORRECTED PVALS = FDR CORRECTION OF SA CORRECTED PVALS
LBCorrectedPvals_Lake = pd.read_csv('SpatialAutocorrelationPermutation_CSV/SAcorrectedSpearmanLB-Schaefer200_1000perm.csv')
FDRCorrectedPvals = pd.DataFrame(stats.false_discovery_control(LBCorrectedPvals_Lake['CorrectedPval'],method = 'BH'))

df_correctedLabelledPvalues_LB = pd.concat([LBCorrectedPvals_Lake,FDRCorrectedPvals],ignore_index=True,axis=1)
df_correctedLabelledPvalues_LB
df_correctedLabelledPvalues_LB = df_correctedLabelledPvalues_LB.set_axis(['GeneLabel',
                                                            'EmpiricalSpearman',
                                                            'SA_CorrectedPval',
                                                            'FDRandSA_Corrected'],
                                                           axis = 1)

print(df_correctedLabelledPvalues_LB[df_correctedLabelledPvalues_LB['FDRandSA_Corrected']<0.00005].sort_values
      ('EmpiricalSpearman',ascending=False).head(20))
print(df_correctedLabelledPvalues_LB[df_correctedLabelledPvalues_LB['FDRandSA_Corrected']<0.00005].sort_values
      ('EmpiricalSpearman',ascending=True).head(20))


# In[17]:


# Save the dataframes of corrected pvalues
df_correctedLabelledPvalues_BP.to_csv('BP_Schaefer200Correlations_n40.csv',index=False)
df_correctedLabelledPvalues_LB.to_csv('LB_Schaefer200Correlations_n40.csv',index=False)
df_correctedLabelledPvalues_LP.to_csv('LP_Schaefer200Correlations_n40.csv',index=False)


# In[18]:


# Plotting most positively correlated significant genes for LP
geneExpressionCSV = 'LH_Schaefer200GeneExpression.csv'
geneData = pd.read_csv(geneExpressionCSV)
brain_phenotype_file = 'C13Data/formatted_LP-Schaefer200.txt'
brainMapFile = np.loadtxt(brain_phenotype_file)

fig, axes = plt.subplots(3,1,figsize=(4,15))

# Find most positive correlations with fdrandSA_corrected pvalues less than 0.00005, get the gene label for that and extract all gene data 
# for each of those correlations
positiveCorrelationData = geneData[
    df_correctedLabelledPvalues_LP[
    df_correctedLabelledPvalues_LP['FDRandSA_Corrected']<0.00005].sort_values('EmpiricalSpearman',ascending=False).head(3)
    ['GeneLabel']]
for i in range(len(positiveCorrelationData.columns)):
    axes[i].scatter(brainMapFile,positiveCorrelationData[positiveCorrelationData.columns[i]]) 
    axes[i].plot(brainMapFile,np.poly1d(np.polyfit(brainMapFile,np.array(
    positiveCorrelationData.iloc[:,i]),1))(brainMapFile))
    axes[i].set_ylabel(positiveCorrelationData.columns[i],fontsize=16)
    axes[2].set_xlabel('Lactate-to-Pyruvate ratio',fontsize=16)
    axes[i].tick_params(axis='x', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12)

plt.tight_layout()
plt.savefig('MostPositiveLP.pdf')
    
# Only 3 significant correlations, so need for both negative and positive correlation plots


# In[19]:


# Plotting most positively correlated significant genes for LB
geneExpressionCSV = 'LH_Schaefer200GeneExpression.csv'
geneData = pd.read_csv(geneExpressionCSV)
brain_phenotype_file = 'C13Data/formatted_LB-Schaefer200.txt'
brainMapFile = np.loadtxt(brain_phenotype_file)

fig, axes = plt.subplots(3,1,figsize=(4,15))

# Find most positive correlations with fdrandSA_corrected pvalues less than 0.00005, get the gene label for that and extract all gene data 
# for each of those correlations
positiveCorrelationData = geneData[
    df_correctedLabelledPvalues_LB[
    df_correctedLabelledPvalues_LB['FDRandSA_Corrected']<0.00005].sort_values('EmpiricalSpearman',ascending=False).head(3)
    ['GeneLabel']]
for i in range(len(positiveCorrelationData.columns)):
    axes[i].scatter(brainMapFile,positiveCorrelationData[positiveCorrelationData.columns[i]]) 
    axes[i].plot(brainMapFile,np.poly1d(np.polyfit(brainMapFile,np.array(
    positiveCorrelationData.iloc[:,i]),1))(brainMapFile))
    axes[i].set_ylabel(positiveCorrelationData.columns[i],fontsize=16)
    axes[2].set_xlabel('Lactate-to-Bicarbonate ratio',fontsize=16)
    axes[i].tick_params(axis='x', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12)

plt.tight_layout()
plt.savefig('MostPositiveLB.pdf')
# Only 3 significant correlations, so need for both negative and positive correlation plots


# In[20]:


# Plotting most negatively correlated significant genes for LP
geneExpressionCSV = 'LH_Schaefer200GeneExpression.csv'
geneData = pd.read_csv(geneExpressionCSV)
brain_phenotype_file = 'C13Data/formatted_LP-Schaefer200.txt'
brainMapFile = np.loadtxt(brain_phenotype_file)

fig, axes = plt.subplots(3,1,figsize=(4,15))

# Find most negative correlations with fdrandSA_corrected pvalues less than 0.00005, get the gene label for that and extract all gene data 
# for each of those correlations
negativeCorrelationData = geneData[
    df_correctedLabelledPvalues_LP[
    df_correctedLabelledPvalues_LP['FDRandSA_Corrected']<0.00005].sort_values('EmpiricalSpearman',ascending=True).head(3)
    ['GeneLabel']]
for i in range(len(negativeCorrelationData.columns)):
    axes[i].scatter(brainMapFile,negativeCorrelationData[negativeCorrelationData.columns[i]]) 
    axes[i].plot(brainMapFile,np.poly1d(np.polyfit(brainMapFile,np.array(
    negativeCorrelationData.iloc[:,i]),1))(brainMapFile))
    axes[i].set_ylabel(negativeCorrelationData.columns[i],fontsize=16)
    axes[2].set_xlabel('Lactate-to-Pyruvate ratio',fontsize=16)
    axes[i].tick_params(axis='x', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12)

plt.tight_layout()
plt.savefig('MostNegativeLP.pdf')
# Only 3 significant correlations, so need for both negative and positive correlation plots


# In[22]:


# Plotting most negatively correlated significant genes for LB
geneExpressionCSV = 'LH_Schaefer200GeneExpression.csv'
geneData = pd.read_csv(geneExpressionCSV)
brain_phenotype_file = 'C13Data/formatted_LB-Schaefer200.txt'
brainMapFile = np.loadtxt(brain_phenotype_file)

fig, axes = plt.subplots(3,1,figsize=(4,15))

# Find most negative correlations with fdrandSA_corrected pvalues less than 0.00005, get the gene label for that and extract all gene data 
# for each of those correlations
negativeCorrelationData = geneData[
    df_correctedLabelledPvalues_LB[
    df_correctedLabelledPvalues_LB['FDRandSA_Corrected']<0.05].sort_values('EmpiricalSpearman',ascending=True).head(2)
    ['GeneLabel']]
negativeCorrelationData_wSLC16A7 = pd.concat([negativeCorrelationData,pd.DataFrame(geneData['SLC16A7'])],axis=1)
for i in range(len(negativeCorrelationData_wSLC16A7.columns)):
    axes[i].scatter(brainMapFile,negativeCorrelationData_wSLC16A7[negativeCorrelationData_wSLC16A7.columns[i]]) 
    axes[i].plot(brainMapFile,np.poly1d(np.polyfit(brainMapFile,np.array(
    negativeCorrelationData_wSLC16A7.iloc[:,i]),1))(brainMapFile))
    axes[i].set_ylabel(negativeCorrelationData_wSLC16A7.columns[i],fontsize=16)
    axes[2].set_xlabel('Lactate-to-Bicarbonate ratio',fontsize=16)
    axes[i].tick_params(axis='x', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12)

plt.tight_layout()
plt.savefig('MostNegativeLB.pdf')    
# Only 3 significant correlations, so need for both negative and positive correlation plots

